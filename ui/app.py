"""
GeneScript Compiler -- Streamlit UI.

Type (or load) a GeneScript program, compile it, and see what every
phase of the compiler produced: tokens, AST, symbol table, three-address
code before/after optimization, the optimization report, LLVM IR,
assembly, and the program's output.

This file never touches the compiler's internals: it just runs the real
`gsc` binary (built with `make`) as a subprocess with its --dump-* flags
and lays out what comes back.

Run from the repo root:
    streamlit run ui/app.py
"""
from __future__ import annotations

import re
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import streamlit as st
from ml.features import extract_features
from ml.predict import predict_sequence
from ml.validation import DNAValidationError, validate_sequence

REPO = Path(__file__).resolve().parent.parent
GSC = REPO / "gsc"

PHASES = ["Lexer", "Parser", "Semantic analysis", "Optimizer", "LLVM codegen", "Run"]

DEFAULT_PROGRAM = """\
SEQUENCE dna = "ATGCGATCGATCG";

PRINT LENGTH(dna);
PRINT GC_CONTENT(dna);

// the optimizer computes this only once
x = GC_CONTENT(dna);
y = GC_CONTENT(dna);
PRINT y;

// ...and rewrites this into one REVERSE_COMPLEMENT call
PRINT REVERSE(COMPLEMENT(dna));

FIND_MOTIF(dna, "ATG");
TRANSLATE(dna);
"""


# --------------------------------------------------------------------------
# Running the compiler
# --------------------------------------------------------------------------

@dataclass
class Result:
    ok: bool = True
    error: str = ""
    failed_phase: str | None = None
    sections: dict[str, str] = field(default_factory=dict)
    llvm_ir: str = ""
    assembly: str = ""
    output: str = ""
    exit_code: int = 0


def split_sections(text: str) -> dict[str, str]:
    """gsc prints each dump under a '== Title ==' header."""
    sections: dict[str, str] = {}
    current = None
    lines: list[str] = []
    for line in text.splitlines():
        m = re.fullmatch(r"== (.+) ==", line)
        if m:
            if current:
                sections[current] = "\n".join(lines).rstrip()
            current, lines = m.group(1), []
        elif current:
            lines.append(line)
    if current:
        sections[current] = "\n".join(lines).rstrip()
    return sections


def phase_of_error(err: str) -> str:
    if "Lex Error" in err:
        return "Lexer"
    if "Parse Error" in err:
        return "Parser"
    if "Semantic Error" in err:
        return "Semantic analysis"
    return "LLVM codegen"


def clean_stderr(err: str) -> str:
    return "\n".join(l for l in err.splitlines() if not l.startswith("[gsc]")).strip()


def compile_program(source: str, frontend: str, optimize: bool) -> Result:
    res = Result()
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "program.gs"
        src.write_text(source)
        common = [str(GSC), str(src), f"--frontend={frontend}"]
        if not optimize:
            common.append("--no-opt")

        # 1) phase dumps only (no codegen output needed)
        dumps = subprocess.run(
            common + ["--dump-tokens", "--dump-ast", "--dump-tac",
                      "--show-opt", "--dump-symtab", "--no-run"],
            cwd=tmp, capture_output=True, text=True, timeout=60,
        )
        res.sections = split_sections(dumps.stdout)
        if dumps.returncode != 0:
            res.ok = False
            res.error = clean_stderr(dumps.stderr) or "Compilation failed."
            res.failed_phase = phase_of_error(res.error)
            return res

        # 2) compile + run for real; gsc leaves program.ll / program.s behind
        run = subprocess.run(common, cwd=tmp, capture_output=True, text=True, timeout=60)
        ll, asm = Path(tmp) / "program.ll", Path(tmp) / "program.s"
        res.llvm_ir = ll.read_text() if ll.exists() else ""
        res.assembly = asm.read_text() if asm.exists() else ""
        res.output = run.stdout.lstrip("\n")
        res.exit_code = run.returncode
        if run.returncode != 0:
            res.ok = False
            res.error = clean_stderr(run.stderr) or f"Program exited with code {run.returncode}."
            res.failed_phase = "Run" if res.llvm_ir else phase_of_error(res.error)
    return res


# --------------------------------------------------------------------------
# Parsing the dumps into tables / numbers
# --------------------------------------------------------------------------

def tokens_table(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines()[1:]:  # skip the header row
        parts = line.split(None, 2)
        if len(parts) == 3:
            rows.append({"Line": int(parts[0]), "Token": parts[1], "Lexeme": parts[2]})
    return pd.DataFrame(rows)


def symtab_table(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines()[1:]:
        parts = line.split(None, 5)
        if len(parts) == 6:
            name, typ, length, gc, ln, src = parts
            rows.append({"Name": name, "Type": typ, "Length": length, "GC %": gc,
                         "Line": int(ln), "Value / source": src})
    return pd.DataFrame(rows)


def opt_summary(report: str) -> tuple[int, int, int]:
    m = re.search(r"summary: (\d+) domain rewrite\(s\), (\d+) reused value\(s\), (\d+) temporary", report)
    return tuple(int(g) for g in m.groups()) if m else (0, 0, 0)


def opt_rows(report: str) -> pd.DataFrame:
    rows = []
    for line in report.splitlines():
        m = re.match(r"\s*line (\d+)\s+\[(\w+)\] (.*)", line)
        if m:
            kind = "Domain rewrite" if m.group(2) == "domain" else "CSE"
            rows.append({"Line": int(m.group(1)), "Kind": kind, "Change": m.group(3)})
    return pd.DataFrame(rows)


def count_instructions(tac: str) -> int:
    return sum(1 for l in tac.splitlines() if l.strip() and l.strip() != "HALT")


def count_calls(tac: str) -> int:
    """Instructions that call a DNA operation in the runtime -- the real work.
    (Copies like `y = x`, PRINT and LOAD_SEQ are cheap and not counted.)"""
    n = 0
    for l in tac.splitlines():
        l = l.strip()
        if l.startswith("FIND_MOTIF") or re.match(r"\S+ = (?!LOAD_SEQ)[A-Z_]+ ", l):
            n += 1
    return n


# --------------------------------------------------------------------------
# Page
# --------------------------------------------------------------------------

st.set_page_config(page_title="GeneScript Compiler", page_icon="🧬", layout="wide")

st.markdown(
    """
    <style>
      .block-container { padding-top: 2rem; }
      .pipeline { display:flex; flex-wrap:wrap; gap:.4rem; align-items:center; margin:.25rem 0 1rem; }
      .step { padding:.3rem .7rem; border-radius:999px; font-size:.85rem; font-weight:600;
              border:1px solid rgba(128,128,128,.35); }
      .step.done { background:rgba(34,160,90,.15); border-color:rgba(34,160,90,.6); }
      .step.fail { background:rgba(220,60,60,.18); border-color:rgba(220,60,60,.8); }
      .step.skip { opacity:.45; }
      .arrow { opacity:.5; }
      .stTextArea textarea { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
                             font-size: .9rem; line-height: 1.5; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧬 GeneScript Compiler")
st.caption("A DSL for DNA sequence analysis, compiled through Flex/Bison → AST → optimizer → LLVM → native code.")

with st.expander("ML Classification: Promoter vs Non-Promoter", expanded=False):
    st.caption("Classify a DNA sequence with the separately trained Random Forest model.")
    ml_sequence = st.text_area("DNA sequence (A/C/G/T)", key="ml_sequence", height=110)
    if st.button("Predict", key="ml_predict", type="secondary"):
        try:
            normalized = validate_sequence(ml_sequence)
            features = extract_features(normalized)
            prediction = predict_sequence(normalized)
            st.session_state.ml_prediction = prediction
            st.session_state.ml_sequence_info = (len(normalized), features[4])
        except DNAValidationError as exc:
            st.session_state.pop("ml_prediction", None)
            st.error(str(exc))
        except FileNotFoundError as exc:
            st.error(str(exc))
    if st.session_state.get("ml_prediction"):
        prediction = st.session_state.ml_prediction
        length, gc_content = st.session_state.ml_sequence_info
        st.success(f"Predicted class: **{prediction['prediction']}**")
        if prediction["confidence"] is not None:
            st.metric("Model probability for predicted class", f"{prediction['confidence']:.1%}")
        c_len, c_gc = st.columns(2)
        c_len.metric("Sequence length", length)
        c_gc.metric("GC content", f"{gc_content:.1%}")

if not GSC.exists():
    st.error(f"Compiler not found at `{GSC}`. Run `make` in the repo root first, then reload this page.")
    st.stop()

examples = sorted(list((REPO / "examples").glob("*.gs")) + list((REPO / "tests" / "cases").glob("*.gs")))

left, right = st.columns([2, 3], gap="large")

with left:
    choice = st.selectbox(
        "Load an example",
        ["(demo program)"] + [str(p.relative_to(REPO)) for p in examples],
    )
    if st.session_state.get("loaded") != choice:
        st.session_state.loaded = choice
        st.session_state.source = (
            DEFAULT_PROGRAM if choice == "(demo program)" else (REPO / choice).read_text()
        )
        st.session_state.pop("result", None)

    source = st.text_area("GeneScript source", key="source", height=380)

    c1, c2 = st.columns(2)
    frontend_label = c1.radio("Parser", ["Bison (LALR)", "Recursive descent"], horizontal=False)
    optimize = c2.toggle("Optimizer", value=True, help="Off = compile with --no-opt")
    frontend = "bison" if frontend_label.startswith("Bison") else "rd"

    if st.button("Compile & run", type="primary", width="stretch"):
        with st.spinner("Compiling…"):
            st.session_state.result = compile_program(source, frontend, optimize)
            st.session_state.result_opt = optimize

res: Result | None = st.session_state.get("result")

with right:
    if res is None:
        st.info("Write a program (or load an example) and press **Compile & run**.")
        st.stop()

    # ---- pipeline status ----
    failed_at = PHASES.index(res.failed_phase) if res.failed_phase else len(PHASES)
    optimized = st.session_state.get("result_opt", True)
    steps = []
    for i, name in enumerate(PHASES):
        if name == "Optimizer" and not optimized:
            cls, mark = "skip", "–"
        elif i < failed_at:
            cls, mark = "done", "✓"
        elif i == failed_at:
            cls, mark = "fail", "✗"
        else:
            cls, mark = "skip", "·"
        steps.append(f'<span class="step {cls}">{mark} {name}</span>')
    st.markdown('<div class="pipeline">' + '<span class="arrow">→</span>'.join(steps) + "</div>",
                unsafe_allow_html=True)

    if not res.ok:
        st.error(f"**{res.failed_phase} failed**\n\n{res.error}")

    # ---- headline numbers ----
    # Output of phases after the failing one is not meaningful: hide it.
    stopped_early = res.failed_phase in ("Lexer", "Parser", "Semantic analysis")
    if stopped_early:
        for later in ("TAC", "TAC (before optimization)", "TAC (after optimization)",
                      "Optimization report", "Symbol table"):
            res.sections.pop(later, None)
    tac_before = res.sections.get("TAC (before optimization)", res.sections.get("TAC", ""))
    tac_after = res.sections.get("TAC (after optimization)", tac_before)
    report = res.sections.get("Optimization report", "")
    domain, reused, temps = opt_summary(report)
    toks = tokens_table(res.sections.get("Tokens", ""))

    n_before, n_after = count_instructions(tac_before), count_instructions(tac_after)
    c_before, c_after = count_calls(tac_before), count_calls(tac_after)
    if res.ok:  # a half-finished compile's numbers would only mislead
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Tokens", len(toks))
        m2.metric("DNA operations run", c_after,
                  delta=(c_after - c_before) if optimized and c_after != c_before else None,
                  delta_color="inverse",
                  help="Calls into the runtime (GC_CONTENT, REVERSE, ...) in the final code. "
                       "The delta is what the optimizer removed.")
        m3.metric("Domain rewrites", domain)
        m4.metric("Values reused (CSE)", reused)

    tabs = st.tabs(["▶ Output", "Tokens", "AST", "Symbol table", "TAC & optimizer", "LLVM IR", "Assembly"])

    with tabs[0]:
        if res.output:
            st.code(res.output, language=None)
        elif res.ok:
            st.caption("The program printed nothing.")
        else:
            st.caption("No output — compilation stopped at the step marked ✗ above.")

    with tabs[1]:
        if len(toks):
            st.dataframe(toks, hide_index=True, width="stretch")
        else:
            st.caption("No tokens (the lexer stopped at an error).")

    with tabs[2]:
        if "AST" in res.sections:
            st.code(res.sections["AST"], language=None)
        else:
            st.caption("No AST — parsing did not finish.")

    with tabs[3]:
        sym = symtab_table(res.sections.get("Symbol table", ""))
        if len(sym):
            st.dataframe(sym, hide_index=True, width="stretch")
            st.caption("Names starting with `$cse` are temporaries the optimizer created.")
        else:
            st.caption("No symbol table.")

    with tabs[4]:
        if not tac_before:
            st.caption("No intermediate code — compilation stopped earlier.")
        elif optimized:
            a, b = st.columns(2)
            a.markdown(f"**Before optimization** · {n_before} instructions, {c_before} DNA operations")
            a.code(tac_before, language=None)
            b.markdown(f"**After optimization** · {n_after} instructions, {c_after} DNA operations")
            b.code(tac_after, language=None)
            rows = opt_rows(report)
            st.markdown("**What the optimizer changed**")
            if len(rows):
                st.dataframe(rows, hide_index=True, width="stretch")
            else:
                st.caption("Nothing to optimize in this program.")
            if temps:
                st.caption(f"{temps} compiler temporar{'y' if temps == 1 else 'ies'} created to hold shared values.")
        else:
            st.markdown(f"**Three-address code** · {n_before} instructions (optimizer off)")
            st.code(tac_before, language=None)

    with tabs[5]:
        if res.llvm_ir:
            st.code(res.llvm_ir, language="llvm")
        else:
            st.caption("No LLVM IR — compilation stopped earlier.")

    with tabs[6]:
        if res.assembly:
            st.code(res.assembly, language="asm")
        else:
            st.caption("No assembly — compilation stopped earlier.")
