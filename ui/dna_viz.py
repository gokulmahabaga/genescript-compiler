"""
dna_viz.py -- charts for the "DNA Visualizer" tab of the GeneScript UI.

Everything shown here comes from the compiler's own results:
  * the sequences come from the symbol table (`gsc --dump-symtab`),
  * motif hits come from the compiled program's FIND_MOTIF output,
  * mutations come from the compiled program's COMPARE output,
matched up with the FIND_MOTIF / COMPARE instructions in the TAC
(`gsc --dump-tac`) so each result is tied to the sequences it used.
Only simple counting (base composition, a sliding GC window, codon
grouping for display) is done here in Python.

Colors: one fixed color per base, validated with the dataviz palette
checker (light mode, all pairs): A blue, T yellow, G aqua, C red.
Every chart also prints the base letter or a legend, so color is never
the only way to tell bases apart.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import altair as alt
import pandas as pd

BASES = ["A", "T", "G", "C"]
BASE_COLORS = ["#2a78d6", "#eda100", "#1baf7a", "#e34948"]
# Ink for the letter printed on each colored cell (chosen for contrast).
BASE_TEXT = {"A": "#ffffff", "T": "#0b0b0b", "G": "#0b0b0b", "C": "#ffffff"}
INK = "#0b0b0b"
MUTED = "#52514e"
SINGLE = "#2a78d6"  # single-series charts use one hue
SURFACE = "#ffffff"

BASE_SCALE = alt.Scale(domain=BASES, range=BASE_COLORS)
ROW_WIDTH = 50  # bases per row in the sequence map

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# --------------------------------------------------------------------------
# Pulling the facts out of the compiler's output
# --------------------------------------------------------------------------

@dataclass
class MotifResult:
    sequence: str          # name of the declared sequence searched
    motif: str
    positions: list[int]   # 1-indexed, as printed by the runtime


@dataclass
class Mutation:
    position: int
    reference: str
    sample: str
    kind: str              # Substitution / Insertion / Deletion


@dataclass
class CompareResult:
    reference: str
    sample: str
    mutations: list[Mutation] = field(default_factory=list)
    count: int = 0
    similarity: float = 100.0


def declared_sequences(symtab: pd.DataFrame) -> dict[str, str]:
    """name -> DNA string, for every SEQUENCE declaration in the symbol table."""
    seqs: dict[str, str] = {}
    if symtab.empty:
        return seqs
    for _, row in symtab.iterrows():
        src = str(row["Value / source"])
        if row["Type"] == "SEQUENCE" and src.startswith('"') and src.endswith('"'):
            seqs[row["Name"]] = src[1:-1].upper()
    return seqs


def _tac_ops(tac: str, op: str) -> list[tuple[str, str]]:
    ops = []
    for line in tac.splitlines():
        m = re.match(rf"\s*{op} (\S+), (.+)$", line)
        if m:
            ops.append((m.group(1), m.group(2).strip()))
    return ops


def motif_results(tac: str, output: str) -> list[MotifResult]:
    """Pair each FIND_MOTIF instruction with its 'Motif:' block in the output."""
    blocks: list[tuple[str, list[int]]] = []
    lines = output.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "Motif:":
            motif, hits = "", []
            i += 1
            while i < len(lines):
                m = re.match(r"(\S+) -> (?:position (\d+)|not found)$", lines[i].strip())
                if not m:
                    break
                motif = m.group(1)
                if m.group(2):
                    hits.append(int(m.group(2)))
                i += 1
            blocks.append((motif, hits))
            continue
        i += 1

    results = []
    for (seq_name, motif_arg), (motif, hits) in zip(_tac_ops(tac, "FIND_MOTIF"), blocks):
        results.append(MotifResult(seq_name, motif or motif_arg.strip('"'), hits))
    return results


def compare_results(tac: str, output: str) -> list[CompareResult]:
    """Pair each COMPARE instruction with its 'Mutation Report' in the output."""
    reports: list[CompareResult] = []
    current: CompareResult | None = None
    pending: dict[str, str] = {}
    for raw in output.splitlines():
        line = raw.strip()
        if line == "Mutation Report":
            current = CompareResult("", "")
            reports.append(current)
            pending = {}
            continue
        if current is None:
            continue
        m = re.match(r"(Position|Reference|Sample|Type)\s+(\S+)$", line)
        if m:
            pending[m.group(1)] = m.group(2)
            if m.group(1) == "Type":
                current.mutations.append(Mutation(int(pending["Position"]), pending["Reference"],
                                                  pending["Sample"], pending["Type"]))
                pending = {}
            continue
        m = re.match(r"Mutations: (\d+)\s+Similarity: ([\d.]+)%", line)
        if m:
            current.count, current.similarity = int(m.group(1)), float(m.group(2))
            current = None

    for (ref, smp), rep in zip(_tac_ops(tac, "COMPARE"), reports):
        rep.reference, rep.sample = ref, smp
    return reports


# --------------------------------------------------------------------------
# Numbers
# --------------------------------------------------------------------------

def composition(seq: str) -> pd.DataFrame:
    n = len(seq) or 1
    rows = [{"Base": b, "Count": seq.count(b), "Percent": 100 * seq.count(b) / n} for b in BASES]
    df = pd.DataFrame(rows)
    df["Label"] = df.apply(lambda r: f"{r.Count}  ({r.Percent:.1f}%)", axis=1)
    return df


def gc_percent(seq: str) -> float:
    return 100 * sum(seq.count(b) for b in "GC") / len(seq) if seq else 0.0


def gc_windows(seq: str, window: int) -> pd.DataFrame:
    rows = []
    for start in range(0, len(seq) - window + 1):
        chunk = seq[start:start + window]
        rows.append({"Start": start + 1, "End": start + window,
                     "Centre": start + (window + 1) / 2, "GC %": gc_percent(chunk)})
    return pd.DataFrame(rows)


def codons(seq: str) -> list[tuple[str, str]]:
    """(codon, amino acid) pairs, stopping after a stop codon -- same rule as the runtime."""
    out = []
    for i in range(0, len(seq) - 2, 3):
        c = seq[i:i + 3]
        aa = CODON_TABLE.get(c, "X")
        out.append((c, aa))
        if aa == "*":
            break
    return out


# --------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------

def _base_legend(title: str = "Base") -> alt.Legend:
    return alt.Legend(title=title, orient="top", direction="horizontal", symbolType="square")


def composition_chart(seq: str) -> alt.Chart:
    df = composition(seq)
    base = alt.Chart(df).encode(
        y=alt.Y("Base:N", sort=BASES, title=None, scale=alt.Scale(paddingInner=0.3),
                axis=alt.Axis(labelFontSize=14, labelFontWeight="bold", ticks=False, domain=False)),
    )
    bars = base.mark_bar(cornerRadiusEnd=4).encode(
        x=alt.X("Count:Q", title="Number of bases", axis=alt.Axis(tickMinStep=1, grid=True)),
        color=alt.Color("Base:N", scale=BASE_SCALE, legend=None),
        tooltip=[alt.Tooltip("Base:N"), alt.Tooltip("Count:Q"),
                 alt.Tooltip("Percent:Q", format=".1f", title="% of sequence")],
    )
    labels = base.mark_text(align="left", dx=6, color=MUTED, fontSize=12).encode(
        x="Count:Q", text="Label:N")
    return (bars + labels).properties(height=alt.Step(36))


def sequence_map(seq: str, highlight: MotifResult | None = None) -> alt.Chart:
    """Every base as a colored cell, wrapped into rows of ROW_WIDTH.
    With a motif selected, its hits stay bright and everything else fades."""
    in_motif = set()
    if highlight:
        for p in highlight.positions:
            in_motif.update(range(p, p + len(highlight.motif)))
    rows = []
    for i, b in enumerate(seq):
        pos = i + 1
        rows.append({"Position": pos, "Base": b, "Row": i // ROW_WIDTH, "Col": i % ROW_WIDTH,
                     "RowLabel": f"{(i // ROW_WIDTH) * ROW_WIDTH + 1}",
                     "Motif": (highlight.motif if pos in in_motif else "") if highlight else "",
                     "Lit": (pos in in_motif) if highlight else True})
    for r in rows:  # faded cells get dark letters so A and C stay readable
        r["Ink"] = BASE_TEXT.get(r["Base"], INK) if r["Lit"] else MUTED
    df = pd.DataFrame(rows)
    show_letters = True

    enc = dict(
        x=alt.X("Col:O", axis=None),
        y=alt.Y("RowLabel:O", sort=None, title=None,
                axis=alt.Axis(ticks=False, domain=False, labelColor=MUTED, labelPadding=6)),
    )
    tooltip = [alt.Tooltip("Position:Q"), alt.Tooltip("Base:N")]
    if highlight:
        tooltip.append(alt.Tooltip("Motif:N", title="In motif"))
    cells = alt.Chart(df).mark_rect(cornerRadius=3, stroke=SURFACE, strokeWidth=2).encode(
        **enc,
        color=alt.Color("Base:N", scale=BASE_SCALE, legend=_base_legend()),
        opacity=alt.condition("datum.Lit", alt.value(1.0), alt.value(0.18)),
        tooltip=tooltip,
    )
    layers = cells
    if show_letters:
        letters = alt.Chart(df).mark_text(fontSize=11, fontWeight="bold").encode(
            **enc, text="Base:N",
            color=alt.Color("Ink:N", scale=None),
            opacity=alt.condition("datum.Lit", alt.value(1.0), alt.value(0.75)),
            tooltip=tooltip,
        )
        layers = cells + letters
    return layers.properties(height=alt.Step(28))


def gc_window_chart(seq: str, window: int) -> alt.Chart:
    df = gc_windows(seq, window)
    overall = gc_percent(seq)
    line = alt.Chart(df).mark_line(strokeWidth=2, color=SINGLE, interpolate="monotone").encode(
        x=alt.X("Centre:Q", title="Position along the sequence",
                scale=alt.Scale(domain=[1, len(seq)]), axis=alt.Axis(tickMinStep=1)),
        y=alt.Y("GC %:Q", title="GC content (%)", scale=alt.Scale(domain=[0, 100])),
    )
    hover = alt.selection_point(fields=["Centre"], nearest=True, on="pointerover", empty=False)
    points = alt.Chart(df).mark_circle(size=70, color=SINGLE, stroke=SURFACE, strokeWidth=2).encode(
        x="Centre:Q", y="GC %:Q",
        opacity=alt.condition(hover, alt.value(1), alt.value(0)),
        tooltip=[alt.Tooltip("Start:Q", title="From base"), alt.Tooltip("End:Q", title="To base"),
                 alt.Tooltip("GC %:Q", format=".1f")],
    ).add_params(hover)
    avg = pd.DataFrame({"y": [overall], "label": [f"whole sequence {overall:.2f}%"]})
    rule = alt.Chart(avg).mark_rule(strokeDash=[4, 4], color=MUTED, strokeWidth=1).encode(y="y:Q")
    rule_label = alt.Chart(avg).mark_text(align="left", dx=4, dy=-7, color=MUTED, fontSize=11).encode(
        y="y:Q", x=alt.value(0), text="label:N")
    return (rule + rule_label + line + points).properties(height=240)


def gc_by_sequence_chart(seqs: dict[str, str]) -> alt.Chart:
    df = pd.DataFrame([{"Sequence": k, "GC %": gc_percent(v), "Length": len(v)} for k, v in seqs.items()])
    df["Label"] = df["GC %"].map(lambda v: f"{v:.1f}%")
    base = alt.Chart(df).encode(y=alt.Y("Sequence:N", sort="-x", title=None, scale=alt.Scale(paddingInner=0.3),
                                        axis=alt.Axis(ticks=False, domain=False)))
    bars = base.mark_bar(cornerRadiusEnd=4, color=SINGLE).encode(
        x=alt.X("GC %:Q", title="GC content (%)", scale=alt.Scale(domain=[0, 100])),
        tooltip=["Sequence:N", alt.Tooltip("GC %:Q", format=".2f"), "Length:Q"],
    )
    labels = base.mark_text(align="left", dx=6, color=MUTED).encode(x="GC %:Q", text="Label:N")
    return (bars + labels).properties(height=alt.Step(36))


def mutation_chart(ref: str, smp: str, result: CompareResult) -> alt.Chart:
    """Reference and sample stacked base by base; only the mutated positions stay bright."""
    mutated = {m.position: m for m in result.mutations}
    length = max(len(ref), len(smp))
    rows = []
    for track, s in (("Reference", ref), ("Sample", smp)):
        for pos in range(1, length + 1):
            b = s[pos - 1] if pos <= len(s) else "-"
            m = mutated.get(pos)
            rows.append({"Track": track, "Position": pos, "Base": b,
                         "Lit": m is not None, "Mutation": m.kind if m else "",
                         "Ink": BASE_TEXT.get(b, MUTED) if m is not None else MUTED})
    df = pd.DataFrame(rows)
    enc = dict(x=alt.X("Position:O", title="Position", axis=alt.Axis(labelAngle=0, ticks=False, domain=False,
                                                                      labelColor=MUTED,
                                                                      values=[p for p in range(1, length + 1)
                                                                              if p == 1 or p % 5 == 0 or p in mutated])),
               y=alt.Y("Track:N", sort=["Reference", "Sample"], title=None,
                       axis=alt.Axis(ticks=False, domain=False, labelFontWeight="bold")))
    tooltip = ["Track:N", "Position:Q", "Base:N", alt.Tooltip("Mutation:N", title="Mutation")]
    real = df[df["Base"] != "-"]
    gaps = df[df["Base"] == "-"]
    cells = alt.Chart(real).mark_rect(cornerRadius=3, stroke=SURFACE, strokeWidth=2).encode(
        **enc, color=alt.Color("Base:N", scale=BASE_SCALE, legend=_base_legend()),
        opacity=alt.condition("datum.Lit", alt.value(1.0), alt.value(0.25)), tooltip=tooltip)
    gap_cells = alt.Chart(gaps).mark_rect(cornerRadius=3, stroke=MUTED, strokeDash=[2, 2], strokeWidth=1,
                                          fill=SURFACE).encode(**enc, tooltip=tooltip)
    letters = alt.Chart(df).mark_text(fontSize=11, fontWeight="bold").encode(
        **enc, text="Base:N", color=alt.Color("Ink:N", scale=None),
        opacity=alt.condition("datum.Lit", alt.value(1.0), alt.value(0.75)), tooltip=tooltip)
    marks = pd.DataFrame([{"Position": p, "Track": "Reference", "Label": "▼"} for p in mutated])
    markers = alt.Chart(marks).mark_text(dy=-22, fontSize=12, color=INK).encode(
        x=enc["x"], y=enc["y"], text="Label:N") if len(marks) else None
    chart = cells + gap_cells + letters
    if markers is not None:
        chart = chart + markers
    return chart.properties(height=alt.Step(32))


def codon_html(seq: str) -> str:
    """Codons in boxes with the amino acid underneath (stop codon marked)."""
    items = []
    for i, (c, aa) in enumerate(codons(seq)):
        stop = aa == "*"
        items.append(
            f'<div class="codon{" stop" if stop else ""}" title="codon {i + 1}: {c} → '
            f'{"STOP" if stop else aa}"><div class="nt">{c}</div><div class="aa">{"STOP" if stop else aa}</div></div>'
        )
    return '<div class="codons">' + "".join(items) + "</div>"
