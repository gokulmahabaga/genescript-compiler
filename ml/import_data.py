"""Import external promoter data into the training set, without committing it.

The repo ships only the UCI dataset. This tool turns data *you* download
(RegulonDB promoters, any FASTA) into `label,name,sequence` rows via
`ml.add_data`, so third-party data -- and its license terms -- stay on your
machine. Don't commit the downloaded files or a grown `promoters.data`
built from them; see the README for the licensing note.

Two sources:

  # Promoters from a RegulonDB-style TSV + negatives from the E. coli genome.
  # Positives are cut from the genome around each TSS (or read from a sequence
  # column if the file has one); negatives are random genome windows of the
  # SAME lengths that stay away from every known TSS -- so length can't leak
  # the label.
  python -m ml.import_data regulondb --tsv PromoterSet.tsv --genome U00096.3.fasta --dry-run
  python -m ml.import_data regulondb --tsv PromoterSet.tsv --genome U00096.3.fasta --retrain

  # Any FASTA, one label for every record
  python -m ml.import_data fasta --fasta extra.fasta --label Non-Promoter

`--dry-run` validates and reports without touching the dataset. Use
`--dataset scratch.data --model scratch.joblib` to experiment on a copy.
"""

import argparse
import bisect
import json
import random
from pathlib import Path

from .add_data import add_sequences
from .dataset import DEFAULT_DATASET, MIN_SEQUENCE_LENGTH, load_dataset, warn_if_skewed
from .validation import DNAValidationError, validate_sequence

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def read_fasta(path):
    """Return [(header, sequence)] for every record in a FASTA file."""
    records, header, chunks = [], None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip() or f"record_{len(records) + 1}", []
        elif header is not None:
            chunks.append(line)
        else:
            raise ValueError(f"{path} is not FASTA: sequence data before the first '>' header")
    if header is not None:
        records.append((header, "".join(chunks)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def _read_table(path, has_header=True):
    """Read a TSV, skipping '#' comment lines. Returns (header_or_None, rows)."""
    lines = [ln.rstrip("\r\n") for ln in Path(path).read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"{path} has no data rows")
    rows = [ln.split("\t") for ln in lines]
    return (rows[0], rows[1:]) if has_header else (None, rows)


def _column(spec, header, what):
    """Resolve a 1-based index or a (case-insensitive substring of a) header name."""
    if spec is None:
        return None
    if str(spec).isdigit():
        return int(spec) - 1
    if header is None:
        raise ValueError(f"--{what}-col {spec!r} is a name but there is no header row")
    matches = [i for i, name in enumerate(header) if str(spec).lower() in name.lower()]
    if not matches:
        raise ValueError(f"No column matching {spec!r} for {what}; header is {header}")
    return matches[0]


def _auto_column(header, candidates):
    if header is None:
        return None
    for cand in candidates:
        for i, name in enumerate(header):
            if cand in name.lower():
                return i
    return None


def _strand(value):
    v = str(value).strip().lower()
    if v in {"+", "forward", "fwd", "plus"}:
        return "+"
    if v in {"-", "reverse", "rev", "minus"}:
        return "-"
    raise ValueError(f"unrecognized strand {value!r}")


def read_promoter_table(path, seq_col=None, tss_col=None, strand_col=None, has_header=True):
    """Parse a promoter TSV into [{"sequence"?, "tss"?, "strand"?, "name"}].

    Columns are given by 1-based index or header-name substring; when omitted
    they're auto-detected from the header ("sequence", "tss"/"position",
    "strand"). A row needs either a sequence, or a TSS + strand (to be cut
    from the genome).
    """
    header, rows = _read_table(path, has_header)
    s = _column(seq_col, header, "seq") if seq_col else _auto_column(header, ["sequence", "seq"])
    t = _column(tss_col, header, "tss") if tss_col else _auto_column(header, ["tss", "position", "pos"])
    d = _column(strand_col, header, "strand") if strand_col else _auto_column(header, ["strand"])
    if s is None and (t is None or d is None):
        raise ValueError(
            f"Can't find promoter data in {path}: need a sequence column, or TSS + strand columns. "
            f"Header: {header}. Pass --seq-col / --tss-col / --strand-col."
        )
    records = []
    for n, row in enumerate(rows, 1):
        def cell(i):
            return row[i].strip() if i is not None and i < len(row) else ""
        rec = {"name": f"regulondb_{n}", "sequence": cell(s), "tss": None, "strand": None}
        try:
            if cell(t):
                rec["tss"] = int(float(cell(t)))
            if cell(d):
                rec["strand"] = _strand(cell(d))
        except ValueError:
            pass  # row without usable coordinates; kept only if it has a sequence
        if rec["sequence"] or (rec["tss"] and rec["strand"]):
            records.append(rec)
    return records


def window_at(genome, tss, strand, upstream, downstream):
    """Bases `upstream` before and `downstream` from (incl.) a 1-based TSS, on
    `strand`, read 5'->3'. Returns None if the window runs off the genome."""
    if strand == "+":
        start, end = tss - 1 - upstream, tss - 1 + downstream
    else:
        start, end = tss - downstream, tss + upstream
    if start < 0 or end > len(genome):
        return None
    window = genome[start:end]
    return window if strand == "+" else reverse_complement(window)


def sample_negatives(genome, lengths, tss_positions, exclude=200, seed=42, max_tries_per=200):
    """One random genome window per entry of `lengths`, each either strand,
    >= `exclude` bp from every TSS. Same lengths as the positives => length
    carries no label information. Returns [(start, strand, sequence)]."""
    rng = random.Random(seed)
    tss_sorted = sorted(tss_positions)

    def near_tss(start, end):
        i = bisect.bisect_left(tss_sorted, start - exclude)
        return i < len(tss_sorted) and tss_sorted[i] <= end + exclude

    out, seen = [], set()
    for length in lengths:
        if length > len(genome):
            continue
        for _ in range(max_tries_per):
            start = rng.randrange(0, len(genome) - length + 1)
            if (start, length) in seen or near_tss(start + 1, start + length):
                continue
            window = genome[start:start + length]
            if set(window) - set("ACGT"):
                continue
            strand = rng.choice("+-")
            seen.add((start, length))
            out.append((start + 1, strand, window if strand == "+" else reverse_complement(window)))
            break
    return out


def build_regulondb_entries(tsv, genome_fasta, seq_col=None, tss_col=None, strand_col=None,
                            has_header=True, upstream=50, downstream=7, neg_ratio=1.0,
                            exclude=200, seed=42):
    """Return (entries, notes) ready for `add_sequences`."""
    genome = read_fasta(genome_fasta)[0][1].upper() if genome_fasta else None
    records = read_promoter_table(tsv, seq_col, tss_col, strand_col, has_header)
    notes, entries, lengths = [], [], []
    for rec in records:
        seq = rec["sequence"]
        if not seq and genome is not None:
            seq = window_at(genome, rec["tss"], rec["strand"], upstream, downstream) or ""
        if not seq:
            continue
        try:
            seq = validate_sequence(seq)
        except DNAValidationError:
            continue
        if len(seq) < MIN_SEQUENCE_LENGTH:
            continue
        entries.append({"label": "Promoter", "sequence": seq, "name": rec["name"]})
        lengths.append(len(seq))
    if not entries:
        raise ValueError("No usable promoter rows found (need a sequence, or TSS+strand and --genome).")
    if neg_ratio > 0:
        if genome is None:
            raise ValueError("--genome is required to generate negatives (or pass --neg-ratio 0).")
        wanted = [lengths[i % len(lengths)] for i in range(round(len(lengths) * neg_ratio))]
        tss = [r["tss"] for r in records if r["tss"]]
        if not tss:
            notes.append("no TSS coordinates in the file: negatives can't avoid known promoters, "
                         "so some may overlap real ones")
        negs = sample_negatives(genome, wanted, tss, exclude=exclude, seed=seed)
        for start, strand, seq in negs:
            entries.append({"label": "Non-Promoter", "sequence": seq,
                            "name": f"genome_neg_{start}{strand}"})
        if len(negs) < len(wanted):
            notes.append(f"only {len(negs)} of {len(wanted)} negatives could be placed")
    return entries, notes


def build_fasta_entries(fasta, label):
    return [{"label": label, "sequence": seq, "name": header.split()[0]}
            for header, seq in read_fasta(fasta)]


def _finish(entries, notes, args):
    from .add_data import _print_report  # shared report format
    for note in notes:
        print(f"note: {note}")
    report = add_sequences(entries, dataset_path=args.dataset, dry_run=args.dry_run)
    _print_report(report)
    if args.dry_run:
        print("\nDry run -- the dataset was not modified.")
        return
    if not report["added"]:
        return
    warn_if_skewed(*load_dataset(args.dataset))
    if args.retrain:
        from .train import MODEL_PATH, train
        model_path = args.model or MODEL_PATH
        print(f"\nRetraining on the updated dataset -> {model_path} ...")
        print(json.dumps(train(dataset_path=args.dataset, model_path=model_path), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--dataset", default=str(DEFAULT_DATASET), help="dataset file to append to")
    common.add_argument("--dry-run", action="store_true", help="validate and report; write nothing")
    common.add_argument("--retrain", action="store_true", help="retrain after importing")
    common.add_argument("--model", help="where --retrain writes the model (default: the production path)")
    sub = parser.add_subparsers(dest="source", required=True)

    reg = sub.add_parser("regulondb", parents=[common], help="promoter TSV + genome-sampled negatives")
    reg.add_argument("--tsv", required=True, help="promoter table (tab-separated, '#' comments allowed)")
    reg.add_argument("--genome", help="genome FASTA (e.g. E. coli K-12 MG1655, U00096.3); "
                                       "needed for negatives and for TSV rows without a sequence")
    reg.add_argument("--seq-col", help="sequence column (1-based index or header text)")
    reg.add_argument("--tss-col", help="TSS position column")
    reg.add_argument("--strand-col", help="strand column")
    reg.add_argument("--no-header", action="store_true", help="the table has no header row")
    reg.add_argument("--upstream", type=int, default=50, help="bases before the TSS when cutting from the genome")
    reg.add_argument("--downstream", type=int, default=7, help="bases from the TSS on (incl.) when cutting from the genome")
    reg.add_argument("--neg-ratio", type=float, default=1.0, help="negatives per positive (0 = none)")
    reg.add_argument("--exclude", type=int, default=200, help="min distance in bp of a negative from any TSS")
    reg.add_argument("--seed", type=int, default=42)

    fa = sub.add_parser("fasta", parents=[common], help="a FASTA file, one label for all records")
    fa.add_argument("--fasta", required=True)
    fa.add_argument("--label", required=True, help="Promoter or Non-Promoter (or +/-)")

    args = parser.parse_args()
    try:
        if args.source == "regulondb":
            entries, notes = build_regulondb_entries(
                args.tsv, args.genome, args.seq_col, args.tss_col, args.strand_col,
                has_header=not args.no_header, upstream=args.upstream, downstream=args.downstream,
                neg_ratio=args.neg_ratio, exclude=args.exclude, seed=args.seed)
        else:
            entries, notes = build_fasta_entries(args.fasta, args.label), []
    except FileNotFoundError as e:
        parser.exit(2, f"error: file not found: {e.filename}\n"
                       "(the file names in the README are examples -- pass the path of a file you have)\n")
    except ValueError as e:
        parser.exit(2, f"error: {e}\n")
    _finish(entries, notes, args)

if __name__ == "__main__":
    main()
