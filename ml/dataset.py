"""Loader and validator for the promoters.data format.

Row format: `<label>,<name>,<sequence>` where label is `+` (promoter) or
`-` (non-promoter). The original 106 UCI rows are all fixed-width 57-base
*E. coli* sequences; that width is a property of that one source dataset,
not a requirement of the model (features.py's extraction is frequency-based
and works on any length). Rows added later via `ml/add_data.py` may be a
different length -- only a sane minimum is enforced here so dinucleotide
frequencies stay meaningful.
"""

from pathlib import Path

from .validation import validate_sequence

DEFAULT_DATASET = Path(__file__).resolve().parent / "data" / "promoters.data"
MIN_SEQUENCE_LENGTH = 10  # below this, dinucleotide frequencies are too noisy to be useful


def parse_row(line, line_number=1, path="<string>"):
    """Parse one `label,name,sequence` row. Returns (sequence, label) or raises ValueError."""
    parts = line.split(",", 2)
    if len(parts) != 3 or parts[0].strip() not in {"+", "-"}:
        raise ValueError(f"Malformed dataset row {line_number} in {path}")
    label, _name, seq_field = parts
    seq = validate_sequence(seq_field)
    if len(seq) < MIN_SEQUENCE_LENGTH:
        raise ValueError(
            f"Dataset row {line_number} in {path} has only {len(seq)} base(s); "
            f"need at least {MIN_SEQUENCE_LENGTH}."
        )
    return seq, ("Promoter" if label.strip() == "+" else "Non-Promoter")


def load_dataset(path=DEFAULT_DATASET):
    """Return (sequences, labels) from every non-blank row in `path`."""
    path = Path(path)
    sequences, labels = [], []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        seq, label = parse_row(line, line_number, path)
        sequences.append(seq)
        labels.append(label)
    if not sequences or set(labels) != {"Promoter", "Non-Promoter"}:
        raise ValueError("Dataset must contain sequences from both classes.")
    return sequences, labels
