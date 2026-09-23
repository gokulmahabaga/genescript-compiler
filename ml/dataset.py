"""Loader and validator for the UCI promoters.data format."""

from pathlib import Path

from .validation import validate_sequence

DEFAULT_DATASET = Path(__file__).resolve().parent / "data" / "promoters.data"


def load_dataset(path=DEFAULT_DATASET):
    """Return (sequences, labels) from UCI rows: +/- , name, 57 bases."""
    path = Path(path)
    sequences, labels = [], []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split(",", 2)
        if len(parts) != 3 or parts[0].strip() not in {"+", "-"}:
            raise ValueError(f"Malformed dataset row {line_number} in {path}")
        label, _name, seq_field = parts
        seq = validate_sequence(seq_field)
        if len(seq) != 57:
            raise ValueError(f"Expected 57 bases on dataset row {line_number}; got {len(seq)}")
        sequences.append(seq)
        labels.append("Promoter" if label.strip() == "+" else "Non-Promoter")
    if not sequences or set(labels) != {"Promoter", "Non-Promoter"}:
        raise ValueError("Dataset must contain sequences from both classes.")
    return sequences, labels
