"""Loader and validator for the promoters.data format.

Row format: `<label>,<name>,<sequence>` where label is `+` (promoter) or
`-` (non-promoter). The original 106 UCI rows are all fixed-width 57-base
*E. coli* sequences; that width is a property of that one source dataset,
not a requirement of the model (features.py's extraction is frequency-based
and works on any length). Rows added later via `ml/add_data.py` may be a
different length -- only a sane minimum is enforced here so dinucleotide
frequencies stay meaningful.
"""

import warnings
from pathlib import Path

from .validation import validate_sequence

DEFAULT_DATASET = Path(__file__).resolve().parent / "data" / "promoters.data"
MIN_SEQUENCE_LENGTH = 10  # below this, dinucleotide frequencies are too noisy to be useful

# Skew-check thresholds (see dataset_skew_issues).
MIN_CLASS_SHARE = 0.25         # warn if the smaller class is under 25% of the rows
LENGTH_PROXY_AUC = 0.85        # warn if length alone separates the classes this well
MIN_PER_CLASS_FOR_LENGTH = 5   # too few rows per class makes the length check meaningless


class LengthSkewWarning(UserWarning):
    """The dataset's labels are imbalanced or confounded with sequence length."""


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


def _length_auc(sequences, labels):
    """Mann-Whitney AUC of sequence length as a score for "Promoter" (0.5 = no
    signal, 1.0/0.0 = length alone perfectly separates the classes). Uses
    average ranks, so ties (e.g. all-57-base data) correctly give 0.5."""
    order = sorted(range(len(sequences)), key=lambda i: len(sequences[i]))
    ranks = [0.0] * len(order)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and len(sequences[order[j + 1]]) == len(sequences[order[i]]):
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    pos = [r for r, label in zip(ranks, labels) if label == "Promoter"]
    n_pos, n_neg = len(pos), len(labels) - len(pos)
    return (sum(pos) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def dataset_skew_issues(sequences, labels):
    """Return human-readable warnings (an empty list if the data looks fine).

    Two problems that quietly inflate metrics on grown datasets:
    - class imbalance: accuracy looks good by always guessing the big class;
    - length confounding: `length` is one of the 22 features, so if
      promoters and non-promoters were collected or built at different
      lengths, the model can score well by reading length instead of
      composition -- and then fail on real sequences.
    """
    issues = []
    n = len(labels)
    n_pos = sum(1 for label in labels if label == "Promoter")
    n_neg = n - n_pos
    if n and min(n_pos, n_neg) / n < MIN_CLASS_SHARE:
        issues.append(
            f"label imbalance: {n_pos} Promoter vs {n_neg} Non-Promoter "
            f"(smaller class is under {MIN_CLASS_SHARE:.0%}). Accuracy will look better than it is; "
            "check precision/recall/F1 instead, or add more of the smaller class."
        )
    if min(n_pos, n_neg) >= MIN_PER_CLASS_FOR_LENGTH:
        auc = _length_auc(sequences, labels)
        if max(auc, 1 - auc) >= LENGTH_PROXY_AUC:
            issues.append(
                f"length skew: sequence length alone separates the classes (AUC {max(auc, 1 - auc):.2f}). "
                "The model can learn length instead of promoter composition and its metrics will be "
                "inflated. Build both classes at the same length (ml.import_data does this)."
            )
    return issues


def warn_if_skewed(sequences, labels, stacklevel=2):
    """Emit a LengthSkewWarning per issue found by dataset_skew_issues; return the list."""
    issues = dataset_skew_issues(sequences, labels)
    for message in issues:
        warnings.warn(message, LengthSkewWarning, stacklevel=stacklevel + 1)
    return issues
