"""Reusable, fixed-order 22-value DNA feature extraction."""

from .validation import validate_sequence

DINUCLEOTIDES = tuple(a + b for a in "ACGT" for b in "ACGT")
FEATURE_NAMES = ("A_frequency", "C_frequency", "G_frequency", "T_frequency",
                 "GC_content", "length", *(f"{d}_frequency" for d in DINUCLEOTIDES))


def extract_features(sequence: str) -> list[float]:
    seq = validate_sequence(sequence)
    n = len(seq)
    values = [seq.count(base) / n for base in "ACGT"]
    values.extend(((seq.count("G") + seq.count("C")) / n, float(n)))
    pairs = len(seq) - 1
    values.extend((sum(seq[i:i + 2] == d for i in range(pairs)) / pairs if pairs else 0.0)
                  for d in DINUCLEOTIDES)
    return values


def extract_feature_matrix(sequences):
    """Extract a 2-D matrix using the same feature logic as prediction."""
    return [extract_features(sequence) for sequence in sequences]
