"""DNA input validation shared by training, prediction, and the UI."""


class DNAValidationError(ValueError):
    """Raised when a sequence is not a non-empty A/C/G/T string."""


def validate_sequence(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise DNAValidationError("DNA sequence must be text.")
    normalized = "".join(sequence.split()).upper()
    if not normalized:
        raise DNAValidationError("DNA sequence cannot be empty.")
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        raise DNAValidationError("Invalid DNA character(s): " + ", ".join(invalid))
    return normalized
