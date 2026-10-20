"""Add new labeled sequences to the training dataset (ml/data/promoters.data
by default, or any dataset file in the same format).

This is the supported way to grow the training set -- hand-editing the
dataset file works too (it's a shared format: `label,name,sequence`), but
this validates everything first (DNA alphabet, minimum length, duplicates)
and reports exactly what was added vs. skipped and why.

CLI usage:

    # add one sequence
    python -m ml.add_data --label Promoter --sequence ACGTACGT... --name my_seq_1

    # bulk import from a CSV with columns: label,name,sequence (name optional)
    # label accepts +/-, Promoter/Non-Promoter, P/N, 1/0 -- case-insensitive
    python -m ml.add_data --csv new_sequences.csv

    # either form, then retrain immediately and print the updated metrics
    python -m ml.add_data --csv new_sequences.csv --retrain

Library usage:

    from ml.add_data import add_sequences
    report = add_sequences([
        {"label": "Promoter", "sequence": "ACGT...", "name": "sample_1"},
        {"label": "-", "sequence": "TTGA..."},
    ])
    print(report["added"], report["skipped"])
"""

import argparse
import csv
import json
from pathlib import Path

from .dataset import DEFAULT_DATASET, MIN_SEQUENCE_LENGTH, load_dataset, warn_if_skewed
from .validation import DNAValidationError, validate_sequence

_LABEL_ALIASES = {
    "+": "+", "-": "-",
    "promoter": "+", "non-promoter": "-", "nonpromoter": "-", "non_promoter": "-",
    "p": "+", "n": "-",
    "1": "+", "0": "-",
}


def normalize_label(raw):
    """Map any accepted spelling of a label to the on-disk '+'/'-' form."""
    key = str(raw).strip().lower()
    if key not in _LABEL_ALIASES:
        raise ValueError(f"Unrecognized label {raw!r}. Use Promoter/Non-Promoter (or +/-, P/N).")
    return _LABEL_ALIASES[key]


def _format_row(label, name, sequence):
    name = (name or "unnamed").strip().replace(",", ";") or "unnamed"
    return f"{label},{name},{sequence}"


def _existing_sequences(dataset_path):
    if not dataset_path.exists() or dataset_path.stat().st_size == 0:
        return set()
    try:
        sequences, _ = load_dataset(dataset_path)
        return set(sequences)
    except ValueError:
        # Let a genuinely malformed dataset file surface when someone next
        # trains or reads it -- don't block adding new, valid rows here.
        return set()


def add_sequences(entries, dataset_path=DEFAULT_DATASET, allow_duplicates=False, dry_run=False):
    """Validate `entries` and append the valid, non-duplicate ones to
    `dataset_path` in the shared `label,name,sequence` format.

    `entries` is an iterable of dicts, each with:
      - "label": one of Promoter/Non-Promoter/+/-/P/N/1/0 (case-insensitive)
      - "sequence": a DNA string (any length >= ml.dataset.MIN_SEQUENCE_LENGTH)
      - "name": optional; auto-generated if omitted

    Returns {"added": [...], "skipped": [...]} -- each entry in "skipped"
    carries a "reason" explaining why (invalid DNA, too short, duplicate,
    bad label). Nothing is written if every entry is invalid, or if
    `dry_run` is true (everything is still validated and reported).
    """
    dataset_path = Path(dataset_path)
    existing = _existing_sequences(dataset_path)

    added, skipped, new_rows = [], [], []
    seen_this_batch = set()
    for i, entry in enumerate(entries, 1):
        name = (entry.get("name") or f"user_added_{i}").strip() or f"user_added_{i}"
        try:
            label = normalize_label(entry.get("label"))
            seq = validate_sequence(entry.get("sequence"))
            if len(seq) < MIN_SEQUENCE_LENGTH:
                raise ValueError(
                    f"sequence has only {len(seq)} base(s); need at least {MIN_SEQUENCE_LENGTH}"
                )
            if not allow_duplicates and (seq in existing or seq in seen_this_batch):
                raise ValueError("duplicate of a sequence already in the dataset")
        except (ValueError, DNAValidationError) as e:
            skipped.append({"name": name, "reason": str(e)})
            continue
        seen_this_batch.add(seq)
        new_rows.append(_format_row(label, name, seq))
        added.append({
            "name": name,
            "label": "Promoter" if label == "+" else "Non-Promoter",
            "length": len(seq),
        })

    if new_rows and not dry_run:
        dataset_path.parent.mkdir(parents=True, exist_ok=True)
        prior_text = dataset_path.read_text(encoding="utf-8") if dataset_path.exists() else ""
        needs_leading_newline = bool(prior_text) and not prior_text.endswith("\n")
        with dataset_path.open("a", encoding="utf-8") as f:
            if needs_leading_newline:
                f.write("\n")
            f.write("\n".join(new_rows) + "\n")

    return {"added": added, "skipped": skipped}


def _rows_from_csv(csv_path):
    """Read label/name/sequence rows from a CSV (or TSV) file. A header row
    is required; recognized columns (case-insensitive) are label, name,
    sequence. 'name' is optional."""
    csv_path = Path(csv_path)
    with csv_path.open(newline="", encoding="utf-8-sig") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t") if sample.strip() else csv.excel
        reader = csv.DictReader(f, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(f"{csv_path} has no header row (expected columns: label, sequence[, name])")
        fields = {name.strip().lower(): name for name in reader.fieldnames}
        if "label" not in fields or "sequence" not in fields:
            raise ValueError(
                f"{csv_path} must have 'label' and 'sequence' columns; found {reader.fieldnames}"
            )
        rows = []
        for row in reader:
            rows.append({
                "label": row.get(fields["label"], ""),
                "sequence": row.get(fields["sequence"], ""),
                "name": row.get(fields.get("name", ""), "") if "name" in fields else "",
            })
        return rows


def _print_report(report):
    for row in report["added"]:
        print(f"  added    {row['name']:<24} {row['label']:<13} ({row['length']} bases)")
    for row in report["skipped"]:
        print(f"  skipped  {row['name']:<24} -- {row['reason']}")
    print(f"\n{len(report['added'])} added, {len(report['skipped'])} skipped.")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--csv", help="CSV/TSV file with label,sequence[,name] columns (header required)")
    parser.add_argument("--label", help="Promoter/Non-Promoter (or +/-, P/N) for a single --sequence add")
    parser.add_argument("--sequence", help="a single DNA sequence to add (used with --label)")
    parser.add_argument("--name", help="optional name for the single sequence being added")
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET), help="dataset file to append to")
    parser.add_argument("--allow-duplicates", action="store_true",
                         help="add a sequence even if it's already in the dataset")
    parser.add_argument("--retrain", action="store_true",
                         help="retrain the model on the updated dataset and print the new metrics")
    parser.add_argument("--model",
                         help="where --retrain writes the trained model "
                              "(defaults to ml/models/promoter_classifier.joblib -- "
                              "the model the UI and ml.predict actually use). "
                              "Pass a scratch path here when testing against a non-default "
                              "--dataset, or the production model gets overwritten.")
    args = parser.parse_args()

    if not args.csv and not (args.label and args.sequence):
        parser.error("provide --csv FILE, or both --label and --sequence")

    entries = _rows_from_csv(args.csv) if args.csv else [
        {"label": args.label, "sequence": args.sequence, "name": args.name}
    ]

    report = add_sequences(entries, dataset_path=args.dataset, allow_duplicates=args.allow_duplicates)
    _print_report(report)
    if report["added"]:
        warn_if_skewed(*load_dataset(args.dataset))

    if args.retrain:
        if not report["added"]:
            print("\nNothing new was added -- skipping retrain.")
            return
        from .train import MODEL_PATH, train  # imported lazily: keeps `--help` fast and dependency-light
        model_path = args.model or MODEL_PATH
        print(f"\nRetraining on the updated dataset -> {model_path} ...")
        metrics = train(dataset_path=args.dataset, model_path=model_path)
        print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
