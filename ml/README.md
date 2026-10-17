# GeneScript DNA classification

The optional Python ML layer classifies a valid DNA string as `Promoter` or
`Non-Promoter`. It does not change the C compiler or its architecture.

## Dataset and preprocessing

Training uses the included UCI Machine Learning Repository **Molecular Biology
(Promoter Gene Sequences)** dataset (ID 67; DOI 10.24432/C5S01D; CC BY 4.0).
It contains 106 *E. coli* records, with 53 `+` promoter and 53 `-`
non-promoter labels, each 57 bases. The file is `ml/data/promoters.data`
(row format: `label,name,sequence`); its source is described in
`ml/data/README.md`. Records are checked for the expected label, A/C/G/T
alphabet, and a sane minimum length (`ml.dataset.MIN_SEQUENCE_LENGTH`, 10
bases) -- the 57-base width is a property of this one source dataset, not a
model requirement, so rows added later (see below) don't have to match it.
User sequences may be lowercase and are normalized to uppercase; empty or
other-alphabet inputs are rejected.

## Adding more training data

106 samples is a small dataset (see Limitations below) -- the direct way to
improve on it is to add more labeled sequences, through `ml/add_data.py`
rather than hand-editing `promoters.data`. It validates everything first
(DNA alphabet, minimum length, duplicates) and reports exactly what was
added vs. skipped and why; nothing is written if every row is invalid.

Add one sequence:

```bash
python -m ml.add_data --label Promoter --sequence ACGTACGT... --name my_seq_1
```

Bulk import from a CSV (or TSV) with `label` and `sequence` columns --
`name` is optional, and label accepts `Promoter`/`Non-Promoter`, `+`/`-`,
or `P`/`N`, case-insensitively:

```bash
python -m ml.add_data --csv new_sequences.csv
```

Add `--retrain` to either form to retrain on the updated dataset immediately
and print the new metrics -- useful for seeing right away whether the added
rows actually helped:

```bash
python -m ml.add_data --csv new_sequences.csv --retrain
```

`--retrain` writes to `ml/models/promoter_classifier.joblib` by default --
the same file the UI and `ml.predict` use. Pass `--model scratch.joblib`
together with `--dataset scratch.data` when trying this out against a copy
rather than the real dataset, so a test run can't overwrite the production
model.

A sequence that's byte-identical to one already in the dataset is skipped
as a duplicate (pass `--allow-duplicates` to add it anyway). Rows are also
usable directly from Python:

```python
from ml.add_data import add_sequences
report = add_sequences([
    {"label": "Promoter", "sequence": "ACGT...", "name": "sample_1"},
    {"label": "-", "sequence": "TTGA..."},
])
print(report["added"], report["skipped"])
```

## Features and model

`ml/features.py` creates 22 features in one shared implementation: A/C/G/T
frequencies, GC content, sequence length, and frequencies of all 16
dinucleotides in lexicographic A/C/G/T order. Dinucleotide frequencies use
the number of adjacent pairs as denominator (single-base input yields zeros).
The model is `RandomForestClassifier(n_estimators=200, random_state=42)`.
Training uses a stratified 80/20 split with `random_state=42`. Evaluation is
separate in `ml/evaluate.py`: accuracy, positive-class precision/recall/F1,
confusion matrix, classification report, and ROC-AUC when probabilities and
both test labels are available. Metrics are emitted by a real training run,
not hard-coded.

## Usage

Install UI dependencies with `pip install -r ui/requirements.txt`. Train and
persist the model using `python -m ml.train`; this writes
`ml/models/promoter_classifier.joblib`. Inference does not retrain:

```python
from ml.predict import predict_sequence
print(predict_sequence("TTGACA..."))
```

Run the optional Streamlit section with the existing command:
`streamlit run ui/app.py`, then select **ML Classification**, enter a sequence,
and click **Predict**.

## Limitations

This is a small historical *E. coli* benchmark (106 sequences), not a broad
or clinically validated promoter detector. A 22-feature composition model
does not retain positional promoter motifs and predictions should be treated
as educational estimates. Performance reported by the fixed random split can
be unstable for this dataset size and is not independent external validation.
