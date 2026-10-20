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

### Importing external data (kept out of the repo)

`ml/import_data.py` turns data you download yourself into dataset rows, so
third-party data never has to be committed. The repo ships only the UCI
dataset.

```bash
# RegulonDB-style promoter TSV + random E. coli genome windows as negatives
python -m ml.import_data regulondb --tsv PromoterSet.tsv --genome U00096.3.fasta --dry-run
python -m ml.import_data regulondb --tsv PromoterSet.tsv --genome U00096.3.fasta --retrain

# any FASTA, one label for every record
python -m ml.import_data fasta --fasta extra.fasta --label Non-Promoter
```

- Positives come from a sequence column if the TSV has one, otherwise they
  are cut from the genome around each TSS (`--upstream 50 --downstream 7`
  gives 57 bp, the UCI width). Columns are auto-detected from the header or
  set with `--seq-col/--tss-col/--strand-col` (1-based index or header text).
- Negatives are random genome windows (either strand) with the **same
  lengths as the positives** and at least `--exclude 200` bp from every
  known TSS, so length can't reveal the label. `--seed` makes them
  reproducible; `--neg-ratio` sets negatives per positive.
- `--dry-run` validates and reports without writing. Use
  `--dataset scratch.data --model scratch.joblib` to experiment on a copy.

**Licensing.** RegulonDB is free for academic/non-commercial use but does
not allow redistributing it inside other databanks without written consent,
and the genome is an NCBI record. Download them yourself, cite RegulonDB
(Salgado et al., v12.0) in anything you present, and **don't commit** the
downloaded files or a `promoters.data` grown from them -- keep a scratch
copy (`--dataset`) instead.

**Comparing numbers.** Metrics on a grown dataset aren't comparable with the
106-row UCI figures: different data, different negatives, bigger test set.
Report them as a new experiment, not an improvement of the old number.

### Skew warnings

`ml.train`, `ml.compare` and `ml.add_data` warn (`LengthSkewWarning`) when
the data is skewed in a way that inflates metrics: the smaller class is
under 25% of the rows, or sequence length alone separates the classes
(AUC >= 0.85) -- `length` is a model feature, so mismatched lengths let a
model "predict" the label without looking at composition. The shipped UCI
data is clean (all 57 bp, 53/53).

## Features and model

`ml/features.py` creates 22 features in one shared implementation: A/C/G/T
frequencies, GC content, sequence length, and frequencies of all 16
dinucleotides in lexicographic A/C/G/T order. Dinucleotide frequencies use
the number of adjacent pairs as denominator (single-base input yields zeros).
Training uses a stratified 80/20 split with `random_state=42`. Evaluation is
separate in `ml/evaluate.py`: accuracy, positive-class precision/recall/F1,
confusion matrix, classification report, and ROC-AUC when probabilities and
both test labels are available. Metrics are emitted by a real training run,
not hard-coded.

`ml/models.py` is a small registry of four classifiers, all trained on the
same 22 features:

| Name (`ml.models` key) | Algorithm | Notes |
|---|---|---|
| `random_forest` (**production default**) | `RandomForestClassifier(n_estimators=200, random_state=42)` | Tree-based, so it's scale-invariant -- used as-is |
| `logistic_regression` | `LogisticRegression(max_iter=1000, random_state=42)` | Behind a `StandardScaler` (see below) |
| `svm` | `SVC(kernel="rbf", probability=True, random_state=42)` | Behind a `StandardScaler` |
| `knn` | `KNeighborsClassifier(n_neighbors=5)` | Behind a `StandardScaler` |

The three non-tree models are distance- or gradient-based, so they're each
wrapped in a `StandardScaler` pipeline first -- without it, the unbounded
`length` feature (roughly 20-150) would dominate the 21 other features,
which are all bounded to `[0, 1]`, and understate what those models can
actually do. `random_forest` doesn't need this (tree splits don't care
about feature scale), so it's used directly, matching the model this
project has used since the ML layer was first added.

**`random_forest` is still the one production model** -- the only one
`ml.predict` and the UI's **Predict** button load from
`ml/models/promoter_classifier.joblib`. The other three exist for
comparison (below); training one doesn't change what **Predict** uses
unless you explicitly retrain to the production path yourself.

## Comparing models

`ml/compare.py` trains every registered model on the *same* stratified
train/test split -- so the comparison is apples-to-apples -- and reports
the same metrics `ml.evaluate_model` produces elsewhere:

```bash
python -m ml.compare                   # table, all 4 registered models
python -m ml.compare --json            # same run, full metrics as JSON
python -m ml.compare --models svm knn  # only compare these two
python -m ml.compare --save-models ml/models/comparison
                                        # also persist each fitted model
                                        # there as <name>.joblib
```

On the full 106-sequence dataset, a typical run looks like:

```
Model                   Accuracy   Precision   Recall   F1      ROC-AUC
------------------------------------------------------------------------
Random Forest           0.818      0.889       0.727    0.800   0.913
Logistic Regression     0.818      0.889       0.727    0.800   0.868
SVM (RBF kernel)        0.864      0.900       0.818    0.857   0.909
K-Nearest Neighbors     0.864      0.900       0.818    0.857   0.855
```

(Exact numbers will shift slightly as the dataset grows via `ml.add_data`.)
This isn't a claim that SVM or KNN are definitively better -- the test set
here is only ~21 sequences, so these numbers move around on a single split.
It's meant to show the comparison exists and how to run it, not to settle
which algorithm is "best" on a dataset this small. The same comparison is
available from the UI: open **ML Classification**, then click **Compare
models**.

Also usable directly from Python:

```python
from ml.compare import compare_models
for result in compare_models():
    print(result["display_name"], result["accuracy"])
```

### Picking a different production model

By default, `ml.train` (and `ml.add_data --retrain`) both train and persist
`random_forest`. To make a different algorithm the production model that
`ml.predict`/the UI actually use, train it straight to that same path:

```bash
python -m ml.train --algorithm svm --model ml/models/promoter_classifier.joblib
```

Leave off `--model` and it still writes there by default -- so to just try
an algorithm without touching the production model, point `--model`
somewhere else instead:

```bash
python -m ml.train --algorithm knn --model /tmp/knn_test.joblib
```

## Unsupervised exploration

Three more views of a sequence, independent of the promoter/non-promoter
classifier above -- none of them need a label, so they say something even
about a sequence you have no ground truth for. Each is fit fresh against
the dataset on every call (never persisted to disk): the training sets
here are small enough that refitting costs nothing, and this way every
result reflects whatever `ml.add_data` has most recently added, instead of
going stale between retrains.

**GC-content banding** (`ml/cluster.py`) groups the dataset into Low/Medium/
High GC bands with KMeans on GC% alone -- the band boundaries are *learned*
from the data's own distribution, not hardcoded thresholds:

```bash
python -m ml.cluster                              # print the current band ranges
python -m ml.cluster --sequence ACGTACGT...        # classify one sequence's band
```

**Similarity search** (`ml/similarity.py`) finds the `k` dataset sequences
closest to a query in the same 22-feature space the classifiers use (scaled
with `StandardScaler` first) -- so "similar" means close in composition
(base/dinucleotide frequencies, GC%, length), not raw string edit distance.
An exact match to itself is never returned as its own neighbor:

```bash
python -m ml.similarity --sequence ACGTACGT... --k 5
```

**Anomaly detection** (`ml/anomaly.py`) fits an Isolation Forest on the
dataset and scores whether a sequence's composition looks typical or
unusual next to it, regardless of whether it turns out to be a promoter:

```bash
python -m ml.anomaly --sequence ACGTACGT...
```

All three are also usable directly from Python (`classify_gc_band`,
`find_similar`, `check_sequence`) and from the UI's **ML Exploration**
expander, each as its own tab.

Current-dataset caveat for anomaly detection: the original 106 UCI rows are
*all* exactly 57 bases, so `length` has zero variance in training. Isolation
Forest's splits are scale-invariant but not variance-invariant against a
constant feature, so right now any sequence that isn't exactly 57 bases
gets isolated in one split and reads as maximally anomalous regardless of
composition. This resolves itself as `ml.add_data` brings in varied-length
sequences, at which point `length` stops being a constant and becomes one
signal among several, same as any other feature here.

## Usage

Install UI dependencies with `pip install -r ui/requirements.txt`. Train and
persist the production model using `python -m ml.train`; this writes
`ml/models/promoter_classifier.joblib` (Random Forest by default --
`--algorithm` picks a different one, see above). Inference does not retrain:

```python
from ml.predict import predict_sequence
print(predict_sequence("TTGACA..."))
```

Run the optional Streamlit section with the existing command:
`streamlit run ui/app.py`. **ML Classification** holds **Predict** (a single
call against the production model) and **Compare models** (all four
registered models on the same held-out split). **ML Exploration** holds the
three unsupervised views above, each its own tab: **GC band**, **Similar
sequences**, and **Anomaly check**.

## Limitations

This is a small historical *E. coli* benchmark (106 sequences), not a broad
or clinically validated promoter detector. A 22-feature composition model
does not retain positional promoter motifs and predictions should be treated
as educational estimates. Performance reported by the fixed random split can
be unstable for this dataset size and is not independent external validation.
