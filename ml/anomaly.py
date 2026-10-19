"""Outlier detection: flags whether a sequence's feature vector looks
atypical next to the training set, using an Isolation Forest -- a model
that isolates points by random recursive splits, so points that split
off quickly (few splits needed to isolate them) are scored as outliers.
Unlike the promoter/non-promoter classifier, this needs no label: it
just asks "does this composition look unusual compared to what we've
seen," independent of whether it turns out to be a promoter.

Fit fresh on every call, same reasoning as ml.cluster/ml.similarity: the
training set is small enough that persisting a model isn't worth the
staleness risk between dataset updates.

Known current-dataset quirk: the original 106 UCI rows are *all* exactly
57 bases, i.e. the `length` feature has zero variance in training. Isolation
Forest's splits are scale-invariant (so a StandardScaler wouldn't help) but
not variance-invariant: against a feature with zero training variance, any
different length gets isolated in a single split and so reads as maximally
anomalous, regardless of composition. In practice this means right now
*any* sequence that isn't exactly 57 bases is flagged, dominated by length
alone -- genuinely anomalous composition only shows through once length is
held equal (see tests/test_ml.py's AnomalyTests, which does exactly that).
This resolves itself as ml.add_data brings in sequences of varied lengths,
at which point length stops being a constant and becomes one signal among
several, same as any other feature here.
"""

import argparse
import json

from sklearn.ensemble import IsolationForest

from .dataset import DEFAULT_DATASET, load_dataset
from .features import extract_feature_matrix, extract_features
from .validation import validate_sequence

DEFAULT_CONTAMINATION = 0.1  # expected fraction of the training set that's "unusual"


def check_sequence(sequence, dataset_path=DEFAULT_DATASET,
                    contamination=DEFAULT_CONTAMINATION, random_state=42):
    """Return {"is_outlier", "anomaly_score"} for `sequence` against a
    fresh Isolation Forest trained on `dataset_path`.

    `anomaly_score` is scikit-learn's score_samples output: higher means
    more typical/normal, lower (more negative) means more anomalous.
    `is_outlier` is the forest's own inlier/outlier call at the
    `contamination` threshold used to train it.
    """
    normalized = validate_sequence(sequence)
    sequences, _ = load_dataset(dataset_path)
    X = extract_feature_matrix(sequences)
    model = IsolationForest(contamination=contamination, random_state=random_state).fit(X)

    x = [extract_features(normalized)]
    is_outlier = bool(model.predict(x)[0] == -1)
    score = float(model.score_samples(x)[0])
    return {"is_outlier": is_outlier, "anomaly_score": score}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET))
    parser.add_argument("--contamination", type=float, default=DEFAULT_CONTAMINATION,
                         help="expected fraction of outliers in the training set (default 0.1)")
    args = parser.parse_args()
    print(json.dumps(check_sequence(args.sequence, args.dataset, args.contamination), indent=2))


if __name__ == "__main__":
    main()
