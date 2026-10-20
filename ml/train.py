"""Train and persist a classifier (default: Random Forest) on the real UCI dataset."""

import argparse
import json
from pathlib import Path

import joblib
from sklearn.model_selection import train_test_split

from .dataset import DEFAULT_DATASET, load_dataset, warn_if_skewed
from .evaluate import evaluate_model
from .features import extract_feature_matrix
from .models import DEFAULT_MODEL, available_models, build_model

MODEL_PATH = Path(__file__).resolve().parent / "models" / "promoter_classifier.joblib"


def train(dataset_path=DEFAULT_DATASET, model_path=MODEL_PATH, model_name=DEFAULT_MODEL):
    sequences, labels = load_dataset(dataset_path)
    warn_if_skewed(sequences, labels)
    X = extract_feature_matrix(sequences)
    X_train, X_test, y_train, y_test = train_test_split(
        X, labels, test_size=0.2, random_state=42, stratify=labels)
    model = build_model(model_name)
    model.fit(X_train, y_train)
    metrics = evaluate_model(model, X_test, y_test)
    model_path = Path(model_path)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, model_path)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--algorithm", default=DEFAULT_MODEL, choices=available_models(),
                         help=f"which classifier to train and persist (default: {DEFAULT_MODEL}). "
                              "See ml.compare to evaluate all of them before picking one.")
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET), help="dataset file to train on")
    # Named --model (the output path), matching ml.add_data's --retrain flag of the same
    # name -- not --algorithm, which is the choice above.
    parser.add_argument("--model", default=str(MODEL_PATH),
                         help="where to write the trained model (defaults to the production path "
                              "ml.predict and the UI load -- pass a scratch path to try an "
                              "algorithm without replacing the production model)")
    args = parser.parse_args()
    metrics = train(dataset_path=args.dataset, model_path=args.model, model_name=args.algorithm)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
