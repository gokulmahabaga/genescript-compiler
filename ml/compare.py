"""Train every registered model on one shared split and compare them.

Unlike `ml.train` (which persists one model as the production estimator
`ml.predict` and the UI load), this is for evaluation only: it fits each
model in `ml.models` on the *same* stratified train/test split, so the
comparison is apples-to-apples, and reports the same metrics
`ml.evaluate_model` produces elsewhere. Nothing here touches
`ml/models/promoter_classifier.joblib` -- use
`python -m ml.train --model NAME` for that, once you've decided which
model you actually want as the production one.

CLI usage:

    python -m ml.compare                  # table, all 4 registered models
    python -m ml.compare --json           # same run, full metrics as JSON
    python -m ml.compare --models svm knn # only compare these two
    python -m ml.compare --save-models ml/models/comparison
                                           # also persist each fitted model
                                           # there as <name>.joblib, so you
                                           # can load a specific one later
                                           # without retraining
"""

import argparse
import json
from pathlib import Path

import joblib
from sklearn.model_selection import train_test_split

from .dataset import DEFAULT_DATASET, load_dataset
from .evaluate import evaluate_model
from .features import extract_feature_matrix
from .models import DISPLAY_NAMES, available_models, build_model


def compare_models(dataset_path=DEFAULT_DATASET, model_names=None, random_state=42, save_dir=None):
    """Fit every model in `model_names` (default: all registered models) on
    one shared stratified split of `dataset_path`, evaluating each with the
    same metrics `ml.evaluate_model` produces elsewhere.

    Returns a list of dicts, each an `evaluate_model()` result plus
    "model" (the registry name) and "display_name". If `save_dir` is
    given, each fitted model is also persisted there as "<name>.joblib".
    """
    names = list(model_names) if model_names else available_models()
    sequences, labels = load_dataset(dataset_path)
    X = extract_feature_matrix(sequences)
    X_train, X_test, y_train, y_test = train_test_split(
        X, labels, test_size=0.2, random_state=random_state, stratify=labels)

    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for name in names:
        model = build_model(name)
        model.fit(X_train, y_train)
        metrics = evaluate_model(model, X_test, y_test)
        results.append({"model": name, "display_name": DISPLAY_NAMES.get(name, name), **metrics})
        if save_dir:
            joblib.dump(model, save_dir / f"{name}.joblib")
    return results


def _fmt(value):
    return f"{value:.3f}" if isinstance(value, float) else str(value if value is not None else "-")


def _print_table(results):
    columns = [("display_name", "Model", 22), ("accuracy", "Accuracy", 9),
               ("precision", "Precision", 10), ("recall", "Recall", 7),
               ("f1", "F1", 6), ("roc_auc", "ROC-AUC", 8)]
    header = "  ".join(title.ljust(width) for _, title, width in columns)
    print(header)
    print("-" * len(header))
    for row in results:
        print("  ".join(_fmt(row.get(key)).ljust(width) for key, _, width in columns))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET), help="dataset file to compare on")
    parser.add_argument("--models", nargs="+", choices=available_models(), default=None,
                         help="which registered models to compare (default: all of them)")
    parser.add_argument("--save-models", metavar="DIR",
                         help="also persist each fitted model to DIR/<name>.joblib")
    parser.add_argument("--json", action="store_true", help="print full results as JSON instead of a table")
    args = parser.parse_args()

    results = compare_models(dataset_path=args.dataset, model_names=args.models, save_dir=args.save_models)
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        _print_table(results)


if __name__ == "__main__":
    main()
