"""Train and persist the requested random forest on the real UCI dataset."""

import json
from pathlib import Path

import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

from .dataset import DEFAULT_DATASET, load_dataset
from .evaluate import evaluate_model
from .features import extract_feature_matrix

MODEL_PATH = Path(__file__).resolve().parent / "models" / "promoter_classifier.joblib"


def train(dataset_path=DEFAULT_DATASET, model_path=MODEL_PATH):
    sequences, labels = load_dataset(dataset_path)
    X = extract_feature_matrix(sequences)
    X_train, X_test, y_train, y_test = train_test_split(
        X, labels, test_size=0.2, random_state=42, stratify=labels)
    model = RandomForestClassifier(n_estimators=200, random_state=42)
    model.fit(X_train, y_train)
    metrics = evaluate_model(model, X_test, y_test)
    model_path = Path(model_path)
    model_path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, model_path)
    return metrics


if __name__ == "__main__":
    print(json.dumps(train(), indent=2))
