"""Prediction API. Inference loads the persisted estimator and never trains."""

from pathlib import Path

import joblib

from .features import extract_features
from .train import MODEL_PATH
from .validation import validate_sequence


def load_model(model_path=MODEL_PATH):
    model_path = Path(model_path)
    if not model_path.exists():
        raise FileNotFoundError(f"Trained model not found at {model_path}. Run python -m ml.train first.")
    return joblib.load(model_path)


def predict_sequence(sequence, model=None):
    normalized = validate_sequence(sequence)
    estimator = model if model is not None else load_model()
    row = [extract_features(normalized)]
    prediction = str(estimator.predict(row)[0])
    probabilities = estimator.predict_proba(row)[0] if hasattr(estimator, "predict_proba") else None
    confidence = (float(probabilities[list(estimator.classes_).index(prediction)])
                  if probabilities is not None else None)
    return {"prediction": prediction, "confidence": confidence}
