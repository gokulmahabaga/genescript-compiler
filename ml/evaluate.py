"""Evaluation helpers kept outside the Streamlit UI."""

from sklearn.metrics import (accuracy_score, classification_report,
                             confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score)


def evaluate_model(model, X_test, y_test):
    predictions = model.predict(X_test)
    result = {
        "accuracy": accuracy_score(y_test, predictions),
        "precision": precision_score(y_test, predictions, pos_label="Promoter", zero_division=0),
        "recall": recall_score(y_test, predictions, pos_label="Promoter", zero_division=0),
        "f1": f1_score(y_test, predictions, pos_label="Promoter", zero_division=0),
        "confusion_matrix": confusion_matrix(y_test, predictions, labels=["Non-Promoter", "Promoter"]).tolist(),
        "classification_report": classification_report(y_test, predictions, zero_division=0),
    }
    if hasattr(model, "predict_proba") and len(set(y_test)) == 2:
        positive_index = list(model.classes_).index("Promoter")
        result["roc_auc"] = roc_auc_score(y_test, model.predict_proba(X_test)[:, positive_index])
    return result
