"""Registry of classifiers available for the promoter task.

Random Forest stays the production default -- it's what `ml.train` writes
to `ml/models/promoter_classifier.joblib` and what `ml.predict`/the UI load
unless told otherwise -- because it's scale-invariant: tree splits don't
care that `length` (roughly 20-150) sits on a totally different scale from
the 21 bounded [0,1] frequency features in `ml.features`.

The other three are distance- or gradient-based, so that same scale gap
would let `length` dominate them and understate what they can actually do.
Each is wrapped in a `StandardScaler` pipeline to correct for that before
comparing them -- see `ml.compare`.
"""

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

DEFAULT_MODEL = "random_forest"

DISPLAY_NAMES = {
    "random_forest": "Random Forest",
    "logistic_regression": "Logistic Regression",
    "svm": "SVM (RBF kernel)",
    "knn": "K-Nearest Neighbors",
}


def _random_forest():
    return RandomForestClassifier(n_estimators=200, random_state=42)


def _logistic_regression():
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000, random_state=42))


def _svm():
    return make_pipeline(StandardScaler(), SVC(kernel="rbf", probability=True, random_state=42))


def _knn():
    return make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=5))


_BUILDERS = {
    "random_forest": _random_forest,
    "logistic_regression": _logistic_regression,
    "svm": _svm,
    "knn": _knn,
}


def available_models():
    """Names of every registered model, in a fixed, stable order."""
    return list(_BUILDERS)


def build_model(name):
    """Return a fresh, unfitted estimator for `name` (see available_models())."""
    try:
        builder = _BUILDERS[name]
    except KeyError:
        raise ValueError(f"Unknown model {name!r}. Available: {', '.join(available_models())}") from None
    return builder()
