"""Nearest-neighbor similarity search: given a sequence, finds the most
similar labeled sequences already in the training set, using the same
22-value feature representation (ml.features) the classifiers are trained
on -- so "similar" means close in composition (base frequencies, GC
content, length, dinucleotide frequencies), not raw string edit distance.

Fit fresh on every call, same reasoning as ml.cluster: the dataset is a
few hundred rows, so there's nothing worth persisting to disk, and the
index always reflects whatever ml.add_data has most recently added.
"""

import argparse
import json

from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from .dataset import DEFAULT_DATASET, load_dataset
from .features import extract_feature_matrix, extract_features
from .validation import validate_sequence

DEFAULT_K = 5


def find_similar(sequence, k=DEFAULT_K, dataset_path=DEFAULT_DATASET):
    """Return the `k` dataset sequences closest to `sequence` in (scaled)
    feature space, nearest first, each as {"sequence", "label", "distance"}.

    An exact match to `sequence` already in the dataset is excluded from
    its own neighbor list (distance 0 to yourself isn't an interesting
    result) -- so up to `k` neighbors are returned, never `sequence` itself.
    """
    normalized = validate_sequence(sequence)
    sequences, labels = load_dataset(dataset_path)

    X = extract_feature_matrix(sequences)
    scaler = StandardScaler().fit(X)
    Xs = scaler.transform(X)
    # +1 so there's still room for k real neighbors after dropping a self-match
    index = NearestNeighbors(n_neighbors=min(k + 1, len(sequences))).fit(Xs)

    x = scaler.transform([extract_features(normalized)])
    distances, indices = index.kneighbors(x)

    results = []
    for dist, idx in zip(distances[0], indices[0]):
        if sequences[idx] == normalized:
            continue
        results.append({"sequence": sequences[idx], "label": labels[idx], "distance": float(dist)})
    return results[:k]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sequence", required=True, help="DNA sequence to find neighbors for")
    parser.add_argument("--k", type=int, default=DEFAULT_K)
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET))
    args = parser.parse_args()
    print(json.dumps(find_similar(args.sequence, args.k, args.dataset), indent=2))


if __name__ == "__main__":
    main()
