"""GC-content clustering: groups sequences into Low/Medium/High GC bands
using KMeans on GC% -- unsupervised, so band boundaries are *learned* from
the data's own distribution rather than hardcoded (e.g. "GC% < 40 is Low").

Fit fresh on every call rather than persisted to disk: the whole dataset
is a few hundred 1-D points, so clustering it costs nothing, and this way
the bands always reflect the current contents of promoters.data (including
anything ml.add_data has appended) instead of going stale between retrains.
"""

import argparse
import json

import numpy as np
from sklearn.cluster import KMeans

from .dataset import DEFAULT_DATASET, load_dataset
from .validation import validate_sequence

N_BANDS = 3
BAND_NAMES = ["Low GC", "Medium GC", "High GC"]  # ascending GC order


def gc_percent(seq: str) -> float:
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def fit_bands(dataset_path=DEFAULT_DATASET, n_bands=N_BANDS, random_state=42):
    """Fit KMeans on the dataset's GC% values.

    Returns (model, band_by_cluster, band_ranges):
      - model: the fitted KMeans (1 feature: GC%)
      - band_by_cluster: {raw cluster index -> band name}, raw cluster
        indices from KMeans are arbitrary, so this remaps them to band
        names ordered by ascending centroid GC%
      - band_ranges: {band name -> (min GC% seen, max GC% seen)} among the
        training sequences assigned to that band
    """
    sequences, _ = load_dataset(dataset_path)
    gc_values = np.array([[gc_percent(s)] for s in sequences])
    model = KMeans(n_clusters=n_bands, random_state=random_state, n_init=10).fit(gc_values)

    order = np.argsort(model.cluster_centers_.ravel())
    names = BAND_NAMES if n_bands == len(BAND_NAMES) else [f"Band {i + 1}" for i in range(n_bands)]
    band_by_cluster = {int(cluster): names[rank] for rank, cluster in enumerate(order)}

    ranges = {}
    assigned = model.labels_
    for cluster, name in band_by_cluster.items():
        vals = gc_values[assigned == cluster].ravel()
        if len(vals):
            ranges[name] = (float(vals.min()), float(vals.max()))
    return model, band_by_cluster, ranges


def classify_gc_band(sequence, dataset_path=DEFAULT_DATASET, n_bands=N_BANDS):
    """Classify `sequence`'s GC band against bands fit fresh on `dataset_path`.

    Returns {"gc_percent", "band", "band_ranges"}.
    """
    normalized = validate_sequence(sequence)
    model, band_by_cluster, ranges = fit_bands(dataset_path, n_bands)
    gc = gc_percent(normalized)
    cluster = int(model.predict([[gc]])[0])
    return {"gc_percent": gc, "band": band_by_cluster[cluster], "band_ranges": ranges}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sequence", help="classify this sequence's GC band; omit to just print band ranges")
    parser.add_argument("--dataset", default=str(DEFAULT_DATASET))
    parser.add_argument("--bands", type=int, default=N_BANDS)
    args = parser.parse_args()

    if args.sequence:
        print(json.dumps(classify_gc_band(args.sequence, args.dataset, args.bands), indent=2))
    else:
        _, _, ranges = fit_bands(args.dataset, args.bands)
        for name in sorted(ranges, key=lambda n: ranges[n][0]):
            lo, hi = ranges[name]
            print(f"{name:<12} GC% range: {lo:.1f}% - {hi:.1f}%")


if __name__ == "__main__":
    main()
