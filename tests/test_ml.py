import random
import tempfile
import unittest
import warnings
from pathlib import Path

from ml.add_data import add_sequences, normalize_label
from ml.anomaly import check_sequence as check_anomaly
from ml.cluster import classify_gc_band, fit_bands, gc_percent
from ml.compare import compare_models
from ml.dataset import LengthSkewWarning, dataset_skew_issues, load_dataset, warn_if_skewed
from ml.import_data import (build_fasta_entries, build_regulondb_entries, read_fasta,
                            reverse_complement, sample_negatives, window_at)
from ml.features import DINUCLEOTIDES, extract_features
from ml.models import DEFAULT_MODEL, available_models, build_model
from ml.predict import load_model, predict_sequence
from ml.similarity import find_similar
from ml.validation import DNAValidationError, validate_sequence


class ValidationTests(unittest.TestCase):
    def test_valid_dna(self):
        self.assertEqual(validate_sequence("ACGT"), "ACGT")

    def test_lowercase_normalization(self):
        self.assertEqual(validate_sequence("acgt"), "ACGT")

    def test_whitespace_is_ignored(self):
        self.assertEqual(validate_sequence("ac gt\n"), "ACGT")

    def test_invalid_dna(self):
        with self.assertRaises(DNAValidationError):
            validate_sequence("ACNX")

    def test_empty_sequence(self):
        with self.assertRaises(DNAValidationError):
            validate_sequence("  \n")


class FeatureTests(unittest.TestCase):
    def test_feature_values_and_count(self):
        features = extract_features("AACG")
        self.assertEqual(len(features), 22)
        self.assertEqual(len(DINUCLEOTIDES), 16)
        self.assertEqual(features[:4], [0.5, 0.25, 0.25, 0.0])
        self.assertEqual(features[4:6], [0.5, 4.0])


class AddDataTests(unittest.TestCase):
    """Exercises ml.add_data against a throwaway dataset file -- never the
    real ml/data/promoters.data, so a test run can't corrupt the model's
    training data."""

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.dataset_path = Path(self._tmpdir.name) / "dataset.data"
        self.dataset_path.write_text(
            "+,seed_1,ACGTACGTACGTACGTACGT\n"
            "-,seed_2,TTGACATTGACATTGACATT\n"
        )

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_label_aliases_normalize_to_on_disk_form(self):
        for spelling in ("Promoter", "promoter", "+", "P", "1"):
            self.assertEqual(normalize_label(spelling), "+")
        for spelling in ("Non-Promoter", "non-promoter", "-", "N", "0"):
            self.assertEqual(normalize_label(spelling), "-")

    def test_unrecognized_label_raises(self):
        with self.assertRaises(ValueError):
            normalize_label("maybe")

    def test_valid_sequence_is_added_and_persisted(self):
        report = add_sequences(
            [{"label": "Promoter", "sequence": "GATCGATCGATCGATCGATC", "name": "new_1"}],
            dataset_path=self.dataset_path,
        )
        self.assertEqual(len(report["added"]), 1)
        self.assertEqual(report["skipped"], [])
        sequences, labels = load_dataset(self.dataset_path)
        self.assertEqual(len(sequences), 3)
        self.assertIn("GATCGATCGATCGATCGATC", sequences)

    def test_invalid_dna_is_skipped_not_raised(self):
        report = add_sequences(
            [{"label": "Promoter", "sequence": "ACGTXYZACGTACGTACGT", "name": "bad"}],
            dataset_path=self.dataset_path,
        )
        self.assertEqual(report["added"], [])
        self.assertEqual(len(report["skipped"]), 1)
        self.assertIn("Invalid DNA", report["skipped"][0]["reason"])

    def test_too_short_sequence_is_skipped(self):
        report = add_sequences(
            [{"label": "Promoter", "sequence": "ACGT", "name": "short"}],
            dataset_path=self.dataset_path,
        )
        self.assertEqual(report["added"], [])
        self.assertIn("base(s)", report["skipped"][0]["reason"])

    def test_duplicate_of_existing_row_is_skipped(self):
        report = add_sequences(
            [{"label": "Promoter", "sequence": "acgtacgtacgtacgtacgt", "name": "dup"}],
            dataset_path=self.dataset_path,
        )
        self.assertEqual(report["added"], [])
        self.assertIn("duplicate", report["skipped"][0]["reason"])

    def test_duplicate_within_same_batch_is_skipped(self):
        entry = {"label": "Promoter", "sequence": "GGGGCCCCAAAATTTTGGGG", "name": "batch"}
        report = add_sequences([entry, dict(entry, name="batch_again")], dataset_path=self.dataset_path)
        self.assertEqual(len(report["added"]), 1)
        self.assertEqual(len(report["skipped"]), 1)

    def test_bad_label_is_skipped_not_raised(self):
        report = add_sequences(
            [{"label": "unsure", "sequence": "ACGTACGTACGTACGTACGT", "name": "bad_label"}],
            dataset_path=self.dataset_path,
        )
        self.assertEqual(report["added"], [])
        self.assertIn("Unrecognized label", report["skipped"][0]["reason"])

    def test_nothing_written_when_everything_invalid(self):
        before = self.dataset_path.read_text()
        add_sequences([{"label": "nope", "sequence": "AC", "name": "x"}], dataset_path=self.dataset_path)
        self.assertEqual(self.dataset_path.read_text(), before)

    def test_retrained_model_still_predicts_both_classes_are_possible(self):
        # train.py's stratified 80/20 split needs >=3 samples per class to find
        # even 1 test example of each -- add enough rows for that, on top of
        # setUp's 2 seed rows. Not a statistical claim about this tiny set,
        # just confirms the add -> retrain pipeline runs end to end.
        add_sequences(
            [
                {"label": "Promoter", "sequence": "GATCGATCGATCGATCGATC", "name": "p1"},
                {"label": "Promoter", "sequence": "AAACCCGGGTTTAAACCCGG", "name": "p2"},
                {"label": "Non-Promoter", "sequence": "CCCCAAAAGGGGTTTTCCCC", "name": "n1"},
                {"label": "Non-Promoter", "sequence": "TGCATGCATGCATGCATGCA", "name": "n2"},
            ],
            dataset_path=self.dataset_path,
        )
        from ml.train import train
        model_path = Path(self._tmpdir.name) / "model.joblib"
        metrics = train(dataset_path=self.dataset_path, model_path=model_path)
        self.assertIn("accuracy", metrics)
        self.assertTrue(model_path.exists())


class ModelRegistryTests(unittest.TestCase):
    """ml.models is the registry ml.train, ml.compare, and the UI's
    "Compare models" button all build estimators from."""

    def test_registry_has_the_four_expected_models(self):
        self.assertEqual(
            set(available_models()),
            {"random_forest", "logistic_regression", "svm", "knn"},
        )

    def test_default_model_is_random_forest(self):
        # ml.train / ml.predict / the UI's single-prediction path all rely
        # on this staying random_forest -- changing it would silently swap
        # the production model's algorithm.
        self.assertEqual(DEFAULT_MODEL, "random_forest")

    def test_unknown_model_name_raises(self):
        with self.assertRaises(ValueError):
            build_model("not_a_real_model")

    def test_every_registered_model_builds_a_fittable_estimator(self):
        # KNN defaults to k=5, so this needs >= 5 rows -- a size that would
        # silently pass for the other three models but not this one.
        X = [[0.1, 0.2, 10.0], [0.3, 0.1, 12.0], [0.2, 0.4, 9.0], [0.4, 0.3, 11.0],
             [0.15, 0.25, 10.5], [0.35, 0.15, 11.5]]
        y = ["Promoter", "Non-Promoter", "Promoter", "Non-Promoter", "Promoter", "Non-Promoter"]
        for name in available_models():
            model = build_model(name)
            model.fit(X, y)  # must not raise for any registered model
            self.assertIn(model.predict([X[0]])[0], {"Promoter", "Non-Promoter"})


class ModelComparisonTests(unittest.TestCase):
    """Exercises ml.compare against a throwaway dataset -- never the real
    ml/data/promoters.data -- large enough (8 per class) that every
    registered model, including KNN (k=5), has enough training rows after
    the 80/20 split to actually fit."""

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.dataset_path = Path(self._tmpdir.name) / "dataset.data"
        promoters = [
            "TTGACAATTGACAATTGACAATTGACAA", "CATGCATGCATGCATGCATGCATGCATG",
            "GGGAAATTTCCCGGGAAATTTCCCGGGA", "ACGTACGTACGTACGTACGTACGTACGT",
            "TATATATATATATATATATATATATATA", "CCCAAATTTGGGCCCAAATTTGGGCCCA",
            "GATCGATCGATCGATCGATCGATCGATC", "AATTCCGGAATTCCGGAATTCCGGAATT",
        ]
        non_promoters = [
            "CCCGGGCCCGGGCCCGGGCCCGGGCCCG", "GCGCGCGCGCGCGCGCGCGCGCGCGCGC",
            "TTTTAAAATTTTAAAATTTTAAAATTTT", "AGCTAGCTAGCTAGCTAGCTAGCTAGCT",
            "GGATCCGGATCCGGATCCGGATCCGGAT", "CTAGCTAGCTAGCTAGCTAGCTAGCTAG",
            "AACCGGTTAACCGGTTAACCGGTTAACC", "TGCATGCATGCATGCATGCATGCATGCA",
        ]
        rows = [f"+,p{i},{seq}" for i, seq in enumerate(promoters)]
        rows += [f"-,n{i},{seq}" for i, seq in enumerate(non_promoters)]
        self.dataset_path.write_text("\n".join(rows) + "\n")

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_compares_all_four_models_by_default(self):
        results = compare_models(dataset_path=self.dataset_path)
        self.assertEqual({r["model"] for r in results}, set(available_models()))
        for r in results:
            self.assertIn("accuracy", r)
            self.assertIn("roc_auc", r)  # binary task with predict_proba on every registered model
            self.assertIn("display_name", r)

    def test_can_compare_a_subset(self):
        results = compare_models(dataset_path=self.dataset_path, model_names=["svm", "knn"])
        self.assertEqual({r["model"] for r in results}, {"svm", "knn"})

    def test_save_dir_persists_one_file_per_compared_model(self):
        save_dir = Path(self._tmpdir.name) / "saved"
        compare_models(dataset_path=self.dataset_path, model_names=["random_forest", "logistic_regression"],
                        save_dir=save_dir)
        self.assertTrue((save_dir / "random_forest.joblib").exists())
        self.assertTrue((save_dir / "logistic_regression.joblib").exists())

    def test_same_split_is_used_for_every_model(self):
        # Two separate calls with the same random_state must evaluate every
        # model against an identical test set -- otherwise "SVM beat Random
        # Forest" wouldn't be a fair comparison.
        first = compare_models(dataset_path=self.dataset_path, random_state=7)
        second = compare_models(dataset_path=self.dataset_path, random_state=7)
        for r1, r2 in zip(sorted(first, key=lambda r: r["model"]), sorted(second, key=lambda r: r["model"])):
            self.assertEqual(r1["confusion_matrix"], r2["confusion_matrix"])

    def test_train_can_persist_a_non_default_model(self):
        from ml.train import train
        model_path = Path(self._tmpdir.name) / "svm_model.joblib"
        metrics = train(dataset_path=self.dataset_path, model_path=model_path, model_name="svm")
        self.assertIn("accuracy", metrics)
        model = load_model(model_path)
        result = predict_sequence("TTGACAATTGACAATTGACAATTGACAA", model=model)
        self.assertIn(result["prediction"], {"Promoter", "Non-Promoter"})


class ClusterTests(unittest.TestCase):
    """Exercises ml.cluster (unsupervised GC-band grouping) against a
    throwaway dataset -- never the real ml/data/promoters.data."""

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.dataset_path = Path(self._tmpdir.name) / "dataset.data"
        # Three well-separated GC clusters so the fit is unambiguous: this
        # tests that KMeans correctly orders clusters by ascending GC% and
        # that classify_gc_band routes a sequence into the right one, not
        # that KMeans can resolve a noisy real distribution.
        low_gc = ["ATATATATATATATATATATATATATAT", "TATATATATAAATATATATATATATATA"]
        mid_gc = ["ATGCATGCATGCATGCATGCATGCATGC", "GCATGCATGCATGCATGCATGCATGCAT"]
        high_gc = ["GCGCGCGCGCGCGCGCGCGCGCGCGCGCG", "CGCGCGCGCGCGCGCGCGCGCGCGCGCGC"]
        rows = ([f"+,low{i},{s}" for i, s in enumerate(low_gc)]
                + [f"-,mid{i},{s}" for i, s in enumerate(mid_gc)]
                + [f"+,high{i},{s}" for i, s in enumerate(high_gc)])
        self.dataset_path.write_text("\n".join(rows) + "\n")

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_gc_percent_helper(self):
        self.assertEqual(gc_percent("GGCC"), 100.0)
        self.assertEqual(gc_percent("AATT"), 0.0)
        self.assertEqual(gc_percent(""), 0.0)  # guards div-by-zero, not just a typical case

    def test_fit_bands_orders_by_ascending_gc(self):
        _, _, ranges = fit_bands(dataset_path=self.dataset_path)
        self.assertEqual(set(ranges), {"Low GC", "Medium GC", "High GC"})
        self.assertLess(ranges["Low GC"][1], ranges["Medium GC"][0])
        self.assertLess(ranges["Medium GC"][1], ranges["High GC"][0])

    def test_classify_gc_band_routes_correctly(self):
        low = classify_gc_band("ATATATATATATATATATATATATATAT", dataset_path=self.dataset_path)
        high = classify_gc_band("GCGCGCGCGCGCGCGCGCGCGCGCGCGCG", dataset_path=self.dataset_path)
        self.assertEqual(low["band"], "Low GC")
        self.assertEqual(high["band"], "High GC")
        self.assertEqual(set(low), {"gc_percent", "band", "band_ranges"})

    def test_custom_band_count(self):
        _, _, ranges = fit_bands(dataset_path=self.dataset_path, n_bands=2)
        self.assertEqual(len(ranges), 2)


class SimilarityTests(unittest.TestCase):
    """Exercises ml.similarity (nearest-neighbor search in feature space)
    against a throwaway dataset -- never the real ml/data/promoters.data."""

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.dataset_path = Path(self._tmpdir.name) / "dataset.data"
        self.sequences = [
            "TTGACAATTGACAATTGACAATTGACAA", "TTGACAATTGACAATTGACAATTGACAT",
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCG", "GCGCGCGCGCGCGCGCGCGCGCGCGCGCC",
            "ATCGATCGATCGATCGATCGATCGATCG",
        ]
        # load_dataset() requires both classes present (a shared invariant
        # for every caller, not specific to similarity search); labels are
        # otherwise irrelevant here since find_similar() only uses feature
        # distance, not the label, to pick neighbors.
        labels = ["+", "-", "+", "-", "+"]
        rows = [f"{labels[i]},seq{i},{s}" for i, s in enumerate(self.sequences)]
        self.dataset_path.write_text("\n".join(rows) + "\n")

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_returns_k_results_nearest_first(self):
        results = find_similar(self.sequences[0], k=3, dataset_path=self.dataset_path)
        self.assertEqual(len(results), 3)
        distances = [r["distance"] for r in results]
        self.assertEqual(distances, sorted(distances))
        for r in results:
            self.assertEqual(set(r), {"sequence", "label", "distance"})

    def test_exact_match_excludes_itself(self):
        # sequences[0] is byte-identical to a dataset row -- it must never
        # appear in its own neighbor list (a distance-0 self-match is not
        # an interesting result).
        results = find_similar(self.sequences[0], k=len(self.sequences), dataset_path=self.dataset_path)
        self.assertNotIn(self.sequences[0], [r["sequence"] for r in results])

    def test_nearly_identical_sequence_ranks_first(self):
        # sequences[0] and sequences[1] differ by one base -- a query equal
        # to sequences[0] should find sequences[1] closer than the
        # GC-inverted sequences[2]/[3].
        results = find_similar(self.sequences[0], k=1, dataset_path=self.dataset_path)
        self.assertEqual(results[0]["sequence"], self.sequences[1])

    def test_k_larger_than_available_neighbors_does_not_crash(self):
        results = find_similar(self.sequences[0], k=100, dataset_path=self.dataset_path)
        self.assertLessEqual(len(results), len(self.sequences) - 1)  # minus itself


class AnomalyTests(unittest.TestCase):
    """Exercises ml.anomaly (Isolation Forest outlier check) against a
    throwaway dataset -- never the real ml/data/promoters.data."""

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.dataset_path = Path(self._tmpdir.name) / "dataset.data"
        # A consistent family of "normal" sequences (same length, similar,
        # non-repetitive composition) so a clearly different one reads as
        # the outlier -- not a claim about the real dataset's distribution.
        normal = [
            "TTGACAATTGACAATTGACAATTGACAA", "TTGACAATTGACAATTGACAATTGACAT",
            "TTGACAATTGACAATTGACAATTGACAG", "TTGACAATTGACAATTGACAATTGACAC",
            "TTGACTATTGACAATTGACAATTGACAA", "TTGACAATTGACAATCGACAATTGACAA",
        ]
        # load_dataset() requires both classes present (a shared invariant
        # for every caller, not specific to anomaly detection); labels are
        # otherwise irrelevant here since check_sequence() never looks at
        # them -- Isolation Forest only sees the feature vectors.
        labels = ["+", "-", "+", "-", "+", "-"]
        rows = [f"{labels[i]},n{i},{s}" for i, s in enumerate(normal)]
        self.dataset_path.write_text("\n".join(rows) + "\n")

    def tearDown(self):
        self._tmpdir.cleanup()

    def test_result_structure(self):
        r = check_anomaly("TTGACAATTGACAATTGACAATTGACAA", dataset_path=self.dataset_path)
        self.assertEqual(set(r), {"is_outlier", "anomaly_score"})
        self.assertIsInstance(r["is_outlier"], bool)
        self.assertIsInstance(r["anomaly_score"], float)

    def test_in_distribution_sequence_scores_higher_than_extreme_outlier(self):
        # Higher anomaly_score = more typical. A sequence drawn from the
        # same family as training should score higher (more normal) than
        # one with a wildly different, repetitive composition.
        normal_score = check_anomaly("TTGACAATTGACAATTGACAATTGACAG", dataset_path=self.dataset_path)
        outlier_score = check_anomaly("GGGGGGGGGGGGGGGGGGGGGGGGGGGGG", dataset_path=self.dataset_path)
        self.assertGreater(normal_score["anomaly_score"], outlier_score["anomaly_score"])

    def test_contamination_accepts_custom_value(self):
        r = check_anomaly("TTGACAATTGACAATTGACAATTGACAA", dataset_path=self.dataset_path, contamination=0.3)
        self.assertEqual(set(r), {"is_outlier", "anomaly_score"})


class LengthSkewTests(unittest.TestCase):
    def test_uniform_length_balanced_data_is_clean(self):
        seqs = ["ACGT" * 15] * 10 + ["TTGA" * 15] * 10
        labels = ["Promoter"] * 10 + ["Non-Promoter"] * 10
        self.assertEqual(dataset_skew_issues(seqs, labels), [])

    def test_length_confounded_with_label_is_flagged(self):
        seqs = ["ACGT" * 15] * 10 + ["ACGT" * 30] * 10
        labels = ["Promoter"] * 10 + ["Non-Promoter"] * 10
        issues = dataset_skew_issues(seqs, labels)
        self.assertEqual(len(issues), 1)
        self.assertIn("length skew", issues[0])

    def test_overlapping_lengths_are_not_flagged(self):
        rng = random.Random(1)
        seqs = ["A" * rng.randint(40, 80) for _ in range(40)]
        labels = ["Promoter"] * 20 + ["Non-Promoter"] * 20
        self.assertEqual(dataset_skew_issues(seqs, labels), [])

    def test_label_imbalance_is_flagged(self):
        seqs = ["ACGT" * 15] * 20
        labels = ["Promoter"] * 3 + ["Non-Promoter"] * 17
        issues = dataset_skew_issues(seqs, labels)
        self.assertTrue(any("imbalance" in i for i in issues))

    def test_warn_if_skewed_emits_warning(self):
        seqs = ["ACGT" * 15] * 10 + ["ACGT" * 30] * 10
        labels = ["Promoter"] * 10 + ["Non-Promoter"] * 10
        with self.assertWarns(LengthSkewWarning):
            warn_if_skewed(seqs, labels)

    def test_shipped_dataset_has_no_skew(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error", LengthSkewWarning)
            warn_if_skewed(*load_dataset())


def _fake_genome(n=5000, seed=7):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class ImportDataTests(unittest.TestCase):
    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.dir = Path(self._tmpdir.name)
        self.genome = _fake_genome()
        (self.dir / "genome.fa").write_text(">chr test\n" + "\n".join(
            self.genome[i:i + 60] for i in range(0, len(self.genome), 60)) + "\n")
        rows = ["# comment line", "Name\tTSS position\tStrand"]
        self.tss = [(300, "+"), (900, "-"), (1500, "forward"), (2400, "reverse"), (3300, "+")]
        for i, (pos, strand) in enumerate(self.tss):
            rows.append(f"p{i}\t{pos}\t{strand}")
        (self.dir / "promoters.tsv").write_text("\n".join(rows) + "\n")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AACGT"), "ACGTT")

    def test_read_fasta_multi_record(self):
        p = self.dir / "x.fa"
        p.write_text(">a desc\nACGT\nAC\n>b\nGGG\n")
        self.assertEqual(read_fasta(p), [("a desc", "ACGTAC"), ("b", "GGG")])

    def test_read_fasta_rejects_non_fasta(self):
        p = self.dir / "x.fa"
        p.write_text("ACGT\n")
        with self.assertRaises(ValueError):
            read_fasta(p)

    def test_window_at_plus_and_minus_strand(self):
        g = self.genome
        self.assertEqual(window_at(g, 300, "+", 50, 7), g[300 - 1 - 50:300 - 1 + 7])
        self.assertEqual(window_at(g, 900, "-", 50, 7), reverse_complement(g[900 - 7:900 + 50]))
        self.assertEqual(len(window_at(g, 900, "-", 50, 7)), 57)
        self.assertIsNone(window_at(g, 10, "+", 50, 7))  # runs off the start

    def test_regulondb_entries_are_balanced_and_length_matched(self):
        entries, _ = build_regulondb_entries(self.dir / "promoters.tsv", self.dir / "genome.fa")
        pos = [e for e in entries if e["label"] == "Promoter"]
        neg = [e for e in entries if e["label"] == "Non-Promoter"]
        self.assertEqual(len(pos), 5)
        self.assertEqual(len(neg), 5)
        self.assertEqual(sorted(len(e["sequence"]) for e in pos),
                         sorted(len(e["sequence"]) for e in neg))
        self.assertEqual(pos[0]["sequence"], self.genome[299 - 50:299 + 7])

    def test_negatives_stay_away_from_tss(self):
        tss = [p for p, _ in self.tss]
        negs = sample_negatives(self.genome, [57] * 20, tss, exclude=200, seed=3)
        self.assertEqual(len(negs), 20)
        for start, _strand, seq in negs:
            self.assertEqual(len(seq), 57)
            for t in tss:
                self.assertFalse(start - 200 <= t <= start + 56 + 200)

    def test_negatives_are_reproducible(self):
        a = sample_negatives(self.genome, [57] * 5, [300], seed=9)
        b = sample_negatives(self.genome, [57] * 5, [300], seed=9)
        self.assertEqual(a, b)

    def test_sequence_column_used_when_present(self):
        (self.dir / "seq.tsv").write_text("name\tsequence\nx\tACGTACGTACGTACGT\ny\tTTTTGGGGCCCCAAAA\n")
        entries, _ = build_regulondb_entries(self.dir / "seq.tsv", self.dir / "genome.fa", neg_ratio=0)
        self.assertEqual([e["sequence"] for e in entries], ["ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA"])

    def test_neg_ratio_zero_needs_no_genome(self):
        (self.dir / "seq.tsv").write_text("name\tsequence\nx\tACGTACGTACGTACGT\n")
        entries, _ = build_regulondb_entries(self.dir / "seq.tsv", None, neg_ratio=0)
        self.assertEqual(len(entries), 1)

    def test_negatives_without_genome_is_an_error(self):
        (self.dir / "seq.tsv").write_text("name\tsequence\nx\tACGTACGTACGTACGT\n")
        with self.assertRaises(ValueError):
            build_regulondb_entries(self.dir / "seq.tsv", None)

    def test_table_without_usable_columns_is_an_error(self):
        (self.dir / "bad.tsv").write_text("a\tb\n1\t2\n")
        with self.assertRaises(ValueError):
            build_regulondb_entries(self.dir / "bad.tsv", self.dir / "genome.fa")

    def test_fasta_entries_and_dry_run(self):
        p = self.dir / "extra.fa"
        p.write_text(">s1 info\nACGTACGTACGTAC\n>s2\nTTTTTTTTTTTTTT\n")
        entries = build_fasta_entries(p, "Non-Promoter")
        self.assertEqual([e["name"] for e in entries], ["s1", "s2"])
        dataset = self.dir / "d.data"
        dataset.write_text("+,a,ACGTACGTACGTACGT\n-,b,TTGATTGATTGATTGA\n")
        before = dataset.read_text()
        report = add_sequences(entries, dataset_path=dataset, dry_run=True)
        self.assertEqual(len(report["added"]), 2)
        self.assertEqual(dataset.read_text(), before)
        add_sequences(entries, dataset_path=dataset)
        self.assertEqual(len(load_dataset(dataset)[0]), 4)


class PredictionTests(unittest.TestCase):
    def test_model_load_and_prediction_structure(self):
        model = load_model()
        result = predict_sequence("acgt" * 14 + "a", model=model)
        self.assertEqual(set(result), {"prediction", "confidence"})
        self.assertIn(result["prediction"], {"Promoter", "Non-Promoter"})
        self.assertIsInstance(result["confidence"], float)
        self.assertGreaterEqual(result["confidence"], 0.0)
        self.assertLessEqual(result["confidence"], 1.0)


if __name__ == "__main__":
    unittest.main()
