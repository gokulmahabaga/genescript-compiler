import tempfile
import unittest
from pathlib import Path

from ml.add_data import add_sequences, normalize_label
from ml.dataset import load_dataset
from ml.features import DINUCLEOTIDES, extract_features
from ml.predict import load_model, predict_sequence
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
