import unittest

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
