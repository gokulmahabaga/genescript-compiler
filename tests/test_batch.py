"""Focused integration tests for the standalone GeneScript batch CLI."""
import json
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT / "genescript-batch"


class BatchCliTests(unittest.TestCase):
    def invoke(self, request: bytes, ok=True):
        run = subprocess.run([str(BATCH)], input=request, capture_output=True, cwd=ROOT)
        self.assertEqual(run.returncode == 0, ok, run.stderr.decode("utf-8", "replace"))
        parsed = json.loads(run.stdout)
        self.assertIsInstance(parsed, dict)
        return parsed, run

    def test_all_operations_and_input_order(self):
        request = {"operations": ["LENGTH", "GC_CONTENT", "REVERSE", "COMPLEMENT",
                                   "REVERSE_COMPLEMENT", "TRANSLATE"],
                   "sequences": ["ATGCGT", "acgt", "TTT"]}
        result, _ = self.invoke(json.dumps(request).encode())
        rows = result["results"]
        self.assertEqual([r["sequence"] for r in rows], request["sequences"])
        self.assertEqual(rows[0]["LENGTH"], 6)
        self.assertEqual(rows[0]["GC_CONTENT"], 50.0)
        self.assertEqual(rows[0]["REVERSE"], "TGCGTA")
        self.assertEqual(rows[0]["COMPLEMENT"], "TACGCA")
        self.assertEqual(rows[0]["REVERSE_COMPLEMENT"], "ACGCAT")
        self.assertEqual(rows[0]["TRANSLATE"], "MR")
        self.assertEqual(rows[1]["sequence"], "acgt")
        self.assertEqual(rows[1]["REVERSE"], "TGCA")

    def test_invalid_sequence_isolated_and_json_escaped(self):
        request = {"operations": ["LENGTH", "REVERSE"], "sequences": ["ATX", "ACGT", "A\n"]}
        result, _ = self.invoke(json.dumps(request).encode())
        rows = result["results"]
        self.assertIn("error", rows[0])
        self.assertEqual(rows[1]["LENGTH"], 4)
        self.assertEqual(rows[2]["sequence"], "A\n")
        self.assertIn("error", rows[2])

    def test_json_escapes_and_unicode_escape(self):
        raw = b'{"operations":["LENGTH"],"sequences":["\\u0041CGT"]}'
        result, _ = self.invoke(raw)
        self.assertEqual(result["results"][0]["LENGTH"], 4)

    def test_malformed_or_unsupported_requests_return_json_error(self):
        cases = [b'{"operations":["NOPE"],"sequences":["ACGT"]}',
                 b'{"operations":["LENGTH"],"sequences":["ACGT"]',
                 b'{"operations":["LENGTH"],"sequences":["ACGT"],"extra":1}',
                 b'{"operations":["LENGTH","LENGTH"],"sequences":["ACGT"]}']
        for raw in cases:
            with self.subTest(raw=raw):
                parsed, run = self.invoke(raw, ok=False)
                self.assertIn("error", parsed)
                self.assertTrue(run.stderr)

    def test_count_and_size_limits(self):
        too_many = {"operations": ["LENGTH"], "sequences": ["A"] * 101}
        parsed, _ = self.invoke(json.dumps(too_many).encode(), ok=False)
        self.assertIn("error", parsed)
        oversized = {"operations": ["LENGTH"], "sequences": ["A" * 100001]}
        parsed, _ = self.invoke(json.dumps(oversized).encode())
        self.assertIn("error", parsed["results"][0])
        parsed, _ = self.invoke(b" " * (1024 * 1024 + 1), ok=False)
        self.assertIn("error", parsed)

    def test_empty_batches_and_invalid_utf8_rejected(self):
        for request in (b'{"operations":[],"sequences":["A"]}',
                        b'{"operations":["LENGTH"],"sequences":[]}'):
            parsed, _ = self.invoke(request, ok=False)
            self.assertIn("error", parsed)
        parsed, _ = self.invoke(b'{"operations":["LENGTH"],"sequences":["\xff"]}', ok=False)
        self.assertIn("error", parsed)


if __name__ == "__main__":
    unittest.main()
