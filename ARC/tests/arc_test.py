import shutil
import sys
import unittest
import unittest.mock
import io

from ARC.classifier import SeqClassifier
from Bio import SeqIO

sys.path.append("..")


class TestClassifier(unittest.TestCase):
	"""Test framework"""

	def test_dependencies(self):
		self.assertIsNotNone(shutil.which("blastp"), "BLAST was not found")
		self.assertIsNotNone(shutil.which("hmmscan"), "HMMER was not found")

	def test_primary(self):
		sc = SeqClassifier()
		# This just tests to make sure nothing odd happens for standard sequences
		sc.classify_seqfile("ARC/tests/crosscheck_s_ids_names.fasta")

	def test_multiprocessing(self):
		sc = SeqClassifier(threads=10)
		sc.classify_seqfile("ARC/tests/test.fa")

	@unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
	def test_empty(self, mock_stdout):
		sc = SeqClassifier()
		sc.classify_seqfile("ARC/tests/empty_test.fa")
		self.assertEqual(mock_stdout.getvalue(), "fake_seq has empty sequence. Skipping sequence.\n")

	@unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
	def test_invalid(self, mock_stdout):
		sc = SeqClassifier()
		sc.classify_seqfile("ARC/tests/invalid_amino_acid.subset")
		self.assertEqual(mock_stdout.getvalue(), "5E94_A_light contains invalid amino acid sequence. Skipping sequence.\n")


if __name__ == '__main__':
    unittest.main()
