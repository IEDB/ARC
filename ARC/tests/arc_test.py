import shutil
import sys
import unittest

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
		# This just tests to make sure nothing odd happens
		sc.classify_seqfile("crosscheck_s_ids_names.fasta")

	def test_empty(self):
		sc = SeqClassifier()
		with self.assertRaises(Exception):
			sc.classify_seqfile("empty_test.fa")


if __name__ == '__main__':
    unittest.main()
