"""
Unit testing skeleton
"""
import sys
sys.path.append("..") 
import unittest
from ARC.classifier import SeqClassifier
from Bio import SeqIO

class TestClassifier(unittest.TestCase):
	def test_primary(self):
		sc = SeqClassifier()
		sc.classify_seqfile("all_mhcBcrTcr_IEDB.fasta")

	def test_light_chain(self):
		light_seqs = list(SeqIO.parse("light_test.fa", "fasta"))

	def test_heavy_chain(self):
		heavy_seqs = list(SeqIO.parse("heavy_test.fa", "fasta")) 

	def test_MHC(self):
		mhc_seqs = list(SeqIO.parse("mhc_test.fa", "fasta"))

	def test_unrelated(self):
		neg_seqs = list(SeqIO.parse("neg_test.fa", "fasta"))

	def test_empty(self):
		sc = SeqClassifier()
		sc.classify_seqfile("empty_test.fa")

if __name__=='__main__':
    unittest.main()
