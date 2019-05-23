"""
Unit testing skeleton
"""

import unittest
from ConstantClassifier import SeqClassifier
from Bio import SeqIO

class TestHMMProfiles(unittest.TestCase):
    def test_light_chain(self):
        light_seqs = list(SeqIO.parse("light_test.fa", "fasta"))

    def test_heavy_chain(self):
        heavy_seqs = list(SeqIO.parse("heavy_test.fa", "fasta")) 

    def test_MHC(self):
        mhc_seqs = list(SeqIO.parse("mhc_test.fa", "fasta"))

    def test_unrelated(self):
        neg_seqs = list(SeqIO.parse("neg_test.fa", "fasta"))

if __name__=='__main__':
    unittest.main()
