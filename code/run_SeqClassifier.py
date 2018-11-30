# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:10:17 2018

@author: swapnil
"""

from SeqClassifier import SeqClassifier 
import argparse

"""
parser = argparse.ArgumentParser()
parser.add_argument('infile', help = '')
parser.add_argument('outfile', help = '')
args = parser.parse_args()

classification=SeqClassifier(args.infile, args.outfile)
classification.classify()
"""

cl=SeqClassifier()
pdbs=cl.classify_pdb_chains_API()