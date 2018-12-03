# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 16:10:17 2018

@author: swapnil
"""

from SeqClassifier import SeqClassifier 
import pandas as pd
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
#pdbs=cl.classify_pdb_chains_API()
#api_res=cl.get_pdb_seq_API()
cl.seqfile='../out/PDB_29Nov2018.csv'
api_res=pd.read_csv('../out/PDB_29Nov2018.csv')
#pdbseq=self.create_SeqRecord(api_res)
#iedb_PDBs=cl.get_IEDB_PDBs()
iedb_PDBs=pd.read_csv('../out/IEDB_PDBs_PubMed_IDs_29Nov2018.csv')
#mro_out=cl.get_MRO_Gdomains(cl.mro_file)
mro_out=pd.read_csv('../out/MRO_Gdomain_29Nov2018.csv')
new_pdbs=cl.get_PDBs_classication(set(iedb_PDBs['pdb_id']), set(api_res['structureId']))
pdb_api=api_res[api_res.structureId.isin(new_pdbs)]
pdbseq=cl.create_SeqRecord(pdb_api)
cl.classify(pdbseq, mro_out, pdb_api)