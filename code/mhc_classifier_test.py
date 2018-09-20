# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 13:48:12 2018

@author: swapnil
"""
from __future__ import print_function
import pandas as pd

# returns 0 if unusual character in sequeunce
def check_seq(seq):
  pattern=re.compile(r'[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]', re.IGNORECASE)
  if len(seq)>0 and not pattern.findall(seq):
    return 1
  else:
    print( 'ERROR: sequence {} has non amino acid characters'.format(seq))
    return 0

def mro_tsvTofasta(mro_file, outfile=None):
  """
  Converts MRO chain_sequence.tsv file to fasta format file.
  """
  outfile='../data/mro_chain_seq.fasta'
  mro_df=pd.read_csv(mro_file, sep='\t', skiprows=[1])
  
  out=open(outfile, 'wt')
  
  for i in range(mro_df.shape[0]):
    if not pd.isnull(mro_df.loc[i,'Sequence']):
      print('>'+mro_df.loc[i,'Label'].replace(' chain','')+'\n'+mro_df.loc[i,'Sequence'].replace('\n','').replace('\r',''), file=out)
  
  out.close()


# Converts MRO chain_sequence.tsv file to fasta format file.
mro_tsvTofasta('../data/MRO-master_09202018/ontology/chain-sequence.tsv')

