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

def bcr_xlsTofasta(bcr_file, outfile=None):
  """
  Converts IEDB 3D benchmark csv file to fasta format file.
  """
  outfile='../data/bcr_chain_seq.fasta'
  bcr_df=pd.read_excel(bcr_file)
  
  out=open(outfile, 'wt')
  
  for i in range(bcr_df.shape[0]):
    if not pd.isnull(bcr_df.loc[i,'chain1_full_seq']):
      print('>'+bcr_df.loc[i,'pdb_id']+'_'+str(bcr_df.loc[i, 'ab_c1_pdb_chain'])+'_'+bcr_df.loc[i, 'chain1_type']+'\n'+bcr_df.loc[i,'chain1_full_seq'].replace('\n','').replace('\r',''), file=out)
    
    if not pd.isnull(bcr_df.loc[i,'chain2_full_seq']):
      print('>'+bcr_df.loc[i,'pdb_id']+'_'+str(bcr_df.loc[i, 'ab_c2_pdb_chain'])+'_'+bcr_df.loc[i, 'chain2_type']+'\n'+bcr_df.loc[i,'chain2_full_seq'].replace('\n','').replace('\r',''), file=out)
  
  out.close()

def tcr_xlsTofasta(tcr_file, outfile=None):
  """
  Converts IEDB 3D benchmark csv file to fasta format file.
  """
  outfile='../data/tcr_chain_seq.fasta'
  tcr_df=pd.read_excel(tcr_file)
  
  out=open(outfile, 'wt')
  
  for i in range(tcr_df.shape[0]):
    if not pd.isnull(tcr_df.loc[i,'chain1_full_seq']):
      print('>'+tcr_df.loc[i,'pdb_id']+'_'+str(tcr_df.loc[i, 'tcr_c1_pdb_chain'])+'_'+tcr_df.loc[i, 'chain1_type']+'\n'+tcr_df.loc[i,'chain1_full_seq'].replace('\n','').replace('\r',''), file=out)
    
    if not pd.isnull(tcr_df.loc[i,'chain2_full_seq']):
      print('>'+tcr_df.loc[i,'pdb_id']+'_'+str(tcr_df.loc[i, 'tcr_c2_pdb_chain'])+'_'+tcr_df.loc[i, 'chain2_type']+'\n'+tcr_df.loc[i,'chain2_full_seq'].replace('\n','').replace('\r',''), file=out)
  
  out.close()

def mhc_xlsTofasta(mhc_file, outfile=None):
  """
  Converts IEDB 3D benchmark csv file to fasta format file.
  """
  outfile='../data/mhc_chain_seq.fasta'
  mhc_df=pd.read_excel(mhc_file)
  
  out=open(outfile, 'wt')
  
  for i in range(mhc_df.shape[0]):
    if not pd.isnull(mhc_df.loc[i,'mhc_chain_i_seq']):
      print('>'+str(mhc_df.loc[i,'pdb_id'])+'_'+str(mhc_df.loc[i, 'mhc_c1_pdb_chain'])+'_'+str(mhc_df.loc[i, 'mhc_chain1'])+';'+str(mhc_df.loc[i,'mhc_class'])+';'+str(mhc_df.loc[i,'mhc_allele_name'])+'\n'+mhc_df.loc[i,'mhc_chain_i_seq'].replace('\n','').replace('\r',''), file=out)
    
    if not pd.isnull(mhc_df.loc[i,'mhc_chain_ii_seq']):
      print('>'+str(mhc_df.loc[i,'pdb_id'])+'_'+str(mhc_df.loc[i, 'mhc_c2_pdb_chain'])+'_'+str(mhc_df.loc[i, 'mhc_chain2'])+';'+str(mhc_df.loc[i,'mhc_class'])+';'+str(mhc_df.loc[i,'mhc_allele_name'])+'\n'+mhc_df.loc[i,'mhc_chain_ii_seq'].replace('\n','').replace('\r',''), file=out)
  
  out.close()
# Converts MRO chain_sequence.tsv file to fasta format file.
#mro_tsvTofasta('../data/MRO-master_09202018/ontology/chain-sequence.tsv')
#bcr_xlsTofasta('../data/BCR_benchmark_r3.5_l50_t0.7_06Jul2018.xlsx')
#tcr_xlsTofasta('../data/TCR_benchmark_r3.5_l5_t0.7_11Jul2018.xlsx')
mhc_xlsTofasta('../data/MHC_benchmark_r3.5_l5_t0.7_11Jul2018.xlsx')