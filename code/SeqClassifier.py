# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 14:25:07 2018

@author: swapnil

Classifies input sequence/s into BCR, TCR or MHC.
"""
from __future__ import print_function
import subprocess
import tempfile 
import re
import os
import sys
import pandas as pd
from Bio import SeqIO
#import random

class SeqClassifier:
  
  def __init__(self,seqfile, outfile='SeqClassifier_output.csv', hmm_score_threshold=25):
    """
    Classifies input sequence/s into BCR, TCR or MHC chains.
    """
    self.seqfile=seqfile
    self.outfile=outfile
    self.hmm_score_threshold=hmm_score_threshold
    
    self.mhc_I_hmm='../data/Pfam_MHC_I.hmm'
    self.mhc_II_alpha_hmm='../data/Pfam_MHC_II_alpha.hmm'
    self.mhc_II_beta_hmm='../data/Pfam_MHC_II_beta.hmm'
  
  # returns 0 if unusual character in sequeunce
  def check_seq(self,seq_record):
    print(seq_record.id)
    pattern=re.compile(r'[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]', re.IGNORECASE)
    if len(str(seq_record.seq))>0 and not pattern.findall(str(seq_record.seq)):
      return 1
    else:
      print( 'ERROR: ID: {} sequence of has non amino acid characters'.format(seq_record.description))
      return 0
  
  def run_cmd(self,cmd, input_string=''):
    """Run the cmd with input_string as stdin and return output."""
    
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
               stderr=subprocess.PIPE, universal_newlines=True, close_fds=True,
               env=dict(os.environ, my_env_prop='value'), shell=True)
    out, stderr = p.communicate(input=input_string)
    if p.returncode:
      raise Exception('Cmd {} failed: {}'.format(cmd[0], stderr))
    
    return out
  
  def run_anarci(self,seq_record, imgtfile):
    """
    Runs ANARCI on BCR/TCR chain sequence and write IMGT numbering file.
    """
    #hmmerpath='/home/swapnil/Softwares/hmmer-3.1b2/binaries'
    if not seq_record.seq:
      print( 'ERROR: ID: {} sequence was not found.'.format(seq_record.description))
    
    # Use chopped sequence for chain 1 for ANARCI
    #args=['ANARCI','-s', 'i', '--assign_germline', '--hmmerpath', hmmerpath,'-i', seqfile,'-o',imgtfile]
    args=['ANARCI','-s', 'i', '-i', str(seq_record.seq),'-o',imgtfile]
    cmd = ' '.join(args)
    #print(cmd)
    out=self.run_cmd(cmd)
    
    if not (os.path.exists(imgtfile)  and os.access(imgtfile, os.R_OK)):
      print( 'ERROR: ID: {} imgtfile is not found or is not readable.'.format(seq_record.description))
      return
    if os.path.getsize(imgtfile) ==0:
      print('ERROR: ID: {} imgtfile is empty. ANARCI may have failed to run.\n Please add hmmerpath in /etc/environment file.'.format(seq_record.description))
      return
    return 1
      
  
  # To get BCR and TCR chain_type
  def get_chain_type(self,seq_record):
    """
    Returns BCR or TCR chain_types
    """
    receptor=None
    chain_type=None
    c_type=[]
    cnt_var_domain=0
    flag_var=0
    pattern=re.compile(r'\#\|species\|chain_type')
    #pattern_germline=re.compile(r'\#\|species\|v_gene')
  
    imgtfile=str(seq_record.id)+'.imgt'
    imgt_out=self.run_anarci(seq_record, imgtfile)
    if not imgt_out:
      return
    with open(imgtfile) as imgt_num:
      for line in imgt_num:
        if pattern.findall(line) and cnt_var_domain==0:
          flag_var=1
          cnt_var_domain+=1
          continue
        if pattern.findall(line) and cnt_var_domain==1:
          flag_var=1
          cnt_var_domain+=1
          continue
        if flag_var==1 and cnt_var_domain==1:
          val=line.split('|')
          # BCR
          if val[2]=='H':
            chain_type='heavy'
            c_type.append('H')
            receptor='BCR'
          if val[2]=='L':
            chain_type='lambda'
            c_type.append('L')
            receptor='BCR'
          if val[2]=='K':
            chain_type='kappa'
            c_type.append('K')
            receptor='BCR'
  
          # TCR        
          if val[2]=='A':
            chain_type='alpha'
            c_type.append('A')
            receptor='TCR'
          if val[2]=='B':
            chain_type='beta'
            c_type.append('B')
            receptor='TCR'
          if val[2]=='G':
            chain_type='gamma'
            c_type.append('G')
            receptor='TCR'
          if val[2]=='D':
            chain_type='delta'
            c_type.append('D')
            receptor='TCR'
          flag_var=0
  
        if flag_var==1 and cnt_var_domain==2:
          val=line.split('|')
          if chain_type=='heavy' and val[2]=='H':
            #chain_type='dual VH construct'
            chain_type='construct'
            c_type.append('H')
          elif (chain_type=='light') and (val[2]=='K' or val[2]=='L'):
            #chain_type='dual VL construct'
            chain_type='construct'
            c_type.append('L')
          elif chain_type=='heavy' and (val[2]=='K' or val[2]=='L'):
            chain_type='scFv'
            c_type.append('L')
          elif (chain_type=='light') and val[2]=='H':
            chain_type='scFv'
            c_type.append('H')
          # TCR
          if chain_type=='alpha' and val[2]=='A':
            #chain_type='dual VA construct'
            chain_type='construct'
            c_type.append('A')
          elif chain_type=='beta' and val[2]=='B':
            #chain_type='dual VB construct'
            chain_type='construct'
            c_type.append('B')
          if chain_type=='gamma' and val[2]=='G':
            #chain_type='dual VG construct'
            chain_type='construct'
            c_type.append('G')
          elif chain_type=='delta' and val[2]=='D':
            #chain_type='dual VD construct'
            chain_type='construct'
            c_type.append('D')
          elif chain_type=='alpha' and val[2]=='B':
            chain_type='TscFv'
            c_type.append('B')
          elif chain_type=='beta' and val[2]=='A':
            chain_type='TscFv'
            c_type.append('A')
          elif chain_type=='gamma' and val[2]=='D':
            chain_type='TscFv'
            c_type.append('D')
          elif chain_type=='delta' and val[2]=='G':
            chain_type='TscFv'
            c_type.append('G')
    os.remove(imgtfile)
    return (receptor,chain_type,c_type)
  
  def is_MHC(self,sequence,hmm):
    """
    Input: protein sequence and HMM file
    Output: HMM bit score
    """
    # create a temporary file
    fp = tempfile.NamedTemporaryFile()
    fp.write('>seq\n')
    fp.write('{}'.format(sequence))
    fp.flush() 
  
    #Find MHC sequences
    args = ['hmmsearch', hmm, fp.name ]
    cmd = ' '.join(args)
    output = self.run_cmd(cmd)
    aln = [ line.split() for line in output.splitlines() ]
  
    # Search for score to see if there is a match
    score = None
    for i, line in enumerate(aln):
      if line[0:3] == ['E-value', 'score', 'bias'] and aln[i+2]:
        try:
          E_value = float(aln[i+2][0])
          score = float(aln[i+2][1])
          break
        except ValueError:
          E_value = float(aln[i+3][0])
          score = float(aln[i+3][1])
          break
    
    # close the file. When the file is closed it will be removed.
    fp.close()  
  
    return score
  
  def classify(self):
    """
    Classifies input sequence/s into BCR, TCR or MHC.
    """
    
    sequences = SeqIO.parse(self.seqfile, 'fasta')
    out=pd.DataFrame()
    cnt=0
    
    for seq in sequences:
      #print(seq.description)
      #print(seq.seq)
      chain_type=None
      c_type=None
      receptor=None
      out.loc[cnt,'ID']="'"+str(seq.description)
      if self.check_seq(seq)==1:
        receptor, chain_type,c_type=self.get_chain_type(seq)
        #print(chain_type, str(c_type))
        if chain_type:
          out.loc[cnt,'class']=receptor
          out.loc[cnt,'chain_type']=str(chain_type)
        else:
          chain_type=None
          c_type=None
          mhc_I_score=None
          
          mhc_I_score=self.is_MHC(str(seq.seq), self.mhc_I_hmm)
          if mhc_I_score >= self.hmm_score_threshold:
            out.loc[cnt,'class']='MHC-I'
            out.loc[cnt,'chain_type']='alpha'
          else:
            mhc_II_alpha_score=self.is_MHC(str(seq.seq), self.mhc_II_alpha_hmm)
            if mhc_II_alpha_score >= self.hmm_score_threshold:
              out.loc[cnt,'class']='MHC-II'
              out.loc[cnt,'chain_type']='alpha'
            else:
              mhc_II_beta_score=self.is_MHC(str(seq.seq), self.mhc_II_beta_hmm)
              if mhc_II_beta_score >= self.hmm_score_threshold:
                out.loc[cnt,'class']='MHC-II'
                out.loc[cnt,'chain_type']='beta'
              else:
                cnt+=1
                continue
      cnt+=1
    out.to_csv(self.outfile, index=False)
      
      