# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 16:16:48 2018

@author: swapnil

Input: Full MHC chain proteins sequence
Output: mhci_G or mhcii_G domain sequence.
"""

from __future__ import print_function, division
import sys, os, subprocess, re, copy
import pandas as pd
import random

class mhc_G_domain:
  def __init__(self, chain1_seq, chain2_seq=None, ch1_id=None, ch2_id=None):
    """
    Input: Full MHC chain proteins sequence
    Output: mhci_G, mhcii_alpha or mhcii_beta domain sequence.
    
    External softwares needed: BLAST
    """
    
    self.chain1_seq= chain1_seq
    self.chain2_seq=chain2_seq
    self.ch1_id= ch1_id
    self.ch2_id= ch2_id
    
    self.e_val=pow(10,-10)
    self.max_out_seq =10
    self.hit_coverage= 0.85
    self.rand= str(random.randint(0,1000000))
  
  def check_seq(self, seq):
    """
    Returns 0 if unusual character in sequeunce.
    """
    pattern=re.compile(r'[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]', re.IGNORECASE)
    if len(seq)>0 and not pattern.findall(seq):
      return 1
    else:
      print('ERROR: {} sequence has non amino acid characters'.format(self.ch1_id))
      return 0
      
  def run_blast(self, seqfile, db, blastout, qcov=None, e_val=None):
    """
    Runs Blast for the give sequence file against the given database and return a output file in CSV format.
    Command:
    blastp -query ../test_seq.fas -db G_ALPHA1.fasta -outfmt '10 std slen' -max_target_seqs 10 -evalue 1e-10 > ./test_alpha1.out
    
    blastp -query ../test_seq.fas -db G_ALPHA1.fasta -outfmt '10 std slen' -max_target_seqs 10 -evalue 1e-10 -qcov_hsp_perc [90] > ./test_alpha1.out
    
    blastout columns:
    'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen'
    """
    if not e_val:
      e_val=self.e_val
    if not qcov:
      args=['blastp','-query', seqfile, '-db', db, '-outfmt', '"10 std slen"',\
          '-max_target_seqs', str(self.max_out_seq), '-evalue' , str(e_val), '>', blastout]
    else:
      args=['blastp','-query', seqfile, '-db', db, '-outfmt', '"10 std slen"',\
          '-max_target_seqs', str(self.max_out_seq), '-evalue' , str(e_val), '-qcov_hsp_perc',qcov,'>', blastout]
    cmd = ' '.join(args)
    #print(cmd)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,\
                            env=dict(os.environ, my_env_prop='value'), shell=True)
    out, err=proc.communicate()
    return blastout
    
  def blast_all(self, seqfile):
    """
    Runs blast for a given seq file against all G domains (alpha1, alpha2, alpha, beta).
    """
    g_alpha1_db= '../data/blastdb/G_ALPHA1.fasta'
    g_alpha1_db=os.path.join(os.path.dirname(__file__),  g_alpha1_db)
    g_alpha2_db= '../data/blastdb/G_ALPHA2.fasta'
    g_alpha2_db=os.path.join(os.path.dirname(__file__),  g_alpha2_db)
    g_alpha_db= '../data/blastdb/G_ALPHA.fasta'
    g_alpha_db=os.path.join(os.path.dirname(__file__),  g_alpha_db)
    g_beta_db= '../data/blastdb/G_BETA.fasta'
    g_beta_db=os.path.join(os.path.dirname(__file__),  g_beta_db)
    
    seq_id= str(seqfile.split('.')[0])
    alpha1_blout= self.run_blast(seqfile, g_alpha1_db, seq_id+'_alpha1'+self.rand+'.out')
    alpha2_blout= self.run_blast(seqfile, g_alpha2_db, seq_id+'_alpha2'+self.rand+'.out')
    alpha_blout= self.run_blast(seqfile, g_alpha_db, seq_id+'_alpha'+self.rand+'.out')
    beta_blout= self.run_blast(seqfile, g_beta_db, seq_id+'_beta'+self.rand+'.out')
    
    return alpha1_blout, alpha2_blout, alpha_blout, beta_blout
  
  def get_domain(self, blastout):
    """
    Identify the domain boundaries from blast results.
    """
    columns='qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen'.split(' ')
    
    df = pd.read_csv(blastout)
    df.columns= columns
    for i in range(df.shape[0]):
      hit_cover = (df.loc[i, 'send'] - df.loc[i, 'sstart'] + 1) / df.loc[i, 'slen']
      if hit_cover >= self.hit_coverage:
        return df.loc[i, 'qstart'], df.loc[i, 'qend']
    print('WARNING: No hits found. May be {} seq is partial.'.format(self.ch1_id))
    return None
        
     
  def get_all_domains(self, alpha1_blout, alpha2_blout, alpha_blout, beta_blout):
    """
    Returns mhc_class and domain boundaries.
    """
    
    if os.path.getsize(alpha1_blout) >0:
      dom1_st_end = self.get_domain(alpha1_blout)
      if not dom1_st_end:
        print('No alpha1 domain was found in seq {}.'.format(self.ch1_id))
        return None
      if os.path.getsize(alpha2_blout) >0:
        dom2_st_end = self.get_domain(alpha2_blout)
        if dom2_st_end:
          return 'I', [dom1_st_end[0], dom1_st_end[1], dom2_st_end[0], dom2_st_end[1]]
        else:
          print('No alpha2 domain was found in seq {}.'.format(self.ch1_id))
          return None
      
    if os.path.getsize(alpha_blout) >0:
      dom1_st_end = self.get_domain(alpha_blout)
      if not dom1_st_end:
        print('No alpha domain was found in seq {}.'.format(self.ch1_id))
        return None
      return 'IIa', dom1_st_end
    if os.path.getsize(beta_blout) >0:
      dom1_st_end = self.get_domain(beta_blout)
      if not dom1_st_end:
        print('No beta domain was found in seq {}.'.format(self.ch1_id))
        return None
      return 'IIb', dom1_st_end
    
    print('ERROR: No G domain was not found in {}.'.format(self.ch1_id))
    return None
  
  def check_b2m(self, seqfile):
    """
    Check if the seq is beta-2-microglobulin. Returns 0 if seq is b2m.
    """
    #g_dom_db= '../data/blastdb/IMGT_G_domain_Species.fasta'
    #g_dom_db=os.path.join(os.path.dirname(__file__),  g_dom_db)
    b2m_db= '../data/blastdb/b2m.fasta'
    b2m_db=os.path.join(os.path.dirname(__file__),  b2m_db)
    
    #blout= self.run_blast(seqfile, g_dom_db, 'b2m_'+self.rand+'.out')
    blout= self.run_blast(seqfile, b2m_db, 'b2m_'+self.rand+'.out', '90', pow(10,-50))
    if os.path.getsize(blout) >0:
      os.remove(blout)
      #return 1
      return 0
    else:
      os.remove(blout)
      #print('WARNING: chain might not be a MHC chain. Ignore this message if chain is beta-2-microglobulin.')
      return 1
      #return 0
    
  def get_subseq(self, seq, dom_st, dom_end):
    """
    Returns a subseq of give seq using start and end positions.
    """
    subseq=''
    subseq=seq[dom_st-1:dom_end]
    return(subseq)
    
  def chain_g_domain(self, seq, seq_id=None):
    """
    Returns G-domain sequence of the given sequence.
    """
    if not seq_id:
      seq_id='tmp_seq'+self.rand
    
    g_domain=''
    inpfile=os.path.join(os.path.dirname(__file__),  str(seq_id)+'.fasta')
    inp=open(inpfile, 'wt')
    print('>'+str(seq_id), file=inp)
    print(seq, file=inp)
    inp.close()
    
    # check if seq is b2m or mhc chain
    b2m=self.check_b2m(inpfile)
    if b2m==0:
      os.remove(inpfile)
      return None
    
    alpha1_blout, alpha2_blout, alpha_blout, beta_blout = self.blast_all(inpfile)
    res = self.get_all_domains(alpha1_blout, alpha2_blout, alpha_blout, beta_blout)
    if res:
      mhc_class, st_end=res
    
      if mhc_class and st_end:
        if mhc_class=='I':
          if st_end[2]-st_end[1] !=1:
            st_end[2]=st_end[1]+1
          dom1= self.get_subseq(seq, st_end[0], st_end[1])
          dom2= self.get_subseq(seq, st_end[2], st_end[3])
          g_domain= dom1+dom2
        else:
          g_domain=self.get_subseq(seq, st_end[0], st_end[1])
  
        os.remove(inpfile)
        os.remove(alpha1_blout)
        os.remove(alpha2_blout)
        os.remove(alpha_blout)
        os.remove(beta_blout)
        return mhc_class, g_domain
    else:
      os.remove(inpfile)
      os.remove(alpha1_blout)
      os.remove(alpha2_blout)
      os.remove(alpha_blout)
      os.remove(beta_blout)
      print('No G domain was found in seq {}.'.format(self.ch1_id))
      return None
  
  def get_g_domain(self):
    """
    """
    mhc_class1=None; dom1_seq=None
    mhc_class2=None; dom2_seq=None
    g_dom1=None; g_dom2=None
    
    if self.chain1_seq:
      if self.check_seq(self.chain1_seq)==0:
        print('ERROR: {} chain1_seq has non standard amino acids.'.format(self.ch1_id))
        return None
      res1= self.chain_g_domain(self.chain1_seq, self.ch1_id)
      if res1:
        mhc_class1, dom1_seq =res1
        if mhc_class1=='I':
          return mhc_class1, dom1_seq
        if mhc_class1=='IIa':
          g_dom1=dom1_seq
          if not self.chain2_seq:
            return mhc_class1, dom1_seq
        if mhc_class1=='IIb':
          g_dom2=dom1_seq
          if not self.chain2_seq:
            return mhc_class1, dom1_seq
    else:
      print('ERROR: {} chain1_seq has non standard amino acids.'.format(self.ch1_id))
      return None
    
    if self.chain2_seq:
      if self.check_seq(self.chain2_seq)==0:
        print('ERROR: {} chain2_seq has non standard amino acids.'.format(self.ch2_id))
        return None
      res2= self.chain_g_domain(self.chain2_seq, self.ch2_id)
      if res2:
        mhc_class2, dom2_seq =res2
        if mhc_class2=='I':
          return mhc_class2, dom2_seq
        if mhc_class2=='IIa':
          g_dom1=dom2_seq
          if not self.chain1_seq:
            return mhc_class2, dom2_seq
        if mhc_class2=='IIb':
          g_dom2=dom2_seq
          if not self.chain1_seq:
            return mhc_class2, dom2_seq
          
    if self.chain1_seq and self.chain2_seq:
      return 'II', g_dom1+g_dom2
      