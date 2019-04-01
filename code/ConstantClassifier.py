"""
Created on Mon Apr 1 2019

@author: Austin Crinklaw

Classifies input sequences into BCR, TCR, or MHC.
Specifies chain type including constant regions
"""
import re
import os
from Bio import SearchIO
from Bio import SeqIO
import datetime
import subprocess

class SeqClassifier:
  def __init__(self, seqfile=None, outfile=None, hmm_score_threshold=100):
    """
    Classifies input sequence/s into BCR, TCR or MHC chains, including
    constant regions.

    Parameters
    ----------
    seqfile: FASTA formatted sequence file

    outfile: Output file in .csv format
    """
    self.seqfile = seqfile
    self.outfile = outfile
    self.hmm_score_threshold = hmm_score_threshold

    self.hmm_file = os.path.join(os.path.dirname(__file__),'../data/constant_sequences/hmms/ALL_with_constant.hmm')
    now = datetime.datetime.now()
    if not self.outfile:
      self.outfile = os.path.join(os.path.dirname(__file__),'../out/SeqClassifier_output_'+now.strftime("%d%b%Y")+'.csv')

  # returns 0 if unusual character in sequeunce
  def check_seq(self, seq_record):
    pattern = re.compile(r'[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]', re.IGNORECASE)
    if len(str(seq_record.seq))>0 and not pattern.findall(str(seq_record.seq)):
      return 1
    else:
      print( 'ERROR: ID: {} sequence of has non amino acid characters'.format(seq_record.description))
      return 0

  def run_cmd(self, cmd, input_string=''):
    """
    Run the cmd with input_string as stdin and return output.
    """
    print("input string is " + input_string)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
               stderr=subprocess.PIPE, universal_newlines=True, close_fds=True,
               env=dict(os.environ, my_env_prop='value'), shell=True)
    out, stderr = p.communicate(input=input_string)
    if p.returncode:
      raise Exception('Cmd {} failed: {}'.format(cmd, stderr))

    return out
  
  def run_hmmscan(self, seq_record):
      hmm_out = seq_record.id + ".txt"
      if not seq_record.seq:
          print('ERROR: ID: {} sequence was not found'.format(seq_record.description))
          return 0 
      args = ['hmmscan','-o', hmm_out, "/home/austin/classifier_tool/data/constant_sequences/hmms/ALL_with_constant.hmm", '-']
      cmd = (' ').join(args)
      print("CMD is " + cmd)
      self.run_cmd(cmd, str(seq_record.seq))

      if not(os.path.exists(hmm_out) and os.access(hmm_out, os.R_OK)):
          print('ERROR: ID {} hmmer out is not found or is not readable.'.format(seq_record.description))
          return 0
      if os.path.getsize(hmm_out) == 0:
          print('ERROR: ID {} hmmer out is empty. Please add path to hmmer to your environment variables'.format(seq_record.description))
          return 0
      return 1

  def get_chain_type(self, seq_record):
    """
    Returns BCR or TCR chain types
    """
    receptor = None
    chain_type = None
    
    hmm_out = self.run_hmmscan(seq_record)
    if not hmm_out:
        return (receptor, chain_type)
    scan_results = list(SearchIO.parse((seq_record.id+'.txt', 'hmmer3-text')))
    sig_hits = set()
    for x in scan_results:
        for hit in x.hits:
            if hit.bitscore > hmm_score_threshold:
                sig_hits.append(hit.description.split("_")[1])

    return sig_hits

  def assign_class(self, seq):
    """
    Returns BCR, TCR or MHC class and chain type for an input sequence
    """
    if self.check_seq(seq) == 1:
        sig_hits = self.get_chain_type(seq)
    
    if chain_type:
        return sig_hits
    else:
        return (None, None)
