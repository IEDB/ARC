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
import tempfile

class SeqClassifier:
  def __init__(self, seqfile=None, outfile=None, hmm_score_threshold=100):
    """
    Classifies input sequence/s into BCR, TCR or MHC chains, including
    constant regions.

    @param seqfile: Input sequence file in FASTA format
    @param outfile: Name of output file
    @param hmm_score_threshold: Minimum score for a hit against HMM to be significant
    """
    self.seqfile = seqfile
    self.outfile = outfile
    self.hmm_score_threshold = hmm_score_threshold

  def check_seq(self, seq_record):
    """
    Checks validity of an amino acid sequence

    @param seq_record: A biopython sequence record object
    """
    pattern = re.compile(r'[^A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y]', re.IGNORECASE)
    if len(str(seq_record.seq))>0 and not pattern.findall(str(seq_record.seq)):
      return 1
    else:
      print( 'ERROR: ID: {} sequence of has non amino acid characters'.format(seq_record.description))
      return 0

  def run_cmd(self, cmd, input_string=''):
    """
    Run the cmd with input_string as stdin and return output.

    @param cmd: The shell command to run
    @param input_string: String to pass in via stdin
    """
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
               stderr=subprocess.PIPE, universal_newlines=True, close_fds=True,
               env=dict(os.environ, my_env_prop='value'), shell=True)
    
    out, stderr = p.communicate(input=input_string)
    if p.returncode:
        raise Exception('Cmd {} failed: {}'.format(cmd, stderr))
    
    return out
  
  def run_hmmscan(self, seq_record, hmm_out):
    """
    Runs hmmscan from the HMMER3 software suite
    Note: writes temporary output files

    @param seq_record: A biopython sequence record object
    @param hmm_out: tempfile object for hmm output
    """
    with tempfile.NamedTemporaryFile(mode="w") as temp_out:
        if not seq_record.seq:
            print('ERROR: ID: {} sequence was not found'.format(seq_record.description))
            return False
        SeqIO.write(seq_record, temp_out.name, "fasta")
        
        args = ['hmmscan','-o', hmm_out.name, "../data/constant_sequences/hmms/ALL_with_constant.hmm", temp_out.name]
        cmd = (' ').join(args)
        self.run_cmd(cmd, str(seq_record.seq))

        if not(os.path.exists(hmm_out.name) and os.access(hmm_out.name, os.R_OK)):
            print('ERROR: ID {} hmmer out is not found or is not readable.'.format(seq_record.id))
            return False
        if os.path.getsize(hmm_out.name) == 0:
            print('ERROR: ID {} hmmer out is empty. Please add path to hmmer to your environment variables'.format(seq_record.description))
            return False
        return True

  def domains_are_same(self, dom1, dom2):
    """
    Check to see if two domains are overlapping.

    @param dom1:
    @param dom2:

    @return: True or False
    """
    dom1, dom2 = sorted([dom1, dom2], key=lambda x: x.query_start)
    if dom2.query_start >= dom1.query_end:
      return False
    return True

  def parse_hmmer_query(self, query, bit_score_threshold=100):
    """
    The function will identify multiple domains if they have been found and provide the details for the best alignment for each domain.
    This allows the ability to identify single chain fvs and engineered antibody sequences as well as the capability in the future for identifying constant domains. 
    
    @param query: hmmer query object from Biopython
    @param bit_score_threshold: the threshold for which to consider a hit a hit. 

    """
    hit_table = [ ['id', 'description', 'evalue', 'bitscore', 'bias', 
                    'query_start', 'query_end' ] ]

    # Find the best hit for each domain in the sequence.

    top_descriptions, domains,state_vectors = [], [], []

    if query.hsps: # We have some hits
        for hsp in sorted(query.hsps, key=lambda x: x.evalue): # Iterate over the matches of the domains in order of their e-value (most significant first)
            new=True
            if hsp.bitscore >= bit_score_threshold: # Only look at those with hits that are over the threshold bit-score.
                for i in range( len(domains) ): # Check to see if we already have seen the domain
                    if self.domains_are_same( domains[i], hsp ):
                        new = False
                        break      
                hit_table.append( [ hsp.hit_id, hsp.hit_description, hsp.evalue, hsp.bitscore, hsp.bias, hsp.query_start, hsp.query_end] )
                if new: # It is a new domain and this is the best hit. Add it for further processing.
                    domains.append( hsp )
                    top_descriptions.append(  dict( zip(hit_table[0], hit_table[-1]) ) ) # Add the last added to the descriptions list. 

        # Reorder the domains according to the order they appear in the sequence.         
        ordering = sorted( range(len(domains)), key=lambda x: domains[x].query_start)
        domains = [ domains[_] for _ in ordering ]
        top_descriptions = [ top_descriptions[_] for _ in ordering ]         
   
    ndomains = len( domains )
    for i in range(ndomains): # If any significant hits were identified parse and align them to the reference state.
        domains[i].order = i
        species, chain = top_descriptions[i]["id"].split("_")
        top_descriptions[i][ "species"] = species # Reparse
        top_descriptions[i][ "chain_type"] = chain

    return hit_table, top_descriptions


  def get_chain_type(self, top_hits):
    """
    Returns the chain type from the list of top hits

    @param top_hits: the highest scoring hits per domain of a query
    """
    if len(top_hits) == 1: #Only one domain present in query
      print("")
    #Check if it is just variable and constant case

    #Check for scFV case


  def assign_class(self, seq_record):
    """
    Returns BCR, TCR or MHC class and chain type for an input sequence

    @param seq_recored: A biopython sequence record object
    """
    if self.check_seq(seq_record) == 1:
        with tempfile.NamedTemporaryFile(mode="w") as hmm_out:
            self.run_hmmscan(seq_record, hmm_out)
            hmmer_query = SearchIO.read(hmm_out.name, 'hmmer3-text')
            hit_table, top_descriptions = self.parse_hmmer_query(hmmer_query)

