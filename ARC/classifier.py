"""
Created on Mon Apr 1 2019

@author: Austin Crinklaw, Swapnil Mahajan

Classifies input sequences into BCR, TCR, or MHC.
Specifies chain type including constant regions.
Contains code from ANARCI and SeqClassifier.
"""
import re
import os
from ARC.mhc_G_domain import mhc_G_domain
import numpy as np
import pandas as pd
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
    self.package_directory = os.path.dirname(os.path.abspath(__file__))
    self.seqfile = seqfile
    self.outfile = outfile
    self.hmm_score_threshold = hmm_score_threshold
    self.mhc_I_hmm = os.path.join(self.package_directory,'data/MHC_HMMs/Pfam_MHC_I.hmm')
    self.mhc_II_alpha_hmm = os.path.join(self.package_directory,'data/MHC_HMMs/Pfam_MHC_II_alpha.hmm')
    self.mhc_II_beta_hmm = os.path.join(self.package_directory, 'data/MHC_HMMs/Pfam_MHC_II_beta.hmm')
    self.mro_file = os.path.join(self.package_directory,'data/MRO/ontology/chain-sequence.tsv')
    self.mro_gdomain_file = os.path.join(self.package_directory,'data/MRO_Gdomain.csv')


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
        
        args = ['hmmscan','-o', hmm_out.name, os.path.join(self.package_directory, "data/HMMs/ALL_AND_C.hmm"), temp_out.name]
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

    ndomains = len(top_hits)
    top_domains = { x["id"].split("_")[1] for x in top_hits }
    #These sets simplify checking for various conditions
    bcr_constant = {"KCC": "IGKC", "LCC": "IGLC", 
                    "HCC": "IGHC", "HC1": "IGHC domain 1", 
                    "HC2": "IGHC domain 2", "HC3":"IGHC domain 3"}
    tcr_constant = {"TRAC": "TRAC", "TRBC": "TRBC", "TRDC": "TRDC", "TRGC": "TRGC"}
    tcr_var = {"A": "TRAV", "B": "TRBV", "G": "TRGV", "D": "TRDV"} 
    bcr_var = {"H": "IGHV", "K": "IGKV", "L": "IGLV"}

    #We have no hits
    if ndomains == 0:
        return None, None
    
    if ndomains == 1:
        #Check for single constant domains
        if top_domains.issubset(bcr_constant):
            #sets don't support indexing so this gets messy
            return ("BCR", bcr_constant[next(iter(top_domains))])
        if top_domains.issubset(tcr_constant):
            return ("TCR", tcr_constant[next(iter(top_domains))])

        #Check for single variable domains
        if top_domains.issubset(tcr_var.keys()):
            return ("TCR", tcr_var[next(iter(top_domains))])
        if top_domains.issubset(bcr_var.keys()):
            return ("BCR", bcr_var[next(iter(top_domains))])

    #Check if the construct is artificial scfv BCR
    if ndomains == 2 and top_domains.issubset(bcr_var.keys()):
        domain_1 = next(iter(top_domains))
        domain_2 = next(iter(top_domains))
        if domain_1 == "H" and domain_2 == "L":
            return ("BCR", "scFv")
        elif domain_1 == "L" and domain_2 == "H":
            return ("BCR", "scFv")
        else:
            return ("BCR", "construct")

    #Check if the construct is artificial scfv TCR 
    if ndomains == 2 and top_domains.issubset(tcr_var.keys()):
        return ("TCR", "TscFv")

    #Check for tandem scfv's and other rare 3+ domain constructs
    if ndomains >= 3 and top_domains.issubset(tcr_var.keys()):
        return ("TCR", "construct")
    if ndomains >= 3 and top_domains.issubset(bcr_var.keys()):
        return ("BCR", "construct")

    #Handle variable with constant
    if any(x in iter(tcr_constant) for x in iter(top_domains)):
        for x in iter(top_domains):
            if x in tcr_var:
                top_domains.remove(x)
                return "TCR", tcr_var[x] + " ," + tcr_constant[next(iter(top_domains))]
    
    if any(x in iter(bcr_constant) for x in iter(top_domains)):
        for x in iter(top_domains):
            if x in bcr_var:
                top_domains.remove(x)
                return "BCR", bcr_var[x] + " ," + bcr_constant[next(iter(top_domains))]

    return None, None

  def assign_Gdomain(self, seq, seq_id=None):
    """
    Returns G domain of a MHC sequence.
    """
    gd = mhc_G_domain(chain1_seq=seq, ch1_id=re.sub(r'[^\w|\.]','',seq_id))
    res = gd.get_g_domain()
    if res:
      return res
    else:
      return None, None
  
  def get_MRO(self):
    """
    Clone or pull MRO GitHub repository.
    """
    mro_path = os.path.join(self.package_directory,'data/MRO')
    if os.path.exists(mro_path):
      print('Updating MRO repository..')
      self.run_cmd('git -C %s pull' % mro_path)
      return
    else:
      print('Getting MRO repository..')
      self.run_cmd('git clone https://github.com/IEDB/MRO.git %s' % mro_path)
      return
  
  def get_MRO_Gdomains(self, mro_TSVfile):
    """
    Returns G doamins of the MRO chain sequences.
    """
    self.get_MRO()
    mro = pd.read_csv(mro_TSVfile, sep='\t', skiprows=[1])
    if os.path.exists(self.mro_gdomain_file)  and os.path.getsize(self.mro_gdomain_file) > 0:
      mro_out=pd.read_csv(self.mro_gdomain_file)
      cnt=mro_out.Label.index[-1]+1
    else:
      mro_out= pd.DataFrame(columns= ['Label','Sequence', 'calc_mhc_class','ch_g_dom'])
      cnt=0
    for i in list(range(mro.shape[0])):
      if pd.isnull(mro.loc[i, 'Sequence']):
        continue
      if mro.loc[i, 'Label'] in list(mro_out['Label']) and mro.loc[i, 'Sequence'] in list(mro_out['Sequence']):
        continue
      mro_out.loc[cnt, 'Label']= mro.loc[i,'Label']
      mro_out.loc[cnt, 'Sequence']= mro.loc[i,'Sequence']
      mro_out.loc[cnt, 'calc_mhc_class'], mro_out.loc[cnt, 'ch_g_dom'] = self.assign_Gdomain(mro.loc[i,'Sequence'], mro.loc[i,'Accession'])
      cnt+=1
    mro_out.to_csv(self.mro_gdomain_file, index=False)
    return mro_out
  
  def get_MRO_allele(self, mro_df, seq, seq_id=None):
    mhc_class, pdb_g_dom=self.assign_Gdomain(seq, seq_id)
    if pdb_g_dom:
      mro_allele= str(list(mro_df[mro_df.fillna({'ch_g_dom' :''}).apply(lambda r : r['ch_g_dom']!='' and (pdb_g_dom in r['ch_g_dom'] or r['ch_g_dom'] in pdb_g_dom) , axis=1)]['Label'])).strip('[]').replace(',','#')
      return mro_allele
    else:
      #print('Unable to assign G domain to the {} chain sequence'.format(seq_id))
      return
  
  def is_MHC(self, sequence, hmm):
    """
    Input: protein sequence and HMM file
    Output: HMM bit score
    """
    score = 0
    # create a temporary file
    fp = tempfile.NamedTemporaryFile(mode="w")
    fp.write('>seq\n')
    fp.write('{}'.format(sequence))
    fp.flush() 
  
    #Find MHC sequences
    args = ['hmmscan', hmm, fp.name ]
    cmd = ' '.join(args)
    output = self.run_cmd(cmd)
    aln = [ line.split() for line in output.splitlines() ]
  
    # Search for score to see if there is a match
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

  def assign_class(self, seq_record):
    """
    Returns BCR, TCR or MHC class and chain type for an input sequence

    @param seq_recored: A biopython sequence record object
    """
#    if self.check_seq(seq_record) == 1:
    with tempfile.NamedTemporaryFile(mode="w") as hmm_out:
        receptor, chain_type = None, None
        self.run_hmmscan(seq_record, hmm_out)
        hmmer_query = SearchIO.read(hmm_out.name, 'hmmer3-text')
        hit_table, top_descriptions = self.parse_hmmer_query(hmmer_query)
        receptor, chain_type = self.get_chain_type(top_descriptions)

        #We have no hits so now we check for MHC, avoid excessive computations this way
        if not receptor or not chain_type:
            mhc_I_score = None
            mhc_I_score = self.is_MHC(str(seq_record.seq), self.mhc_I_hmm)
            if mhc_I_score >= self.hmm_score_threshold:
              return('MHC-I', 'alpha')
            else:
              mhc_II_alpha_score = None
              mhc_II_alpha_score = self.is_MHC(str(seq_record.seq), self.mhc_II_alpha_hmm)
              if mhc_II_alpha_score and mhc_II_alpha_score >= self.hmm_score_threshold:
                return('MHC-II', 'alpha')
              else:
                mhc_II_beta_score = None
                mhc_II_beta_score = self.is_MHC(str(seq_record.seq), self.mhc_II_beta_hmm)
                if mhc_II_beta_score and mhc_II_beta_score >= self.hmm_score_threshold:
                  return('MHC-II', 'beta')
                else:
                  return(None, None)
        else:
          return(receptor, chain_type)

  def classify(self, seq_record, mro_df = None):
    """
    Returns BCR, TCR or MHC class and chain type for an input sequence.
    If sequence is MHC, finds its g-domain and returns its corresponding
    allele

    @param seq_record: a biopython sequence record object
    @param mro_df: dataframe containing the MRO data for allele assignment from g-domain
    """
    g_domain = ""
    calc_mhc_allele = ""
    receptor, chain_type = self.assign_class(seq_record)
    if receptor == "MHC-I" or receptor == "MHC-II":
        g_domain = self.assign_Gdomain(str(seq_record.seq), seq_record.id)
        if mro_df.empty:
            mro_df = self.get_MRO_Gdomains(self.mro_file)
        calc_mhc_allele= self.get_MRO_allele(mro_df, str(seq_record.seq), str(seq_record.description))

    return receptor, chain_type, calc_mhc_allele

  def classify_seqfile(self, seq_file):
    """
    Takes a file of sequences in FASTA format and writes a CSV with their
    receptor type and chain classification

    @param seq_file: the name of a FASTA file of sequences
    """
    seq_records = list(SeqIO.parse(seq_file, "fasta"))
    out = pd.DataFrame(columns=["id", "class", "chain_type", "calc_mhc_allele"])
    cnt = 0
    mro_df = self.get_MRO_Gdomains(self.mro_file)
    for seq in seq_records:
        receptor, chain_type, calc_mhc_allele = self.classify(seq, mro_df)
        out.loc[cnt,'id'] = seq.description
        out.loc[cnt, 'class'] = receptor
        out.loc[cnt, 'chain_type'] = chain_type
        out.loc[cnt, 'calc_mhc_allele'] = calc_mhc_allele
        cnt += 1

    out.to_csv(self.outfile, index=False)