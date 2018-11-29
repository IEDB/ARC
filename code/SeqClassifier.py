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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import urllib2
import urllib
import gzip
import datetime
#import random

class SeqClassifier:
  
  def __init__(self,seqfile=None, outfile=None, hmm_score_threshold=25, length_threshold=50):
    """
    Classifies input sequence/s into BCR, TCR or MHC chains.
    
    Parameters
    ----------
    
    seqfile: Fasta formatted sequence file (Optional: Name is not needed for PDB sequence classification).
    
    """
    self.seqfile=seqfile
    self.outfile=outfile
    self.hmm_score_threshold=hmm_score_threshold
    self.length_threshold=length_threshold
    self.mro_file=os.path.abspath('../data/MRO/ontology/chain-sequence.tsv')
    self.mhc_I_hmm='../data/Pfam_MHC_I.hmm'
    self.mhc_II_alpha_hmm='../data/Pfam_MHC_II_alpha.hmm'
    self.mhc_II_beta_hmm='../data/Pfam_MHC_II_beta.hmm'
    
    now=datetime.datetime.now()
    if not self.seqfile:
      self.seqfile=os.path.abspath('../out/PDB_'+now.strftime("%d%b%Y")+'.csv')
    if not self.outfile:
      self.outfile=os.path.abspath('../out/SeqClassifier_output_'+now.strftime("%d%b%Y")+'.csv')
    self.sql_results_file= os.path.abspath('../out/IEDB_PDBs_PubMed_IDs_'+now.strftime("%d%b%Y")+'.csv')
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
  
  def geturl(self, url, outpath):
    """
    Retrieves file database from any url text file given.
    Will save in output designated by second variable.  
    """
    try:
      urllib.urlretrieve(url, outpath)
      print("Data retrieved at: {}".format(outpath))
    except:
      print("Error in retrieving url {}, please check paths.".format(url))
    urllib.urlcleanup()
    return
  
  def create_SeqRecord(self, pdb_api_file):
    """
    Create sequences in biopython format from PDB API CSV file.
    """
    df=pd.read_csv(pdb_api_file)
    sequences=[]
    for i in range(df.shape[0]):
      record=None
      if df.loc[i,'sequence'] and not pd.isnull(df.loc[i,'sequence']) \
      and len(df.loc[i,'sequence']) >=self.length_threshold:
        record=SeqRecord(Seq(str(df.loc[i,'sequence']), IUPAC.protein), 
                         id= str(df.loc[i,'structureId'])+'_'+str(df.loc[i,'entityId']),
                         dbxref=df.loc[i,'pubmedId'])
        sequences.append(record)
    return sequences
  
  def get_pdb_seq_API(self):
    """
    PDB Restful API request to get a customized report
    --------------------------------------------------
    
    Check for query fields: http://www.rcsb.org/pdb/results/reportField.do
    
    Get all the current PDB IDs in JSON format: http://www.rcsb.org/pdb/json/getCurrent
    """
    #import urllib2
    #import datetime
    #seqfile=None
    import requests
    print('Getting sequences from PDB API..')
    
    pdburl="""http://www.rcsb.org/pdb/rest/customReport"""
    query="""?pdbids=*&customReportColumns=structureId,releaseDate,pubmedId,publicationYear,entityMacromoleculeType,sequence&primaryOnly=1&format=csv&service=wsfile"""
    #query="""?pdbids=5JZI,1MCN,1FBI &customReportColumns=structureId,pubmedId,releaseDate,publicationYear,revisionDate,entityMacromoleculeType,sequence&primaryOnly=1&format=csv&service=wsfile"""
    result= requests.get(pdburl, data=query)
    #f = urllib2.urlopen(req)
    #result = f.read().decode('utf-8')
    if result:
      """
      if os.path.exists(self.seqfile) and os.path.getsize(self.seqfile) >0:
        raise Exception('The {} file already exists. Please delete or rename this file and run the script again.'.format(self.seqfile))
      else:
        with open(self.seqfile, 'wt') as out:
          print(result, file=out)
      """
      with open(self.seqfile, 'wt') as out:
        print(result.text, file=out)
      res_df = pd.read_csv(self.seqfile)
      return (res_df)
    else:
      raise Exception("Failed to retrieve results from PDB API.")
    #return result
  
  def get_revised_PDBs(self):
    """
    Returns a list of recently revised PDB structures. 
    """
    
    url = 'http://www.rcsb.org/pdb/rest/search'
    queryText = """
    <orgPdbQuery>
    <queryType>org.pdb.query.simple.ModifiedStructuresQuery</queryType>
    </orgPdbQuery>
    """
    #print("query:", queryText)
    print("Getting revised PDB structures...")
    req = urllib2.Request(url, data=queryText)
    f = urllib2.urlopen(req)
    result = f.read()
    revised_pdbs=result.split('\n')[:-1]
    return revised_pdbs

  def get_pdb_seq_ftp(self):
    """
    Get all the PDB chain sequences using FTP.
    Return a SeqIO object.
    """
    print('Getting sequences from PDB FTP..')
    pdburl = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
    outpdb = "../data/pdb_seqres.txt.gz"
    self.geturl(pdburl, outpdb)
    with gzip.open(outpdb, 'rb') as pdbfile:
      pdbseq=SeqIO.read(pdbfile, 'fasta')
    return pdbseq
  
  def assign_Gdomain(self, seq, seq_id=None):
    """
    Returns G domain of a MHC sequence.
    """
    from mhc_G_domain import mhc_G_domain
    gd=mhc_G_domain(chain1_seq=seq, ch1_id=seq_id)
    res= gd.get_g_domain()
    if res:
      return res
    else:
      #print('ERROR: G domain cannot be assigned.')
      return None,None
  
  def get_MRO(self):
    """
    Clone or pull MRO GitHub repository.
    """
    mro_path=os.path.abspath('../data/MRO')
    wd=os.path.dirname(os.path.realpath(__file__))
    if os.path.exists(mro_path):
      print('Updating MRO repository..')
      os.chdir(mro_path)
      self.run_cmd('git pull')
      os.chdir(wd)
      return
    else:
      print('Getting MRO repository..')
      os.chdir(os.path.abspath('../data'))
      self.run_cmd('git clone https://github.com/IEDB/MRO.git')
      os.chdir(wd)
      return
  
  def get_MRO_Gdomains(self, mro_TSVfile):
    """
    Returns G doamins of the MRO chain sequences.
    """
    print('Assigning G domains to the MRO chain sequences..')
    mro = pd.read_csv(mro_TSVfile, sep='\t', skiprows=[1])
    #mro_header= list(mro.columns)
    mro_out= pd.DataFrame(columns= ['Label','Sequence', 'calc_mhc_class','ch_g_dom'])
    cnt=0
    for i in list(range(mro.shape[0])):
      if pd.isnull(mro.loc[i, 'Sequence']):
        continue
      mro_out.loc[cnt, 'Label']= mro.loc[i,'Label']
      mro_out.loc[cnt, 'Sequence']= mro.loc[i,'Sequence']
      mro_out.loc[cnt, 'calc_mhc_class'], mro_out.loc[cnt, 'ch_g_dom'] = self.assign_Gdomain(mro.loc[i,'Sequence'], mro.loc[i,'Accession'])
      cnt+=1
    #mro_out.to_csv(mro_gdom_file, index=False)
    return mro_out
  
  def get_MRO_allele(self, mro_df, seq, seq_id=None):
    mhc_class, pdb_g_dom=self.assign_Gdomain(seq, seq_id)
    if pdb_g_dom:
      mro_allele= str(list(mro_df[mro_df.fillna({'ch_g_dom' :''}).apply(lambda r : r['ch_g_dom']!='' and (pdb_g_dom in r['ch_g_dom'] or r['ch_g_dom'] in pdb_g_dom) , axis=1)]['Label'])).strip('[]')
      return mro_allele
    else:
      #print('Unable to assign G domain to the {} chain sequence'.format(seq_id))
      return
  
  def run_sql_query(self, sql_command):
    """
    Get PDB_IDs and their PubMed_IDs from iedb_curation SQL database.
    
    SQLconfig.py file should contain the following dict:
      sql={'host':'', 'port':3306, 'user':'', 'passwd':''}
    """
    #import MySQLdb as mdb
    import pymysql as mdb
    import SQLconfig as config
    
    try:                 
      #connection = mdb.connect(host='', port=306, user='', passwd='', db='iedb_production', charset='utf8', use_unicode=False)
      connection = mdb.connect(host= config.sql['host'],
                               port= config.sql['port'], 
                               user= config.sql['user'],
                               passwd=config.sql['passwd'],
                               db='iedb_curation',
                               charset='utf8')
    except (mdb.OperationalError, mdb.ProgrammingError) as msg:
      raise Exception(msg)
    except mdb.InternalError as msg:
      raise Exception(msg)
    #except mdb.Error:
    #  print("ERROR: %d: %s" % (mdb.Error.args[0],mdb.Error.args[1]))
    #  sys.exit(1)
  
    cursor = connection.cursor()
    cursor.execute(sql_command)
    header = [i[0] for i in cursor.description]
    header=list(header)
    result = cursor.fetchall()
    cursor.close()
    connection.close()
    return result, header
  
  def get_IEDB_PDBs(self):
    # BCR
    print('Running BCR SQL query..')
    bcell_query="""select b.bcell_id as assay_id, c.pdb_id, (select pubmed_id from article where reference_id=b.reference_id) as pubmed_id from bcell b, complex c where b.complex_id=c.complex_id and c.e_viewer_status='Y';
    """
    bcell_results, header=self.run_sql_query(bcell_query)
    # TCR
    print('Running TCR SQL query..')
    tcell_query="""select t.tcell_id as assay_id, c.pdb_id, (select pubmed_id from article where reference_id=t.reference_id) as pubmed_id from tcell t, complex c where t.complex_id=c.complex_id and c.e_viewer_status='Y';
    """
    tcell_results, header=self.run_sql_query(tcell_query)
    # MHC
    print('Running MHC SQL query..')
    mhc_query="""select m.mhc_bind_id as assay_id, c.pdb_id, (select pubmed_id from article where reference_id=m.reference_id) as pubmed_id from mhc_bind m, complex c where m.complex_id=c.complex_id and c.e_viewer_status='Y';
    """
    mhc_results, header=self.run_sql_query(mhc_query)
    print('Writting SQL query results to {}'.format(self.sql_results_file))
    sql_res=[]
    sql_res.extend(list(bcell_results))
    sql_res.extend(list(tcell_results))
    sql_res.extend(list(mhc_results))
    inp = pd.DataFrame(sql_res, columns=header)
    inp.to_csv(self.sql_results_file, index=False, encoding='utf-8')
    return inp
 
  def get_PDBs_classication(self, iedb_pdbs, revised_pdbs, all_pdbs):
    """
    Get a list of new or revised PDB IDs for classification. 
    All inputs are sets of PDB IDs.
    """
    new_pdbs=[]
    all_nonIEDB_pdbs= all_pdbs-iedb_pdbs
    revised_iedb_pdbs= iedb_pdbs & revised_pdbs
    # Check for PDBs in the revised_iedb_pdbs with major revisions (version: 2, 3, ...)
    # Add these PDBs with major revisions to the all_nonIEDB_pdbs for further classification
    
    
    return new_pdbs
    
    
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
    self.run_cmd(cmd)
    
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
  
  def assign_class(self, seq):
    """
    Returns BCR, TCR or MHC class and a chain type for an input sequence.
    """
    if self.check_seq(seq)==1:
      receptor, chain_type,c_type=self.get_chain_type(seq)
      #print(chain_type, str(c_type))
      if chain_type:
        return (receptor, str(chain_type))
      else:
        chain_type=None
        mhc_I_score=None
        mhc_I_score=self.is_MHC(str(seq.seq), self.mhc_I_hmm)
        if mhc_I_score >= self.hmm_score_threshold:
          return('MHC-I', 'alpha')
        else:
          mhc_II_alpha_score=None
          mhc_II_alpha_score=self.is_MHC(str(seq.seq), self.mhc_II_alpha_hmm)
          if mhc_II_alpha_score and mhc_II_alpha_score >= self.hmm_score_threshold:
            return('MHC-II', 'alpha')
          else:
            mhc_II_beta_score=None
            mhc_II_beta_score=self.is_MHC(str(seq.seq), self.mhc_II_beta_hmm)
            if mhc_II_beta_score and mhc_II_beta_score >= self.hmm_score_threshold:
              return('MHC-II', 'beta')
            else:
              return(None,None)
  
  def classify(self, sequences=None, mro_df=None):
    """
    Returns a csv file with BCR, TCR or MHC class and chain type assignments to the 
    sequences in the provided input fasta sequence file.
    """
    if not sequences:
      if not self.seqfile:
        raise Exception('Please provide a fasta formatted protein sequence file.')
      else:
        sequences = SeqIO.parse(self.seqfile, 'fasta')
    out=pd.DataFrame()
    cnt=0
    
    for seq in sequences:
      #print(seq.description)
      #print(seq.seq)
      chain_type=None
      receptor=None
      out.loc[cnt,'ID']=str(seq.description)
      receptor, chain_type=self.assign_class(seq)
      out.loc[cnt,'class']=receptor
      out.loc[cnt,'chain_type']=str(chain_type)
      if mro_df and receptor in ('MHC-I', 'MHC-II'):
        out.loc[cnt,'calc_mhc_allele']= self.get_MRO_allele(mro_df, str(seq.seq), str(seq.description))
      #out.loc[cnt,'ID']='"'+str(seq.description)+'"'
      cnt+=1
    out.to_csv(self.outfile, index=False)
      
    
  def classify_pdb_chains_FTP(self):
    pdbseq=self.get_pdb_seq_ftp()
    iedb_PDBs=self.get_IEDB_PDBs()
    
    self.get_MRO()
    mro_out=self.get_MRO_Gdomains(self.mro_file)
    self.classify(pdbseq, mro_out)
  
  def classify_pdb_chains_API(self):
    api_res=self.get_pdb_seq_API()
    pdbseq=self.create_SeqRecord(self.seqfile)
    iedb_PDBs=self.get_IEDB_PDBs()
    
    self.get_MRO()
    mro_out=self.get_MRO_Gdomains(self.mro_file)
    self.classify(pdbseq, mro_out)