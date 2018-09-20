import subprocess
import argparse
import tempfile 
from Bio import SeqIO

"""
How to run:
# python DTU_find_MHC_dataset.py ../data/mro_chain_seq.fasta ../out/DTU_MHC-I_out.csv ../data/DTU_MHC_I_complete.hmm
"""

score_threshold=25

##################################
##          Functions           ##
##################################

def run_cmd(cmd, input_string=''):
  """Run the cmd with input_string as stdin and return output."""
  
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True)
  out, stderr = p.communicate(input=input_string)
  if p.returncode:
    raise Exception('Cmd {} failed: {}'.format(cmd[0], stderr))
  
  return out

def is_MHC(seq,hmm):
  # create a temporary file

  fp = tempfile.NamedTemporaryFile()
  fp.write('>seq\n'.format(seq))
  fp.write('{}'.format(seq))
  fp.flush() 

  #Find MHC sequences
  cmd = ['hmmsearch', hmm, fp.name ]
  
  output = run_cmd(cmd)
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

##################################
##            MAIN              ##
##################################

parser = argparse.ArgumentParser()
parser.add_argument('infile', help = '')
parser.add_argument('outfile', help = '')
parser.add_argument('hmm', help = '')

args = parser.parse_args()

outfile = open(args.outfile,'w')
outfile.write('PDBID_ChainID,HMMscore\n')

sequences = SeqIO.parse(args.infile, 'fasta')

for record in sequences:
  # Return None if the HMM does not align to the sequence 
  score = is_MHC(record.seq,args.hmm) 
  if score != None and score >= score_threshold:
    # Only write out positive HMM scores 
    outfile.write('{},{}\n'.format(record.description, score))

outfile.close()








