from Bio import SeqIO
from ConstantClassifier import SeqClassifier
import time

start = time.time()
records = list(SeqIO.parse("../data/all_mhcBcrTcr_IEDB.fasta", "fasta"))
sc = SeqClassifier()
sc.classify_seqfile("../data/all_mhcBcrTcr_IEDB.fasta")
set_size = str(len(records))
end = time.time()

print("Script processed " + set_size + " sequences in " + str(end - start) + " seconds")
