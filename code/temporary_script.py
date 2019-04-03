from Bio import SeqIO
from ConstantClassifier import SeqClassifier

records = list(SeqIO.parse("../data/all_mhcBcrTcr_IEDB.fasta", "fasta"))
seq = records[2812]

sc = SeqClassifier()
sc.assign_class(seq)
