
# ClassifierTool: Classifies protein sequences into BCR/TCR/MHC chains.
### @author: Swapnil Mahajan

## Requirements:
- Linux OS
- HMMER3: http://hmmer.org/
- ANARCI v1.1 : http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php#download
	Dunbar J and Deane CM. ANARCI: Antigen receptor numbering and receptor classification. Bioinformatics (2016)
- Python 2.7
	Python packages: Pandas, BioPython

## How to use:

```shell
cd dir_path/ClassifierTool
cd code
python run_SeqClassifier.py ../data/all_mhcBcrTcr_IEDB.fasta ../out/all_mhcBcrTcr_IEDB.csv
```

- ../data/all_mhcBcrTcr_IEDB.fasta file has one or more protien sequences in fasta format.
- ../out/all_mhcBcrTcr_IEDB.csv file is a output file name.

## How it works:
- BCR and TCR chains are identified using ANARCI. ANARCI searches a given protein sequence against HMMs built using BCR and TCR chain sequences from IMGT. HMMER is used to align an input sequence to the HMMs.
- MHC class I (alpha1-alpha2 domains) and MHC class I alpha and beta chain HMMs were downloaded from Pfam website. An input protein sequence is searched against these HMMs. A HMMER bit score threshold of 25 was used to identify MHC chain sequences. DTU uses 250 as a score cutoff.