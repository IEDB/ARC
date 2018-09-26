# ClassifierTool: Classifies protein sequences into BCR/TCR/MHC chains.
### @author: @Swapnil

## Requirements:
- Linux OS
- HMMER3: http://hmmer.org/
- ANARCI v1.1 : http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php#download  
  Dunbar J and Deane CM. ANARCI: Antigen receptor numbering and receptor classification. Bioinformatics (2016).  

- Python 2.7
	Python packages: Pandas, BioPython

## How to use:
### Input  
  A fasta format file with one or more protein sequences.  
  e.g. ../data/all_mhcBcrTcr_IEDB.fasta file has one or more protien sequences in the fasta format.  
  ```
  >1WBZ_A_alpha I H2-Kb
MVPCTLLLLLAAALAPTQTRAGPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTDSPKAHVTHHSRPEDKVTLRCWALGFYPADITLTWQLNGEELIQDMELVETRPAGDGTFQKWASVVVPLGKEQYYTCHVYHQGLPEPLTLRWEPPPSTVSNMATVAVLVVLGAAIVTGAVVAFVMKMRRRNTGGKGGDYALAPGSQTSDLSLPDCKVMVHDPHSLA
>1WBZ_B_b2m I H2-Kb
MARSVTLVFLVLVSLTGLYAIQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHASMAEPKTVYWDRDM
  ```
  

### Command  
```shell
cd dir_path/ClassifierTool
cd code
python run_SeqClassifier.py ../data/all_mhcBcrTcr_IEDB.fasta ../out/all_mhcBcrTcr_IEDB.csv
```

### Output  
  CSV file with 3 columns.  
  e.g. ../out/all_mhcBcrTcr_IEDB.csv file is a output file name.

|ID	| class	| chain_type |
|--- |--- |--- |
|1WBY_A_alpha I H2-Db |	MHC-I|	alpha|
|1WBY_B_b2m I H2-Db	|	|
|1HQR_A_alpha II HLA-DRA*01:01/DRB5*01:01|	MHC-II|	alpha|
|1HQR_B_beta II HLA-DRA*01:01/DRB5*01:01|	MHC-II|	beta|
|2CMR_H_heavy|	BCR	|heavy|
|2CMR_L_light|	BCR	|kappa|
|4RFO_L_light|	BCR	|lambda|
|3UZE_A_heavy|	BCR	|scFv|
|1FYT_D_alpha|	TCR	|alpha|
|1FYT_E_beta	|TCR	|beta|
|3TF7_C_alpha|	TCR|	TscFv|

## How it works:
- BCR and TCR chains are identified using ANARCI. ANARCI searches a given protein sequence against HMMs built using BCR and TCR chain sequences from IMGT. HMMER is used to align an input sequence to the HMMs.
- MHC class I (alpha1-alpha2 domains) and MHC class I alpha and beta chain HMMs were downloaded from Pfam website. An input protein sequence is searched against these HMMs. A HMMER bit score threshold of 25 was used to identify MHC chain sequences. DTU uses 250 as a score cutoff.