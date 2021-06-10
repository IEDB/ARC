
# ARC (Antigen Receptor Classifier)
### Authors: Austin Crinklaw, Swapnil Mahajan

## Requirements:
- Linux OS
- [HMMER3](http://hmmer.org/)
- NCBI Blast+
- Python 3+
  - Python packages: Pandas, BioPython

## Installation:
We provide a Dockerfile for ease of use.

ARC can also be downloaded through PyPI using the following pip command.
```shell
pip install bio-arc
```

### Testing Installation:
A quick check for proper dependencies and successful installation can be performed by navigating to your pip package install directory (which can be located by executing ```pip show bio-arc```) and running the following command:
```shell
python3 -m arc_test
```
Passing all unit-tests means that your system is configured properly and ready to classify some protein sequences.

## Usage:
### Input  
-  A fasta format file with one or more protein sequences.  
  ```
  >1WBZ_A_alpha I H2-Kb
MVPCTLLLLLAAALAPTQTRAGPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTDSPKAHVTHHSRPEDKVTLRCWALGFYPADITLTWQLNGEELIQDMELVETRPAGDGTFQKWASVVVPLGKEQYYTCHVYHQGLPEPLTLRWEPPPSTVSNMATVAVLVVLGAAIVTGAVVAFVMKMRRRNTGGKGGDYALAPGSQTSDLSLPDCKVMVHDPHSLA
>1WBZ_B_b2m I H2-Kb
MARSVTLVFLVLVSLTGLYAIQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHASMAEPKTVYWDRDM
  ```

  

### Commands
-  Using Fasta file as an input:
```shell
python -m ARC classify -i /path/to/input.fasta -o /path/to/output.csv
```
### Output  
-  Output file has 4 columns in CSV format. 
-  First column named 'ID' is the description provoded in the fasta for each sequence.  
-  Second column named 'class' is the assigned molecule class for each sequence.
   -  e.g. MHC-I, MHC-II, BCR or TCR.  
-  The third column named 'chain_type' is the assigned chain type for each sequence.
   -  e.g. alpha, beta, heavy, lambda, kappa, scFv, TscFv or construct. These will also be labelled as V for variable domain or C for constant domain.
-  The fourth column named 'calc_mhc_allele' is the MHC allele identified using groove domain similarity to MRO alleles.

| ID	                                  | class  | chain_type | calc_mhc_allele|
|---------------------------------------- |------- |----------- |---------------|
| 1WBY_A_alpha I H2-Db                    |	MHC-I  | alpha V     | |
| 1WBY_B_b2m I H2-Db	                  |	       |            | |
| 1HQR_A_alpha II HLA-DRA*01:01/DRB5*01:01|	MHC-II | alpha C     | HLA-DRA*01:01 |
| 1HQR_B_beta II HLA-DRA*01:01/DRB5*01:01 |	MHC-II | beta C     | HLA-DRB5*01:01 |
| 2CMR_H_heavy                            |	BCR	   | heavy V      | |
| 2CMR_L_light                            |	BCR	   | kappa C     | |
| 4RFO_L_light                            |	BCR	   | lambda V    | |
| 3UZE_A_heavy                            |	BCR	   | scFv       | |
| 1FYT_D_alpha                            |	TCR	   | alpha V     | |
| 1FYT_E_beta                             | TCR	   | beta C      | |
| 3TF7_C_alpha                            |	TCR    | TscFv      | |

## How it works:
- BCR and TCR chains are identified using HMMs. A given protein sequence is searched against HMMs built using BCR and TCR chain sequences from IMGT. HMMER is used to align an input sequence to the HMMs.
- MHC class I (alpha1-alpha2 domains) and MHC class I alpha and beta chain HMMs are downloaded from Pfam website. An input protein sequence is searched against these HMMs. A HMMER bit score threshold of 25 was used to identify MHC chain sequences.
- To identify MHC alleles, groove domains (G-domains) are assigned based on the MRO repository. 
- IgNAR sequences are identified through querying against a custom blast database.

## References:
Several methods for HMMER result parsing were sourced from ANARCI.

[Dunbar J and Deane CM. ANARCI: Antigen receptor numbering and receptor classification. Bioinformatics (2016)](https://academic.oup.com/bioinformatics/article/32/2/298/1743894)
