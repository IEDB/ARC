# ARC (Antigen Receptor Classifier)
### @authors: Swapnil Mahajan, Austin Crinklaw

## Requirements:
- Linux OS
- [HMMER3](http://hmmer.org/)
- NCBI Blast+
- Python 3+
  - Python packages: Pandas, BioPython
- Git

## How to use:

### Installation:
The easiest way to use the software is to download and utilize the Dockerfile.

ARC can also be downloaded through PyPI using the following pip command.
```shell
pip install bio-arc
```

### Input  
-  A fasta format file with one or more protein sequences.  
  ```
  >1WBZ_A_alpha I H2-Kb
MVPCTLLLLLAAALAPTQTRAGPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTDSPKAHVTHHSRPEDKVTLRCWALGFYPADITLTWQLNGEELIQDMELVETRPAGDGTFQKWASVVVPLGKEQYYTCHVYHQGLPEPLTLRWEPPPSTVSNMATVAVLVVLGAAIVTGAVVAFVMKMRRRNTGGKGGDYALAPGSQTSDLSLPDCKVMVHDPHSLA
>1WBZ_B_b2m I H2-Kb
MARSVTLVFLVLVSLTGLYAIQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHASMAEPKTVYWDRDM
  ```

  

### Commands
- A list of commands can be found via the -h flag
```shell
python -m ARC -h
```
- Specific commands are explained in a similar manner
```shell
python -m ARC <command> -h
```
- HMMs come by default but can be updated. If they are missing for some they can be added via the install command
```shell
python -m ARC update -archive
```
Note: Archive is an optional parameter that will create an archive folder with HMMs stored by date. You will find this in the data folder wherever you have installed ARC.

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
   -  e.g. alpha, beta, heavy, lambda, kappa, scFv, TscFv or construct.
-  The fourth column named 'calc_mhc_allele' is the MHC allele identified using groove domain similarity to MRO alleles.

| ID	                                  | class  | chain_type | calc_mhc_allele|
|---------------------------------------- |------- |----------- |---------------|
| 1WBY_A_alpha I H2-Db                    |	MHC-I  | alpha      | |
| 1WBY_B_b2m I H2-Db	                  |	       |            | |
| 1HQR_A_alpha II HLA-DRA*01:01/DRB5*01:01|	MHC-II | alpha      | HLA-DRA*01:01 |
| 1HQR_B_beta II HLA-DRA*01:01/DRB5*01:01 |	MHC-II | beta       | HLA-DRB5*01:01 |
| 2CMR_H_heavy                            |	BCR	   | heavy      | |
| 2CMR_L_light                            |	BCR	   | kappa      | |
| 4RFO_L_light                            |	BCR	   | lambda     | |
| 3UZE_A_heavy                            |	BCR	   | scFv       | |
| 1FYT_D_alpha                            |	TCR	   | alpha      | |
| 1FYT_E_beta                             | TCR	   | beta       | |
| 3TF7_C_alpha                            |	TCR    | TscFv      | |

## How it works:
- BCR and TCR chains are identified using HMMs. A given protein sequence is searched against HMMs built using BCR and TCR chain sequences from IMGT. HMMER is used to align an input sequence to the HMMs.
- MHC class I (alpha1-alpha2 domains) and MHC class I alpha and beta chain HMMs are downloaded from Pfam website. An input protein sequence is searched against these HMMs. A HMMER bit score threshold of 250 was used to identify MHC chain sequences. DTU uses 250 as a score cutoff which can exclude MHC like molecules such as Human and Mouse CD1d molecules.
-To identify MHC alleles, MRO repository is downloaded every time the script is run. Groove domains (G-domains) are assigned to new MRO allles and stored in a CSV file. If this file does not exist then G-domains are assigned to all the MRO alleles (which may slow down the script).

## References:
Several pieces of code, including the IMGT ripping / HMM generation, was sourced from ANARCI.

[Dunbar J and Deane CM. ANARCI: Antigen receptor numbering and receptor classification. Bioinformatics (2016)](https://academic.oup.com/bioinformatics/article/32/2/298/1743894)
