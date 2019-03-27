# ClassifierTool: Classifies protein sequences into BCR/TCR/MHC chains.
### @author: @Swapnil

## Requirements:
- Linux OS
- HMMER3: http://hmmer.org/
- ANARCI v1.3 : http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php#download  
  - Dunbar J and Deane CM. ANARCI: Antigen receptor numbering and receptor classification. Bioinformatics (2016).  
  - Modified setup.py to create ANARCI1.3 executable instead of ANARCI to keep different versions.
  - Update the ANARCI_LOC variable in setup.py for correct installation.
- Python 2.7
  - Python packages: Pandas, BioPython
  - 
- Git

## How to use:
### Input  
-  A fasta format file with one or more protein sequences.  
  ```
  >1WBZ_A_alpha I H2-Kb
MVPCTLLLLLAAALAPTQTRAGPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTDSPKAHVTHHSRPEDKVTLRCWALGFYPADITLTWQLNGEELIQDMELVETRPAGDGTFQKWASVVVPLGKEQYYTCHVYHQGLPEPLTLRWEPPPSTVSNMATVAVLVVLGAAIVTGAVVAFVMKMRRRNTGGKGGDYALAPGSQTSDLSLPDCKVMVHDPHSLA
>1WBZ_B_b2m I H2-Kb
MARSVTLVFLVLVSLTGLYAIQKTPQIQVYSRHPPENGKPNILNCYVTQFHPPHIEIQMLKNGKKIPKVEMSDMSFSKDWSFYILAHTEFTPTETDTYACRVKHASMAEPKTVYWDRDM
  ```

-  e.g. ../data/all_mhcBcrTcr_IEDB.fasta file has multiple protein sequences in the fasta format.  

  

### Command  
-  Using Fasta file as an input:
```shell
cd dir_path/ClassifierTool/code
python run_FastaSeqClassifier.py ../data/all_mhcBcrTcr_IEDB.fasta ../out/all_mhcBcrTcr_IEDB.xlsx
```
-  When a fasta file is used e.g. `../data/all_mhcBcrTcr_IEDB.fasta` then the output excel file name should be provided e.g. `../out/all_mhcBcrTcr_IEDB.xlsx`.
Or
-  Classify the new released PDB sequences (weekly update)
```shell
cd dir_path/ClassifierTool/code
python run_RecentPDBClassifier.py > ../log/recentPDB_log.txt
```
Or
-  Classify all current PDB sequences (Run this command before running regular weekly updates)
```shell
cd dir_path/ClassifierTool/code
python run_PDBClassifier.py > ../log/All_PDB_log.txt
```
-  When PDB sequence classifier scripts are used then the output excel file is written in `ClassifierTool/out` dir with names such as `SeqClassifier_output_07Jan2019.xlsx`. The date in the file name represents the day on which the script was run.
### Output  
-  Output excel file has 4 columns. 
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
- BCR and TCR chains are identified using ANARCI. ANARCI searches a given protein sequence against HMMs built using BCR and TCR chain sequences from IMGT. HMMER is used to align an input sequence to the HMMs.
- MHC class I (alpha1-alpha2 domains) and MHC class I alpha and beta chain HMMs are downloaded from Pfam website. An input protein sequence is searched against these HMMs. A HMMER bit score threshold of 25 was used to identify MHC chain sequences. DTU uses 250 as a score cutoff which can exclude MHC like molecules such as Human and Mouse CD1d molecules.
-To identify MHC alleles, MRO repository is downloaded every time the script is run. Groove domains (G-domains) are assigned to new MRO allles and stored in the file `ClassifierTool/out/MRO_Gdomain.csv`. If this file is not in the out directory then G-domains are assigned to all the MRO alleles (which may slow down the script).
- PDB classifier script tracks the PDB chains which were not immune receptors in the file `ClassifierTool/out/previous_ClassifiedPDBs_woImmuneReceptors.csv`. To speed up, these PDB chains are not classified next time the script is run.
- PDB classifier also tracks the PDBs which had immune receptor sequences but did not have PubMedID in the file `ClassifierTool/out/previous_ClassifiedPDBs_woPubMedIDs.csv`. The script checks if any of these PDBs were published every time it is run.
