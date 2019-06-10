# Run the pipeline to get the sequences, format them and build the databases

DIR=`python -c "import ARC as _; print(_.__path__[0])"`

rm -rf $DIR/data/HMMs

# Rip the sequences from the imgt website. HTML may change in the future. 
mkdir -p $DIR/build_pipeline/IMGT_sequence_files/htmlfiles
mkdir -p $DIR/build_pipeline/IMGT_sequence_files/fastafiles
python3 $DIR/build_pipeline/RipIMGT.py

# Format the alignments and handle imgt oddities to put into a consistent alignment format
mkdir -p $DIR/build_pipeline/curated_alignments
mkdir -p $DIR/build_pipeline/muscle_alignments
python3 $DIR/build_pipeline/FormatAlignments.py

# Build the hmms for each species and chain.
# --hand option required otherwise it will delete columns that are mainly gaps. We want 128 columns otherwise ARNACI will fall over.
mkdir -p $DIR/build_pipeline/HMMs
hmmbuild --hand $DIR/build_pipeline/HMMs/ALL.hmm $DIR/build_pipeline/curated_alignments/ALL.stockholm
hmmbuild --hand $DIR/build_pipeline/HMMs/ALL_AND_C.hmm $DIR/build_pipeline/curated_alignments/ALL_AND_C.stockholm

# Turn the output HMMs file into a binary form. This is required for hmmscan that is used in ARNACI.
hmmpress -f $DIR/build_pipeline/HMMs/ALL.hmm 
hmmpress -f $DIR/build_pipeline/HMMs/ALL_AND_C.hmm

mv $DIR/build_pipeline/HMMs $DIR/data/
rm -rf $DIR/build_pipeline/curated_alignments
rm -rf $DIR/build_pipeline/IMGT_sequence_files
rm -rf $DIR/build_pipeline/muscle_alignments
