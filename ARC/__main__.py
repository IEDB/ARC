import argparse
import sys

#from ARC.build_pipeline import build //Currently deprecated
from ARC.classifier import SeqClassifier

tasks = ["classify"]
if len(sys.argv) < 2 or (sys.argv[1] not in tasks):
    usage = """
    ___    ____  ______
   /   |  / __ \/ ____/
  / /| | / /_/ / /     
 / ___ |/ _, _/ /___   
/_/  |_/_/ |_|\____/   (Antigen Receptor Classifier)
                       
	Usage: python arc <task> [options]

	Available tasks are:
		{:s}

	For help:
	> python ARC <task> -h
	""".format(*tasks)
    print(usage)
    sys.exit()

elif sys.argv[1] == "classify":
    prsr = argparse.ArgumentParser(
        prog='classify', description='Classify protein sequences using HMMs')
    prsr.add_argument(
        '-p',
        help="the number of threads to use (default=1)",
        type=int,
        metavar='num_threads',
        default=1
    )
    prsr.add_argument(
        '-i',
        help="input file containing protein sequences in FASTA sequence format",
        type=str,
        metavar='infile_name',
        required=True)
    prsr.add_argument(
        '-o',
        help="output file name and location, ex: data/my_outfile.csv",
        type=str,
        metavar='outfile_name',
        required=True)
    prsr.add_argument(
        '-hmmer',
        help="Path to your local HMMer installation",
        type=str,
        metavar='hmmer_path',
        required=False)
    prsr.add_argument(
        '-blast',
        help="Path to your local blast installation",
        type=str,
        metavar='blast_path',
        required=False)
    args = prsr.parse_args(sys.argv[2:])
    classifier = SeqClassifier(args.i, args.o, args.p, args.hmmer, args.blast)
    classifier.classify_seqfile(args.i)
