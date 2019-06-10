import argparse
import sys
from ARC.build_pipeline import build
from ARC.classifier import SeqClassifier

tasks = ["classify", "install", "update"]
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
		{:s}
		{:s}

	For help:
	> python ARC <task> -h (--help)
	""".format(*tasks)
	print(usage)
	sys.exit()

elif sys.argv[1]  == "classify":
	prsr = argparse.ArgumentParser(prog='classify',
		description='Classify protein sequences using HMMs')
	prsr.add_argument('-i', help="input file containing protein sequences in FASTA sequence format",
					type=str, metavar='infile_name', required=True)
	prsr.add_argument('-o', help="output file name and location, ex: data/my_outfile.csv",
					type=str, metavar='outfile_name', required=True)
	args = prsr.parse_args(sys.argv[2:])

	classifier = SeqClassifier(args.i, args.o)
	classifier.classify_seqfile(args.i)

elif sys.argv[1] == "install":
	prsr = argparse.ArgumentParser(prog='install',
		description='Build HMMs used for classification by downloading sequences from IMGT')
	prsr.add_argument('-quiet', help="Suppress installation output",
					action='store_true', default=False)
	args = prsr.parse_args(sys.argv[2:])

	build.install(args.quiet)

elif sys.argv[1] == "update":
	prsr = argparse.ArgumentParser(prog='update',
		description='Update HMMs used for classification')
	prsr.add_argument('-archive', help='Choose to archive old HMMs',
					action='store_true', default=False)
	args = prsr.parse_args(sys.argv[2:])

	build.update(args.archive)