"""
Created on Mon July 15 2019

@author: Austin Crinklaw

Utilizes BLAST to identify IgNAR protein sequences.
"""

from __future__ import print_function, division
import sys
import os
import subprocess
import re
import copy
import pandas as pd
import random


class IgNAR:
    def __init__(self, seq):
        """
        Input: Full MHC chain proteins sequence
        Output: mhci_G, mhcii_alpha or mhcii_beta domain sequence.

        External softwares needed: BLAST
        """

        self.e_val = pow(10, -4)
        self.max_out_seq = 10
        self.hit_coverage = 0.75
        self.hit_perc_id = 0.50
        self.rand = str(random.randint(0, 1000000))

    def run_blast(self, seqfile, db, blastout, qcov=None, e_val=None):
        """
        Runs Blast for the give sequence file against the given database and return a output file in CSV format.
        Command:
        blastp -query ../test_seq.fas -db G_ALPHA1.fasta -outfmt '10 std slen' -max_target_seqs 10 -evalue 1e-10 > ./test_alpha1.out

        blastp -query ../test_seq.fas -db G_ALPHA1.fasta -outfmt '10 std slen' -max_target_seqs 10 -evalue 1e-10 -qcov_hsp_perc [90] > ./test_alpha1.out

        blastout columns:
        'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen'
        """
        if not e_val:
            e_val = self.e_val
        if not qcov:
            args = ['blastp', '-query', seqfile, '-db', db, '-outfmt', '"10 std slen"',
                    '-max_target_seqs', str(self.max_out_seq), '-evalue', str(e_val), '>', blastout]
        else:
            args = ['blastp', '-query', seqfile, '-db', db, '-outfmt', '"10 std slen"',
                    '-max_target_seqs', str(self.max_out_seq), '-evalue', str(e_val), '-qcov_hsp_perc', qcov, '>', blastout]
        cmd = ' '.join(args)
        # print(cmd)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                env=dict(os.environ, my_env_prop='value'), shell=True)
        out, err = proc.communicate()
        return blastout

    def blast_all(self, seqfile):
        """
        Runs blast for a given seq file against all G domains (alpha1, alpha2, alpha, beta).
        """
        ignar_db = '../data/blastdb/G_ALPHA1.fasta'
        ignar_db = os.path.join(os.path.dirname(__file__),  ignar_db)

        seq_id = str(seqfile.split('.')[0])
        ignar_blout = self.run_blast(
            seqfile, ignar_db, seq_id+'_ignar'+self.rand+'.out')

        return ignar_blout

