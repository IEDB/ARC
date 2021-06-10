"""
Obtains the G-domain (peptide binding groove domain)
of an MHC sequence. Contains methods related to this purpose.

Authors: Swapnil Mahajan
"""
from __future__ import print_function, division

import copy
import pandas as pd
import random
import re
import subprocess
import sys
import os
import tempfile

class mhc_G_domain:
    """Assigns mhci_G, mhcii_alpha or mhcii_beta domains to a sequence

    Requires BLAST as an external dependency.
    
    Atrributes:
        chain1_seq: full MHC chain sequence
        chain2_seq: full MHC chain sequence
        ch1_id: ID of chain1_seq
        ch2_id: ID of chain2_seq
    """
    def __init__(self, chain1_seq, chain2_seq=None, ch1_id=None, ch2_id=None):
        """Inits the mhc_G_domain class with required attributes"""
        self.chain1_seq = chain1_seq
        self.chain2_seq = chain2_seq
        self.ch1_id = ch1_id
        self.ch2_id = ch2_id

        self.e_val = pow(10, -10)
        self.max_out_seq = 10
        self.hit_coverage = 0.85
        self.rand = str(random.randint(0, 1000000))

    def run_blast(self, seqfile, db, blastout, qcov=None, e_val=None):
        """Runs Blast for the give sequence file against the given database
        
        Full command formats:
            blastp -query ../test_seq.fas -db G_ALPHA1.fasta /
                -outfmt '10 std slen' -max_target_seqs 10 /
                -evalue 1e-10 > ./test_alpha1.out

            blastp -query ../test_seq.fas -db G_ALPHA1.fasta /
                -outfmt '10 std slen' -max_target_seqs 10 /
                -evalue 1e-10 -qcov_hsp_perc [90] > ./test_alpha1.out

        blastout columns:
        'qseqid sseqid pident length mismatch gapopen  /
            qstart qend sstart send evalue bitscore slen'
        
        Args:
            seqfile: FASTA file containing the sequence to process
            db: database to use with BLAST command
            blastout: output file for the BLAST command
            qcov: the query coverage threshold (should be a float i.e. 75% = .75)
            e_val: e-value threshold (should be in the format 10e-10)
                See e_val default in initialization for example

        Returns:
            The output BLAST file
        """
        if not e_val:
            e_val = self.e_val
        if not qcov:
            args = [
                'blastp', '-query', seqfile, '-db', db, '-outfmt',
                '"10 std slen"', '-max_target_seqs',
                str(self.max_out_seq), '-evalue',
                str(e_val), '>', blastout
            ]
        else:
            args = [
                'blastp', '-query', seqfile, '-db', db, '-outfmt',
                '"10 std slen"', '-max_target_seqs',
                str(self.max_out_seq), '-evalue',
                str(e_val), '-qcov_hsp_perc', qcov, '>', blastout
            ]
        cmd = ' '.join(args)
        proc = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                env=dict(os.environ, my_env_prop='value'),
                                shell=True)
        proc.communicate()

        return blastout

    def blast_all(self, seqfile):
        """Runs BLAST for a given seq file against all G domains 
            (alpha1, alpha2, alpha, beta).

        Args:
            seqfile: the sequence file in FASTA format
        
        Returns:
            BLAST results from alpha1, alpha2, alpha, beta chain BLAST dbs
        """
        g_alpha1_db = 'data/blastdb/G_ALPHA1.fasta'
        g_alpha1_db = os.path.join(os.path.dirname(__file__), g_alpha1_db)
        g_alpha2_db = 'data/blastdb/G_ALPHA2.fasta'
        g_alpha2_db = os.path.join(os.path.dirname(__file__), g_alpha2_db)
        g_alpha_db = 'data/blastdb/G_ALPHA.fasta'
        g_alpha_db = os.path.join(os.path.dirname(__file__), g_alpha_db)
        g_beta_db = 'data/blastdb/G_BETA.fasta'
        g_beta_db = os.path.join(os.path.dirname(__file__), g_beta_db)

        alpha1_file = tempfile.NamedTemporaryFile(delete=False)
        alpha2_file = tempfile.NamedTemporaryFile(delete=False)
        alpha_file = tempfile.NamedTemporaryFile(delete=False)
        beta_file = tempfile.NamedTemporaryFile(delete=False)
        alpha1_blout = self.run_blast(seqfile, g_alpha1_db, alpha1_file.name)
        alpha2_blout = self.run_blast(seqfile, g_alpha2_db, alpha2_file.name)
        alpha_blout = self.run_blast(seqfile, g_alpha_db, alpha_file.name)
        beta_blout = self.run_blast(seqfile, g_beta_db, beta_file.name)

        return alpha1_blout, alpha2_blout, alpha_blout, beta_blout

    def get_domain(self, blastout):
        """Identify the domain boundaries from blast results.

        Args:
            blastout: The BLAST output file to parse

        Returns:
            The protein sequence of the matched domain, None if no domain identified
        """
        columns = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen'.split(
            ' ')

        df = pd.read_csv(blastout)
        df.columns = columns
        for i in range(df.shape[0]):
            hit_cover = (df.loc[i, 'send'] - df.loc[i, 'sstart'] +
                         1) / df.loc[i, 'slen']
            if hit_cover >= self.hit_coverage:
                return df.loc[i, 'qstart'], df.loc[i, 'qend']
        return None

    def get_all_domains(self, alpha1_blout, alpha2_blout, alpha_blout,
                        beta_blout):
        """Retrieves all mhc_class and domain boundaries.

        Args:
            alpha1_blout: BLAST output from alpha1 query
            alpha2_blout: BLAST output from alpha2 query
            alpha_blout: BLAST output from full alpha query
            beta_blout: BLAST output from full beta query
        
        Returns:
            Relevant domain boundaries from the BLAST queries
        """
        if os.path.getsize(alpha1_blout) > 0:
            dom1_st_end = self.get_domain(alpha1_blout)
            if not dom1_st_end:
                return None
            if os.path.getsize(alpha2_blout) > 0:
                dom2_st_end = self.get_domain(alpha2_blout)
                if dom2_st_end:
                    return 'I', [
                        dom1_st_end[0], dom1_st_end[1], dom2_st_end[0],
                        dom2_st_end[1]
                    ]
                else:
                    return None

        if os.path.getsize(alpha_blout) > 0:
            dom1_st_end = self.get_domain(alpha_blout)
            if not dom1_st_end:
                return None
            return 'IIa', dom1_st_end
        if os.path.getsize(beta_blout) > 0:
            dom1_st_end = self.get_domain(beta_blout)
            if not dom1_st_end:
                return None
            return 'IIb', dom1_st_end

        return None

    def check_b2m(self, seqfile):
        """Check if the seq is beta-2-microglobulin.
        
        Args:
            seqfile: sequence file in FASTA format to query
        
        Returns:
            0 if sequence is b2m, 1 if sequence is NOT b2m
        """
        #g_dom_db= '../data/blastdb/IMGT_G_domain_Species.fasta'
        #g_dom_db=os.path.join(os.path.dirname(__file__),  g_dom_db)
        b2m_db = 'data/blastdb/b2m.fasta'
        b2m_db = os.path.join(os.path.dirname(__file__), b2m_db)
        blout_file = tempfile.NamedTemporaryFile(delete=False)
        #blout= self.run_blast(seqfile, g_dom_db, 'b2m_'+self.rand+'.out')
        blout = self.run_blast(seqfile, b2m_db, blout_file.name,
                               '90', pow(10, -50))
        if os.path.getsize(blout) > 0:
            os.unlink(blout)
            # return 1
            return 0
        else:
            os.unlink(blout)
            return 1

    def get_subseq(self, seq, dom_st, dom_end):
        """Returns a subsequence of a protein given the start and end positions

        Args:
            seq: sequence to parse
            dom_st: start position of the domain
            dom_end: end position of the domain
        """
        subseq = ''
        subseq = seq[dom_st - 1:dom_end]
        return (subseq)

    def chain_g_domain(self, seq, seq_id=None):
        """Returns G-domain sequence of the given sequence.
        
        Args:
            seq: protein sequence to retrieve G domain for
            seq_id: the sequence identifier (>"identifier" in FASTA)
        """
        inpfile = tempfile.NamedTemporaryFile(delete=False) 
        if not seq_id:
            seq_id = 'tmp_seq' + self.rand

        g_domain = ''
        #inpfile = os.path.join(os.path.dirname(__file__),
        #                       str(seq_id) + '.fasta')
        inp = open(inpfile.name, 'wt')
        print('>' + str(seq_id), file=inp)
        print(seq, file=inp)
        inp.close()

        # check if seq is b2m or mhc chain
        b2m = self.check_b2m(inpfile.name)
        if b2m == 0:
            return None

        alpha1_blout, alpha2_blout, alpha_blout, beta_blout = self.blast_all(
            inpfile.name)
        res = self.get_all_domains(alpha1_blout, alpha2_blout, alpha_blout,
                                   beta_blout)
        if res:
            mhc_class, st_end = res

            if mhc_class and st_end:
                if mhc_class == 'I':
                    if st_end[2] - st_end[1] != 1:
                        st_end[2] = st_end[1] + 1
                    dom1 = self.get_subseq(seq, st_end[0], st_end[1])
                    dom2 = self.get_subseq(seq, st_end[2], st_end[3])
                    g_domain = dom1 + dom2
                else:
                    g_domain = self.get_subseq(seq, st_end[0], st_end[1])
                os.unlink(inpfile.name)
                os.unlink(alpha1_blout)
                os.unlink(alpha2_blout)
                os.unlink(alpha_blout)
                os.unlink(beta_blout)
                return mhc_class, g_domain
        else:
            os.unlink(inpfile.name)
            os.unlink(alpha1_blout)
            os.unlink(alpha2_blout)
            os.unlink(alpha_blout)
            os.unlink(beta_blout)
            return None

    def get_g_domain(self):
        """Retrieves the G domain of a given sequence

        Relies on the attributes defined upon instantiation of mhc_G_domain class
        (Ex: self.chain1_seq, self.ch1_id, self.chain2_seq, self.ch2_id)

        If these are not initialized, the script will return "None"
        """
        mhc_class1 = None
        dom1_seq = None
        mhc_class2 = None
        dom2_seq = None
        g_dom1 = None
        g_dom2 = None

        if self.chain1_seq:
            res1 = self.chain_g_domain(self.chain1_seq, self.ch1_id)
            if res1:
                mhc_class1, dom1_seq = res1
                if mhc_class1 == 'I':
                    return mhc_class1, dom1_seq
                if mhc_class1 == 'IIa':
                    g_dom1 = dom1_seq
                    if not self.chain2_seq:
                        return mhc_class1, dom1_seq
                if mhc_class1 == 'IIb':
                    g_dom2 = dom1_seq
                    if not self.chain2_seq:
                        return mhc_class1, dom1_seq
        else:
            return None

        if self.chain2_seq:
            res2 = self.chain_g_domain(self.chain2_seq, self.ch2_id)
            if res2:
                mhc_class2, dom2_seq = res2
                if mhc_class2 == 'I':
                    return mhc_class2, dom2_seq
                if mhc_class2 == 'IIa':
                    g_dom1 = dom2_seq
                    if not self.chain1_seq:
                        return mhc_class2, dom2_seq
                if mhc_class2 == 'IIb':
                    g_dom2 = dom2_seq
                    if not self.chain1_seq:
                        return mhc_class2, dom2_seq

        if self.chain1_seq and self.chain2_seq:
            return 'II', g_dom1 + g_dom2
