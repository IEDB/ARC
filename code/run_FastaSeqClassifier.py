# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 09:44:50 2019

@author: swapnil
"""

from SeqClassifier import SeqClassifier
from datetime import datetime
#import pandas as pd
import sys

startTime = datetime.now()
classification=SeqClassifier(sys.argv[1], sys.argv[2])

# =============================================================================
# To classify all the protein sequences in the input fasta file
# =============================================================================
classification.classify_seqfile()

print(datetime.now() - startTime)