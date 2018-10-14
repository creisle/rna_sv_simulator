import pandas as pd
import numpy as np
import os
import sys
import re

# Read metadata (add later)

# Read cosmic data 
cosmic_GRCh37 = pd.read_table("../raw_data/CosmicCompleteGeneExpression_GRCh37.tsv")

# Read GTF file 
ens69_gtf = pd.read_table("../raw_data/Homo_sapiens.GRCh37.69.gtf", header = None, sep = '\t')

