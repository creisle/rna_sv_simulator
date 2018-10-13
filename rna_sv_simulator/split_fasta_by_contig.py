# -*- coding: utf-8 -*-
"""Divide a fasta file such that each sequence gets its own file with header matching name."""

import os
import sys

FASTA_FILE=sys.argv[1]

if __name__ == "__main__":
    current_filename = None
    output_file = open(os.devnull, 'w')

    with open(FASTA_FILE, 'r') as fasta:
        for row in fasta:
            if row[0] == ">": # Case where row is a header. Create new file.
                current_filename = row.split(" ")[0] # Get contig name.
                current_filename = current_filename[1:] # Remove `>`.
                current_filename = current_filename.strip() # Strip trailling/leading spaces.
                current_filename = current_filename + ".fa" # Add fasta file extension.
        
                output_file.close()
                output_file = open(current_filename, 'w')
                
            # Dump row to file.
            output_file.write(row)
