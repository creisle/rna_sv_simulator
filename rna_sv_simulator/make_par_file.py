# -*- coding: utf-8 -*-
"""Helper script to generate Flux Simulator .PAR files"""
# Spec: http://confluence.sammeth.net/display/SIM/.PAR+Simulation+Parameters

import os
import sys

REQUIRED_FIELDS = ["REF_FILE_NAME",
                   "GEN_DIR",
                   "NB_MOLECULES",
                   "READ_NUMBER",
                   "ERR_FILE",
                   "FASTA"
                   ]

BUNDLE_PATH = "flux-bundle/"
PAR_FILE = "test.par"

def make_par(ref_file_name="annotation.gtf",
             gen_dir="genome",
             nb_molecules=100000,
             read_number=10000,
             err_file=76,
             fasta='YES',
             **kwargs):
    """Write given parameters to PAR file."""

    properties = {**locals(), **kwargs}
    properties.pop("kwargs")

    # Clear stale file
    if os.path.isfile(PAR_FILE):
        os.remove(PAR_FILE)

    for key, value in properties.items():
        print_key_value(key=key.upper(), value=value, output_file=PAR_FILE)


def print_row(elements, output_filename, delimiter='\t'):
    """Print elements to output_file, delimited by delimiter."""
    with open(os.path.join(BUNDLE_PATH, output_filename), 'a') as output_file:
        if elements is not None:
            output_file.write(delimiter.join([str(x) for x in elements]))  # Data row
        output_file.write('\n')  # Empty line


def print_key_value(key=None, value=None, output_file=None):
    """Helper function to generate key-value rows."""
    try:
        print_row([key, value], output_file)
    except IOError as error:
        print(error.errno)
        print(error)
        sys.exit(-1)

if __name__ == "__main__":
    # Should fail
    # make_par(apropertythatdoesntexist="12345")

    # Default, should work.
    # make_par()

    # Now adding a specific parameter
    make_par(READ_LENGTH=50)
