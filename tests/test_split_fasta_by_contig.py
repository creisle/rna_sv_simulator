#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `make_par.py` script."""

import pytest
import os
import sys

from rna_sv_simulator import split_fasta_by_contig

@pytest.fixture()
def par_resources():
    print('\nAppending `flux-bundle` to sys.path.')
    os.symlink("rna_sv_simulator/flux-bundle", "flux-bundle")
    
def test_main(par_resources):
    """Test the .PAR file maker."""
    make_par.make_par()
    #assert os.path.isfile(os.path.join(make_par.BUNDLE_PATH, make_par.PAR_FILE)) == True
