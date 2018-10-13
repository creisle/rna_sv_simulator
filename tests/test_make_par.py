#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `make_par.py` script."""

import pytest
import os
import sys

from rna_sv_simulator import make_par

BUNDLE_PATH="flux-bundle"

@pytest.fixture()
def par_resources():
    print('\nAppending `'+BUNDLE_PATH+'` to sys.path.')
    
    if os.path.islink(BUNDLE_PATH):
        # Clear if exists.
        os.unlink(BUNDLE_PATH)
        
    os.symlink("../rna_sv_simulator/flux-bundle", BUNDLE_PATH)
    
def test_make_par(par_resources):
    """Test the .PAR file maker."""

    # Make file
    make_par.make_par()

    # Check if it exists
    PAR_FILE = os.path.join(make_par.BUNDLE_PATH, make_par.PAR_FILE)
    assert os.path.isfile(PAR_FILE) == True

    # Check if the required fields are there.
    with open(PAR_FILE, 'r') as myfile:
        data = myfile.read()
        
    assert all( [x in data for x in make_par.REQUIRED_FIELDS ]  )

    
    # Make custom PAR file and test for custom attributes
    make_par.make_par(foo1="fooA",
                      bar2="barB",
                      foobarX="foobarY")

    with open(PAR_FILE, 'r') as myfile:
        data = myfile.read()
    assert all( [x in data for x in ["foo1".upper(), "fooA", "bar2".upper(), "barB", "foobarX".upper(), "foobarY"]  ] ) 

    # Count number of expected strings (i.e. test for duplicated properties/failure to delete stale file."
    EXPECTED_KVS = make_par.REQUIRED_FIELDS + ["foo1".upper(), "fooA", "bar2".upper(), "barB", "foobarX".upper(), "foobarY"] 
    assert sum([data.count(x) for x in EXPECTED_KVS ])  == len(EXPECTED_KVS)
    
