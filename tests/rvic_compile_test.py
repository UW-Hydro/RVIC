#!/opt/local/bin/python
"""
rvic_compile_test.py

Set to run with pytest

Usage: py.test (from RVIC or test directory)
"""

import traceback
import sys
import py_compile as pyc
sys.path.append("../")


# -------------------------------------------------------------------- #
# Check that all the top level modules compile

def test_pylint_uhs2paramfile():
    try:
        pyc.compile('../uhs2paramfile.py')
    except:
        assert traceback.format_exc() == None

def test_pylint_make_parameters():
    try:
        pyc.compile('../make_parameters.py')
    except:
        assert traceback.format_exc() == None

def test_pylint_rvic_model():
    try:
        pyc.compile('../rvic_model.py')
    except:
        assert traceback.format_exc() == None
# -------------------------------------------------------------------- #
