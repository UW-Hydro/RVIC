#!/opt/local/bin/python
"""
rvic_compile_test.py

Set to run with pytest

Usage: py.test (from RVIC or test directory)
"""

import traceback
import sys
sys.path.append("../")


# -------------------------------------------------------------------- #
# Check that all the top level modules compile

def test_pylint_uhs2paramfile():
    try:
        import uhs2paramfile
    except:
        assert traceback.format_exc() == None

def test_pylint_make_parameters():
    try:
        import make_parameters
    except:
        assert traceback.format_exc() == None

def test_pylint_rvic_model():
    try:
        import rvic_model
    except:
        assert traceback.format_exc() == None
# -------------------------------------------------------------------- #
