#!/opt/local/bin/python
"""
rvic_unit_test.py

Set to run with pytest

Usage: py.test (from RVIC or test directory)
"""
import numpy as np
import sys
sys.path.append("../")

# -------------------------------------------------------------------- #
# Unit tests for utilities.py
from rvic.utilities import *
from rvic.config import *

def test_config_type_int():
    assert config_type('1') == 1

def test_config_type_float():
    assert config_type('1.75') == 1.75

def test_config_type_bool():
    assert config_type('True') == True

def test_find_nearest():
    assert find_nearest(np.array([8, 19, 39, 100, 399]), 20) == 1

def test_find_nearest_max():
    assert find_nearest(np.array([8, 19, 39, 100, 399]), 500) == 4

def test_find_nearest_min():
    assert find_nearest(np.array([8, 19, 39, 100, 399]), 1) == 0
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Unit tests for make_uh.py
from rvic.make_uh import *

def test_find_ts():
    assert find_ts(np.array([0, 86400, 172800])) == 86400
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Unit tests for make_uh.py
from rvic.time_utility import *
from rvic.share import TIMEUNITS
from netCDF4 import date2num
from datetime import datetime

def test_ord_to_datetime():
    # Independence day
    date = datetime(1776, 7, 4, 12, 0, 0, 0)
    ord_time = date2num(date, TIMEUNITS)
    # Independence day (note that this fails if date has microseconds != 0)
    assert ord_to_datetime(ord_time, TIMEUNITS) == date
# -------------------------------------------------------------------- #
