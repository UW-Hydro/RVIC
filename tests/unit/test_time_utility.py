# -------------------------------------------------------------------- #
# Unit tests for make_uh.py
from rvic.core.time_utility import ord_to_datetime, Dtime
from rvic.core.share import TIMEUNITS
from netCDF4 import date2num
from datetime import datetime


def test_ord_to_datetime():
    # Independence day
    date = datetime(1776, 7, 4, 12, 0, 0, 0)
    ord_time = date2num(date, TIMEUNITS)
    # Independence day (note that this fails if date has microseconds != 0)
    assert ord_to_datetime(ord_time, TIMEUNITS) == date


def test_dtime():
    dt = Dtime('2014-12-01-00', 'ndays', 5, None, 'ndays', 5, None, 'noleap',
               3600.00001)
    assert dt.timestamp.year == 2014
    assert dt.timestamp.month == 12
    assert dt.timestamp.day == 1
    assert dt.timestamp.hour == 0
    dt.advance_timestep()
    assert dt.timestamp.year == 2014
    assert dt.timestamp.month == 12
    assert dt.timestamp.day == 1
    assert dt.timestamp.hour == 1
