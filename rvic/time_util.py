"""
time_util.py
"""

from log import log_name
from netCDF4 import num2date, date2num
from datetime import datetime
from dateutil.relativedelta import relativedelta
import logging
from share import timeUnits, secsPerDay, minsPerDay, hoursPerDay, timeStampForm

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(log_name)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# RVIC Time Class
class Dtime(object):
    """ A Time Module for handling flags and timesteps """

    # ---------------------------------------------------------------- #
    # Setup Dtim object
    def __init__(self, start_date, stop_option, stop_n, stop_date,
                 rest_option, rest_n, rest_date, calendar, dt):

        self.start_date = datetime.strptime(start_date, timeStampForm)
        self.calendar = calendar
        self.end = False
        self.dt = float(dt) / float(secsPerDay)  # In days

        # Setup Current Time
        self.timestamp = self.start_date
        self.time_ord = date2num(self.timestamp, timeUnits, calendar=self.calendar)

        # Stop option
        self.stop_option = stop_option
        self.stop_n = stop_n
        if self.stop_option == 'date':
            date = map(int, stop_date.split('-'))
            self.stop_date = datetime(*date)
        else:
            self.stop_date = False

        # Rest option
        self.rest_option = rest_option
        self.rest_n = rest_n
        if self.rest_option == 'date':
            date = map(int, rest_date.split('-'))
            self.rest_date = datetime(*date)
        else:
            self.rest_date = False

        # Counters
        self.timesteps = 0
        self.statefiles = 0

        # Flags
        self.stop_flag = False
        self.rest_flag = False

    def advance_timestep(self):
        self.time_ord += self.dt
        self.timestamp = num2date(self.time_ord, self.time_units, calendar=self.calendar)
        self.timesteps += 1
        self.stop_flag = self.stop()
        self.rest_flag = self.rest()
        return self.timestamp
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Time to stop run
    def stop(self):
        if self.stop_option == 'nsteps':
            if self.timesteps >= self.stop_n:
                return True
        elif self.stop_option == 'nseconds':
            if (self.timesteps * self.dt / secsPerDay) >= self.stop_n:
                return True
        elif self.stop_option == 'nminutes':
            if (self.timesteps * self.dt / minsPerDay) >= self.stop_n:
                return True
        elif self.stop_option == 'nhours':
            if (self.timesteps * self.dt / hoursPerDay) >= self.stop_n:
                return True
        elif self.stop_option == 'ndays':
            if (self.timesteps * self.dt) >= self.stop_n:
                return True
        elif self.stop_option == 'nmonths':
            if relativedelta(self.timestamp, self.run_startdate).months >= self.stop_n:
                return True
        elif self.stop_option == 'nyears':
            if relativedelta(self.timestamp, self.run_startdate).years >= self.stop_n:
                return True
        elif self.stop_option == 'date':
            if self.timestamp >= self.stop_date:
                return True
        elif self.end:
            return True
        else:
            return False
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Time to write restart?
    def rest(self):
        if self.rest_option == 'nsteps':
            if self.timesteps >= self.rest_n:
                return True
        elif self.rest_option == 'nseconds':
            if (self.timesteps * self.dt / secsPerDay) >= self.rest_n:
                return True
        elif self.rest_option == 'nminutes':
            if (self.timesteps * self.dt / minsPerDay) >= self.rest_n:
                return True
        elif self.rest_option == 'nhours':
            if (self.timesteps * self.dt / hoursPerDay) >= self.rest_n:
                return True
        elif self.rest_option == 'ndays':
            if (self.timesteps * self.dt) >= self.rest_n:
                return True
        elif self.rest_option == 'nmonths':
            if relativedelta(self.timestamp, self.run_startdate).months >= self.rest_n:
                return True
        elif self.rest_option == 'nyears':
            if relativedelta(self.timestamp, self.run_startdate).years >= self.rest_n:
                return True
        elif self.rest_option == 'date':
            if self.timestamp >= self.rest_date:
                return True
        elif self.end:
            return True
        else:
            return False
    # ---------------------------------------------------------------- #

# -------------------------------------------------------------------- #
