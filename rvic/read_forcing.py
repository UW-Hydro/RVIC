"""
read_forcings.py
"""

from netCDF4 import Dataset, num2date, date2num, date2index
from log import log_name
import calendar
import bisect
import logging

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(log_name)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Data Model
class DataModel(object):
    """ RVIC Forcing Data Model Class"""

    # ---------------------------------------------------------------- #
    # Initialize
    def __init__(self, path, file_str, time_fld, liq_flds,
                 start, end, timestamp):

        self.path = path
        self.time_fld = time_fld
        self.liq_flds = liq_flds
        self.timestamp = timestamp
        self.files = []
        self.start_dates = []
        self.end_dates = []
        self.current_file = 'Null'
        self.current_tint = 0

        if type(start) == float:
            start = [int(start)]
            end = [int(end)]
        else:
            start = map(int, start.split('-'))
            end = map(int, end.split('-'))

        # yearly files
        if len(start) == 1:
            for year in xrange(start[0], end[0]):
                self.files.append(file_str.replace('$YYYY', str(year)))

        # Monthly
        elif len(start) == 2:
            year = start[0]
            month = start[1]
            while True:
                self.files.append(file_str.replace('$YYYY', str(year)).replace('$MM', str(month)))
                if year == end[0] and month == end[1]:
                    break
                else:
                    if month == 12:
                        year += 1
                        month = 1
                    else:
                        month += 1

        # Daily files
        elif len(start) == 3:
            year = start[0]
            month = start[1]
            day = start[2]
            while True:
                self.files.append(file_str.replace('$YYYY', str(year)).replace('$MM', str(month)).replace('$DD', str(day)))
                if year == end[0] and month == end[1] and day == end[2]:
                    break
                else:
                    if day == calendar.monthrange(year, month)[1]:
                        day = 1
                        if month == 12:
                            year += 1
                            month = 1
                        else:
                            month += 1
                    else:
                        day += 1
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # find time bounds for input files
        for i, fname in enumerate(self.files):
            f = Dataset(fname, 'r+')
            self.start_dates.append(f.variables[self.time_fld][0])
            self.end_dates.append(f.variables[self.time_fld][-1])
            if i == 0:
                # get the calendar and units information
                self.calendar = f.variables[self.time_fld].calendar
                self.time_units = f.variables[self.time_fld].units
            else:
                # check that the units match the first file
                if f.variables[self.time_fld].units != self.units:
                    log.error('Units do not match in input files')
                    raise ValueError
                if f.variables[self.time_fld].calendar != self.calendar:
                    log.error('Calendars do not match in input files')
                    raise ValueError
            f.close()
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # find and open first file
        self.ordtime = date2num(timestamp, self.units, calendar=self.calendar)

        self.current_filenum = bisect.bisect_right(self.start_dates, self.ordtime)
        self.current_file = self.files[self.current_filenum]
        self.current_fhdl = Dataset(self.current_file, 'r+')

        self.current_tind = date2index(timestamp, self.current_fhdl.variables[self.time_fld],
                                       select='exact')
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def read(self, timestamp):
        """ Read the current timestamp from the data stream """

        # ------------------------------------------------------------ #
        # check that the timestamp is what was expected
        if timestamp != num2date(self.current_fhdl.variables[self.time_fld][self.current_tind],
                                 self.time_units, calendar=self.calendar):
            log.warning('Timestamp is not what was expected')

            # Attempt to find the correct timestamp
            self.ordtime = date2num(timestamp, self.units, calendar=self.calendar)

            # check to make sure the timestamp exists in input dataset
            if self.ordtime > max(self.end_dates):
                log.error('Timestamp is exceeds date range in input files')
                raise

            new_filenum = bisect.bisect_right(self.start_dates, self.ordtime)
            if new_filenum != self.current_filenum:
                # close the current file and open a new one
                self.current_fhdl.close()
                self.current_file = self.files[self.current_filenum]
                self.current_fhdl = Dataset(self.current_file, 'r+')

            self.current_tind = date2index(timestamp, self.current_fhdl.variables[self.time_fld],
                                           select='exact')
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get the liquid fluxes
        for i, fld in enumerate(self.liq_flds):
            if i == 0:
                forcing = self.current_fhdl.variables[fld][self.current_tind, :, :]
            else:
                forcing += self.current_fhdl.variables[fld][self.current_tind, :, :]

        if self.current_tind == len(self.current_fhdl[self.time_fld]):
            # close file and open next
            self.current_fhdl.close()
            self.current_filenum += 1
            self.current_file = self.files[self.current_filenum]
            self.current_fhdl = Dataset(self.current_file, 'r+')
            self.current_tind = 0
        else:
            # move forward one step
            self.current_tind += 1

        return forcing
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
