"""
read_forcings.py
"""

import os
from calendar import monthrange
from bisect import bisect_left
from netCDF4 import Dataset, date2num, date2index
from logging import getLogger
from log import LOG_NAME
from time_utility import ord_to_datetime

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
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
        if isinstance(liq_flds, list):
            self.liq_flds = liq_flds
        else:
            self.liq_flds = [liq_flds]
        self.ice_flds = []
        self.timestamp = timestamp
        self.files = []
        self.start_dates = []
        self.end_dates = []
        self.current_file = 'Null'
        self.current_tint = 0

        if start:
            if type(start) == float:
                start = [int(start)]
                end = [int(end)]
            else:
                start = map(int, start.split('-'))
                end = map(int, end.split('-'))

        # single files without year in them
        if not start:
            self.files.append(os.path.join(self.path, file_str))

        # yearly files
        elif len(start) == 1:
            for year in xrange(start[0], end[0]):
                self.files.append(os.path.join(self.path,
                                  file_str.replace('$YYYY', "{0:04d}".format(year))))

        # Monthly
        elif len(start) == 2:
            year = start[0]
            month = start[1]
            while True:
                self.files.append(os.path.join(self.path,
                                  file_str.replace('$YYYY', "{0:04d}".format(year)).replace('$MM', "{0:02d}".format(month))))
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
                self.files.append(os.path.join(self.path,
                                  file_str.replace('$YYYY', "{0:04d}".format(year)).replace('$MM', "{0:02d}".format(month)).replace('$DD', "{0:02d}".format(day))))
                if year == end[0] and month == end[1] and day == end[2]:
                    break
                else:
                    if day == monthrange(year, month)[1]:
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
            log.info('reading forcing file: %s' % fname)
            f = Dataset(fname, 'r+')
            self.start_dates.append(f.variables[self.time_fld][0])
            self.end_dates.append(f.variables[self.time_fld][-1])
            if i == 0:
                # get the calendar and units information
                self.calendar = f.variables[self.time_fld].calendar
                self.time_units = f.variables[self.time_fld].units
            else:
                # check that the units match the first file
                if f.variables[self.time_fld].units != self.time_units:
                    raise ValueError('Units do not match in input files')
                if f.variables[self.time_fld].calendar != self.calendar:
                    raise ValueError('Calendars do not match in input files')

            f.close()

        self.end = ord_to_datetime(self.end_dates[-1], self.time_units, calendar=self.calendar)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # find and open first file
        self.ordtime = date2num(timestamp, self.time_units, calendar=self.calendar)

        if len(self.files) == 1:
            self.current_filenum = 0
        else:
            self.current_filenum = bisect_left(self.start_dates, self.ordtime)
        print self.current_filenum
        self.current_file = self.files[self.current_filenum]
        self.current_fhdl = Dataset(self.current_file, 'r+')

        self.current_tind = date2index(timestamp, self.current_fhdl.variables[self.time_fld])
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def read(self, timestamp):
        """ Read the current timestamp from the data stream """

        # ------------------------------------------------------------ #
        # Get the current data index location
        if self.current_tind == len(self.current_fhdl.variables[self.time_fld])-1:
            # close file and open next
            self.current_fhdl.close()
            self.current_filenum += 1
            self.current_file = self.files[self.current_filenum]
            self.current_fhdl = Dataset(self.current_file, 'r+')
            self.current_tind = 0
        else:
            # move forward one step
            self.current_tind += 1
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # check that the timestamp is what was expected
        expected_timestamp = ord_to_datetime(self.current_fhdl.variables[self.time_fld][self.current_tind],
                                             self.time_units, calendar=self.calendar)
        if timestamp != expected_timestamp:
            log.warning('Timestamp is not what was expected')
            log.warning('Got timestamp %s' %timestamp)
            log.warning('Expected timestamp %s' %expected_timestamp)

            # Attempt to find the correct timestamp
            self.ordtime = date2num(timestamp, self.current_fhdl.variables[self.time_fld].units, calendar=self.calendar)

            # check to make sure the timestamp exists in input dataset
            if self.ordtime > max(self.end_dates):
                raise ValueError('Timestamp is exceeds date range in input files')

            new_filenum = bisect_left(self.start_dates, self.ordtime)
            if new_filenum != self.current_filenum:
                # close the current file and open a new one
                self.current_fhdl.close()
                self.current_file = self.files[new_filenum]
                self.current_fhdl = Dataset(self.current_file, 'r+')

            self.current_tind = date2index(timestamp, self.current_fhdl.variables[self.time_fld])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get the liquid fluxes
        log.info('getting fluxes for %s (%s)' %(timestamp, self.current_tind))
        forcing = {}
        for i, fld in enumerate(self.liq_flds):
            if i == 0:
                forcing['LIQ'] = self.current_fhdl.variables[fld][self.current_tind, :, :]
            else:
                forcing['LIQ'] += self.current_fhdl.variables[fld][self.current_tind, :, :]

        for i, fld in enumerate(self.ice_flds):
            if i == 0:
                forcing['ICE'] = self.current_fhdl.variables[fld][self.current_tind, :, :]
            else:
                forcing['ICE'] += self.current_fhdl.variables[fld][self.current_tind, :, :]

        return forcing
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
