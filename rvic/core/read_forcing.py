# -*- coding: utf-8 -*-
'''
read_forcings.py
'''

import os
import numpy as np
from calendar import monthrange
from bisect import bisect_right
from netCDF4 import Dataset, date2num, date2index, num2date
from logging import getLogger
from .log import LOG_NAME
from .time_utility import ord_to_datetime
from .share import MMPERMETER, CMPERMETER, WATERDENSITY, TIMEUNITS, SECSPERDAY
from .pycompat import pyrange

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Data Model
class DataModel(object):
    '''RVIC Forcing Data Model Class'''

    # ---------------------------------------------------------------- #
    # Initialize
    def __init__(self, path, file_str, time_fld, lat_fld, liq_flds,
                 start, end):

        self.path = path
        self.time_fld = time_fld
        self.lat_fld = lat_fld
        if isinstance(liq_flds, list):
            self.liq_flds = liq_flds
        else:
            self.liq_flds = [liq_flds]
        self.ice_flds = []
        self.files = []
        self.start_dates = []
        self.end_ords = []
        self.current_file = 'Null'
        self.current_tint = 0
        self.fld_mult = {}

        if start:
            if type(start) in [float, int]:
                start = [int(start)]
                end = [int(end)]
            else:
                start = list(map(int, start.split('-')))
                end = list(map(int, end.split('-')))

        # single files without year in them
        if not start:
            self.files.append(os.path.join(self.path, file_str))

        # yearly files
        elif len(start) == 1:
            for year in pyrange(start[0], end[0] + 1):
                self.files.append(os.path.join(self.path,
                                  file_str.replace('$YYYY',
                                                   '{0:04d}'.format(year))))

        # Monthly
        elif len(start) == 2:
            year = start[0]
            month = start[1]
            while True:
                self.files.append(
                    os.path.join(
                        self.path,
                        file_str.replace(
                            '$YYYY',
                            '{0:04d}'.format(year)).replace(
                                '$MM', '{0:02d}'.format(month))))
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
                self.files.append(
                    os.path.join(
                        self.path,
                        file_str.replace(
                            '$YYYY',
                            '{0:04d}'.format(year)).replace(
                                '$MM',
                                '{0:02d}'.format(month)).replace(
                                    '$DD',
                                    '{0:02d}'.format(day))))
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
        # find time bounds and timestep for input files
        for i, fname in enumerate(self.files):
            log.info('reading forcing file: %s', fname)
            f = Dataset(fname, 'r')
            self.start_dates.append(f.variables[self.time_fld][0])
            self.end_ords.append(f.variables[self.time_fld][-1])

            if i == 0:
                # get the calendar and units information
                self.calendar = f.variables[self.time_fld].calendar
                self.time_units = f.variables[self.time_fld].units
                time_series = f.variables[self.time_fld][:]

                # determine the latitude order
                lats = f.variables[self.lat_fld][:]
                if lats.ndim == 1:
                    if lats[-1] > lats[0]:
                        log.debug('Input fluxes came in upside down, flipping '
                                  'params and maybe domain.')
                        self.lat0_is_min = True
                    else:
                        self.lat0_is_min = False
                else:
                    self.lat0_is_min = False
            else:
                # check that the units match the first file
                if f.variables[self.time_fld].units != self.time_units:
                    raise ValueError('Units do not match in input files')
                if f.variables[self.time_fld].calendar != self.calendar:
                    raise ValueError('Calendars do not match in input files')

                time_series = np.append(time_series,
                                        f.variables[self.time_fld][:])

            f.close()

        self.end = ord_to_datetime(self.end_ords[-1], self.time_units,
                                   calendar=self.calendar)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # find timestep information
        if len(time_series) > 1:
            t0 = date2num(num2date(time_series[0], self.time_units,
                                   calendar=self.calendar),
                          TIMEUNITS, calendar=self.calendar)
            t1 = date2num(num2date(time_series[1], self.time_units,
                                   calendar=self.calendar),
                          TIMEUNITS, calendar=self.calendar)
            self.secs_per_step = (t1 - t0) * SECSPERDAY
        else:
            raise ValueError('Need more than 1 forcing timestep in order to '
                             'calculate timestep')
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def start(self, timestamp, rout_var):
        '''Initialize the first files inputs'''
        # ------------------------------------------------------------ #
        # find and open first file
        self.ordtime = date2num(timestamp, self.time_units,
                                calendar=self.calendar)

        if len(self.files) == 1:
            self.current_filenum = 0
        else:
            self.current_filenum = bisect_right(self.start_dates,
                                                self.ordtime) - 1

        log.debug('Filenum %s', self.current_filenum)

        self.current_file = self.files[self.current_filenum]
        self.current_fhdl = Dataset(self.current_file, 'r')

        try:
            self.current_tind = date2index(
                timestamp, self.current_fhdl.variables[self.time_fld])
        except Exception as e:
            log.error(num2date(self.start_dates, self.time_units,
                               calendar=self.calendar))
            log.error('timestamp: %s', timestamp)
            log.exception(e)
            raise
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # find multiplier for units all fields will be in liquid flux
        # (kg m-2 s-1)
        # Todo: use udunits to convert units
        for fld in self.liq_flds:
            units = self.current_fhdl.variables[fld].units

            if units in ['kg/m2*s', 'kg m-2 s-1', 'kg m^-2 s^-1',
                         'kg*m-2*s-1', 'kg s-1 m-2']:
                self.fld_mult[fld] = 1.0
            elif units in ['mm', 'MM', 'milimeters', 'Milimeters']:
                self.fld_mult[fld] = WATERDENSITY / MMPERMETER / self.secs_per_step
            elif units in ['m', 'M', 'meters', 'Meters']:
                self.fld_mult[fld] = WATERDENSITY / self.secs_per_step
            elif units in ['cm', 'CM', 'centimeters', 'Centimeters']:
                self.fld_mult[fld] = (WATERDENSITY / CMPERMETER /
                                      self.secs_per_step)
            else:
                raise ValueError('unknown forcing units: %s' % units)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Compare rout_var mask to forcing mask
        for fld in self.liq_flds:
            self.current_fhdl.variables[fld].set_auto_maskandscale(False)
            try:
                fill_val = getattr(self.current_fhdl.variables[fld],
                                   '_FillValue')
                fmask = np.squeeze(
                    self.current_fhdl.variables[fld][0]) == fill_val

                if np.any(fmask[rout_var.source_y_ind, rout_var.source_x_ind]):
                    log.error('There are fill values in the routing domain')
                    log.error('Can not continue...')
                    raise ValueError('Exiting now due to fill values being '
                                     'inside the domain defined by the RVIC '
                                     'parameter file')

            except AttributeError:
                pass
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def read(self, timestamp):
        '''Read the current timestamp from the data stream '''

        # ------------------------------------------------------------ #
        # Get the current data index location
        if self.current_tind == len(
                self.current_fhdl.variables[self.time_fld]):
            # close file and open next
            self.current_fhdl.close()
            self.current_filenum += 1
            self.current_file = self.files[self.current_filenum]
            self.current_fhdl = Dataset(self.current_file, 'r')
            self.current_tind = 0
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # check that the timestamp is what was expected
        expected_timestamp = ord_to_datetime(
            self.current_fhdl.variables[self.time_fld][self.current_tind],
            self.time_units,
            calendar=self.calendar)
        if timestamp != expected_timestamp:
            log.warning('Timestamp is not what was expected')
            log.warning('Got timestamp %s', timestamp)
            log.warning('Expected timestamp %s', expected_timestamp)

            # Attempt to find the correct timestamp
            self.ordtime = date2num(
                timestamp,
                self.current_fhdl.variables[self.time_fld].units,
                calendar=self.calendar)

            # check to make sure the timestamp exists in input dataset
            if self.ordtime > max(self.end_ords):
                raise ValueError('Timestamp is exceeds date range in '
                                 'input files')

            new_filenum = bisect_right(self.start_dates, self.ordtime) - 1
            if new_filenum != self.current_filenum:
                # close the current file and open a new one
                self.current_fhdl.close()
                self.current_file = self.files[new_filenum]
                self.current_fhdl = Dataset(self.current_file, 'r')

            self.current_tind = date2index(
                timestamp, self.current_fhdl.variables[self.time_fld])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get the liquid fluxes
        log.info('getting fluxes for %s (%s)', timestamp, self.current_tind)
        forcing = {}
        for flux, fields in (('LIQ', self.liq_flds), ('ICE', self.ice_flds)):
            for i, fld in enumerate(fields):
                # set auto_maskandscale to False to avoid slowdown from numpy
                # masked array
                self.current_fhdl.variables[fld].set_auto_maskandscale(False)
                temp = self.current_fhdl.variables[fld][self.current_tind]
                fill_val = getattr(self.current_fhdl.variables[fld],
                                   '_FillValue', None)
                # replace fill values with numpy.nan
                if fill_val is not None:
                    temp[temp == fill_val] = np.nan

                if i == 0:
                    forcing[flux] = self.fld_mult[fld] * temp
                else:
                    forcing[flux] += self.fld_mult[fld] * temp

        # move forward one step
        self.current_tind += 1

        for field in forcing:
            forcing[field] = forcing[field].astype(np.float64)

        return forcing
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
