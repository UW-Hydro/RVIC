# -*- coding: utf-8 -*-
'''
history.py

Summary:
    This is the core history file module for the rvic model.
    The core of the module is the Tape class.  The basic procedure is as
    follows:
        - initialization:  sets tape options, determines filenames, etc.
        - update: method that incorporates new fluxes into the history tape.
        - __next_update_out_data: method to determine when to update the
        out_data container
        - __next_write_out_data: method to determine when to write the out_data
        container
        - finish: method to close all remaining history tapes.
'''

import os
import numpy as np
from netCDF4 import Dataset, date2num, num2date, stringtochar
from datetime import datetime
from .time_utility import ord_to_datetime
from logging import getLogger
from .log import LOG_NAME
from .share import SECSPERDAY, HOURSPERDAY, TIMEUNITS
from .share import NC_INT, NC_FLOAT, NC_CHAR
from .share import NC_DOUBLE, WATERDENSITY, MONTHSPERYEAR
from .pycompat import iteritems
from . import share


# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# RVIC History File Object
class Tape(object):
    ''' History Tape Object'''

    # ---------------------------------------------------------------- #
    # Init
    def __init__(self, time_ord, caseid, rvar, tape_num=0,
                 fincl=('streamflow'), mfilt=1, ndens=2, nhtfrq=0,
                 avgflag='A', units='kg m-2 s-1',
                 file_format='NETCDF4_CLASSIC', outtype='grid',
                 grid_lons=False, grid_lats=False, grid_area=None, out_dir='.',
                 calendar='standard', glob_ats=None, zlib=True, complevel=4,
                 least_significant_digit=None):
        self._tape_num = tape_num
        self._time_ord = time_ord        # Days since basetime
        self._caseid = caseid            # Case ID and prefix for outfiles
        self._fincl = list(fincl)        # Fields to include in history file
        self._mfilt = mfilt              # Maximum number of time samples
        self._ndens = ndens
        if self._ndens == 1:             # Output file precision
            self._ncprec = NC_FLOAT
        else:
            self._ncprec = NC_DOUBLE
        self._nhtfrq = nhtfrq            # Write frequency
        self._avgflag = avgflag          # Average Flag (A,I,X,M)
        self._outtype = outtype          # Outfile type (grid, array)
        self._count = 0
        self.files_count = 0
        self._file_format = file_format
        self._calendar = calendar
        self._out_dir = out_dir
        self._glob_ats = glob_ats

        self.__get_rvar(rvar)           # Get the initial rvar fields
        self._grid_shape = grid_area.shape
        self._out_data = {}

        # ------------------------------------------------------------ #
        # calculate the step size for each out_data timestep (units=days)
        if self._nhtfrq > 0:
            # If some number of timesteps
            self._out_data_stepsize = self._nhtfrq * self._dt / SECSPERDAY
        elif self._nhtfrq < 0:
            # If some number hours
            self._out_data_stepsize = -1 * self._nhtfrq / HOURSPERDAY
        else:
            # If monthly
            self._out_data_stepsize = None  # varies by month
        log.debug('_out_data_stepsize: %s', self._out_data_stepsize)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get Grid Lons/Lats if outtype is grid
        if outtype.lower() == 'grid':
            self._out_data_shape = self._grid_shape
            if type(grid_lons) == np.ndarray and type(grid_lats) == np.ndarray:
                self._grid_lons = grid_lons
                self._grid_lats = grid_lats
            else:
                raise ValueError('Must include grid lons / lats if '
                                 'outtype == grid')
        elif outtype.lower() == 'array':
            self._out_data_shape = (self._num_outlets, )
        else:
            raise ValueError('Unknown value for outtype: {0}'.format(outtype))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get units multiplier (size of noutlets)
        self._units = units
        if units in ['kg/m2/s', 'kg m-2 s-1', 'kg m^-2 s^-1',
                     'kg*m-2*s-1', 'kg s-1 m-2']:
            self._units_mult = np.ones_like(self._outlet_y_ind,
                                            dtype=np.float64)
        elif units in ['m3/s', 'm^3/s', 'm3 s-1']:
            # kg/m2/s --> m3/s
            self._units_mult = grid_area[self._outlet_y_ind,
                                         self._outlet_x_ind]
            self._units_mult /= WATERDENSITY
        elif units in ['mm/day', 'mm d-1', 'mm d^-1', 'mm/day']:
            # kg/m2/s --> mm/day over basin area
            self._units_mult = grid_area[self._outlet_y_ind,
                                         self._outlet_x_ind]
            self._units_mult *= SECSPERDAY
            self._units_mult /= WATERDENSITY
            self._units_mult /= self._outlet_upstream_area
        elif units in ['gal/day', 'gpd', 'gal d-1']:
            self._units_mult = grid_area[self._outlet_y_ind,
                                         self._outlet_x_ind]
            self._units_mult /= WATERDENSITY
            self._units_mult *= 2.28E7
        elif units in ['cfs', 'ft^3 s-1', 'f3/s']:
            self._units_mult = grid_area[self._outlet_y_ind,
                                         self._outlet_x_ind]
            self._units_mult /= WATERDENSITY
            self._units_mult *= 35.3
        elif units in ['acre-ft/d']:
            self._units_mult = grid_area[self._outlet_y_ind,
                                         self._outlet_x_ind]
            self._units_mult /= WATERDENSITY
            self._units_mult *= 70.0
        else:
            raise ValueError('{0} is not a valid units string'.format(units))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # netCDF variable options
        self.ncvaropts = {'zlib': zlib,
                          'complevel': complevel,
                          'least_significant_digit': least_significant_digit}
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # get current timestamp
        self._timestamp = ord_to_datetime(self._time_ord, TIMEUNITS,
                                          self._calendar)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Initialize the temporary history fields
        self._temp_data = {}
        for field in self._fincl:
            self._temp_data[field] = np.zeros(self._num_outlets,
                                              dtype=np.float64)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Determine the format of the output filename
        if self._avgflag == 'I':
            self._fname_format = os.path.join(
                out_dir, '%s.rvic.h%s%s.%%Y-%%m-%%d-%%H-%%M-%%S.nc' %
                (self._caseid, self._tape_num, self._avgflag.lower()))
        else:
            if self._nhtfrq == 0:
                self._fname_format = os.path.join(
                    out_dir,
                    '%s.rvic.h%s%s.%%Y-%%m.nc' %
                    (self._caseid, self._tape_num, self._avgflag.lower()))
            elif (self._nhtfrq == -24) or (nhtfrq * self._dt == SECSPERDAY):
                self._fname_format = os.path.join(
                    out_dir,
                    '%s.rvic.h%s%s.%%Y-%%m-%%d.nc' %
                    (self._caseid, self._tape_num, self._avgflag.lower()))
            else:
                self._fname_format = os.path.join(
                    out_dir,
                    '%s.rvic.h%s%s.%%Y-%%m-%%d-%%H.nc' %
                    (self._caseid, self._tape_num, self._avgflag.lower()))
        self._rest_fname_format = os.path.join(
            out_dir,
            '%s.rvic.rh%s.%%Y-%%m-%%d-%%H-%%M-%%S.nc' %
            (self._caseid, self._tape_num))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Determine when the next write should be and initialize out_data
        self.__next_write_out_data()
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Determine when the update of out_data should be
        self.__next_update_out_data()
        # ------------------------------------------------------------ #

        log.debug(self.__repr__())
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write a summary
    def __str__(self):
        return 'History Tape - {0}'.format(self.filename)

    def __repr__(self):
        parts = ['------- Summary of History Tape Settings -------',
                 '\t# caseid:       {0}'.format(self._caseid),
                 '\t# fincl:        {0}'.format(','.join(self._fincl)),
                 '\t# nhtfrq:       {0}'.format(self._nhtfrq),
                 '\t# mfilt:        {0}'.format(self._mfilt),
                 '\t# ncprec:       {0}'.format(self._ncprec),
                 '\t# avgflag:      {0}'.format(self._avgflag),
                 '\t# fname_format: {0}'.format(self._fname_format),
                 '\t# file_format:  {0}'.format(self._file_format),
                 '\t# outtype:      {0}'.format(self._outtype),
                 '\t# out_dir:      {0}'.format(self._out_dir),
                 '\t# calendar:     {0}'.format(self._calendar),
                 '\t# units:        {0}'.format(self._units),
                 '  ------- End of History Tape Settings -------']
        return '\n'.join(parts)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Update the history tapes with new fluxes
    def update(self, data2tape, time_ord):
        ''' Update the tape with new data'''

        # ------------------------------------------------------------ #
        # Check that the time_ord is in sync
        if self._time_ord != time_ord:
            raise ValueError('rout_var.time_ord does not match the time_ord '
                             'passed in by the convolution call')
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get the next timestamp
        self._time_ord += self._dt / SECSPERDAY
        self._timestamp = ord_to_datetime(self._time_ord, TIMEUNITS,
                                          calendar=self._calendar)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Advance the Count
        self._count += 1
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Update the fields
        for field in self._fincl:
            tracer = 'LIQ'
            log.debug('updating %s', field)
            fdata = data2tape[field][tracer]
            if self._avgflag == 'A':
                self._temp_data[field] += fdata
            elif self._avgflag == 'I':
                if self._count == self._update_count:
                    self._temp_data[field] = fdata[:]
            elif self._avgflag == 'X':
                self._temp_data[field] = np.maximum(self._temp_data[field],
                                                    fdata)
            elif self._avgflag == 'M':
                self._temp_data[field] = np.minimum(self._temp_data[field],
                                                    fdata)
            else:
                raise ValueError('Average flag ({0}) does not match any of'
                                 ' (A,I,X,M)'.format(self._avgflag))

        # ------------------------------------------------------------ #
        # If count == _update_count, add to _out_data
        # Average first, if necessary
        if (self._avgflag == 'A' and self._count == self._update_count):
            self.__average()

        if self._count == self._update_count:
            # move the data to the out_data structure
            self.__update_out_data()
            # Determine next update
            self.__next_update_out_data()

            # zero out temp_data
            for field in self._fincl:
                self._temp_data[field][:] = 0.0
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def write_initial(self):
        pass
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def __next_write_out_data(self):
        '''determine the maximum size of out_data'''

        log.debug('determining size of out_data')

        self._out_data_i = 0            # position counter for out_data array

        # ------------------------------------------------------------ #
        # b0 is first timestep of next period
        # b1 is end of last timestep of next period

        # time when out_data will start (now)
        b0 = self._time_ord
        self._begtime = b0

        # determine time when out_data will be full
        if self._mfilt == 'year':
            if self._nhtfrq == 0:
                mfilt = MONTHSPERYEAR
            else:
                t1 = datetime(self._timestamp.year + 1, 1, 1)
                b1 = date2num(t1, TIMEUNITS, calendar=self._calendar)

                # calculate the mfilt value
                mfilt = int(round((b1 - b0) / self._out_data_stepsize))
        elif self._mfilt == 'month':
            if self._nhtfrq == 0:
                mfilt = 1
            else:
                if self._timestamp.month == 12:
                    t1 = datetime(self._timestamp.year + 1, 2, 1)
                else:
                    t1 = datetime(self._timestamp.year,
                                  self._timestamp.month + 1, 1)
                b1 = date2num(t1, TIMEUNITS, calendar=self._calendar)

                # calculate the mfilt value
                mfilt = int(round((b1 - b0) / self._out_data_stepsize))
        elif self._mfilt == 'day':
            if self._nhtfrq != 0:
                b1 = b0 + 1.0
            else:
                raise ValueError('Incompatable values for NHTFRQ and MFILT')

            # calculate the mfilt value
            mfilt = int(round((b1 - b0) / self._out_data_stepsize))
        else:
            mfilt = int(self._mfilt)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        if mfilt < 1:
            mfilt = 1

        self._out_data_write = mfilt - 1
        self._out_times = np.empty(mfilt, dtype=np.float64)
        if self._avgflag != 'I':
            self._out_time_bnds = np.empty((mfilt, 2), dtype=np.float64)

        shape = (mfilt, ) + self._out_data_shape

        log.debug('out_data shape: %s', shape)
        log.debug('_out_data_write: %s', self._out_data_write)

        for field in self._fincl:
            self._out_data[field] = np.zeros(shape, dtype=np.float64)

        self._out_data_has_values = False
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # fill in out_data
    def __update_out_data(self):

        self._out_data_has_values = True

        # ------------------------------------------------------------ #
        # Update the _out_data fields
        for field in self._fincl:
            if self._outtype == 'grid':
                # ---------------------------------------------------- #
                # Grid the fields
                self._out_data[field][self._out_data_i,
                                      self._outlet_y_ind,
                                      self._outlet_x_ind] = \
                    self._temp_data[field][:] * self._units_mult
                # ---------------------------------------------------- #
            else:
                self._out_data[field][self._out_data_i, :] = \
                    self._temp_data[field] * self._units_mult

            # Check that all values are valid, if not, exit gracefully
            if np.isnan(self._out_data[field][self._out_data_i].sum()):
                raise ValueError('nan found in output field: {0}, most likely '
                                 'there is a nan/missing/fill value in the'
                                 'input forcings')
        # ------------------------------------------------------------ #

        self._out_times[self._out_data_i] = self._write_ord

        if self._avgflag != 'I':
            self._out_time_bnds[self._out_data_i, :] = self._time_bnds

        # ------------------------------------------------------------ #
        # if out_data is full, write
        if self._out_data_i == self._out_data_write:
            self.finish()
            self._out_data_i = 0

            # Determine when the next write should be and initialize out_data
            self.__next_write_out_data()
        else:
            self._out_data_i += 1
            log.debug('out_data counter is %s of %s', self._out_data_i,
                      self._out_data_write)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def finish(self):
        '''write out_data'''
        log.debug('finishing tape %s', self._tape_num)
        if self._out_data_has_values:
            if self._outtype == 'grid':
                self.__write_grid()
            else:
                self.__write_array()

            self.files_count += 1
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get import rvar fields
    def __get_rvar(self, rvar):
        ''' Get the rvar Fields that are useful for writing output '''
        self._dt = rvar.unit_hydrograph_dt
        self._num_outlets = rvar.n_outlets
        self._outlet_decomp_ind = rvar.outlet_decomp_ind
        self._outlet_x_ind = rvar.outlet_x_ind
        self._outlet_y_ind = rvar.outlet_y_ind
        self._outlet_lon = rvar.outlet_lon
        self._outlet_lat = rvar.outlet_lat
        self._outlet_name = rvar.outlet_name
        self._outlet_upstream_area = rvar.outlet_upstream_area
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine next write time
    def __next_update_out_data(self):
        ''' Determine the count for when the next write should occur '''
        # ------------------------------------------------------------ #
        # If monthly, write at (YYYY,MM,1,0,0)
        # b0 is first timestep of next period
        # b1 is end of last timestep of next period

        b0 = self._time_ord
        self._begtime = b0

        if self._nhtfrq == 0:
            if self._timestamp.month == 12:
                b1 = date2num(datetime(self._timestamp.year + 1, 2, 1),
                              TIMEUNITS, calendar=self._calendar)
            else:
                b1 = date2num(datetime(self._timestamp.year,
                                       self._timestamp.month + 1, 1),
                              TIMEUNITS, calendar=self._calendar)

        else:
            b1 = b0 + self._out_data_stepsize
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get the number of timesteps and datestamp for the next write
        # next_ord is the ord_time when the write will happen
        self._update_count = int(round((b1 - b0) / (self._dt / SECSPERDAY)))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get next file names and timeord
        if self._avgflag == 'I':
            self._write_ord = b1
            self.filename = num2date(
                b1, TIMEUNITS,
                calendar=self._calendar).strftime(self._fname_format)
        else:
            self._time_bnds = np.array([[b0, b1]])
            self._write_ord = np.average(self._time_bnds)
            self.filename = num2date(
                b0, TIMEUNITS,
                calendar=self._calendar).strftime(self._fname_format)
        self.rest_filename = num2date(
            b1, TIMEUNITS,
            calendar=self._calendar).strftime(self._fname_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Set the count to zero
        self._count = 0
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Average fields
    def __average(self):
        ''' Take the average based on the number of accumulated timesteps '''
        for field in self._fincl:
            self._temp_data[field] /= self._count
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write grid style history file
    def __write_grid(self):
        ''' Write history file '''

        # ------------------------------------------------------------ #
        # Open file
        f = Dataset(self.filename, 'w', self._file_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Time Variable
        f.createDimension('time', None)

        time = f.createVariable('time', self._ncprec, ('time',))
        time[:] = self._out_times[:self._out_data_i + 1]
        for key, val in iteritems(share.time):
            if val:
                setattr(time, key, val.encode())
        time.calendar = self._calendar.encode()

        if self._avgflag != 'I':
            f.createDimension('nv', 2)

            time.bounds = 'time_bnds'

            time_bnds = f.createVariable('time_bnds', self._ncprec,
                                         ('time', 'nv',), **self.ncvaropts)
            time_bnds[:, :] = self._out_time_bnds[:self._out_data_i + 1]
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup Coordinate Variables
        if self._grid_lons.ndim > 1:
            coords = ('yc', 'xc',)

            # Grid is not regular
            xc = f.createDimension('xc', self._grid_lons.shape[1])
            yc = f.createDimension('yc', self._grid_lons.shape[0])

            xc = f.createVariable('xc', self._ncprec, coords, **self.ncvaropts)
            yc = f.createVariable('yc', self._ncprec, coords, **self.ncvaropts)
            xc[:, :] = self._grid_lons
            yc[:, :] = self._grid_lats

            for key, val in iteritems(share.xc):
                if val:
                    setattr(xc, key, val.encode())

            for key, val in iteritems(share.yc):
                if val:
                    setattr(yc, key, val.encode())

        else:
            coords = ('lat', 'lon',)

            lon = f.createDimension('lon', len(self._grid_lons))
            lat = f.createDimension('lat', len(self._grid_lats))

            lon = f.createVariable('lon', self._ncprec, ('lon',),
                                   **self.ncvaropts)
            lat = f.createVariable('lat', self._ncprec, ('lat',),
                                   **self.ncvaropts)
            lon[:] = self._grid_lons
            lat[:] = self._grid_lats

            for key, val in iteritems(share.lon):
                if val:
                    setattr(lon, key, val.encode())

            for key, val in iteritems(share.lat):
                if val:
                    setattr(lat, key, val.encode())
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        tcoords = ('time',) + coords

        for field in self._fincl:
            var = f.createVariable(field, self._ncprec, tcoords,
                                   **self.ncvaropts)
            var[:, :] = self._out_data[field][:self._out_data_i + 1]

            for key, val in iteritems(getattr(share, field)):
                if val:
                    setattr(var, key, val.encode())
            var.units = self._units.encode()
            if self._grid_lons.ndim > 1:
                var.coordinates = ' '.join(coords).encode()
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global attributes
        self._glob_ats.update()
        for key, val in iteritems(self._glob_ats.atts):
            if val:
                setattr(f, key, val.encode())
        # ------------------------------------------------------------ #
        f.close()
        log.info('Finished writing %s', self.filename)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write array style history file
    def __write_array(self):
        ''' Write history file '''

        # ------------------------------------------------------------ #
        # Open file
        f = Dataset(self.filename, 'w', self._file_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Time Variable
        f.createDimension('time', None)

        time = f.createVariable('time', self._ncprec, ('time',),
                                **self.ncvaropts)
        time[:] = self._out_times[:self._out_data_i]
        for key, val in iteritems(share.time):
            if val:
                setattr(time, key, val.encode())
        time.calendar = self._calendar.encode()

        if self._avgflag != 'I':
            f.createDimension('nv', 2)

            time.bounds = 'time_bnds'

            time_bnds = f.createVariable('time_bnds', self._ncprec,
                                         ('time', 'nv',), **self.ncvaropts)
            time_bnds[:, :] = self._out_time_bnds[:self._out_data_i]
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup Coordinate Variables
        coords = ('outlets',)

        f.createDimension('outlets', self._num_outlets)

        nocoords = coords + ('nc_chars',)
        char_names = stringtochar(self._outlet_name)
        f.createDimension(nocoords[1], char_names.shape[1])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Variables
        outlet_lon = f.createVariable('lon', self._ncprec, coords,
                                      **self.ncvaropts)
        outlet_lat = f.createVariable('lat', self._ncprec, coords,
                                      **self.ncvaropts)
        outlet_x_ind = f.createVariable('outlet_x_ind', NC_INT, coords,
                                        **self.ncvaropts)
        outlet_y_ind = f.createVariable('outlet_y_ind', NC_INT, coords,
                                        **self.ncvaropts)
        outlet_decomp_ind = f.createVariable('outlet_decomp_ind', NC_INT,
                                             coords, **self.ncvaropts)
        onm = f.createVariable('outlet_name', NC_CHAR, nocoords,
                               **self.ncvaropts)

        outlet_lon[:] = self._outlet_lon
        outlet_lat[:] = self._outlet_lat
        outlet_x_ind[:] = self._outlet_x_ind
        outlet_y_ind[:] = self._outlet_y_ind
        outlet_decomp_ind[:] = self._outlet_decomp_ind
        onm[:, :] = char_names

        for key, val in iteritems(share.outlet_lon):
            if val:
                setattr(outlet_lon, key, val.encode())

        for key, val in iteritems(share.outlet_lat):
            if val:
                setattr(outlet_lat, key, val.encode())

        for key, val in iteritems(share.outlet_y_ind):
            if val:
                setattr(outlet_y_ind, key, val.encode())

        for key, val in iteritems(share.outlet_x_ind):
            if val:
                setattr(outlet_x_ind, key, val.encode())

        for key, val in iteritems(share.outlet_decomp_ind):
            if val:
                setattr(outlet_decomp_ind, key, val.encode())

        for key, val in iteritems(share.outlet_name):
            if val:
                setattr(onm, key, val.encode())
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        tcoords = ('time',) + coords

        for field in self._fincl:
            var = f.createVariable(field, self._ncprec, tcoords,
                                   **self.ncvaropts)
            var[:, :] = self._out_data[field][:self._out_data_i]

            for key, val in iteritems(getattr(share, field)):
                if val:
                    setattr(var, key, val.encode())
            var.units = self._units.encode()
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global attributes
        self._glob_ats.update()
        for key, val in iteritems(self._glob_ats.atts):
            if val:
                setattr(f, key, val.encode())
        f.featureType = 'timeSeries'.encode()
        # ------------------------------------------------------------ #
        f.close()
        log.info('Finished writing %s', self.filename)
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
