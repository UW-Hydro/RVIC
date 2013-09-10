"""
history.py
"""
import os
import numpy as np
from netCDF4 import Dataset, date2num, num2date
from datetime import datetime
from time_utility import ord_to_datetime
from logging import getLogger
from log import LOG_NAME
from share import SECSPERDAY, HOURSPERDAY, TIMEUNITS, NC_INT, NC_FLOAT, NC_DOUBLE, WATERDENSITY
import share

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# RVIC History File Object
class Tape(object):
    """ History Tape Object"""

    # ---------------------------------------------------------------- #
    # Init
    def __init__(self, time_ord, caseid, Rvar, fincl=['streamflow'],
                 mfilt=1, ndens=2, nhtfrq=0, avgflag='A', units='kg m-2 s-1',
                 file_format='NETCDF4_CLASSIC', outtype='grid', grid_lons=False,
                 grid_lats=False, grid_area=None, out_dir='.', calendar=None, glob_ats=None):
        self.time_ord = time_ord        # Days since basetime
        self.caseid = caseid            # Case ID and prefix for outfiles
        self.fincl = fincl              # Fields to include in history file
        self.mfilt = mfilt              # Maximum number of time samples
        if ndens == 1:                  # Output file precision
            self.ndens = NC_FLOAT
        else:
            self.ndens = NC_DOUBLE
        self.nhtfrq = nhtfrq            # Write frequency
        self.avgflag = avgflag          # Average Flag (A,I,X,M)
        self.outtype = outtype          # Outfile type (grid, array)
        self.count = 0
        self.files_count = 0
        self.file_format = file_format
        self.calendar = calendar
        self.out_dir = out_dir
        self.__get_rvar(Rvar)           # Get the initial Rvar fields
        self.glob_ats = glob_ats

        # ------------------------------------------------------------ #
        # Get Grid Lons/Lats if outtype is grid
        if outtype == 'grid':
            if type(grid_lons) == np.ndarray and type(grid_lats) == np.ndarray:
                self.grid_lons = grid_lons
                self.grid_lats = grid_lats
                self._gdata = {}
                for field in self.fincl:
                    if grid_lons.ndim == 1:
                        shape = (self.mfilt,) + self.grid_lats.shape + self.grid_lons.shape
                    else:
                        shape = (self.mfilt,) + self.grid_lons.shape
                    self._gdata[field] = np.zeros(shape)
            else:
                raise ValueError('Must include grid lons / lats if outtype == grid')

        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get units multiplier
        self.units = units
        if units in ['kg/m2*s', 'kg m-2 s-1', 'kg m^-2 s^-1', 'kg*m-2*s-1', 'kg s-1 m-2']:
            self.units_mult = 1.0
        elif units in ['m3/s', 'm^3/s', 'm-3 s-1']:
            self.units_mult = grid_area / WATERDENSITY
        else:
            raise ValueError('%s is not a valid units string' %units)

        # ------------------------------------------------------------ #
        # get current timestamp
        self.timestamp = ord_to_datetime(self.time_ord, TIMEUNITS, self.calendar)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Initialize the history fields
        self._data = {}
        for field in self.fincl:
            self._data[field] = np.zeros((self.mfilt, self.num_outlets))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Determine the format of the output filename
        if self.avgflag == 'I':
            self.fname_format = os.path.join(out_dir,
                                             "%s.rvic.h%s.%%Y-%%m-%%d-%%H-%%M-%%S.nc" % (self.caseid, self.avgflag.lower()))
        else:
            if nhtfrq == 0:
                self.fname_format = os.path.join(out_dir,
                                                 "%s.rvic.h%s.%%Y-%%m.nc" % (self.caseid, self.avgflag.lower()))
            elif (nhtfrq == -24) or (nhtfrq*self.dt == SECSPERDAY):
                self.fname_format = os.path.join(out_dir,
                                                 "%s.rvic.h%s.%%Y-%%m-%%d.nc" % (self.caseid, self.avgflag.lower()))
            else:
                self.fname_format = os.path.join(out_dir,
                                                 "%s.rvic.h%s.%%Y-%%m-%%d-%%H.nc" % (self.caseid, self.avgflag.lower()))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Determine when the next write should be
        self.__next_write()
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write a summary
    def __str__(self):
        parts = ['------- Summary of History Tape Settings -------',
        '\t# caseid:       %s' %self.caseid,
        '\t# fincl:        %s' %','.join(self.fincl),
        '\t# nhtfrq:       %s' %self.nhtfrq,
        '\t# mfilt:        %s' %self.mfilt,
        '\t# ndens:        %s' %self.ndens,
        '\t# avgflag:      %s' %self.avgflag,
        '\t# fname_format: %s' %self.fname_format,
        '\t# file_format:  %s' %self.file_format,
        '\t# outtype:      %s' %self.outtype,
        '\t# out_dir:      %s' %self.out_dir,
        '\t# calendar:     %s' %self.calendar,
        '\t# units:        %s' %self.units,
        '  ------- End of History Tape Settings -------']
        return '\n'.join(parts)

    def __repr__(self):
        return '__repr__'
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Update the history tapes with new fluxes
    def update(self, data2tape, time_ord):
        """ Update the tape with new data"""

        # ------------------------------------------------------------ #
        # Check that the time_ord is in sync
        if self.time_ord != time_ord:
            raise ValueError('rout_var.time_ord does not match the time_ord passed in by the convolution call')
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get the next timestamp
        self.time_ord += self.dt / SECSPERDAY
        self.timestamp = ord_to_datetime(self.time_ord, TIMEUNITS, calendar=self.calendar)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Advance the Count
        self.count += 1
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Update the fields
        for field in self.fincl:
            log.debug('updating %s' %field)
            fdata = data2tape[field]
            if self.avgflag == 'A':
                self._data[field][:, :] += fdata
            elif self.avgflag == 'I':
                if self.count == self.write_count:
                    self._data[field][:, :] = fdata
            elif self.avgflag == 'X':
                self._data[field][:, :] = np.maximum(self._data[field][:, :], fdata)
            elif self.avgflag == 'M':
                self._data[field][:, :] = np.minimum(self._data[field][:, :], fdata)
            else:
                raise ValueError('Average flag (%s) does not match any of (A,I,X,M)' %self.avgflag)

        # ------------------------------------------------------------ #
        # If count == write_count, write
        # Average first, if necessary
        if (self.avgflag == 'A' and self.count == self.write_count):
            self.__average()

        if self.count == self.write_count:
            if self.outtype == 'grid':
                self.__write_grid()
            else:
                self.__write_array()
            self.files_count += 1
            self.__next_write()
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write initial flux
    def write_initial(self):
        pass
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get import rvar fields
    def __get_rvar(self, rvar):
        """ Get the rvar Fields that are useful in writing output """
        self.dt = rvar.unit_hydrograph_dt
        self.num_outlets = rvar.n_outlets
        self.outlet_decomp_ind = rvar.outlet_decomp_ind
        self.outlet_x_ind = rvar.outlet_x_ind
        self.outlet_y_ind = rvar.outlet_y_ind
        self.outlet_lon = rvar.outlet_lon
        self.outlet_lat = rvar.outlet_lat
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine next write time
    def __next_write(self):
        """ Determine the count for when the next write should occur """
        # ------------------------------------------------------------ #
        # If monthly, write at (YYYY,MM,1,0,0)
        # b0 is first timestep of next period
        # b1 is end of last timestep of next period

        b0 = self.time_ord

        if self.nhtfrq == 0:
            if self.timestamp.month == 12:
                b1 = date2num(datetime(self.timestamp.year + 1, 2, 1),
                              TIMEUNITS, calendar=self.calendar)
            else:
                b1 = date2num(datetime(self.timestamp.year, self.timestamp.month + 1, 1),
                                       TIMEUNITS, calendar=self.calendar)

        # If some hours in the future
        elif self.nhtfrq < 0:
            b1 = b0 - (self.nhtfrq / HOURSPERDAY)

        # If some dts in the future
        else:
            b1 = b0 + (self.nhtfrq * self.dt / SECSPERDAY)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get the number of timesteps and datestamp for the next write
        # next_ord is the ord_time when the write will happen
        self.write_count = int(round((b1 - b0) / (self.dt / SECSPERDAY)))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Get next file name and timeord
        if self.avgflag == 'I':
            self.write_ord = b1
            self.filename = num2date(b1, TIMEUNITS,
                                     calendar=self.calendar).strftime(self.fname_format)
        else:
            self.time_bnds = np.array([[b0, b1]])
            self.write_ord = np.average(self.time_bnds)
            self.filename = num2date(b0, TIMEUNITS,
                                     calendar=self.calendar).strftime(self.fname_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Zero out field(s)
        for field in self.fincl:
            self._data[field][:] = 0.0
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Set the count to zero
        self.count = 0
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Average fields
    def __average(self):
        """ Take the average based on the number of accumulated timesteps """
        for field in self.fincl:
            self._data[field] /= self.count
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write grid style history file
    def __write_grid(self):
        """ Write history file """

        # ------------------------------------------------------------ #
        # Grid the fields
        for field in self.fincl:
            self._gdata[field][:, self.outlet_y_ind, self.outlet_x_ind] = self._data[field][:, :]
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Open file
        f = Dataset(self.filename, 'w', self.file_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Time Variable
        time = f.createDimension('time', None)

        time = f.createVariable('time', self.ndens, ('time',))
        time[:] = [self.write_ord]
        for key, val in share.time.__dict__.iteritems():
            if val:
                setattr(time, key, val)
        time.calendar = self.calendar

        if self.avgflag != 'I':
            nv = f.createDimension('nv', 2)

            time.bounds = 'time_bnds'

            time_bnds = f.createVariable('time_bnds', self.ndens, ('time', 'nv',))
            time_bnds[:, :] = self.time_bnds
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup Coordinate Variables
        if self.grid_lons.ndim > 1:
            coords = ('yc', 'xc',)

            # Grid is not regular
            xc = f.createDimension('xc', self.grid_lons.shape[1])
            yc = f.createDimension('yc', self.grid_lons.shape[0])

            xc = f.createVariable('xc', self.ndens, coords)
            yc = f.createVariable('yc', self.ndens, coords)
            xc[:, :] = self.grid_lons
            yc[:, :] = self.grid_lats

            for key, val in share.xc.__dict__.iteritems():
                if val:
                    setattr(xc, key, val)

            for key, val in share.yc.__dict__.iteritems():
                if val:
                    setattr(yc, key, val)

        else:
            coords = ('lat', 'lon',)

            lon = f.createDimension('lon', len(self.grid_lons))
            lat = f.createDimension('lat', len(self.grid_lats))

            lon = f.createVariable('lon', self.ndens, ('lon',))
            lat = f.createVariable('lat', self.ndens, ('lat',))
            lon[:] = self.grid_lons
            lat[:] = self.grid_lats

            for key, val in share.lon.__dict__.iteritems():
                if val:
                    setattr(lon, key, val)

            for key, val in share.lat.__dict__.iteritems():
                if val:
                    setattr(lat, key, val)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        tcoords = ('time',) + coords

        for field in self.fincl:
            var = f.createVariable(field, self.ndens, tcoords)
            var[:, :] = self._gdata[field] * self.units_mult

            for key, val in getattr(share, field).__dict__.iteritems():
                if val:
                    setattr(var, key, val)
            var.units = self.units
            if self.grid_lons.ndim > 1:
                var.coordinates = " ".join(coords)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global attributes
        self.glob_ats.update()
        for key, val in self.glob_ats.__dict__.iteritems():
            if val:
                setattr(f, key, val)
        # ------------------------------------------------------------ #
        f.close()
        log.info('Finished writing %s' % self.filename)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write array style history file
    def __write_array(self):
        """ Write history file """

        # ------------------------------------------------------------ #
        # Open file
        f = Dataset(self.filename, 'w', self.file_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Time Variable
        time = f.createDimension('time', None)

        time = f.createVariable('time', self.ndens, ('time',))
        time[:] = [self.write_ord]
        for key, val in share.time.__dict__.iteritems():
            if val:
                setattr(time, key, val)
        time.calendar = self.calendar

        if self.avgflag != 'I':
            nv = f.createDimension('nv', 2)

            time.bounds = 'time_bnds'

            time_bnds = f.createVariable('time_bnds', self.ndens, ('time', 'nv',))
            time_bnds[:, :] = self.time_bnds
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup Coordinate Variables
        coords = ('outlets',)

        outlets = f.createDimension('outlets', self.num_outlets)

        outlet_lon = f.createVariable('lon', self.ndens, coords)
        outlet_lat = f.createVariable('lat', self.ndens, coords)
        outlet_x_ind = f.createVariable('outlet_x_ind', NC_INT, coords)
        outlet_y_ind = f.createVariable('outlet_y_ind', NC_INT, coords)
        outlet_decomp_ind = f.createVariable('outlet_decomp_ind', NC_INT, coords)
        outlet_lon[:] = self.outlet_lon
        outlet_lat[:] = self.outlet_lat
        outlet_x_ind[:] = self.outlet_x_ind
        outlet_y_ind[:] = self.outlet_y_ind
        outlet_decomp_ind[:] = self.outlet_decomp_ind

        for key, val in share.outlet_lon.__dict__.iteritems():
            if val:
                setattr(outlet_lon, key, val)

        for key, val in share.outlet_lat.__dict__.iteritems():
            if val:
                setattr(outlet_lat, key, val)

        for key, val in share.outlet_y_ind.__dict__.iteritems():
            if val:
                setattr(outlet_y_ind, key, val)

        for key, val in share.outlet_x_ind.__dict__.iteritems():
            if val:
                setattr(outlet_x_ind, key, val)

        for key, val in share.outlet_decomp_ind.__dict__.iteritems():
            if val:
                setattr(outlet_decomp_ind, key, val)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        tcoords = ('time',) + coords

        for field in self.fincl:
            var = f.createVariable(field, self.ndens, tcoords)
            var[:, :] = self._data[field] * self.units_mult[self.outlet_y_ind, self.outlet_x_ind]

            for key, val in getattr(share, field).__dict__.iteritems():
                if val:
                    setattr(var, key, val)
            var.units = self.units
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global attributes
        self.glob_ats.update()
        for key, val in self.glob_ats.__dict__.iteritems():
            if val:
                setattr(f, key, val)
        # ------------------------------------------------------------ #
        f.close()
        log.info('Finished writing %s' % self.filename)
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
