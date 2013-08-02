"""
variables.py
"""
import os
import numpy as np
from netCDF4 import Dataset, date2num
from logging import getLogger
from log import LOG_NAME
from share import TIMEUNITS, NC_INT, NC_DOUBLE, NcGlobals
import share

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Point object
class Point(object):
    '''Creates a point class for intellegently storing coordinate information'''

    def __init__(self, lat='', lon='', x='', y=''):
        '''Defines x and y variables'''
        self.lat = lat
        self.lon = lon
        self.x = x
        self.y = y

    def __str__(self):
        return "Point(%s,%s)" % (self.lon, self.lat)

    def __repr__(self):
        return '__repr__'

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Rvar Object
class Rvar(object):
    """ Creates a RVIC structure """

    # ---------------------------------------------------------------- #
    # Initialize
    def __init__(self, ParamFile, CaseName, timestamp, calendar, GlobAts, out_dir):
        self.ParamFile = ParamFile
        f = Dataset(ParamFile, 'r+')
        self.n_sources = len(f.dimensions['sources'])
        self.n_outlets = len(f.dimensions['outlets'])
        self.subset_length = f.variables['subset_length'][:]
        self.full_time_length = f.variables['full_time_length'][:]
        self.unit_hydrogaph_dt = f.variables['unit_hydrogaph_dt'][:]
        self.x_ind_source = f.variables['x_ind_source'][:]
        self.y_ind_source = f.variables['y_ind_source'][:]
        self.t_offset_source = f.variables['t_offset_source'][:]
        self.source2outlet_index = f.variables['source2outlet_index'][:]
        self.outlet_decomp_id = f.variables['outlet_decomp_id '][:]
        self.lon_outlet = f.variables['lon_outlet '][:]
        self.lat_outlet = f.variables['lat_outlet'][:]
        self.x_ind_outlet = f.variables['x_ind_outlet '][:]
        self.y_ind_outlet = f.variables['y_ind_outlet'][:]
        self.unit_hydrographs = f.variables['unit_hydrographs'][:]
        self.RvicPourPointsFile = f.RvicPourPointsFile
        self.RvicUHFile = f.RvicUHFile
        self.RvicFdrFile = f.RvicFdrFile
        self.RvicDomainFile = f.RvicDomainFile
        f.close()

        # ------------------------------------------------------------ #
        # Initialize state variables
        self.ring = np.zeros((self.full_time_length, self.n_outlets))
        self.ring_time = np.zeros(self.n_outlets)
        # ------------------------------------------------------------ #

        self.caseName = CaseName
        self.calendar = calendar
        self.fname_format = os.path.join(out_dir, "%s.r.%%Y-%%m-%%d-%%H-%%M-%%S.nc" % (self.caseName))
        self.timestamp = timestamp
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initilize State
    def init_state(self, StateFile):
        f = Dataset(StateFile, 'r+')
        self.ring = f.variables['ring'][:]
        self.ring_time = f.variables['time'][:]
        self.StateFile = StateFile

        # Check that timestep and outlet_decomp_ids match ParamFile
        if f.variables['timestep'][:] != self.unit_hydrogaph_dt:
            log.error('Timestep in Statefile does not match timestep in ParamFile')
            raise
        if f.variables['s_outlet_decomp_id'][:] != self.outlet_decomp_id:
            log.error('outlet_decomp_id in Statefile does not match ParamFile')
            raise
        if f.variables['RvicDomainFile'][:] != self.RvicDomainFile:
            log.error('RvicDomainFile in Statefile does not match ParamFile')
            raise
        f.close()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Convolve
    def convolve(self, aggrunin, timestamp):
        """
        This convoluition funciton works by looping over all points and doing the
        convolution one timestep at a time.  This is accomplished by creating an
        convolution ring.  Contributing flow from each timestep is added to the
        convolution ring.  The convolution ring is saved as the state.  The first
        column of values in the ring are the current runoff.
        """
        if timestamp > self.timestamp:
            self.timestamp = timestamp
        else:
            log.error('Timestep did not advance between convolution calls')
            raise

        # First update the ring
        self.ring[0, :] = 0                            # Zero out current ring
        self.ring = self.__cshift(self.ring, 1)

        # this matches the fortran implementation, it may be faster to use np.convolve but testing
        # can be done later
        for s, outlet in enumerate(self.source2outlet_index):   # loop over all source points
            y = self.y_ind_source[s]
            x = self.x_ind_source[s]
            for i in xrange(self.subset_length):
                j = i + self.t_offset_source[s]
                self.ring[j, outlet] = self.ring[j, outlet] + (self.unit_hydrogaphs[i, s] * aggrunin[y, x])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Fortran 90 cshift function
    def __cshift(l, offset):
        """
        see F90  cshift with offset=-offset
        """
        offset %= len(l)
        return np.concatenate((l[offset:], l[:offset]))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Extract the current rof
    def get_rof(self):
        return self.ring[0, :]     # Current timestep flux (units=kg m-2 s-1)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Extract the current storage
    def get_storage(self):
        return self.ring.sum(axis=1)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write the current state
    def write_state(self):

        # ------------------------------------------------------------ #
        # Open file
        filename = self.timestamp.strftime(self.fname_format)
        f = Dataset(filename, 'w', self.file_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Time Variables

        # Current time
        time = f.createDimension('time', 1)
        time = f.createVariable('time', NC_DOUBLE, ('time',))
        time[:] = date2num(self.timestamp, TIMEUNITS, calendar=self.calendar)

        for key, val in share.time.__dict__.iteritems():
            if val:
                setattr(time, key, val)

        # Timesteps
        timesteps = f.createDimension('timesteps', self.full_time_length)
        timesteps = f.createVariable('timesteps', NC_DOUBLE, ('timesteps',))
        timesteps[:] = np.arange(self.full_time_length)

        for key, val in share.timesteps.__dict__.iteritems():
            if val:
                setattr(timesteps, key, val)
        timesteps.timestep_length = 'timestep'

        # UH timestep
        timestep = f.createVariable('timestep', NC_DOUBLE, ())
        timestep[:] = self.unit_hydrogaph_dt
        for key, val in share.timestep.__dict__.iteritems():
            if val:
                setattr(timestep, key, val)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup Coordinate Variables
        coords = ('outlets',)

        outlets = f.createDimension('outlets', self.num_outlets)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        lon_outlet = f.createVariable('lon_outlet', NC_DOUBLE, coords)
        lon_outlet[:] = self.lon_outlet
        for key, val in share.lon_outlet.__dict__.iteritems():
            if val:
                setattr(lon_outlet, key, val)

        lat_outlet = f.createVariable('lat_outlet', NC_DOUBLE, coords)
        lat_outlet[:] = self.lat_outlet
        for key, val in share.lat_outlet.__dict__.iteritems():
            if val:
                setattr(lat_outlet, key, val)

        y_ind_outlet = f.createVariable('y_ind_outlet', NC_INT, coords)
        y_ind_outlet[:] = self.y_ind_outlet
        for key, val in share.y_ind_outlet.__dict__.iteritems():
            if val:
                setattr(y_ind_outlet, key, val)

        x_ind_outlet = f.createVariable('x_ind_outlet', NC_INT, coords)
        x_ind_outlet[:] = self.x_ind_outlet
        for key, val in share.x_ind_outlet.__dict__.iteritems():
            if val:
                setattr(x_ind_outlet, key, val)

        s_outlet_decomp_id = f.createVariable('s_outlet_decomp_id', NC_INT, coords)
        s_outlet_decomp_id[:] = self.outlet_decomp_id
        for key, val in share.outlet_decomp_id.__dict__.iteritems():
            if val:
                setattr(s_outlet_decomp_id, key, val)

        tcoords = ('timesteps',) + coords

        ring = f.createVariable('ring', NC_DOUBLE, tcoords)
        ring[:, :] = self.ring

        for key, val in share.ring.__dict__.iteritems():
            if val:
                setattr(ring, key, val)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global attributes
        if self.GlobAts:
            self.GlobAts.update()
        else:
            self.GlobAts = NcGlobals(title='RVIC restart file',
                                     RvicPourPointsFile=self.RvicPourPointsFile,
                                     RvicUHFile=self.RvicUHFile,
                                     RvicFdrFile=self.RvicFdrFile,
                                     RvicDomainFile=self.RvicDomainFile)
        for key, val in self.GlobAts.__dict__.iteritems():
            if val:
                setattr(f, key, val)
        # ------------------------------------------------------------ #

        f.close()
        log.info('Finished writing %s' % filename)

        return filename
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
