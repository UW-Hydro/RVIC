"""
variables.py
"""
import os
import numpy as np
from netCDF4 import Dataset, date2num
from logging import getLogger
from log import LOG_NAME
from time_utility import ord_to_datetime
from share import TIMEUNITS, NC_INT, NC_DOUBLE, NC_CHAR, REFERENCE_DATE, REFERENCE_TIME
from share import CALENDAR_KEYS, RVIC_TRACERS, NcGlobals, SECSPERDAY, MAX_NC_CHARS
import share

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Point object
class Point(object):
    '''Creates a point class for intellegently storing coordinate information'''

    def __init__(self, lat='', lon='', gridx='', gridy='', routx='',
                 routy='', name=None, cell_id=None):
        '''Defines x and y variables'''
        self.lat = lat
        self.lon = lon
        self.gridx = gridx
        self.gridy = gridy
        self.routx = routx
        self.routy = routy
        self.cell_id = cell_id

        if name:
            self.name = name
        else:
            self.name = 'outlet_%s_%s' %(self.lat, self.lon)

    def __str__(self):
        return "Point(%s,%s,%s,%s)" % (self.lat, self.lon, self.gridy, self.gridx)

    def __repr__(self):
        return '__repr__'

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Rvar Object
class Rvar(object):
    """ Creates a RVIC structure """

    # ---------------------------------------------------------------- #
    # Initialize
    def __init__(self, param_file, case_name, calendar, out_dir, file_format):
        self.param_file = param_file
        f = Dataset(param_file, 'r+')
        self.n_sources = len(f.dimensions['sources'])
        self.n_outlets = len(f.dimensions['outlets'])
        self.subset_length = f.variables['subset_length'][0]
        self.full_time_length = f.variables['full_time_length'][0]
        self.unit_hydrograph_dt = f.variables['unit_hydrograph_dt'][0]
        self.source_lon = f.variables['source_lon'][:]
        self.source_lat = f.variables['source_lat'][:]
        self.source_x_ind = f.variables['source_x_ind'][:]
        self.source_y_ind = f.variables['source_y_ind'][:]
        self.source_time_offset = f.variables['source_time_offset'][:]
        self.source2outlet_ind = f.variables['source2outlet_ind'][:]
        self.outlet_x_ind = f.variables['outlet_x_ind'][:]
        self.outlet_y_ind = f.variables['outlet_y_ind'][:]
        self.outlet_lon = f.variables['outlet_lon'][:]
        self.outlet_lat = f.variables['outlet_lat'][:]
        self.outlet_mask = f.variables['outlet_mask'][:]
        self.outlet_decomp_ind = f.variables['outlet_decomp_ind'][:]
        self.unit_hydrograph = f.variables['unit_hydrograph'][:]
        self.RvicDomainFile = f.RvicDomainFile
        self.RvicPourPointsFile = f.RvicPourPointsFile
        self.RvicUHFile = f.RvicUHFile
        self.RvicFdrFile = f.RvicFdrFile
        self.file_format = file_format
        self.glob_atts = NcGlobals(title='RVIC restart file',
                                   RvicPourPointsFile=f.RvicPourPointsFile,
                                   RvicUHFile=f.RvicUHFile,
                                   RvicFdrFile=f.RvicFdrFile,
                                   RvicDomainFile=f.RvicDomainFile,
                                   casename=case_name)
        f.close()

        # ------------------------------------------------------------ #
        # Initialize state variables
        self.ring = np.zeros((self.full_time_length, self.n_outlets, len(RVIC_TRACERS)))
        # ------------------------------------------------------------ #


        self._calendar = calendar
        self.__fname_format = os.path.join(out_dir, "%s.r.%%Y-%%m-%%d-%%H-%%M-%%S.nc" % (case_name))

        # CESM calendar key (only NO_LEAP_C, GREGORIAN are supported in CESM)
        self._calendar_key = 0
        for key, cals in CALENDAR_KEYS.iteritems():
            if self._calendar in cals:
                self._calendar_key = key
                break
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Check that grid file matches
    def check_grid_file(self, domain_file):
        """Confirm that the grid files match in the parameter and domain files"""
        input_file = os.path.split(domain_file)[1]
        log.info('domain_file: %s' % input_file)
        log.info('Parameter RvicDomainFile: %s' % self.RvicDomainFile)

        if input_file == self.RvicDomainFile:
            log.info('Grid files match in parameter and domain file')
        else:
            raise  ValueError('Grid files do not match in parameter and domain file')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initilize State
    def init_state(self, state_file, run_type, timestamp):
        if run_type in ['startup', 'restart']:
            log.info('reading state_file: %s' %state_file)
            f = Dataset(state_file, 'r+')
            self.ring = f.variables['ring'][:]
            file_timestamp = ord_to_datetime(f.variables['time'][:], f.variables['time'].units, calendar=f.variables['time'].calendar)

            if run_type == 'restart':
                self.timestamp = file_timestamp

            elif run_type == 'startup':
                self.timestamp = timestamp
                if timestamp != file_timestamp:
                    log.warning('restart timestamps do not match (%s, %s', file_timestamp, self.timestamp)
                    log.warning('Runtype is startup so model will continue')
            else:
                raise ValueError('unknown run_type: %s' %run_type)

            # Check that timestep and outlet_decomp_ids match ParamFile
            if f.variables['unit_hydrograph_dt'][:] != self.unit_hydrograph_dt:
                raise ValueError('Timestep in Statefile does not match timestep in ParamFile')

            if not np.array_equal(f.variables['outlet_decomp_ind'][:], self.outlet_decomp_ind):
                raise ValueError('outlet_decomp_ind in Statefile does not match ParamFile')

            if f.RvicDomainFile != self.RvicDomainFile:
                raise ValueError('RvicDomainFile in StateFile does not match ParamFile')

            f.close()

        elif run_type == 'drystart':
            log.info('run_type is drystart so no state_file will be read')
            self.timestamp = timestamp

        self.time_ord = date2num(self.timestamp, TIMEUNITS, calendar=self._calendar)

        self._start_date = self.timestamp
        self._start_ord = self.time_ord
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Convolve
    def convolve(self, aggrunin, time_ord):
        """
        This convoluition funciton works by looping over all points and doing the
        convolution one timestep at a time.  This is accomplished by creating an
        convolution ring.  Contributing flow from each timestep is added to the
        convolution ring.  The convolution ring is saved as the state.  The first
        column of values in the ring are the current runoff.
        """
        # ------------------------------------------------------------ #
        # Check that the time_ord is in sync
        # This is the time at the start of the current step (end of last step)
        if self.time_ord != time_ord:
            log.error('rout_var.time_ord = %s, time_ord = %s' %(self.time_ord, time_ord))
            raise ValueError('rout_var.time_ord does not match the time_ord passed in by the convolution call')

        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # First update the ring
        log.debug('rolling the ring')
        self.ring[0, :, 0] = 0                      # Zero out current ring
        self.ring = np.roll(self.ring, 1, axis=0)   # Equivalent to Fortran 90 cshift function
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Do the convolution
        log.debug('convolving')
        # this matches the fortran implementation, it may be faster to use np.convolve but testing
        # can be done later
        # also this is where the parallelization will happen
        for nt, tracer in enumerate(RVIC_TRACERS):
            for s, outlet in enumerate(self.source2outlet_ind):   # loop over all source points
                y = self.source_y_ind[s]
                x = self.source_x_ind[s]
                for i in xrange(self.subset_length):
                    j = i + self.source_time_offset[s]
                    self.ring[j, outlet, nt] = self.ring[j, outlet, nt] + (self.unit_hydrograph[i, s, nt] * aggrunin[tracer][y, x])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # move the time_ord forward
        self.time_ord += self.unit_hydrograph_dt / SECSPERDAY
        self.timestamp = ord_to_datetime(self.time_ord, TIMEUNITS, calendar=self._calendar)

        return self.timestamp
        # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Extract the current rof
    def get_rof(self):
        return self.ring[0, :, 0]     # Current timestep flux (units=kg m-2 s-1)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Extract the current storage
    def get_storage(self):
        return self.ring.sum(axis=1)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write the current state
    def write_restart(self, current_history_files, history_restart_files):

        # ------------------------------------------------------------ #
        # Open file
        filename = self.timestamp.strftime(self.__fname_format)
        f = Dataset(filename, 'w', self.file_format)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Time Variables

        # Current time
        time = f.createDimension('time', 1)
        time = f.createVariable('time', NC_DOUBLE, ('time',))
        time[:] = date2num(self.timestamp, TIMEUNITS, calendar=self._calendar)

        for key, val in share.time.__dict__.iteritems():
            if val:
                setattr(time, key, val)
        time.calendar = self._calendar

        # Timesteps
        timesteps = f.createDimension('timesteps', self.full_time_length)
        timesteps = f.createVariable('timesteps', NC_DOUBLE, ('timesteps',))
        timesteps[:] = np.arange(self.full_time_length)

        for key, val in share.timesteps.__dict__.iteritems():
            if val:
                setattr(timesteps, key, val)
        timesteps.timestep_length = 'unit_hydrograph_dt'

        # UH timestep
        unit_hydrograph_dt = f.createVariable('unit_hydrograph_dt', NC_DOUBLE, ())
        unit_hydrograph_dt[:] = self.unit_hydrograph_dt
        for key, val in share.unit_hydrograph_dt.__dict__.iteritems():
            if val:
                setattr(unit_hydrograph_dt, key, val)

        # timemgr_rst_type = f.createVariable('timemgr_rst_type', NC_DOUBLE, ())
        # timemgr_rst_type[:] = self._calendar_key
        # for key, val in share.timemgr_rst_type.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_type, key, val)

        # timemgr_rst_step_sec = f.createVariable('timemgr_rst_step_sec', NC_DOUBLE, ())
        # timemgr_rst_step_sec[:] = unit_hydrograph_dt
        # for key, val in share.timemgr_rst_step_sec.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_step_sec, key, val)

        # timemgr_rst_start_ymd = f.createVariable('timemgr_rst_start_ymd', NC_DOUBLE, ())
        # timemgr_rst_start_ymd[:] = self._start_date.year*10000+self._start_date.month*100+self._start_date.day
        # for key, val in share.timemgr_rst_start_ymd.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_start_ymd, key, val)

        # timemgr_rst_start_tod = f.createVariable('timemgr_rst_start_tod', NC_DOUBLE, ())
        # timemgr_rst_start_tod[:] = (self._start_ord%1)*SECSPERDAY
        # for key, val in share.timemgr_rst_start_tod.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_start_tod, key, val)

        # timemgr_rst_ref_ymd = f.createVariable('timemgr_rst_ref_ymd', NC_DOUBLE, ())
        # timemgr_rst_ref_ymd[:] = REFERENCE_DATE
        # for key, val in share.timemgr_rst_ref_ymd.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_ref_ymd, key, val)

        # timemgr_rst_ref_tod = f.createVariable('timemgr_rst_ref_tod', NC_DOUBLE, ())
        # timemgr_rst_ref_tod[:] = REFERENCE_TIME
        # for key, val in share.timemgr_rst_ref_tod.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_ref_tod, key, val)

        # timemgr_rst_curr_ymd = f.createVariable('timemgr_rst_curr_ymd', NC_DOUBLE, ())
        # timemgr_rst_curr_ymd[:] = self.timestamp.year*10000+self.timestamp.month*100+self.timestamp.day
        # for key, val in share.timemgr_rst_curr_ymd.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_curr_ymd, key, val)

        # timemgr_rst_curr_tod = f.createVariable('timemgr_rst_curr_tod', NC_DOUBLE, ())
        # timemgr_rst_curr_tod[:] = (self.time_ord%1)*SECSPERDAY
        # for key, val in share.timemgr_rst_curr_tod.__dict__.iteritems():
        #     if val:
        #         setattr(timemgr_rst_curr_tod, key, val)


        # ------------------------------------------------------------ #
        # Setup Tape Dimensions
        # coords = ('tapes', 'mak_chars')
        # ntapes = f.createDimension(coords[0], len(history_restart_files))
        # ntapes = f.createDimension(coords[1], MAX_NC_CHARS)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        # locfnh = f.createVariable('locfnh', NC_CHAR, ('ntapes', 'max_chars',))
        # locfnh[:] = current_history_files
        # locfnh.long_name = 'History filename'
        # locfnh.comment = 'This variable NOT needed for startup or branch simulations'

        # locfnhr = f.createVariable('locfnhr', NC_CHAR, ('ntapes', 'max_chars',))
        # locfnhr[:] = history_restart_files
        # locfnhr.long_name = 'Restart history filename'
        # locfnhr.comment = 'This variable NOT needed for startup or branch simulations'

        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup Point Dimensions
        coords = ('outlets', 'tracers')

        outlets = f.createDimension(coords[0], self.n_outlets)
        tracers = f.createDimension(coords[1], len(RVIC_TRACERS))
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        oyi = f.createVariable('outlet_y_ind', NC_INT, coords[0])
        oyi[:] = self.outlet_y_ind
        for key, val in share.outlet_y_ind.__dict__.iteritems():
            if val:
                setattr(oyi, key, val)

        oxi = f.createVariable('outlet_x_ind', NC_INT, coords[0])
        oxi[:] = self.outlet_x_ind
        for key, val in share.outlet_x_ind.__dict__.iteritems():
            if val:
                setattr(oxi, key, val)

        odi = f.createVariable('outlet_decomp_ind', NC_INT, coords[0])
        odi[:] = self.outlet_decomp_ind
        for key, val in share.outlet_decomp_ind.__dict__.iteritems():
            if val:
                setattr(odi, key, val)

        tcoords = ('timesteps',) + coords

        ring = f.createVariable('ring', NC_DOUBLE, tcoords)
        ring[:, :, :] = self.ring

        for key, val in share.ring.__dict__.iteritems():
            if val:
                setattr(ring, key, val)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global attributes
        self.glob_atts.update()

        for key, val in self.glob_atts.__dict__.iteritems():
            if val:
                setattr(f, key, val)
        # ------------------------------------------------------------ #

        f.close()
        log.info('Finished writing %s' % filename)

        return filename
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
