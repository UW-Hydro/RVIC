# -*- coding: utf-8 -*-
'''
variables.py
'''
import os
import numpy as np
from netCDF4 import Dataset, date2num, stringtochar
from logging import getLogger
from .log import LOG_NAME
from .time_utility import ord_to_datetime
from .share import TIMEUNITS, NC_INT, NC_DOUBLE, NC_CHAR
from .share import RVIC_TRACERS, NcGlobals, SECSPERDAY, MAX_NC_CHARS
from .share import CALENDAR_KEYS, REFERENCE_DATE, REFERENCE_TIME
from .convolution_wrapper import rvic_convolve
from .pycompat import iteritems
from . import share

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Point object
class Point(object):
    '''
    Creates a point class for intellegently storing coordinate information
    '''

    def __init__(self, lat=-9999.9, lon=-9999.9, domx=-9999, domy=-9999,
                 routx=-9999, routy=-9999, name=None, cell_id=-9999):
        '''Defines x and y variables'''
        self.lat = lat
        self.lon = lon
        self.domx = domx
        self.domy = domy
        self.routx = routx
        self.routy = routy
        self.cell_id = cell_id

        if name:
            self.name = name
        else:
            self.name = 'outlet_{0:3.4f}_{1:3.4f}'.format(self.lat, self.lon)

        return

    def __str__(self):
        return ('Point({0}, lat:{1:3.4f}, lon:{2:3.4f}, y:{3:d}, '
                'x:{4:d})'.format(self.name, self.lat, self.lon, self.domy,
                                  self.domx))

    def __repr__(self):
        return ('    -- Point --    \n'
                'name:\t{0}\n'
                'lat:\t{1:3.4f}\n'
                'lon:\t{2:3.4f}\n'
                'domy:\t{3:d}\n'
                'domx:\t{4:d}\n'
                'routy:\t{5:d}\n'
                'routx:\t{6:d}\n'.format(self.name, self.lat, self.lon,
                                         self.domy, self.domx,
                                         self.routy, self.routx))
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Rvar Object
class Rvar(object):
    ''' Creates a RVIC structure '''

    # ---------------------------------------------------------------- #
    # Initialize
    def __init__(self, param_file, case_name, calendar, out_dir, file_format,
                 zlib=True, complevel=4, least_significant_digit=None):
        self.param_file = param_file
        f = Dataset(param_file, 'r')
        self.n_sources = len(f.dimensions['sources'])
        self.n_outlets = len(f.dimensions['outlets'])
        self.subset_length = f.variables['subset_length'][:]
        self.full_time_length = f.variables['full_time_length'][:]
        self.unit_hydrograph_dt = f.variables['unit_hydrograph_dt'][:]
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
        self.unit_hydrograph = {}
        for tracer in RVIC_TRACERS:
            tname = 'unit_hydrograph_{0}'.format(tracer)
            try:
                self.unit_hydrograph[tracer] = f.variables[tname][:]
            except KeyError:
                log.warning('Could not find unit hydrograph var %s', tname)
                log.warning('trying var name unit_hydrograph')
                self.unit_hydrograph[tracer] = \
                    f.variables['unit_hydrograph'][:]
            except:
                raise ValueError('Cant find Unit Hydrograph Variable')
        self.outlet_name = f.variables['outlet_name'][:]
        self.RvicDomainFile = f.RvicDomainFile
        self.RvicPourPointsFile = f.RvicPourPointsFile
        self.RvicUHFile = f.RvicUHFile
        self.RvicFdrFile = f.RvicFdrFile
        self.file_format = file_format
        try:
            self.outlet_upstream_area = f.variables['outlet_upstream_area'][:]
        except KeyError:
            self.outlet_upstream_area = None
        self.glob_atts = NcGlobals(title='RVIC restart file',
                                   RvicPourPointsFile=f.RvicPourPointsFile,
                                   RvicUHFile=f.RvicUHFile,
                                   RvicFdrFile=f.RvicFdrFile,
                                   RvicDomainFile=f.RvicDomainFile,
                                   casename=case_name)

        f.close()

        # ------------------------------------------------------------ #
        # Initialize state variables
        self.ring = {}
        for tracer in RVIC_TRACERS:

            self.ring[tracer] = np.zeros((self.full_time_length,
                                          self.n_outlets,),
                                         dtype=np.float64)
        # ------------------------------------------------------------ #

        self._calendar = calendar
        self.__fname_format = os.path.join(
            out_dir, '%s.r.%%Y-%%m-%%d-%%H-%%M-%%S.nc' % (case_name))

        # ------------------------------------------------------------ #
        # CESM calendar key (only NO_LEAP_C, GREGORIAN are supported in CESM)
        self._calendar_key = 0
        for key, cals in iteritems(CALENDAR_KEYS):
            if self._calendar in cals:
                self._calendar_key = key
                break
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # netCDF variable options
        self.ncvaropts = {'zlib': zlib,
                          'complevel': complevel,
                          'least_significant_digit': least_significant_digit}
        # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Check that dom file matches
    def _check_domain_file(self, domain_file):
        '''
        Confirm that the dom files match in the parameter and domain files
        '''
        input_file = os.path.split(domain_file)[1]
        log.info('domain_file: %s', input_file)
        log.info('Parameter RvicDomainFile: %s', self.RvicDomainFile)

        if input_file == self.RvicDomainFile:
            log.info('dom files match in parameter and domain file')
        else:
            raise ValueError('dom files do not match in parameter and '
                             'domain file')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def set_domain(self, dom_data, domain, lat0_is_min):
        ''' Set the domain size '''
        self._check_domain_file(domain['FILE_NAME'])

        self.domain_shape = dom_data[domain['LAND_MASK_VAR']].shape

        self.ysize = self.domain_shape[0]
        self.xsize = self.domain_shape[1]

        if self.source_y_ind.max() >= self.ysize:
            raise ValueError('source_y_ind.max() ({0}) > domain ysize'
                             ' ({1})'.format(self.source_y_ind, self.ysize))
        if self.source_x_ind.max() >= self.xsize:
            raise ValueError('source_x_ind.max() ({0}) > domain xsize'
                             ' ({1})'.format(self.source_x_ind, self.xsize))
        log.info('set domain')

        if lat0_is_min:
            log.info('Flipping Parameter File Y inds...')
            self._flip_y_inds()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Flip the y index order
    def _flip_y_inds(self):
        '''
        Flip the y index order
        '''
        self.source_y_ind = self.ysize - self.source_y_ind - 1
        self.outlet_y_ind = self.ysize - self.outlet_y_ind - 1
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initilize State
    def init_state(self, state_file, run_type, timestamp):
        if run_type in ['startup', 'restart']:
            log.info('reading state_file: %s', state_file)
            f = Dataset(state_file, 'r')
            for tracer in RVIC_TRACERS:
                self.ring[tracer] = f.variables['{0}_ring'.format(tracer)][:]

            file_timestamp = ord_to_datetime(
                f.variables['time'][:],
                f.variables['time'].units,
                calendar=f.variables['time'].calendar)

            if run_type == 'restart':
                self.timestamp = file_timestamp

            elif run_type == 'startup':
                self.timestamp = timestamp
                if timestamp != file_timestamp:
                    log.warning('restart timestamps do not match (%s, %s',
                                file_timestamp, self.timestamp)
                    log.warning('Runtype is startup so model will continue')
            else:
                raise ValueError('unknown run_type: {0}'.format(run_type))

            # Check that timestep and outlet_decomp_ids match ParamFile
            if f.variables['unit_hydrograph_dt'][:] != self.unit_hydrograph_dt:
                raise ValueError('Timestep in Statefile does not match '
                                 'timestep in ParamFile')

            if not np.array_equal(f.variables['outlet_decomp_ind'][:],
                                  self.outlet_decomp_ind):
                raise ValueError('outlet_decomp_ind in Statefile does not '
                                 'match ParamFile')

            if f.RvicDomainFile != self.RvicDomainFile:
                raise ValueError('RvicDomainFile in StateFile does not match '
                                 'ParamFile')

            f.close()

        elif run_type == 'drystart':
            log.info('run_type is drystart so no state_file will be read')
            self.timestamp = timestamp

        self.time_ord = date2num(self.timestamp, TIMEUNITS,
                                 calendar=self._calendar)

        self._start_date = self.timestamp
        self._start_ord = self.time_ord
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Convolve
    def convolve(self, aggrunin, time_ord):
        '''
        This convoluition funciton works by looping over all points and doing
        the convolution one timestep at a time.  This is accomplished by
        creating a convolution ring.  Contributing flow from each timestep is
        added to the convolution ring.  The convolution ring is saved as the
        state.  The first column of values in the ring are the current runoff.
        '''
        # ------------------------------------------------------------ #
        # Check that the time_ord is in sync
        # This is the time at the start of the current step (end of last step)
        if self.time_ord != time_ord:
            log.error('rout_var.time_ord = %s, time_ord = %s', self.time_ord,
                      time_ord)
            raise ValueError('rout_var.time_ord does not match the time_ord '
                             'passed in by the convolution call')

        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Do the convolution
        log.debug('convolving')

        for tracer in RVIC_TRACERS:
            # -------------------------------------------------------- #
            # First update the ring
            log.debug('rolling the ring')

            # Zero out current ring
            self.ring[tracer][0, :] = 0.

            # Equivalent to Fortran 90 cshift function
            self.ring[tracer] = np.roll(self.ring[tracer], -1, axis=0)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # C convolution call
            rvic_convolve(self.n_sources,
                          self.n_outlets,
                          self.subset_length,
                          self.xsize,
                          self.source2outlet_ind,
                          self.source_y_ind,
                          self.source_x_ind,
                          self.source_time_offset,
                          self.unit_hydrograph[tracer],
                          aggrunin[tracer],
                          self.ring[tracer])
            # -------------------------------------------------------- #
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # move the time_ord forward
        self.time_ord += self.unit_hydrograph_dt / SECSPERDAY
        self.timestamp = ord_to_datetime(self.time_ord, TIMEUNITS,
                                         calendar=self._calendar)

        return self.timestamp
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def get_time_mode(self, cpl_secs_per_step):
        '''
        Determine the relationship between the coupling period and the unit-
        hydrograph period.  In cases where they do not match, the model will
        aggregate the appropriate quantities before/after the confolution step.
        '''

        log.info('Coupling Timestep is (seconds): %s', cpl_secs_per_step)
        log.info('RVIC Timestep is (seconds): %s', self.unit_hydrograph_dt)

        if (self.unit_hydrograph_dt % cpl_secs_per_step == 0) and \
           (self.unit_hydrograph_dt >= cpl_secs_per_step):
            self.agg_tsteps = self.unit_hydrograph_dt / cpl_secs_per_step
        else:
            log.error('unit_hydrograph_dt must be a multiple of the '
                      'cpl_secs_per_step')
            raise ValueError('Stopped due to error in determining agg_tsteps')

        log.info('RVIC will run 1 time for every %i coupling periods',
                 self.agg_tsteps)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def get_rof(self):
        '''Extract the current rof'''
        rof = {}
        for tracer in RVIC_TRACERS:
            rof[tracer] = self.ring[tracer][0, :]
        return rof  # Current timestep flux (units=kg m-2 s-1)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def get_storage(self):
        '''Extract the current storage'''
        storage = {}
        for tracer in RVIC_TRACERS:
            storage[tracer] = self.ring[tracer].sum(axis=1)
        return storage
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    def write_initial(self):
        '''write initial flux'''
        pass
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
        time = f.createVariable('time', NC_DOUBLE, ('time',), **self.ncvaropts)
        time[:] = date2num(self.timestamp, TIMEUNITS, calendar=self._calendar)

        for key, val in iteritems(share.time):
            if val:
                setattr(time, key, val)
        time.calendar = self._calendar

        # Timesteps
        timesteps = f.createDimension('timesteps', self.full_time_length)
        timesteps = f.createVariable('timesteps', NC_DOUBLE, ('timesteps',),
                                     **self.ncvaropts)
        timesteps[:] = np.arange(self.full_time_length)

        for key, val in iteritems(share.timesteps):
            if val:
                setattr(timesteps, key, val)
        timesteps.timestep_length = 'unit_hydrograph_dt'

        # UH timestep
        unit_hydrograph_dt = f.createVariable('unit_hydrograph_dt',
                                              NC_DOUBLE, (), **self.ncvaropts)
        unit_hydrograph_dt[:] = self.unit_hydrograph_dt
        for key, val in iteritems(share.unit_hydrograph_dt):
            if val:
                setattr(unit_hydrograph_dt, key, val)

        timemgr_rst_type = f.createVariable('timemgr_rst_type', NC_DOUBLE, (),
                                            **self.ncvaropts)
        timemgr_rst_type[:] = self._calendar_key
        for key, val in iteritems(share.timemgr_rst_type):
            if val:
                setattr(timemgr_rst_type, key, val)

        timemgr_rst_step_sec = f.createVariable('timemgr_rst_step_sec',
                                                NC_DOUBLE, (),
                                                **self.ncvaropts)
        timemgr_rst_step_sec[:] = self.unit_hydrograph_dt
        for key, val in iteritems(share.timemgr_rst_step_sec):
            if val:
                setattr(timemgr_rst_step_sec, key, val)

        timemgr_rst_start_ymd = f.createVariable('timemgr_rst_start_ymd',
                                                 NC_DOUBLE, (),
                                                 **self.ncvaropts)
        timemgr_rst_start_ymd[:] = self._start_date.year * 10000 \
            + self._start_date.month * 100 + self._start_date.day
        for key, val in iteritems(share.timemgr_rst_start_ymd):
            if val:
                setattr(timemgr_rst_start_ymd, key, val)

        timemgr_rst_start_tod = f.createVariable('timemgr_rst_start_tod',
                                                 NC_DOUBLE, (),
                                                 **self.ncvaropts)
        timemgr_rst_start_tod[:] = (self._start_ord % 1) * SECSPERDAY
        for key, val in iteritems(share.timemgr_rst_start_tod):
            if val:
                setattr(timemgr_rst_start_tod, key, val)

        timemgr_rst_ref_ymd = f.createVariable('timemgr_rst_ref_ymd',
                                               NC_DOUBLE, (),
                                               **self.ncvaropts)
        timemgr_rst_ref_ymd[:] = REFERENCE_DATE
        for key, val in iteritems(share.timemgr_rst_ref_ymd):
            if val:
                setattr(timemgr_rst_ref_ymd, key, val)

        timemgr_rst_ref_tod = f.createVariable('timemgr_rst_ref_tod',
                                               NC_DOUBLE, (),
                                               **self.ncvaropts)
        timemgr_rst_ref_tod[:] = REFERENCE_TIME
        for key, val in iteritems(share.timemgr_rst_ref_tod):
            if val:
                setattr(timemgr_rst_ref_tod, key, val)

        timemgr_rst_curr_ymd = f.createVariable('timemgr_rst_curr_ymd',
                                                NC_DOUBLE, (),
                                                **self.ncvaropts)
        timemgr_rst_curr_ymd[:] = self.timestamp.year * 10000 + \
            self.timestamp.month * 100 + self.timestamp.day
        for key, val in iteritems(share.timemgr_rst_curr_ymd):
            if val:
                setattr(timemgr_rst_curr_ymd, key, val)

        timemgr_rst_curr_tod = f.createVariable('timemgr_rst_curr_tod',
                                                NC_DOUBLE, (),
                                                **self.ncvaropts)
        timemgr_rst_curr_tod[:] = (self.time_ord % 1) * SECSPERDAY
        for key, val in iteritems(share.timemgr_rst_curr_tod):
            if val:
                setattr(timemgr_rst_curr_tod, key, val)

        # ------------------------------------------------------------ #
        # Setup Tape Dimensions
        coords = ('tapes', 'max_chars')
        f.createDimension(coords[0], len(history_restart_files))
        f.createDimension(coords[1], MAX_NC_CHARS)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        locfnh = f.createVariable('locfnh', NC_CHAR, coords, **self.ncvaropts)
        print('locfnh.shape', locfnh.shape)
        for i, string in enumerate(current_history_files):
            b_string = string.encode()
            assert type(b_string) == bytes
            locfnh[i, :] = stringtochar(np.array(b_string.ljust(MAX_NC_CHARS)))
        locfnh.long_name = 'History filename'
        locfnh.comment = 'This variable is NOT needed for startup or branch '\
                         'simulations'

        locfnhr = f.createVariable('locfnhr', NC_CHAR, coords,
                                   **self.ncvaropts)
        for i, string in enumerate(history_restart_files):
            b_string = string.encode()
            assert type(b_string) == bytes
            locfnh[i, :] = stringtochar(np.array(b_string.ljust(MAX_NC_CHARS)))
        locfnhr.long_name = 'History restart filename'
        locfnhr.comment = 'This variable NOT needed for startup or branch '\
                          'simulations'
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Setup Point Dimensions
        coords = ('outlets', )

        f.createDimension(coords[0], self.n_outlets)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write Fields
        oyi = f.createVariable('outlet_y_ind', NC_INT, coords[0],
                               **self.ncvaropts)
        oyi[:] = self.outlet_y_ind
        for key, val in iteritems(share.outlet_y_ind):
            if val:
                setattr(oyi, key, val)

        oxi = f.createVariable('outlet_x_ind', NC_INT, coords[0],
                               **self.ncvaropts)
        oxi[:] = self.outlet_x_ind
        for key, val in iteritems(share.outlet_x_ind):
            if val:
                setattr(oxi, key, val)

        odi = f.createVariable('outlet_decomp_ind', NC_INT, coords[0],
                               **self.ncvaropts)
        odi[:] = self.outlet_decomp_ind
        for key, val in iteritems(share.outlet_decomp_ind):
            if val:
                setattr(odi, key, val)

        tcoords = ('timesteps', ) + coords

        for tracer in RVIC_TRACERS:
            ring = f.createVariable('{0}_ring'.format(tracer),
                                    NC_DOUBLE, tcoords, **self.ncvaropts)
            ring[:, :] = self.ring[tracer][:, :]

            for key, val in iteritems(share.ring):
                if val:
                    setattr(ring, key, val)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # write global attributes
        self.glob_atts.update()

        for key, val in iteritems(self.glob_atts.atts):
            if val:
                setattr(f, key, val)
        # ------------------------------------------------------------ #

        f.close()
        log.info('Finished writing %s', filename)

        return filename
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
