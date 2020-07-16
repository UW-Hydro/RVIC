# -*- coding: utf-8 -*-
'''
share.py
'''
import sys
import socket
import string
import time as time_mod
from .pycompat import OrderedDict, iteritems
from netCDF4 import default_fillvals
from getpass import getuser


# ----------------------- CONSTANTS --------------------------------- #
EARTHRADIUS = 6.37122e6  # meters
WATERDENSITY = 1000.

# area
METERSPERKM = 1000.
METERSPERMILE = 1609.34
METERS2PERACRE = 4046.856

# time
# reference time
REFERENCE_STRING = '0001-1-1 0:0:0'
REFERENCE_DATE = 10101  # i.e. REFERENCE_STRING
REFERENCE_TIME = 0  # i.e. REFERENCE_STRING
TIMEUNITS = 'days since ' + REFERENCE_STRING  # do not change (MUST BE DAYS)!
TIMESTAMPFORM = '%Y-%m-%d-%H'
CALENDAR = 'noleap'
HOURSPERDAY = 24.
SECSPERHOUR = 3600.
MINSPERHOUR = 60.
MINSPERDAY = HOURSPERDAY * MINSPERHOUR
SECSPERDAY = HOURSPERDAY * SECSPERHOUR
MONTHSPERYEAR = 12

# length
MMPERMETER = 1000.
CMPERMETER = 100.

# precision
PRECISION = 1.0e-30
NC_DOUBLE = 'f8'
NC_FLOAT = 'f4'
NC_INT = 'i4'
NC_CHAR = 'S1'
MAX_NC_CHARS = 256

# fill values
FILLVALUE_F = default_fillvals[NC_DOUBLE]
FILLVALUE_I = default_fillvals[NC_INT]

# filenames
RPOINTER = 'rpointer'

# tracers
RVIC_TRACERS = ('LIQ',)  # Before changing, update history module

# Calendar key number for linking with CESM
CALENDAR_KEYS = {0: ['None'],
                 1: ['noleap', '365_day'],
                 2: ['gregorian', 'standard'],
                 3: ['proleptic_gregorian'],
                 4: ['all_leap', '366_day'],
                 5: ['360_day'],
                 6: ['julian']}

VALID_CHARS = '-_. %s%s' % (string.ascii_letters, string.digits)


# ----------------------- NETCDF VARIABLES --------------------------------- #
class NcGlobals(object):
    def __init__(self,
                 title=None,
                 casename=None,
                 casestr=None,
                 history='Created: {}'.format(time_mod.ctime(time_mod.time())),
                 institution='University of Washington',
                 source=sys.argv[0],
                 references='Based on the initial model of Lohmann, et al., '
                            '1996, Tellus, 48(A), 708-721',
                 comment='Output from the RVIC Streamflow Routing Model.',
                 Conventions='CF-1.6',
                 RvicPourPointsFile=None,
                 RvicFdrFile=None,
                 RvicUHFile=None,
                 RvicDomainFile=None,
                 version=None,
                 hostname=None,
                 username=None):

        self.atts = OrderedDict()

        if title is not None:
            self.atts['title'] = title

        if comment is not None:
            self.atts['comment'] = comment

        if Conventions is not None:
            self.atts['Conventions'] = Conventions

        if history is not None:
            self.atts['history'] = history

        if source is not None:
            self.atts['source'] = source

        if institution is not None:
            self.atts['institution'] = institution

        if hostname is not None:
            self.atts['hostname'] = hostname
        else:
            self.atts['hostname'] = socket.gethostname()

        if username is not None:
            self.atts['username'] = username
        else:
            self.atts['username'] = getuser()

        if casename is not None:
            self.atts['casename'] = casename

        if casestr is not None:
            self.atts['casestr'] = casestr

        if references is not None:
            self.atts['references'] = references

        if version is not None:
            self.atts['version'] = version
        else:
            try:
                from rvic import version
                self.atts['version'] = version.short_version
            except ImportError:
                self.atts['version'] = 'unknown'

        if RvicPourPointsFile is not None:
            self.atts['RvicPourPointsFile'] = RvicPourPointsFile

        if RvicUHFile is not None:
            self.atts['RvicUHFile'] = RvicUHFile

        if RvicFdrFile is not None:
            self.atts['RvicFdrFile'] = RvicFdrFile

        if RvicDomainFile is not None:
            self.atts['RvicDomainFile'] = RvicDomainFile

    def update(self):
        self.atts['history'] = \
            'Created: {0}'.format(time_mod.ctime(time_mod.time()))


# Coordinate Variables
time = dict(long_name='time',
            units=TIMEUNITS)

time_bnds = dict()

timesteps = dict(long_name='Series of timesteps',
                 nits='unitless')

lon = dict(long_name='longitude',
           nits='degrees_east')

lat = dict(long_name='latitude',
           units='degrees_north')

xc = dict(long_name='longitude',
          units='degrees_east')

yc = dict(long_name='latitude',
          units='degrees_north')

# Data Variables
fraction = dict(long_name='fraction of grid cell that is active',
                units='unitless')

unit_hydrograph = dict(long_name='Unit Hydrograph',
                       units='unitless')

avg_velocity = dict(long_name='Flow Velocity Parameter',
                    units='m s-1')

avg_diffusion = dict(long_name='Diffusion Parameter',
                     units='m2 s-1')

global_basin_id = dict(long_name='Global Basin ID from RvicFdrFile',
                       units='unitless')

full_time_length = dict(long_name='Length of original unit hydrograph',
                        units='timesteps')

subset_length = dict(long_name='Shortened length of the unit hydrograph',
                     units='timesteps')

unit_hydrograph_dt = dict(long_name='Unit hydrograph timestep',
                          units='seconds')

outlet_x_ind = dict(long_name='x grid coordinate of outlet grid cell',
                    units='unitless')

outlet_y_ind = dict(long_name='y grid coordinate of outlet grid cell',
                    units='unitless')

outlet_lon = dict(long_name='Longitude coordinate of outlet grid cell',
                  units='degrees_east')

outlet_lat = dict(long_name='Latitude coordinate of outlet grid cell',
                  units='degrees_north')

outlet_decomp_ind = dict(long_name='1d grid location of outlet grid cell',
                         units='unitless')

outlet_number = dict(long_name='outlet number',
                     units='unitless')

outlet_mask = dict(long_name='type of outlet point',
                   units='0-ocean, 1-land, 2-guage, 3-none')

outlet_name = dict(long_name='Outlet guage name',
                   units='unitless')

outlet_upstream_area = dict(long_name='Upstream catchment area contributing '
                                      'to outlet',
                            units='m2')

outlet_upstream_gridcells = dict(long_name='Number of upstream grid cells '
                                           'contributing to outlet',
                                 units='number of grid cells')

source_x_ind = dict(long_name='x grid coordinate of source grid cell',
                    units='unitless')

source_y_ind = dict(long_name='y grid coordinate of source grid cell',
                    units='unitless')

source_lon = dict(long_name='Longitude coordinate of source grid cell',
                  units='degrees_east')

source_lat = dict(long_name='Latitude coordinate of source grid cell',
                  units='degrees_north')

source_decomp_ind = dict(long_name='1d grid location of source grid cell',
                         units='unitless')
source_time_offset = dict(long_name='Number of leading timesteps ommited',
                          units='timesteps')

source2outlet_ind = dict(long_name='source to outlet index mapping',
                         units='unitless')

ring = dict(long_name='Convolution Ring',
            units='kg m-2 s-1')

streamflow = dict(long_name='Streamflow at outlet grid cell',
                  units='kg m-2 s-1')

storage = dict(long_name='Mass storage in stream upstream of outlet '
                         'grid cell',
               units='kg m-2 s-1')

# valid values http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/\
# cf-conventions.html#calendar
timemgr_rst_type = dict(long_name='calendar type',
                        units='unitless',
                        flag_values='0, 1, 2, 3, 4, 5, 6',
                        flag_meanings='NONE, NO_LEAP_C, GREGORIAN, '
                                      'PROLEPTIC_GREGORIAN, ALL_LEAP, '
                                      '360_DAY, JULIAN')

timemgr_rst_step_sec = dict(long_name='seconds component of timestep size',
                            units='sec',
                            valid_range='0, 86400')

timemgr_rst_start_ymd = dict(long_name='start date',
                             units='YYYYMMDD')

timemgr_rst_start_tod = dict(long_name='start time of day',
                             units='sec',
                             valid_range='0, 86400')

timemgr_rst_ref_ymd = dict(long_name='reference date',
                           units='YYYYMMDD')

timemgr_rst_ref_tod = dict(long_name='reference time of day',
                           units='sec',
                           valid_range='0, 86400')

timemgr_rst_curr_ymd = dict(long_name='current date',
                            units='YYYYMMDD')

timemgr_rst_curr_tod = dict(long_name='current time of day',
                            units='sec',
                            valid_range='0, 86400')

nhtfrq = dict(long_name='Frequency of history writes',
              units='absolute value of negative is in hours, 0=monthly, '
                    'positive is time-steps',
              comment='Namelist item')

mfilt = dict(long_name='Number of history time samples on a file',
             units='initless',
             comment='Namelist item')

ncprec = dict(long_name='Flag for data precision',
              flag_values='1, 2',
              flag_meanings='single-precision double-precision',
              comment='Namelist item',
              valid_range='1, 2')

fincl = dict(long_name='Fieldnames to include',
             comment='Namelist item')

fexcl = dict(long_name='Fieldnames to exclude',
             comment='Namelist item')

nflds = dict(long_name='Number of fields on file',
             units='unitless')

ntimes = dict(long_name='Number of time steps on file',
              units='time-step')
is_endhist = dict(long_name='End of history file',
                  flag_values='0, 1',
                  flag_meanings='FALSE TRUE',
                  comment='Namelist item',
                  valid_range='0, 1')

begtime = dict(long_name='Beginning time',
               units='time units')

hpindex = dict(long_name='History pointer index',
               units='units')

avgflag = dict(long_name='Averaging flag',
               units='A=Average, X=Maximum, M=Minimum, I=Instantaneous')

name = dict(long_name='Fieldnames')

long_name = dict(long_name='Long descriptive names for fields')

units = dict(long_name='Units for each history field output')
