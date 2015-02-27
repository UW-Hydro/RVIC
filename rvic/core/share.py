"""
share.py
"""
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
REFERENCE_DATE = 10101                        # i.e. REFERENCE_STRING
REFERENCE_TIME = 0                            # i.e. REFERENCE_STRING
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

VALID_CHARS = "-_. %s%s" % (string.ascii_letters, string.digits)


# ----------------------- NETCDF VARIABLES --------------------------------- #
class NcGlobals:
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

        if title:
            self.atts['title'] = title

        if comment:
            self.atts['comment'] = comment

        if Conventions:
            self.atts['Conventions'] = Conventions

        if history:
            self.atts['history'] = history

        if source:
            self.atts['source'] = source

        if institution:
            self.atts['institution'] = institution

        if hostname:
            self.atts['hostname'] = hostname
        else:
            self.atts['hostname'] = socket.gethostname()

        if username:
            self.atts['username'] = username
        else:
            self.atts['username'] = getuser()

        if casename:
            self.atts['casename'] = casename

        if references:
            self.atts['references'] = references

        if version:
            self.atts['version'] = version
        else:
            try:
                from rvic import version
                self.atts['version'] = version.short_version
            except:
                self.atts['version'] = 'unknown'

        if RvicPourPointsFile:
            self.atts['RvicPourPointsFile'] = RvicPourPointsFile

        if RvicUHFile:
            self.atts['RvicUHFile'] = RvicUHFile

        if RvicFdrFile:
            self.atts['RvicFdrFile'] = RvicFdrFile

        if RvicDomainFile:
            self.atts['RvicDomainFile'] = RvicDomainFile

    def update(self):
        self.atts['history'] = \
            'Created: {0}'.format(time_mod.ctime(time_mod.time()))


class NcVar:
    def __init__(self, **kwargs):
        for key, val in iteritems(kwargs):
            setattr(self, key, val)

# Coordinate Variables
time = NcVar(long_name='time',
             units=TIMEUNITS)

time_bnds = NcVar()

timesteps = NcVar(long_name='Series of timesteps',
                  units='unitless')

lon = NcVar(long_name='longitude',
            units='degrees_east')

lat = NcVar(long_name='latitude',
            units='degrees_north')

xc = NcVar(long_name='longitude',
           units='degrees_east')

yc = NcVar(long_name='latitude',
           units='degrees_north')

# Data Variables
fraction = NcVar(long_name='fraction of grid cell that is active',
                 units='unitless')

unit_hydrograph = NcVar(long_name='Unit Hydrograph',
                        units='unitless')

avg_velocity = NcVar(long_name='Flow Velocity Parameter',
                     units='m s-1')

avg_diffusion = NcVar(long_name='Diffusion Parameter',
                      units='m2 s-1')

global_basin_id = NcVar(long_name='Global Basin ID from RvicFdrFile',
                        units='unitless')

full_time_length = NcVar(long_name='Length of original unit hydrograph',
                         units='timesteps')

subset_length = NcVar(long_name='Shortened length of the unit hydrograph',
                      units='timesteps')

unit_hydrograph_dt = NcVar(long_name='Unit hydrograph timestep',
                           units='seconds')

outlet_x_ind = NcVar(long_name='x grid coordinate of outlet grid cell',
                     units='unitless')

outlet_y_ind = NcVar(long_name='y grid coordinate of outlet grid cell',
                     units='unitless')

outlet_lon = NcVar(long_name='Longitude coordinate of outlet grid cell',
                   units='degrees_east')

outlet_lat = NcVar(long_name='Latitude coordinate of outlet grid cell',
                   units='degrees_north')

outlet_decomp_ind = NcVar(long_name='1d grid location of outlet grid cell',
                          units='unitless')

outlet_number = NcVar(long_name='outlet number',
                      units='unitless')

outlet_mask = NcVar(long_name='type of outlet point',
                    units='0-ocean, 1-land, 2-guage, 3-none')

outlet_name = NcVar(long_name='Outlet guage name',
                    units='unitless')

outlet_upstream_area = NcVar(long_name='Upstream catchment area contributing '
                                       'to outlet',
                             units='m2')

outlet_upstream_gridcells = NcVar(long_name='Number of upstream grid cells '
                                            'contributing to outlet',
                                  units='number of grid cells')

source_x_ind = NcVar(long_name='x grid coordinate of source grid cell',
                     units='unitless')

source_y_ind = NcVar(long_name='y grid coordinate of source grid cell',
                     units='unitless')

source_lon = NcVar(long_name='Longitude coordinate of source grid cell',
                   units='degrees_east')

source_lat = NcVar(long_name='Latitude coordinate of source grid cell',
                   units='degrees_north')

source_decomp_ind = NcVar(long_name='1d grid location of source grid cell',
                          units='unitless')
source_time_offset = NcVar(long_name='Number of leading timesteps ommited',
                           units='timesteps')

source2outlet_ind = NcVar(long_name='source to outlet index mapping',
                          units='unitless')

ring = NcVar(long_name='Convolution Ring',
             units='kg m-2 s-1')

streamflow = NcVar(long_name='Streamflow at outlet grid cell',
                   units='kg m-2 s-1')

storage = NcVar(long_name='Mass storage in stream upstream of outlet '
                          'grid cell',
                units='kg m-2 s-1')

# valid values http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html#calendar
timemgr_rst_type = NcVar(long_name='calendar type',
                         units='unitless',
                         flag_values='0, 1, 2, 3, 4, 5, 6',
                         flag_meanings='NONE, NO_LEAP_C, GREGORIAN, '
                                       'PROLEPTIC_GREGORIAN, ALL_LEAP, '
                                       '360_DAY, JULIAN')

timemgr_rst_step_sec = NcVar(long_name='seconds component of timestep size',
                             units='sec',
                             valid_range='0, 86400')

timemgr_rst_start_ymd = NcVar(long_name='start date',
                              units='YYYYMMDD')

timemgr_rst_start_tod = NcVar(long_name='start time of day',
                              units='sec',
                              valid_range='0, 86400')

timemgr_rst_ref_ymd = NcVar(long_name='reference date',
                            units='YYYYMMDD')

timemgr_rst_ref_tod = NcVar(long_name='reference time of day',
                            units='sec',
                            valid_range='0, 86400')

timemgr_rst_curr_ymd = NcVar(long_name='current date',
                             units='YYYYMMDD')

timemgr_rst_curr_tod = NcVar(long_name='current time of day',
                             units='sec',
                             valid_range='0, 86400')

nhtfrq = NcVar(long_name='Frequency of history writes',
               units='absolute value of negative is in hours, 0=monthly, '
                     'positive is time-steps',
               comment='Namelist item')

mfilt = NcVar(long_name='Number of history time samples on a file',
              units='initless',
              comment='Namelist item')

ncprec = NcVar(long_name='Flag for data precision',
               flag_values='1, 2',
               flag_meanings='single-precision double-precision',
               comment='Namelist item',
               valid_range='1, 2')

fincl = NcVar(long_name='Fieldnames to include',
              comment='Namelist item')

fexcl = NcVar(long_name='Fieldnames to exclude',
              comment='Namelist item')

nflds = NcVar(long_name='Number of fields on file',
              units='unitless')

ntimes = NcVar(long_name='Number of time steps on file',
               units='time-step')
is_endhist = NcVar(long_name='End of history file',
                   flag_values='0, 1',
                   flag_meanings='FALSE TRUE',
                   comment='Namelist item',
                   valid_range='0, 1')

begtime = NcVar(long_name='Beginning time',
                units='time units')

hpindex = NcVar(long_name='History pointer index',
                units='units')

avgflag = NcVar(long_name='Averaging flag',
                units='"A=Average, X=Maximum, M=Minimum, I=Instantaneous')

name = NcVar(long_name='Fieldnames')

long_name = NcVar(long_name='Long descriptive names for fields')

units = NcVar(long_name='Units for each history field output')
