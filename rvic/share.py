"""
share.py
"""
import sys
from netCDF4 import default_fillvals
import time as time_mod

# ----------------------- CONSTANTS --------------------------------- #
EARTHRADIUS = 6.37122e6  # meters
WATERDENSITY = 1000.

# area
METERSPERKM = 1000.
METERSPERMILE = 1609.34
METERS2PERACRE = 4046.856

# time
BASETIME = '0000-1-1 0:0:0'
TIMEUNITS = 'days since ' + BASETIME
TIMESTAMPFORM = '%Y-%m-%d-%H'
CALENDAR = 'noleap'
HOURSPERDAY = 24.
SECSPERHOUR = 3600.
MINSPERHOUR = 60.
MINSPERDAY = HOURSPERDAY * MINSPERHOUR
SECSPERDAY = HOURSPERDAY * SECSPERHOUR

# length
MMPERMETER = 1000.
CMPERMETER = 100.

# precision
PRECISION = 1.0e-30
NC_DOUBLE = 'f8'
NC_FLOAT = 'f4'
NC_INT = 'i4'

# fill values
FILLVALUE_F = default_fillvals[NC_DOUBLE]
FILLVALUE_I = default_fillvals[NC_INT]

# filenames
RPOINTER = 'rpointer'


# ----------------------- NETCDF VARIABLES --------------------------------- #
class NcGlobals:
    def __init__(self, title='',
                 casename='',
                 casestr='',
                 history='Created: {}'.format(time_mod.ctime(time_mod.time())),
                 institution='Univeristy of Washington', source=sys.argv[0],
                 references='Based on the initial model of Lohmann, et al., 1996, Tellus, 48(A), 708-721',
                 comment='Output from the RVIC Streamflow Routing Model.',
                 Conventions='CF-1.6',
                 RvicPourPointsFile='',
                 RvicFdrFile='',
                 RvicUHFile='',
                 RvicDomainFile=''):
        self.title = title
        self.casename = casename
        self.history = history
        self.institution = institution
        self.source = source
        self.references = references
        self.comment = comment
        self.Conventions = Conventions
        self.RvicPourPointsFile = RvicPourPointsFile
        self.RvicUHFile = RvicUHFile
        self.RvicFdrFile = RvicFdrFile
        self.RvicDomainFile = RvicDomainFile

    def update(self):
        self.history = 'Created: {}'.format(time_mod.ctime(time_mod.time()))


class NcVar:
    def __init__(self, long_name='', units=''):
        self.long_name = long_name
        self.units = units

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
fractions = NcVar(long_name='fraction of grid cell that is active',
                  units='unitless')

unit_hydrographs = NcVar(long_name='Unit Hydrograph',
                         units='unitless')

avg_velocity = NcVar(long_name='Flow Velocity Parameter',
                     units='m s-1')

avg_diffusion = NcVar(long_name='Diffusion Parameter',
                      units='m2 s-1')

lon_outlet = NcVar(long_name='Longitude coordinate of outlet grid cell',
                   units='degrees_east')

lat_outlet = NcVar(long_name='Latitude coordinate of outlet grid cell',
                   units='degrees_north')

global_basin_id = NcVar(long_name='Global Basin ID from RvicFdrFile',
                        units='unitless')

y_ind_outlet = NcVar(long_name='y grid coordinate of outlet grid cell',
                     units='unitless')

x_ind_outlet = NcVar(long_name='x grid coordinate of outlet grid cell',
                     units='unitless')

outlet_decomp_id = NcVar(long_name='1d grid location of outlet grid cell',
                        units='unitless')

subset_length = NcVar(long_name='Shortened length of the unit hydrograph',
                      units='timesteps')

full_time_length = NcVar(long_name='Length of original unit hydrograph',
                         units='timesteps')

unit_hydrogaph_dt = NcVar(long_name='Unit hydrograph timestep',
                          units='seconds')

source_decomp_id = NcVar(long_name='1d grid location of source grid cell',
                         units='unitless')

y_ind_source = NcVar(long_name='y grid coordinate of source grid cell',
                     units='unitless')

x_ind_source = NcVar(long_name='x grid coordinate of source grid cell',
                     units='unitless')

lon_source = NcVar(long_name='Longitude coordinate of source grid cell',
                   units='degrees_east')

lat_source = NcVar(long_name='Latitude coordinate of source grid cell',
                   units='degrees_north')

time_offset_source = NcVar(long_name='Number of leading timesteps ommited',
                           units='timesteps')

source2outlet_index = NcVar(long_name='source to outlet index mapping',
                            units='unitless')

outlet_num = NcVar(long_name='outlet number',
                   units='unitless')

streamflow = NcVar(long_name='Streamflow at outlet grid cell',
                   units='kg m-2 s-1')

storage = NcVar(long_name='Mass storage in stream upstream of outlet grid cell',
                units='kg m-2 s-1')
