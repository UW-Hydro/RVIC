"""
share.py
"""
from netCDF4 import default_fillvals
import sys
import time as time_mod

# ----------------------- CONSTANTS --------------------------------- #
earthRadius = 6.37122e6  # meters
waterDensity = 1000.

# area
metersPerKm = 1000.
metersPerMile = 1609.34
meters2PerAcre = 4046.856

# time
baseTime = '0000-1-1 0:0:0'
timeUnits = 'days since ' + baseTime
timeStampForm = '%Y-%m-%d-%H'
calendar = 'noleap'
hoursPerDay = 24.
secsPerHour = 3600.
minsPerHour = 60.
minsPerDay = hoursPerDay * minsPerHour
secsPerDay = hoursPerDay * secsPerHour

# length
mmPerMeter = 1000.
cmPerMeter = 100.

# fill values
fillValue_f = default_fillvals['f8']
fillValue_i = default_fillvals['i8']

# precision
precision = 1.0e-30
nc_double = 'f8'
nc_float = 'f4'
nc_int = 'i4'

# filenames
rpointer = 'rpointer'


# ----------------------- NETCDF VARIABLES --------------------------------- #
class ncGlobals:
    def __init__(self, title='',
                 casename='',
                 casestr='',
                 history='Created: {}'.format(time_mod.ctime(time_mod.time())),
                 institution='Univeristy of Washington', source=sys.argv[0],
                 references='Based on the initial model of Lohmann, et al., 1996, Tellus, 48(A), 708-721',
                 comment='Output from the RVIC Streamflow Routing Model.',
                 convention='CF-1.6',
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
        self.convention = convention
        self.RvicPourPointsFile = RvicPourPointsFile
        self.RvicUHFile = RvicUHFile
        self.RvicFdrFile = RvicFdrFile
        self.RvicDomainFile = RvicDomainFile

    def update(self):
        self.history = 'Created: {}'.format(time_mod.ctime(time_mod.time()))


class ncVar:
    def __init__(self, long_name='', units=''):
        self.long_name = long_name
        self.units = units

# Coordinate Variables
time = ncVar(long_name='time',
             units=timeUnits)

time_bnds = ncVar()

timesteps = ncVar(long_name='Series of timesteps',
                  units='unitless')

lon = ncVar(long_name='longitude',
            units='degrees_east')

lat = ncVar(long_name='latitude',
            units='degrees_north')

xc = ncVar(long_name='longitude',
           units='degrees_east')

yc = ncVar(long_name='latitude',
           units='degrees_north')

# Data Variables
fractions = ncVar(long_name='fraction of grid cell that is active',
                  units='unitless')

unit_hydrographs = ncVar(long_name='Unit Hydrograph',
                         units='unitless')

avg_velocity = ncVar(long_name='Flow Velocity Parameter',
                     units='m s-1')

avg_diffusion = ncVar(long_name='Diffusion Parameter',
                      units='m2 s-1')

lon_outlet = ncVar(long_name='Longitude coordinate of outlet grid cell',
                   units='degrees_east')

lat_outlet = ncVar(long_name='Latitude coordinate of outlet grid cell',
                   units='degrees_north')

global_basin_id = ncVar(long_name='Global Basin ID from RvicFdrFile',
                        units='unitless')

y_ind_outlet = ncVar(long_name='y grid coordinate of outlet grid cell',
                     units='unitless')

x_ind_outlet = ncVar(long_name='x grid coordinate of outlet grid cell',
                     units='unitless')

outlet_decomp_id = ncVar(long_name='1d grid location of outlet grid cell',
                        units='unitless')

subset_length = ncVar(long_name='Shortened length of the unit hydrograph',
                      units='timesteps')

full_time_length = ncVar(long_name='Length of original unit hydrograph',
                         units='timesteps')

unit_hydrogaph_dt = ncVar(long_name='Unit hydrograph timestep',
                          units='seconds')

source_decomp_id = ncVar(long_name='1d grid location of source grid cell',
                         units='unitless')

y_ind_source = ncVar(long_name='y grid coordinate of source grid cell',
                     units='unitless')

x_ind_source = ncVar(long_name='x grid coordinate of source grid cell',
                     units='unitless')

lon_source = ncVar(long_name='Longitude coordinate of source grid cell',
                   units='degrees_east')

lat_source = ncVar(long_name='Latitude coordinate of source grid cell',
                   units='degrees_north')

time_offset_source = ncVar(long_name='Number of leading timesteps ommited',
                           units='timesteps')

source2outlet_index = ncVar(long_name='source to outlet index mapping',
                            units='unitless')

outlet_num = ncVar(long_name='outlet number',
                   units='unitless')

streamflow = ncVar(long_name='Streamflow at outlet grid cell',
                   units='kg m-2 s-1')

storage = ncVar(long_name='Mass storage in stream upstream of outlet grid cell',
                units='kg m-2 s-1')
