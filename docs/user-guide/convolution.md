# Flux Convolution

Once the impulse response function parameters have been generated, the flux convolution step may be done.

#### Doing the Streamflow Routing
  1.  Setup the `convolution` configuration file ([see below](user-guide/convolution/#rvic-convolution-configuration-file))
  2.  run `rvic convolution your_config_file` from the command line

Examples of all configuration files can be found in `[RVIC/config/](https://github.com/UW-Hydro/RVIC/tree/master/config)` directory.

## Flux File(s)

The flux file(s) must be in netCDF format and have a `time` dimension as well as matching spatial dimensions (`y` and `x`) as the domain file.

## RVIC Convolution Configuration File

*Note:  configuration file is parsed by the [python ConfigParser module](https://docs.python.org/2/library/configparser.html#module-ConfigParser).  %(Interploation) is supported inside [sections] only.*

### OPTIONS
1.  **LOG_LEVEL**
    - Description: Level to log output at
    - Type: char
    - valid values: DEBUG, INFO, WARNING, ERROR, CRITICAL
2.  **VERBOSE**
    - Description: Print output to console in addition to the log file
    - Type: bool
    - valid values: True, False
3.  **CASE_DIR**
    - Description: case run directory
    - Type: char
5.  **CASEID**
    - Description: Case ID
    - Type: char
6.  **CASESTR**
    - Description: Case description
    - Type: char
7.  **CALENDAR**
    - Description: Calendar
    - Type: char
    - Valid Values: standard, gregorian, proleptic_gregorian noleap, 365_day, 360_day, julian, all_leap, 366_day
8.  **RUN_TYPE**
    - Description: Run initialization type
    - Type: char
    - Valid Values: startup, drystart, restart
9.  **RUN_STARTDATE**
    - Description: Run start date (yyyy-mm-dd-hh). Only used for startup and drystart runs.
    - Type: char
10.  **STOP_OPTION**
    - Description:  Run stop condition
    - Type: char
    - Valid Values: none, never, nsteps, nseconds, nminutes, nhours, ndays, nmonths, nyears, date, end
11.  **STOP_N**
    - Description:  Run length based on STOP_OPTION
    - Type: int
12.  **STOP_DATE**
    - Description:  Run stop date based on STOP_OPTION
    - Type: char
13.  **REST_OPTION**
    - Description:  Frequency of model restart writes
    - Type: char
    - Valid Values: none, never, nsteps, nseconds, nminutes, nhours, ndays, nmonths, nyears, date, end
14.  **REST_N**
    - Description:  Write restart frequency based on REST_OPTION
    - Type: int
15.  **STOP_DATE**
    - Description:  Write restart date based on REST_OPTION
    - Type: char
16.  **REST_NCFORM**
    - Description: Restart file format
    - Type: char
    - Valid values: NETCDF3_CLASSIC, NETCDF3_64BIT, NETCDF4_CLASSIC, and NETCDF4
17.  **Output file compression options**

    Descriptions of these options can be found in the [netCDF4-Python package](http://unidata.github.io/netcdf4-python/netCDF4.Dataset-class.html#createVariable)

    - NETCDF_ZLIB: False (bool)
    - NETCDF_COMPLEVEL: 4 (int)
    - NETCDF_SIGFIGS: None (bool or int)

### HISTORY
1.  **RVICHIST_NTAPES**
    - Description:  Number of history file output streams (a.k.a. history tapes).
    - Type: int
2.  **RVICHIST_MFILT**
    - Description:  Per tape series maximum number of time samples per output file.
    - Type: int
3.  **RVICHIST_NDENS**
    - Description:  Per tape series history file density (i.e. output precision)
        - 1=double precision
        - 2=single precision
    - Type: int
    - Valid Values: 1, 2
4.  **RVICHIST_NHTFRQ**
    - Description:  Per tape series history write frequency.
        - positive means in time steps
        - 0 = monthly
        - negative means hours
    - Type: int
5.  **RVICHIST_AVGFLAG**
    - Description:  Per tape series history output type.
        - A - Average, over the output interval.
        - I - Instantaneous, output the value at the output interval.
        - X - Maximum, over the output interval.
        - M - Minimum, over the output interval.
    - Type: char
    - Valid Values: A, I, X, M
6.  **RVICHIST_OUTTYPE**
    - Description: History file output shape
        - grid - shape is (time, y, x)
        - array - shape is (time, outlets)
    - Type: char
    - Valid values: grid, array
7.  **RVICHIST_NCFORM**
    - Description: Restart file format
    - Type: char
    - Valid values: NETCDF3_CLASSIC, NETCDF3_64BIT, NETCDF4_CLASSIC, and NETCDF4
8.  **RVICHIST_UNITS**
    - Description: Per tape series output units
    - Type: char
    - Valid values: kg m-2 s-1, m3/s

### DOMAIN
1.  **FILE_NAME**
    - Description: Path to CESM complaint domain file
    - Type: char
2.  **LONGITUDE_VAR**
    - Description: Longitude variable name
    - Type: char
3.  **LATITUDE_VAR**
    - Description: Latitude variable name
    - Type: char
4.  **LAND_MASK_VAR**
    - Description: Land Mask variable name
    - Type: char
5.  **FRACTION_VAR**
    - Description: Land fraction of Grid Cell
    - Type: char
6.  **AREA_VAR**
    - Description: Grid Cell Area
    - Type: char

### INITIAL_STATE
1.  **FILE_NAME**
    - Description: RVIC state file
    - Type: char

### PARAM_FILE
1.  **FILE_NAME**
    - Description: rvic parameter file file
    - Type: char

### INPUT_FORCINGS
1.  **DATL_PATH**
    - Description: Path to directory with land data netCDF forcings
    - Type: char
2.  **DATL_FILE**
    - Description: format of land data files (prfix.$YYYY[-$MM-[$DD[-$HH]]].nc)
    - Type: char
3.  **TIME_VAR**
    - Description: Time variable name
    - Type: char
4.  **LATITUDE_VAR**
    - Description: Latitude variable name
    - Type: char
5.  **DATL_LIQ_FLDS**
    - Description: Liquid variable names (e.g. runoff, baseflow)
    - Type: char
6.  **START**
    - Description: start date, date format YYYY[-MM[-DD]]
    - Type: char
7.  **END**
    - Description: start date, date format YYYY[-MM[-DD]]
    - Type: char
