# Developing Parameters

This documentation outlines how to configure RVIC to generate the impulse response function paramters.

  1.  Begin by developing the necessary inputs (see below)
  2.  Setup the `parameters` configuration file](https://github.com/UW-Hydro/RVIC/wiki/parameters-configuration-file)
  3.  run `rvic parameters your_config_file -np=$NPROCESSORS`

    where `$NPROCESSORS` is the number of processors to use.

# Input Files

## Flow Direction File
####1.  Obtain Flow Direction Raster
This file is typically a ArcGIS raster that describes the flow direction from a grid cell to its downstream neighbor.
source: [Wu et al. 2012](http://www.ntsg.umt.edu/project/drt) or elsewhere
reference: [ArcGIS Help](http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//009z00000063000000.htm)

*Note: If you application is global and you want to use the Wu et al. (2012) dataset, use the Hydro1K product instead of Hydro1K-HydroSHED combined product since there were issues with the 1/16th degree level Hydro1K_HydroSHED data.*

####2.  Calculate grid box Flow Distance
The flow distance is the distance in the direction of the flow direction either zonal, meridional, or diagonal.  This is a per grid cell value and is *NOT* the total flow distance between a point and its basin outlet.

Use `create_xmask.c` in [rout_prep](ftp://ftp.hydro.washington.edu/pub/HYDRO/models/VIC/Utility_Programs/rout_prep.tgz).

*Note 1: the ArcGIS flow directions are different from the VIC flow directions used in this C program, so you'd need to make some minor modifications to the C program for it the work correctly*

*Note 2:  Future version of the RVIC model `parameters.py` will include a function to calculate this field on the fly.*

####3.  Calculate basin_id and land mask
Using the Flow Direction raster and ArcGIS, develop two new raster datasets:

- basin_id (Spatial Analyst - Hydrology - Basin)
- land_mask (Raster Calculator: NoData = 0, everything else = 1)

####4. Calculate Source_Area
You can use the "upstream_drainagearea" Ascii files that are part of Wu's data (units = m^2), or obtain them using ArcGIS (Spatial Analyst - Hydrology - Flow Accumulation) (units = # of grid cells).

####5.  Combine all into netcdf format file
Use gdal translate to convert all Ascii files to netcdf, then use ncks to combine all variables into 1 file (append). Use ncatted to modify the long_names and units for future reference.

`gdal_translate -of netCDF <ascii input file> <netCDF file>`

*Note: gdal and ArcGIS can both produce netCDF formatted files.  Do not use them both!  They use opposite conventions for the y-axis array order.*

#### Optional input fields:
You may provide the variables velocity and/or diffusion in the netCDF file.

## Domain File
RVIC has the capacity to remap from the flow direction grid to the land grid.  The domain file describes the land model grid information.

The easiest way to develop a domain file is from an existing routing setup fraction file.  Use [fraction2domain.bash](https://github.com/UW-Hydro/RVIC/blob/master/tools/fraction2domain.bash).

For irregular grids, you need to come up with a netCDF file that contains these variables:

- land_mask (0 or 1)
- land_fractions (between zero and 1)
- area (in m2 or rad2 or similar)
- coordinate vars (lon and lat or similar)

## Pour Points File
A pour points CSV is needed to develop the parameter file.  This is simply a csv file containing a list of the pour point locations with 2 or 3 columns (required: lons, lats; optional: names/contributing area/etc).  One header row is required.

***Table:*** *Example Pour Points File Format*

| lons          |  lats         |names                          |
| ------------- |:-------------:| -----:                        |
| -113.94       |  49.48        |PINCHER CREEK AT PINCHER CREEK |
| -114.12       |  49.65        |TODD CREEK AT ELTON'S RANCH    |
| -114.4        |  49.59        |CROWSNEST RIVER AT FRANK       |

## UH BOX File
A csv file that describes the routing of flow to the edge of the origin grid cell. This is essentially a unit hydrograph for a basin with area = average size of the grid cells, and characteristic length = average flow distance of the grid cells. One way to calculate this is to use the SCS dimensionless unit hydrograph approach (also need to assume an overland flow velocity to make the calculations).

***Table:*** *Example UHBOX File Format*

| time          |  UHBOX      |
| ------------- |-------------|
| 0             |  0.00000006 |
| 3600          |  0.00048497 |
| 7200          |  0.00714538 |
| 10800         |  0.02295375 |
| ...           |  ...        |
| 169200        |  0.00068766 |

![UH BOX](./images/uh-box.png)
***Figure:*** *Example UH Box.  Values should sum to 1 (blue).*

## RVIC Parameter Configuration File

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
4.  **RVIC_TAG**
    - Description: RVIC Tag
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
