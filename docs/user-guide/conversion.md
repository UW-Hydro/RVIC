# Parameter Conversion

**If you are migrating from an existing routing model setup based on the C version of the original [VIC routing model](https://github.com/UW-Hydro/VIC_Routing):**

  1.  Setup the `convert` config file: [see below](conversion/#conversion-configuration-file)
  2.  run `rvic convert your_config_file`

*Note, only the C version of the routing model is currently able to be converted to RVIC parameters.  If you need to convert a routing setup from the Fortran version of the routing model, a converter will need to be written in [`rvic/core/convert.py`](https://github.com/UW-Hydro/RVIC/blob/feature/output_file_chunks/rvic/core/convert.py).  A stub for this is already present in the code.*

## Conversion Configuration File

### OPTIONS
1.  **LOG_LEVEL**
    - Description: Level to log output at
    - Type: char
    - valid values: DEBUG, INFO, WARNING, ERROR, CRITICAL
2.  **VERBOSE**
    - Description: Print output to console in addition to the log file
    - Type: bool
    - valid values: True, False
4.  **CASEID**
    - Description: Case ID
    - Type: char
5.  **GRIDID**
    - Description: routing domain grid shortname
    - Type: char
6.  **CASE_DIR**
    - Description: case run directory
    - Type: char
7.  **NETCDF_FORMAT**

    _Note: For use with CESM, NETCDF3_CLASSIC is recommended._

    - Description: Output parameter file format
    - Type: char
    - Valid values: NETCDF3_CLASSIC, NETCDF3_64BIT, NETCDF4_CLASSIC, and NETCDF4

8.  **Output parameter file compression options**

    Descriptions of these options can be found in the [netCDF4-Python package](http://unidata.github.io/netcdf4-python/netCDF4.Dataset-class.html#createVariable)

    - NETCDF_ZLIB: False (bool)
    - NETCDF_COMPLEVEL: 4 (int)
    - NETCDF_SIGFIGS: None (bool or int)

9.  **SUBSET_DAYS**
    - Description: Length of unit hydrograph subset in days
    - Type: int
10.  **CONSTRAIN_FRACTIONS**

    _Note: True when routing to coastal grid cells, else False_

    - Description: Constrain the final unit hydrographs sum to be less than or equal to the domain fractions
    - Type: bool
    - Valid values: True, False

### UHS_FILES
1.  **ROUT_PROGRAM**
    - Description: Routing program used to create UHS files
    - Type: char
    - Valid Values: C, Fortran
2.  **ROUT_DIR**
    - Description: Location of UHS files
    - Type: char
3.  **STATION_FILE**
    - Description: Path to stations file
    - Type: char

### ROUTING
1.  **OUTPUT_INTERVAL**
    - Description: Timestep of output unit hydrographs.  Must be a multiple of the timestep in the UH_BOX
    - Type: int

### DOMAIN

Domain file describing the grid that the UHS files were developed on.

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

###NEW DOMAIN
Domain file describing the grid that the routing will be done on.  (Optional)

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
    - Description: Latitude variable name
    - Type: char
6.  **AREA_VAR**
    - Description: Longitude variable name
    - Type: char
