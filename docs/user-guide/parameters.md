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

*Note: If `REMAP==False`, the resolution of the routing data grid and domain grid must be the same.  Furthermore, the domain grid must be a subset of the routing grid.  RVIC will raise an error if either of these conditions is not met.*

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
3.  **CLEAN**
    - Description: Delete temporary files, only used if REMAP=True
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
7.  **TEMP_DIR**
    - Description: Directory to use for temporary read/write operations, only used if REMAP=True.
    - Type: char
8.  **REMAP**
    - Description: Remap Unit Hydrographs from [ROUTING] grid to [DOMAIN] grid
    - Type: bool
    - Valid values: True, False
9.  **AGGREGATE**
    - Description: Aggregate all [POUR_POINTS] inside each [DOMAIN] grid cell
    - Type: bool
    - Valid values: True, False

    Note: *This should only be used when routing to coastal grid cells for CESM.*
10.  **AGG_PAD**
    - Description: Size of pad to add to aggregated files prior to remapping
    - Type: int
11.  **NETCDF_FORMAT**
    - Description: Output parameter file format
    - Type: char
    - Valid values: NETCDF3_CLASSIC, NETCDF3_64BIT, NETCDF4_CLASSIC, and NETCDF4

    Note: *For use with CESM, NETCDF3_CLASSIC is recommended.*
12.  **SUBSET_DAYS**
    - Description: Length of unit hydrograph subset in days
    - Type: int

13.  **CONSTRAIN_FRACTIONS**

    - Description: Constrain the final unit hydrographs sum to be less than or equal to the domain fractions
    - Type: bool
    - Valid values: True, False

    Note: *True when routing to coastal grid cells, else False.*

###POUR_POINTS
1.  **FILE_NAME**
    - Description: Path to Pour Points File, A comma separated file of outlets to route to [lons, lats] - one coordinate pair per line (order not important).  May optionally include a column [names] - which will (if not aggregating) be included in param file.
    - Type: char

###UH_BOX
1.  **FILE_NAME**
    - Description: Path to UH Box File.  This defines the unit hydrograph to rout flow to the edge of each grid cell.  A comma separated file of [time in seconds, unit hydrograph ordinate] - one timestep per line.  The timestep should be 1hr (3600 sec) or less.
    - Type: char
2.  **HEADER_LINES**
    - Description: Number of Header lines to ignore in [UH_BOX]FILE_NAME
    - Type: int

###ROUTING
1.  **FILE_NAME**
    - Description: Path to routing inputs netCDF.
    - Type: char
2.  **LONGITUDE_VAR**
    - Description: Longitude variable name
    - Type: char
3.  **LATITUDE_VAR**
    - Description: Latitude variable name
    - Type: char
4.  **FLOW_DISTANCE_VAR**
    - Description: Flow Distance variable name
    - Type: char
5.  **FLOW_DIRECTION_VAR**
    - Description: Flow Direction variable name
    - Type: char
6.  **BASIN_ID_VAR**
    - Description: Basin ID variable name
    - Type: char
7.  **VELOCITY**
    - Description: Velocity variable name or value
    - Type: char, float
8.  **DIFFUSION**
    - Description: Diffusion variable name or value
    - Type: char, float
9.  **OUTPUT_INTERVAL**
    - Description: Timestep of output unit hydrographs.  Must be a multiple of the timestep in the UH_BOX
    - Type: int
10.  **BASIN_FLOWDAYS**
    - Description: Maximum time for runoff to reach outlet in days
    - Type: int
11. **CELL_FLOWDAYS**
    - Description: Maximum time for runoff to pass through a grid cell in days
    - Type: int

###DOMAIN
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
