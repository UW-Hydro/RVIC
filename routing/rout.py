#!/usr/local/bin/python
"""
PROGRAM rout, Python-Version, written by Joe Hamman winter 2012/2013
Routing algorithm developed by D. Lohmann.
"""

##############################################################################
import sys
import numpy as np
import argparse
import ConfigParser
import time as tm
from netCDF4 import Dataset

##############################################################################
###############################  MAIN PROGRAM ################################
##############################################################################
def main(config_file = None):

    (infile, UHfile, Plons, Plats, 
     velocity, diffusion, verbose,
     NODATA, CELL_FLOWTIME, BASIN_FLOWTIME, 
     PREC, OUTPUT_INTERVAL, DAY_SECONDS) = process_command_line(config_file = config_file)

    for i, (basin_y, basin_x) in enumerate(zip(Plats, Plons)):
        out_file = rout(infile, UHfile, basin_y, basin_x, velocity,diffusion, 
                        verbose, NODATA, CELL_FLOWTIME, BASIN_FLOWTIME, PREC, 
                        OUTPUT_INTERVAL, DAY_SECONDS)
        if verbose:
            print 'Finished routing to point %i of %i (%f, %f)' \
                    % (i+1, len(Plons), basin_y, basin_x)
            print 'Wrote %s' % out_file
    if verbose:
        print 'Routing Program Finished.'
    return
        
def rout(infile, UHfile, basin_y, basin_x, velocity, diffusion, verbose,
    NODATA, CELL_FLOWTIME, BASIN_FLOWTIME, PREC, OUTPUT_INTERVAL, DAY_SECONDS):

    Basin, dy, dx, basin_id = init(infile, basin_y, basin_x, velocity, diffusion, verbose)
       
    # Load UH_BOX input
    (uh_t,UH_Box) = load_uh(UHfile, verbose)
    
    # Find timestep (timestep is determined from UH_BOX input file)
    INPUT_INTERVAL = find_TS(uh_t,verbose)
    T_Cell = CELL_FLOWTIME*DAY_SECONDS/INPUT_INTERVAL
    T_UH = BASIN_FLOWTIME*DAY_SECONDS/INPUT_INTERVAL
    
    # Read direction grid and find to_col (to_x) and to_row (to_y)
    to_y, to_x = read_direction(Basin['Flow_Direction'], Basin['Basin_ID'],
                                dy, dx, basin_id, NODATA, verbose)
    
    # Find row/column indicies of lat/lon inputs
    x_ind = find_nearest(Basin['lon'], basin_x)
    y_ind = find_nearest(Basin['lat'], basin_y)
    
    # Find all grid cells upstream of pour point
    Catchment, fractions = search_catchment(to_y, to_x, y_ind, x_ind,
                                            Basin['Basin_ID'], basin_id, 
                                            verbose)
    
    # Make UH for each grid cell upstream of basin pour point 
    # (linear routing model - Saint-Venant equation)
    UH = make_UH(INPUT_INTERVAL, T_Cell,Catchment['y_inds'], 
                Catchment['x_inds'], Basin['Velocity'], Basin['Diffusion'], 
                Basin['Flow_Distance'], PREC,verbose)
    
    # Make UH_RIVER by incrementally moving upstream comining UH functions
    UH_RIVER = make_grid_UH_river(T_UH, T_Cell, UH, to_y, to_x, y_ind, x_ind,
                                  Catchment['y_inds'], Catchment['x_inds'], 
                                  Catchment['count_ds'], PREC, verbose)
    
    # Make UH_S for each grid cell upstream of basin pour point 
    # (combine IRFs for all grid cells in flow path)
    UH_S = make_grid_UH(T_UH, T_Cell, UH_RIVER, UH_Box, to_y, to_x,
                        Catchment['y_inds'], Catchment['x_inds'], 
                        Catchment['count_ds'], PREC,NODATA,verbose)
    
    # Agregate to output timestep
    UH_out = aggregate(UH_S, T_UH, INPUT_INTERVAL, OUTPUT_INTERVAL,
                       Catchment['x_inds'], Catchment['y_inds'], NODATA, verbose)
    
    #Write to output netcdf
    time_steps = np.arange(UH_out.shape[0])
    times = np.linspace(0, OUTPUT_INTERVAL*UH_out.shape[0], UH_out.shape[0],
                        endpoint=False)
    time_res = "%s seconds" % OUTPUT_INTERVAL
    out_file = write_netcdf(basin_x, basin_y, Basin['lon'], Basin['lat'],
                            times, time_steps, time_res, UH_out, fractions, 
                            velocity, diffusion, basin_id, NODATA, verbose)

    return out_file
    
##############################################################################
############### Routines #####################################################
##############################################################################
def init(infile, basin_y, basin_x, velocity, diffusion, verbose):
    """
    Handle the initial reading/clipping of input grids
    """
    Inputs = {}
    vars = ['Basin_ID', 'lat', 'lon']
    f = Dataset(infile, 'r')
    for var in vars:
        Inputs[var] = f.variables[var][:]
    
    # Find Basin Dims and ID
    # Reads input lons/lats/basins_ids and returns basin bounds.
    if verbose:
        print 'Reading Global Inputs'
    
    outlet_loc = (find_nearest(Inputs['lat'], basin_y), find_nearest(Inputs['lon'], basin_x))
    basin_id = Inputs['Basin_ID'][outlet_loc]
    
    if verbose:
        print 'Input Latitude:', basin_y
        print 'Input Longitude:', basin_x
        print 'Input Basid ID:', basin_id
    inds = np.nonzero(Inputs['Basin_ID'] == basin_id)
    x,y = np.meshgrid(np.arange(len(Inputs['lon'])), np.arange(len(Inputs['lat'])))
    x_min = min(x[inds])
    x_max = max(x[inds])+1
    y_min = min(y[inds])
    y_max = max(y[inds])+1

    # Load input arrays, store in python dictionary.  (Format -Basin['var'])
    vars = ['Basin_ID', 'Flow_Direction', 'Flow_Distance', 'lon', 'lat']
    if not velocity:
        vars.append('velocity')
    if not diffusion:
        vars.append('diffusion')

    #  Clip netCDF Inputs
    if verbose:
        print 'Reading input data vars: %s' % ", ".join(vars)
    Basin={}
    for var in vars:
        try:
            temp = f.variables[var]
        except NameError:
            print 'Unable to clip %s. Confirm that the var exists in %s' \
                    % (var, nc_str)
            raise
        if np.rank(temp) > 1:
            Basin[var] = temp[y_min:y_max, x_min:x_max]
        elif var == 'lon':
            Basin['lon'] = temp[x_min:x_max]
        elif var == 'lat':
            Basin['lat'] = temp[y_min:y_max]
    
    if velocity:
        Basin['Velocity'] = np.zeros((Basin['Flow_Direction'].shape))+velocity
    if diffusion:
        Basin['Diffusion'] = np.zeros((Basin['Flow_Direction'].shape))+diffusion
    
    if 'VIC' in f.variables['Flow_Direction'].units:
        # VIC Directions: http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/Routing/FlowDirection.shtml
        dy = {1:-1, 2:-1, 3:0, 4:1, 5:1, 6:1, 7:0, 8:-1}
        dx = {1:0, 2:1, 3:1, 4:1, 5:0, 6:-1, 7:-1, 8:-1}
        if verbose:
            print 'Using VIC flow directions (1-8).'
    else:
        # ARCMAP Directions: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=flow_direction
        dy = {1:0, 2:1, 4:1, 8:1, 16:0, 32:-1, 64:-1, 128:-1}
        dx = {1:1, 2:1, 4:0, 8:-1, 16:-1, 32:-1, 64:0, 128:1}
        if verbose:
            print 'Using ARCMAP flow directions (1-128).'
    
    f.close()
    if verbose:
        print 'grid cells in subset: %i' % Basin['Velocity'].size

    # Check latitude order, flip if necessary.
    if Basin['lat'][-1]>Basin['lat'][0]:
        if verbose:
            print 'Inputs came in upside down, flipping everything now.'
        vars.remove('lon')
        for var in vars:
            Basin[var] = np.flipud(Basin[var])
        
    return Basin, dy, dx, basin_id

def process_command_line(config_file = None):
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h

    A configuration file may be provided with -C configFile
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-C", "--configFile", type = str, 
        help = "Input configuration file", default = False)
    parser.add_argument("-i", "--infile", type = str, 
        help = "Input netCDF containing all input grids")
    parser.add_argument("-UH", "--UHfile", type = str, 
        help = "Input UH_BOX hydrograph")
    parser.add_argument("-lon", "--longitude", type = float, nargs = '*',
        help = "Input longitude of point to route to")
    parser.add_argument("-lat", "--latitude", type = float,nargs = '*',
        help = "Input latitude of point to route to")
    parser.add_argument("-vel", "--velocity", type = float, 
        help = "Input wave velocity, may be variable in infile or scalar")
    parser.add_argument("-diff", "--diffusion", type = float, 
        help = "Input diffusion")
    parser.add_argument("-v", "--verbose", action = "store_true", 
        help = "Increase output verbosity")
    parser.add_argument("--NODATA", type = float, default = 9.96920996839e+36,
        help = "Input integer to bedefined as NODATA value")
    parser.add_argument("--CELL_FLOWTIME", type = int, default = 2,
        help = "Input integer to be defined as max number of days for flow to"
        "pass through cell")
    parser.add_argument("--BASIN_FLOWTIME", type = int, default = 50,
        help = "Input integer to be defined as max number of days for flow to "
        "pass through basin")
    parser.add_argument("-PREC", "--Precision", type = int, default = 1e-30,
        help = "Input integer to be defined asminimum precision for calculations")
    parser.add_argument("-TS","--OUTPUT_INTERVAL", type = int, default = 86400,
        help = "Output timestep in seconds for Unit Hydrographs")
    parser.add_argument("--DAY_SECONDS", type = int,default = 86400,
        help = "Seconds per day")
    args = parser.parse_args()

    # Assign values
    if config_file:
        process_config(config_file)
    elif args.configFile:
        file_paths,inputs = process_config(args.configFile)
    if args.infile:
        infile = args.infile
    else:
        try:
            infile = file_paths['infile']
        except:
            raise IOError('Need input file from command line or ' 
                'configuration file')

    if args.UHfile:
        UHfile = args.UHfile
    else:
        try:
            UHfile = file_paths['uhfile']
        except:
            raise IOError('Need input Unit Hydrograph File from command line ' 
                'or configuration file')

    if args.longitude:
        Plons = args.longitude
    else:
        try:
            Plons = map(float,inputs['longitude'])
        except:
            raise IOError('Need logitude(s) from command line or '
                'configuration file')
        
    if args.latitude:
        Plats = args.latitude
    else:
        try:
            Plats = map(float,inputs['latitude'])
        except:
            raise IOError('Need latitude(s) from command line or '
                'configuration file')        
    
    try:
        velocity = float(inputs['velocity'])
    except:
        if not args.velocity:
            velocity = None
        else:
            velocity = args.velocity
    try:
        diffusion = float(inputs['diffusion'])
    except:
        if not args.diffusion:
            diffusion = None
        else:
            diffusion = args.diffusion

    try:
        if inputs['verbose'] == 'True':
            verbose = True
        else:
            verbose = False
    except:
        verbose = args.verbose

    try:
        NODATA = float(inputs['nodata'])
    except:
        NODATA = args.NODATA

    try:
        CELL_FLOWTIME = int(inputs['cell_flowtime'])
    except:
        CELL_FLOWTIME = args.CELL_FLOWTIME

    try:
        BASIN_FLOWTIME = int(inputs['basin_flowtime'])
    except:
        BASIN_FLOWTIME = args.BASIN_FLOWTIME

    try:
        PREC = float(inputs['precision'])
    except:
        PREC = args.Precision

    try:
        OUTPUT_INTERVAL = int(inputs['output_interval'])
    except:
        OUTPUT_INTERVAL = args.OUTPUT_INTERVAL

    try:
        DAY_SECONDS = int(inputs['day_seconds'])
    except:
        DAY_SECONDS = args.DAY_SECONDS
        
    return (infile, UHfile, Plons, Plats, velocity, diffusion, verbose, 
            NODATA, CELL_FLOWTIME, BASIN_FLOWTIME, PREC, OUTPUT_INTERVAL, 
            DAY_SECONDS)

##############################################################################
##  Process Configuration File
##  Values that aren't provided will be given a False bool and filled in by
##  the command line processor
##############################################################################
def process_config(configFile):
    """
    Parse arguments and assign flags for further loading of files.  
    Configuration flag must be raised on command line (-C) and configuration 
    file must be provided.
    Usage:  rout.py -C rout.cfg
    There is no set of variables that must be included, any variable not given 
    in configuration file will be filled in during the rest of the command 
    line parsing.
    """
    config = ConfigParser.ConfigParser()
    config.read(configFile)
    file_paths = ConfigSectionMap(config,'file_paths')
    inputs = ConfigSectionMap(config,'inputs')
    return file_paths,inputs

def ConfigSectionMap(config,section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        dict1[option] = config.get(section, option)
    return dict1

##############################################################################
##  Read UH_BOX
##  Read the UH_BOX timeseries.  Save both the UH timeseries and the timestamps
##############################################################################
def load_uh(infile, verbose):
    """ 
    Loads UH from (infile) and returns a timeseries Unit Hydrograph (uh_t, uh)
    """ 
    if verbose:
        print 'Reading UH_Box from file: %s' % infile
    (uh_t,uh) = np.genfromtxt(infile, delimiter = ',', skip_header = 1, 
                              unpack = True)
    uh_t = uh_t.astype(int)
    return (uh_t, uh)

##############################################################################
##  Read Timestep from UH_BOX
##  In this version of the code, the timestep of the UH_BOX must match DELTA_T
##############################################################################
def find_TS(uh_t, verbose):
    """
    Determines the (INPUT_INTERVAL) based on the timestep given in UHfile
    """
    INPUT_INTERVAL = uh_t[1]-uh_t[0]
    if verbose:
        print 'Input Timestep = %i seconds' % INPUT_INTERVAL
    return INPUT_INTERVAL



##############################################################################
##  Read Direction
##  Determines downstream row/col numbers from flow direction input
##############################################################################
def read_direction(fdr, basin_ids, dy, dx, basin_id, NODATA, verbose):
    """
    Reads the direction file and makes two grids (to_x) and (to_y).
    The input grids follow the 1-8 or 1-128 grid directions as shown below.
    val = direction  [to_y][to_x]
    """
    if verbose:
        print 'Reading direction input and finding target row/columns'
        
    to_y = np.zeros(fdr.shape, dtype=int)
    to_x = np.zeros(fdr.shape, dtype=int)

    for (y, x), d in np.ndenumerate(fdr):
        try:
            to_y[y, x] = y+dy[d]
            to_x[y, x] = x+dx[d]
        except:
           to_y[y, x] = -9999
           to_x[y, x] = -9999
           
    return to_y, to_x

##############################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##############################################################################
def find_nearest(array, value):
    """ Find the index location in (array) with value nearest to (value)"""
    idx = (np.abs(array-value)).argmin()
    return idx

##############################################################################
## Seach Catchment
## Find all cells upstream of pour point.  Retrun a dictionary with x_inds, yinds,
## and #of cell to downstream pour point.  All are sorted the by the latter.
##############################################################################
def search_catchment(to_y, to_x, y_ind, x_ind, basin_ids, basin_id, verbose):
    """
    Find all cells upstream of pour point.  Retrun a dictionary with x_inds, 
    yinds, and #of cell to downstream pour point.  All are sorted the by the 
    latter. For each x,y pair, the flow path is followed until either the 
    catchment outlet is encountered (if (yy==y_ind and xx==x_ind):) or the 
    flowpath leads outside of grid.  
    Does not handle wrapped coordinates.  
    """
    if verbose:
        print 'Searching Catchment'
    
    COUNT = 0
    (len_y,len_x) = to_x.shape
    CATCH = {}
    CATCH['x_inds'] = np.array([], dtype=int)
    CATCH['y_inds'] = np.array([], dtype=int)
    CATCH['count_ds'] = np.array([], dtype=int)
    fractions = np.zeros((len_y, len_x))
    for (y, x), basin_num in np.ndenumerate(basin_ids):
        if basin_num == basin_id:
            yy, xx = y, x
            op, cells = 0,0
            while op == 0:
                if (yy == y_ind and xx == x_ind):
                    op = 1
                    CATCH['x_inds'] = np.append(CATCH['x_inds'], x)
                    CATCH['y_inds'] = np.append(CATCH['y_inds'], y)
                    CATCH['count_ds'] = np.append(CATCH['count_ds'], cells)
                    COUNT += 1
                    fractions[y, x] = 1.
                else:
                    yy, xx = to_y[yy, xx], to_x[yy, xx]
                    cells += 1
                    if (xx>len_x-1 or xx<0 or yy>len_y-1 or yy<0):
                        op =-1
    if verbose:
        print "Upstream grid cells from present station: %i" % COUNT
    
    # sort CATCH
    ii = np.argsort(CATCH['count_ds'])
    CATCH['count_ds'] = CATCH['count_ds'][ii]
    CATCH['x_inds'] = CATCH['x_inds'][ii]
    CATCH['y_inds'] = CATCH['y_inds'][ii]

    return (CATCH, fractions)

##############################################################################
##  MakeUH
##  Calculate impulse response function for grid cells using equation (15) 
##  from Lohmann, et al. (1996) Tellus article.  Return 3d UH grid.
##############################################################################
def make_UH(DELTA_T, T_Cell, y_inds, x_inds, velocity, diffusion, xmask, 
            PREC, verbose):
    """
    Calculate the impulse response function for grid cells using equation 15 
    from Lohmann, et al. (1996) Tellus article.  Return 3d UH grid.
    """
    if verbose:
        print 'Making UH for each cell'
    UH = np.zeros((T_Cell, xmask.shape[0], xmask.shape[1]))
    for (y, x) in zip(y_inds, x_inds):
        time = DELTA_T
        flag = 0
        t = 0
        green = np.zeros(T_Cell)
        while (t<T_Cell and flag==0):
            exponent = -1*np.power(velocity[y, x]*time-xmask[y, x], 2)/(4*diffusion[y, x]*time)
            if exponent > np.log(PREC):
                green[t] = xmask[y, x]/(2*time*np.sqrt(np.pi*time*diffusion[y, x]))*np.exp(exponent)
                t += 1
                time = time+DELTA_T
            else:
                flag = 1
        sum = np.sum(green)
        if sum>0.:
            UH[:, y, x] = green[:]/sum
    return UH

##############################################################################
##  Make RIVER UH
##  Calculate impulse response function for river routing
##  Steps upstream combining unit hydrographs
##############################################################################
def make_grid_UH_river(T_UH, T_Cell, UH,to_y, to_x, y_ind, x_ind, y_inds,
                       x_inds, count_ds, PREC, verbose):
    """
    Calculate impulse response function for river routing.  Starts at 
    downstream point incrementally moves upstream.
    """
    if verbose:
        print "Making UH_RIVER grid.... It takes a while..."
    UH_RIVER = np.zeros((T_UH, UH.shape[1], UH.shape[2]))
    for (y, x, d) in zip(y_inds, x_inds, count_ds):
        if d>0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            IRF_temp = np.zeros(T_UH+T_Cell)
            active_timesteps = np.nonzero(UH_RIVER[:, yy, xx]>PREC)[0]
            for t in active_timesteps:
                for l in xrange(T_Cell):
                    IRF_temp[t+l] = IRF_temp[t+l] + UH[l, y, x]*UH_RIVER[t, yy, xx]
            sum = np.sum(IRF_temp[:T_UH])
            if sum>0:
                UH_RIVER[:, y, x] = IRF_temp[:T_UH]/sum
        elif d==0:
            UH_RIVER[:T_Cell, y_ind, x_ind] = UH[:, y_ind, x_ind]
       
    return UH_RIVER

##############################################################################
## Make Grid UH
## Combines the UH_BOX with downstream cell UH_River IRF.
## Cell [0] is given the UH_Box without river routing
##############################################################################
def make_grid_UH(T_UH, T_Cell, UH_RIVER, UH_BOX, to_y, to_x, y_inds, x_inds, 
    count_ds, PREC, NODATA, verbose):
    """
    Combines the UH_BOX with downstream cell UH_RIVER.  Cell [0] is given the
    UH_Box without river routing
    """
    if verbose:
        print "Making UH_S grid"
    UH_S = np.zeros((T_UH,UH_RIVER.shape[1],UH_RIVER.shape[2]))+NODATA
    for (y, x, d) in zip(y_inds, x_inds, count_ds):
        IRF_temp = np.zeros(T_UH+T_Cell)
        if d > 0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            active_timesteps = np.nonzero(UH_RIVER[:, yy, xx]>PREC)[0]
            for t in active_timesteps:
                for l in xrange(len(UH_BOX)):
                    IRF_temp[t+l] = IRF_temp[t+l] + UH_BOX[l]*UH_RIVER[t, yy, xx]
            sum = np.sum(IRF_temp[:T_UH])
            if sum>0:
                UH_S[:, y, x] = IRF_temp[:T_UH]/sum
        else:
            IRF_temp[:len(UH_BOX)] = UH_BOX[:]
            UH_S[:, y, x] = IRF_temp[:T_UH]
    return UH_S

##############################################################################
## Aggregate to larger timestep
##############################################################################
def aggregate(UH_S, T_UH, INPUT_INTERVAL, OUTPUT_INTERVAL, x_inds, y_inds, 
              NODATA, verbose):
    """
    Aggregates to timestep (OUTPUT_INTERVAL).  OUTPUT_INTERVAL must be a 
    multiple of INPUT_INTERVAL.  This function is not setup to disaggregate 
    the Unit Hydrographs to a OUTPUT_INTERVAL<INPUT_INTERVAL.
    """
    if OUTPUT_INTERVAL == INPUT_INTERVAL:
        if verbose:
            print 'No need to aggregate (OUTPUT_INTERVAL = INPUT_INTERVAL)'
            print 'Skipping the aggregate step'
        UH_out = UH_S
    elif np.remainder(OUTPUT_INTERVAL,INPUT_INTERVAL) == 0:
        if verbose:
            print 'Aggregating to %i from %i seconds' \
            % (OUTPUT_INTERVAL, INPUT_INTERVAL)
        fac = int(OUTPUT_INTERVAL/INPUT_INTERVAL)
        T_UH_out = int(T_UH/fac)
        UH_out = np.zeros((T_UH_out,UH_S.shape[1], UH_S.shape[2]))+NODATA
        for (y, x) in zip(y_inds, x_inds):
            for t in xrange(T_UH_out):
                UH_out[t, y, x] = np.sum(UH_S[t*fac:t*fac+fac, y, x])
    else:
        if verbose:
            raise NameError("Not setup to disaggregate. \
                OUTPUT_INTERVAL must be a multiple of INPUT_INTERVAL")
    return UH_out

##############################################################################
##  Write output to netCDF
##  Writes out a netCDF3-64BIT data file containing the UH_S and fractions
##############################################################################
def write_netcdf(basin_x, basin_y, lons, lats, times, time_steps, time_res, UH_S, 
                 fractions, velocity, diffusion, basin_id, NODATA, verbose):
    """
    Write output to netCDF.  Writes out a netCDF4 data file containing the
    UH_S and fractions.
    """
    string = 'UH_%.3f_%.3f.nc' % (basin_x, basin_y)
    f = Dataset(string,'w', format = 'NETCDF4')

    # set dimensions
    time = f.createDimension('time', None)
    lon = f.createDimension('lon', (len(lons)))
    lat = f.createDimension('lat', (len(lats)))

    # initialize variables
    time = f.createVariable('time', 'f8', ('time'))
    time_step = f.createVariable('time_step', 'i8', ('time', ))
    lon = f.createVariable('lon', 'f8', ('lon', ))
    lat = f.createVariable('lat', 'f8', ('lat', ))
    fraction = f.createVariable('fraction', 'f8', ('lat', 'lon', ),
                                fill_value = NODATA)
    UHS = f.createVariable('unit_hydrograph', 'f8', ('time', 'lat', 'lon', ),
                            fill_value = NODATA)

    # write attributes for netcdf
    f.description = 'UH_S grid'
    f.created = tm.ctime(tm.time())
    f.history = ' '.join(sys.argv)
    f.source = sys.argv[0] # returns the name of script used
    f.velocity = velocity
    f.diffusion = diffusion
    f.outlet_lon = ('%.8f' % basin_x)
    f.outlet_lat = ('%.8f' % basin_y)
    f.global_basin_id = basin_id

    lat.long_name = 'latitude coordinate'
    lat.standard_name = 'latitude'
    lat.units = 'degrees_north'

    lon.long_name = 'longitude coordinate'
    lon.standard_name = 'longitude'
    lon.units = 'degrees_east'

    time.units = 'seconds since 0001-1-1 0:0:0'
    time.calendar = 'noleap'
    time.longname = 'time'
    time.type_prefered = 'int'
    time.description = 'Seconds since initial impulse'

    time_step.longname = 'timestep'
    time_step.type_prefered = 'int'
    time_step.description = 'timestep number'
    time_step.resolution = time_res

    UHS.units = 'unitless'
    UHS.description = 'unit hydrograph'
    
    fraction.description = 'fraction of grid cell contributing to outlet location'

    # write data to variables initialized above
    time_step[:]= time_steps
    time[:] = times
    lon[:] = lons
    lat[:] = lats
    UHS[:, :, :] = UH_S
    fraction[:, :]= fractions
    f.close()

    return string

##############################################################################
# Run Program
##############################################################################
if __name__ == "__main__":
    main()
