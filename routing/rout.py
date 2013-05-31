#!/usr/local/bin/python
"""
PROGRAM rout, Python-Version, written by Joe Hamman winter 2012/2013
Routing algorithm developed by D. Lohmann.
"""

####################################################################################
import sys
import numpy as np
import argparse
import ConfigParser
import time as tm
from netCDF4 import Dataset

##################################################################################
###############################  MAIN PROGRAM ####################################
##################################################################################
def main():
    (infile,UHfile,Plons,Plats,load_velocity,velocity,
     load_diffusion,diffusion,verbose,NODATA,CELL_FLOWTIME,
        BASIN_FLOWTIME,PREC,OUTPUT_INTERVAL,DAY_SECONDS)=process_command_line()

    for i,basin_x,basin_y in enumerate(zip(lons,lats)):
        rout(infile,UHfile,basin_x,basin_y,load_velocity,velocity,
             load_diffusion,diffusion,verbose,NODATA,CELL_FLOWTIME,
            BASIN_FLOWTIME,PREC,OUTPUT_INTERVAL,DAY_SECONDS)
        if verbose:
            print 'finished point #',i,'(',basin_x,',',basin_y,')'
    if verbose:
        print 'routing program finished.'
    return
        
def rout(infile,UHfile,basin_x,basin_y,load_velocity,velocity,
     load_diffusion,diffusion,verbose,NODATA,CELL_FLOWTIME,
        BASIN_FLOWTIME,PREC,OUTPUT_INTERVAL,DAY_SECONDS):

    # Load input arrays, store in dictionary.  (Format - Inputs['var'])
    Inputs = read_netcdf(infile,('Basin_ID','lon','lat'),verbose)

    # Find bounds of basin and basin_id
    (basin_id,x_min,x_max,y_min,y_max) = read_global_inputs(basin_x,basin_y,Inputs['lon'],Inputs['lat'],
                                                            Inputs['Basin_ID'],verbose)
    
    # Load input arrays, store in python list.  (Format -Basin['var'])
    vars =['Basin_ID','Flow_Direction','Flow_Distance','lon','lat']
    if load_velocity:
        vars.append('velocity')
    if load_diffusion:
        vars.append('diffusion')
    Basin = clip_netcdf(infile,vars,x_min,x_max,y_min,y_max,
                        load_velocity,velocity,load_diffusion,diffusion,verbose)
    
    # Load UH_BOX input
    (uh_t,UH_Box) = load_uh(UHfile,verbose)
    
    # Find timestep (timestep is determined from UH_BOX input file)
    INPUT_INTERVAL = find_TS(uh_t,verbose)
    T_Cell = CELL_FLOWTIME*DAY_SECONDS/INPUT_INTERVAL
    T_UH = BASIN_FLOWTIME*DAY_SECONDS/INPUT_INTERVAL
    
    # Read direction grid and find to_col (to_x) and to_row (to_y)
    (to_x,to_y) = read_direction(Basin['Flow_Direction'],Basin['Basin_ID'],
                                 basin_id,NODATA,verbose)
    
    # Find row/column indicies of lat/lon inputs
    x_ind = find_nearest(Basin['lon'],basin_x)
    y_ind = find_nearest(Basin['lat'],basin_y)
    
    # Find all grid cells upstream of pour point
    (Catchment, fractions) = search_catchment(to_x,to_y,x_ind,y_ind,
                                              Basin['Basin_ID'],basin_id,verbose)
    
    # Make UH for each grid cell upstream of basin pour point (linear routing model - Saint-Venant equation)
    UH = make_UH(INPUT_INTERVAL,T_Cell,
                 Catchment['x_inds'],Catchment['y_inds'],Basin['Velocity'],
        Basin['Diffusion'],Basin['Flow_Distance'],PREC,verbose)
    
    # Make UH_RIVER by incrementally moving upstream comining UH functions
    UH_RIVER = make_grid_UH_river(T_UH,T_Cell,UH,to_x,to_y,x_ind,y_ind,
                                  Catchment['x_inds'], Catchment['y_inds'], Catchment['count_ds'],PREC,verbose)
    
    # Make UH_S for each grid cell upstream of basin pour point (combine IRFs for all grid cells in flow path)
    UH_S = make_grid_UH(T_UH,T_Cell,UH_RIVER,UH_Box,to_x,to_y,
                        Catchment['x_inds'], Catchment['y_inds'], Catchment['count_ds'],PREC,NODATA,verbose)
    
    # Agregate to output timestep
    UH_out = aggregate(UH_S,T_UH,INPUT_INTERVAL,OUTPUT_INTERVAL,
                       Catchment['x_inds'], Catchment['y_inds'],NODATA,verbose)
    
    #Write to output netcdf
    times = np.linspace(0,OUTPUT_INTERVAL*UH_out.shape[0],UH_out.shape[0],endpoint=False)
    write_netcdf(basin_x,basin_y,Basin['lon'],Basin['lat'],
                 times,UH_out,fractions,velocity,diffusion,basin_id,NODATA,verbose)

    return
    
##################################################################################
############### Routines #########################################################
##################################################################################
def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h

    A configuration file may be provided with -C configFile
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-C","--configFile", type=str, help="Input configuration file",default=False)
    parser.add_argument("-i","--infile", type=str, help="Input netCDF containing all input grids")
    parser.add_argument("-UH","--UHfile", type=str, help="Input UH_BOX hydrograph")
    parser.add_argument("-lon","--longitude", type=float, nargs='*',help="Input longitude of point to route to")
    parser.add_argument("-lat","--latitude", type=float,nargs='*',help="Input latitude of point to route to")
    parser.add_argument("-vel","--velocity", type=float, help="Input wave velocity, may be variable in infile or scalar")
    parser.add_argument("-diff","--diffusion", type=float, help="Input diffusion")
    parser.add_argument("-v", "--verbose", help="Increase output verbosity", action="store_true")
    parser.add_argument("--NODATA",type=float, help="Input integer to bedefined as NODATA value",default=9.96920996839e+36)
    parser.add_argument("--CELL_FLOWTIME",type=int, help="Input integer to be defined as max number of days for flow to pass through cell",default=2)
    parser.add_argument("--BASIN_FLOWTIME",type=int, help="Input integer to be defined as max number of days for flow to pass through basin",default=50)
    parser.add_argument("-PREC","--Precision",type=int, help="Input integer to be defined asminimum precision for calculations",default =1e-30)
    parser.add_argument("-TS","--OUTPUT_INTERVAL",type=int, help="Output timestep in seconds for Unit Hydrographs",default=86400)
    parser.add_argument("--DAY_SECONDS",type=int, help="Seconds per day",default=86400)
    args = parser.parse_args()

    # Assign values
    if args.configFile:
        file_paths,inputs = process_config(args.configFile)
    if args.infile:
        infile = args.infile
    else:
        try:
            infile = file_paths['infile']
        except:
            raise IOError('Need input file from command line or configuration file')

    if args.UHfile:
        UHfile = args.UHfile
    else:
        try:
            UHfile = file_paths['uhfile']
        except:
            raise IOError('Need input Unit Hydrograph File from command line or configuration file')

    if args.longitude:
        Plons = args.longitude
    else:
        try:
            Plons = map(float,inputs['longitude'])
        except:
            raise IOError('Need logitude(s) from command line or configuration file')
        
    if args.latitude:
        Plats = args.latitude
    else:
        try:
            Plats = map(float,inputs['latitude'])
        except:
            raise IOError('Need latitude(s) from command line or configuration file')        
    
    try:
        velocity = float(inputs['velocity'])
        load_velocity = False
    except:
        if not args.velocitye:
            load_velocity=True
            velocity = 'From Grid Inputs'
        else:
            load_velocity=False
            velocity = args.velocity
    try:
        diffusion = float(inputs['diffusion'])
        load_diffusion = False
    except:
        if not args.diffusion:
            load_diffusion = True
            diffusion = 'From Grid Inputs'
        else:
            load_diffusion = False
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
        
    return (infile,UHfile,Plons,Plats,load_velocity,velocity,load_diffusion,
            diffusion,verbose,NODATA,CELL_FLOWTIME,BASIN_FLOWTIME,PREC,OUTPUT_INTERVAL,DAY_SECONDS)

##################################################################################
##  Process Configuration File
##  Values that aren't provided will be given a False bool and filled in by the
##  command line processor
##################################################################################
def process_config(configFile):
    """
    Parse arguments and assign flags for further loading of files.  Configuration
    flag must be raised on command line (-C) to configuration file must be provided.
    Usage:  rout.py -C rout.cfg
    There is no set of variables that must be included, any variable not given in
    configuration file will be filled in during the rest of the command line parsing.
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

##################################################################################
##  Read netCDF Inputs
##  Read data from input netCDF.  Input netCDF should include the following vars:
##  ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat')
##################################################################################
def read_netcdf(nc_str,vars,verbose):
    """
    Read data from input netCDF.  Input netCDF should include the following vars:
    ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat') and
    may also include ('Velocity') and/or ('Diffusion')
    """
    if verbose:
        print 'Reading input data vars:', vars, 'from file:',nc_str
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
    f.close()
    return d

##################################################################################
## Find Basin Dims and ID
##################################################################################
def read_global_inputs(basin_x,basin_y,lons,lats,basins,verbose):
    """Reads input lons/lats/basins_ids and returns basin bounds."""
    if verbose:
        print 'reading global inputs'
    basin_id = basins[find_nearest(lats,basin_y),find_nearest(lons,basin_x)]
    if verbose:
        print 'input lon:', basin_x
        print 'input lat:', basin_y
        print 'basid id:',basin_id
    inds = np.nonzero(basins==basin_id)
    x,y = np.meshgrid(np.arange(len(lons)), np.arange(len(lats)))
    x_inds = x[inds]
    y_inds = y[inds]
    x_min = min(x_inds)
    x_max = max(x_inds)
    y_min = min(y_inds)
    y_max = max(y_inds)
    return (basin_id,x_min,x_max,y_min,y_max)

##################################################################################
##  Clip netCDF Inputs
##  Read data from input netCDF.  Input netCDF should include the following vars:
##  ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat')
##  Velocity and Diffusion are assumed to be constant in this version
##################################################################################
def clip_netcdf(nc_str,vars,x_min,x_max,y_min,y_max,
                load_velocity,velocity,load_diffusion,diffusion,verbose):
    """
    Read data from input netCDF.  Clip to basin sizde based on outputs from
    read_global_inputs.  Input netCDF should include the following vars:
    ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat') and
    may also include ('Velocity') and/or ('Diffusion')
    """
    if verbose:
        print 'Reading input data vars:', vars
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        temp = f.variables[var][:]
        try:
            if np.rank(temp) > 1:
                d[var] = temp[y_min:y_max+1,x_min:x_max+1]
            elif var == 'lon':
                d['lon'] = temp[x_min:x_max+1]
            elif var == 'lat':
                d['lat'] = temp[y_min:y_max+1]
        except NameError:
            print 'Unable to clip ', var, '. Confirm that the var exists in ', nc_str
            raise
    if not load_velocity:
        d['Velocity']=np.zeros((d['Flow_Direction'].shape))+velocity
    if not load_diffusion:
        d['Diffusion']=np.zeros((d['Flow_Direction'].shape))+diffusion
    f.close()
    if verbose:
        print 'grid cells in subset: ', len(d['lon'])*len(d['lat'])
    return d

##################################################################################
##  Read UH_BOX
##  Read the UH_BOX timeseries.  Save both the UH timeseries and the timestamps
##################################################################################
def load_uh(infile,verbose):
    """ Loads UH from (infile) and returns a timeseries Unit Hydrograph (uh_t, uh)""" 
    if verbose:
        print 'Reading UH_Box from file: ', infile
    (uh_t,uh) = np.genfromtxt(infile, delimiter=',', skip_header=1, unpack=True)
    uh_t = uh_t.astype(int)
    return (uh_t, uh)

##################################################################################
##  Read Timestep from UH_BOX
##  In this version of the code, the timestep of the UH_BOX must match DELTA_T
##################################################################################
def find_TS(uh_t,verbose):
    """Determines the (INPUT_INTERVAL) based on the timestep given in UHfile"""
    INPUT_INTERVAL = uh_t[1]-uh_t[0]
    if verbose:
        print 'Input Timestep = '+ str(INPUT_INTERVAL) + ' seconds'
    return INPUT_INTERVAL

##################################################################################
##  Read Direction
##  Determines downstream row/col numbers from flow direction input
##################################################################################
def read_direction(fdr,basin_ids,basin_id,NODATA,verbose):
    """
    Reads the direction file and makes two grids (to_x) and (to_y).
    The input grids follow the 1-8 grid directions as shown below.
    val = direction  [to_y][to_x]
      1 = north      [-1][0]
      2 = northeast  [-1][+1]
      3 = east       [0][+1]
      4 = southeast  [+1][+1]
      5 = south      [+1][0]
      6 = southwest  [+1][-1]
      7 = west       [0][-1]
      8 = northwest  [-1][-1] 
    """
    if verbose:
        print 'Reading Direction input and finding target row/columns'
    to_x = np.zeros((fdr.shape[0],fdr.shape[1]),dtype=int)
    to_y = np.zeros((fdr.shape[0],fdr.shape[1]),dtype=int)
    (len_y,len_x) = to_x.shape
    for y in xrange(len_y):
        for x in xrange(len_x):
            if basin_ids[y,x]==basin_id:
                if fdr[y,x] == 1:
                    to_x[y,x] = x
                    to_y[y,x] = y-1
                elif fdr[y,x] == 2:
                    to_x[y,x] = x+1
                    to_y[y,x] = y-1
                elif fdr[y,x] == 3:
                    to_x[y,x] = x+1
                    to_y[y,x] = y
                elif fdr[y,x] == 4:
                    to_x[y,x] = x+1
                    to_y[y,x] = y+1
                elif fdr[y,x] == 5:
                    to_x[y,x] = x
                    to_y[y,x] = y+1
                elif fdr[y,x] == 6:
                    to_x[y,x] = x-1
                    to_y[y,x] = y+1
                elif fdr[y,x] == 7:
                    to_x[y,x] = x-1
                    to_y[y,x] = y
                elif fdr[y,x] == 8:
                    to_x[y,x] = x-1
                    to_y[y,x] = y-1
                else:
                    to_x[y,x] = -9999
                    to_y[y,x] = -9999
    return (to_x,to_y)

##################################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##################################################################################
def find_nearest(array,value):
    """ Find the index location in (array) with value nearest to (value)"""
    idx = (np.abs(array-value)).argmin()
    return idx
##################################################################################
## Seach Catchment
## Find all cells upstream of pour point.  Retrun a dictionary with x_inds, yinds,
## and #of cell to downstream pour point.  All are sorted the by the latter.
##################################################################################
def search_catchment(to_x,to_y,x_ind,y_ind,basin_ids,basin_id,verbose):
    """
    Find all cells upstream of pour point.  Retrun a dictionary with x_inds, yinds,
    and #of cell to downstream pour point.  All are sorted the by the latter.
    For each x,y pair, the flow path is followed until either the catchment outlet
    is encountered (if (yy==y_ind and xx==x_ind):) or the flowpath leads outside of
    grid (if (xx>len(lons)-1 or xx<0 or yy>len(lats)-1 or yy<0):).
    """
    if verbose:
        print 'Searching Catchment'
    COUNT = 0
    x_inds = np.zeros(to_x.shape[0]*to_x.shape[1],dtype=int)
    y_inds = np.zeros(to_x.shape[0]*to_x.shape[1],dtype=int)
    count_ds = np.zeros(to_x.shape[0]*to_x.shape[1],dtype=int)
    fractions = np.zeros((to_x.shape[0],to_x.shape[1]))
    (len_y,len_x) = to_x.shape
    for y in xrange(len_y):
        for x in xrange(len_x):
            if basin_ids[y,x]==basin_id:
                yy,xx = y,x
                op,cells = 0,0
                while op == 0:
                    if (yy==y_ind and xx==x_ind):
                        op = 1
                        x_inds[COUNT] = x
                        y_inds[COUNT] = y
                        count_ds[COUNT] = cells
                        COUNT += 1
                        fractions[y,x] = 1.
                    elif (to_x[yy,xx] > -1 and to_y[yy,xx] > -1):
                        (xx, yy) = (to_x[yy][xx], to_y[yy,xx])
                        cells += 1
                        if (xx>len_x-1 or xx<0 or yy>len_y-1 or yy<0):
                            op =-1
                    else:
                        op = -1
    if verbose:
        print "Upstream grid cells from present station: ", COUNT
    # Save to sorted dictionary
    CATCH = {}
    ii = np.argsort(count_ds[:COUNT])
    x_inds = x_inds[:COUNT]
    y_inds = y_inds[:COUNT]
    CATCH['x_inds'] = x_inds[ii]
    CATCH['y_inds'] = y_inds[ii]
    CATCH['count_ds'] = count_ds[ii]
    return (CATCH, fractions)

##################################################################################
##  MakeUH
##  Calculate impulse response function for grid cells using equation (15) from
##  Lohmann, et al. (1996) Tellus article.  Return 3d UH grid.
###################################################################################
def make_UH(DELTA_T,T_Cell, x_inds,y_inds,velocity,diffusion,xmask,PREC,verbose):
    """
    Calculate the impulse response function for grid cells using equation 15 from
    Lohmann, et al. (1996) Tellus article.  Return 3d UH grid.
    """
    if verbose:
        print 'Making UH for each cell'
    UH = np.zeros((T_Cell, xmask.shape[0],xmask.shape[1]))
    for i in xrange(len(x_inds)):
        x = x_inds[i]
        y = y_inds[i]
        time = DELTA_T
        flag = t = 0
        green = np.zeros(T_Cell)
        while (t<T_Cell and flag==0):
            exponent = -1*np.power(velocity[y,x]*time-xmask[y,x],2)/(4*diffusion[y,x]*time)
            if exponent > np.log(PREC):
                green[t] = xmask[y,x]/(2*time*np.sqrt(np.pi*time*diffusion[y,x]))*np.exp(exponent)
                t += 1
                time = time+DELTA_T
            else:
                flag=1
        sum = np.sum(green)
        if sum>0.:
            UH[:,y,x] = green[:]/sum
    return UH

##################################################################################
##  Make RIVER UH
##  Calculate impulse response function for river routing
##  Steps upstream combining unit hydrographs
###################################################################################
def make_grid_UH_river(T_UH,T_Cell,UH,to_x,to_y,x_ind,y_ind,
                       x_inds,y_inds,count_ds,PREC,verbose):
    """
    Calculate impulse response function for river routing.  Starts at downstream
    point incrementally moves upstream.
    """
    if verbose:
        print "Making UH_RIVER grid.... It takes a while..."
    UH_RIVER = np.zeros((T_UH,UH.shape[1],UH.shape[2]))
    for i in xrange(len(x_inds)):
        d = count_ds[i]
        if d>0:
            x = x_inds[i]
            y = y_inds[i]
            yy = to_y[y,x]
            xx = to_x[y,x]
            IRF_temp = np.zeros(T_UH+T_Cell)
            active_timesteps = np.nonzero(UH_RIVER[:,yy,xx]>PREC)[0]
            for t in active_timesteps:
                for l in xrange(T_Cell):
                    IRF_temp[t+l] = IRF_temp[t+l] + UH[l,y,x]*UH_RIVER[t,yy,xx]
            sum = np.sum(IRF_temp[:T_UH])
            if sum>0:
                UH_RIVER[:,y,x] = IRF_temp[:T_UH]/sum
        elif d==0:
            UH_RIVER[:T_Cell,y_ind,x_ind] = UH[:,y_ind,x_ind]
        active_timesteps=[]
       
    return UH_RIVER

##################################################################################
## Make Grid UH
## Combines the UH_BOX with downstream cell UH_River IRF.
## Cell [0] is given the UH_Box without river routing
##################################################################################
def make_grid_UH(T_UH,T_Cell,UH_RIVER,UH_BOX,to_x,to_y,x_inds,y_inds,count_ds,PREC,NODATA,verbose):
    """
    Combines the UH_BOX with downstream cell UH_RIVER.  Cell [0] is given the
    UH_Box without river routing
    """
    if verbose:
        print "Making UH_S grid"
    UH_S = np.zeros((T_UH,UH_RIVER.shape[1],UH_RIVER.shape[2]))+NODATA
    for i in xrange(len(x_inds)):
        x = x_inds[i]
        y = y_inds[i]
        d = count_ds[i]
        IRF_temp = np.zeros(T_UH+T_Cell)
        if d > 0:
            yy = to_y[y,x]
            xx = to_x[y,x]
            active_timesteps = np.nonzero(UH_RIVER[:,yy,xx]>PREC)[0]
            for t in active_timesteps:
                for l in xrange(len(UH_BOX)):
                    IRF_temp[t+l] = IRF_temp[t+l] + UH_BOX[l]*UH_RIVER[t,yy,xx]
            sum = np.sum(IRF_temp[:T_UH])
            if sum>0:
                UH_S[:,y,x] = IRF_temp[:T_UH]/sum
        elif d==0:
            IRF_temp[:len(UH_BOX)] = UH_BOX[:]
            UH_S[:,y,x] = IRF_temp[:T_UH]
        active_timesteps=[]
    return UH_S

##################################################################################
## Aggregate to larger timestep
##################################################################################
def aggregate(UH_S,T_UH, INPUT_INTERVAL,OUTPUT_INTERVAL,x_inds,y_inds,NODATA,verbose):
    """
    Aggregates to timestep (OUTPUT_INTERVAL).  OUTPUT_INTERVAL must be a multiple of
    INPUT_INTERVAL.  This function is not setup to disaggregate the Unit Hydrographs
    to a OUTPUT_INTERVAL<INPUT_INTERVAL.
    """
    if OUTPUT_INTERVAL==INPUT_INTERVAL:
        if verbose:
            print 'no need to aggregate, skipping this step (OUTPUT_INTERVAL = INPUT_INTERVAL)'
        UH_out = UH_S
    elif np.remainder(OUTPUT_INTERVAL,INPUT_INTERVAL)==0:
        if verbose:
            print 'aggregating to ', OUTPUT_INTERVAL, 'from ', INPUT_INTERVAL, 'seconds'
        fac = int(OUTPUT_INTERVAL/INPUT_INTERVAL)
        T_UH_out = int(T_UH/fac)
        UH_out = np.zeros((T_UH_out,UH_S.shape[1],UH_S.shape[2]))+NODATA
        for i in xrange(len(x_inds)):
            x = x_inds[i]
            y = y_inds[i]
            for t in xrange(T_UH_out):
                UH_out[t,y,x] = np.sum(UH_S[t*fac:t*fac+fac,y,x])
    else:
        if verbose:
            raise NameError("Not setup to disaggregate. OUTPUT_INTERVAL must be a multiple of INPUT_INTERVAL")
    return UH_out

##################################################################################
##  Write output to netCDF
##  Writes out a netCDF3-64BIT data file containing the UH_S and fractions
##################################################################################
def write_netcdf(basin_x,basin_y,lons,lats,times,UH_S,fractions,
                 velocity,diffusion,basin_id,NODATA,verbose):
    """
    Write output to netCDF.  Writes out a netCDF4-64BIT data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    string = 'UH_'+('%.3f' % basin_x)+'_'+('%.3f' % basin_y)+'.nc'
    f = Dataset(string,'w', format='NETCDF4_CLASSIC')

    # set dimensions
    time = f.createDimension('time', None)
    lon = f.createDimension('lon', (len(lons)))
    lat = f.createDimension('lat', (len(lats)))

    # initialize variables
    time = f.createVariable('time','f8',('time',))
    lon = f.createVariable('lon','f8',('lon',))
    lat = f.createVariable('lat','f8',('lat',))
    fraction = f.createVariable('fraction','f8',('lat','lon',),fill_value=NODATA)
    UHS = f.createVariable('unit_hydrograph','f8',('time','lat','lon',),fill_value=NODATA)

    # write attributes for netcdf
    f.description = 'UH_S for point '+('%.8f' % basin_x)+', '+('%.8f' % basin_y)
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0] # prints the name of script used
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

    time.units = 'seconds since 0000-1-1 0:0:0'
    time.calendar = 'noleap'
    time.longname = 'time'
    time.type_prefered = 'int'
    time.description = 'Seconds since initial impulse'

    UHS.units = 'unitless'
    UHS.description = 'unit hydrograph for each grid cell with respect to downstream grid location'
    
    fraction.description = 'fraction of grid cell contributing to guage location'

    # write data to variables initialized above
    time[:]= times
    lon[:] = lons
    lat[:] = lats
    UHS[:,:,:] = UH_S
    fraction[:,:]= fractions
    f.close()
    if verbose:
        print 'wrote UH_S netCDF for ' + ('%.3f' % basin_x)+'_'+('%.3f' % basin_y)

##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
