#!/usr/bin/python
"""
Things left to do...
0.  Handle timestamps from flux/uh file 
1.  Check timesteps
2.  Check that lat/lon matches in statefile (could add method to find the point even if its not in order)
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import time as tm
import argparse
import os
import glob

def main():

    gridFile,fluxFiles,uhFiles,stateFiles,outFile,var_names,options = process_command_line()

    convolve(gridFile,fluxFiles,uhFiles,stateFiles,outFile,var_names,options)

def convolve(gridFile,fluxFiles,uhFiles,stateFiles,outFile,var_names,options,re = 6.37122e6):
    
    #get grid data
    if options['verbose']:
        print 'opening grid data'
    GridData,GridAttrs,GridGlobs = read_netcdf(gridFile,verbose=options['verbose'])
    area = GridData[var_names['area_var']]*re*re

    if stateFiles:
        # Load State Files
        stateDict = read_state_files(stateFiles)
    else:
        # Make stateDict from uhFiles (grids)
        stateDict = make_state_dict()
    
    lons = np.zeros(len(stateDict))
    lats = np.zeros(len(stateDict))
    data = {}
    times = {}
    inds = {}
    count = {}
    
    # Loop over all flux files 
    for i,fluxFile in enumerate(fluxFiles):
        read_vars = var_names['flux_vars'][:]
        read_vars.append(var_names['flux_time_var'])
        FluxData,FluxAttrs,FluxGlobs = read_netcdf(fluxFile,vars=read_vars,verbose=options['verbose'])
        time = FluxData[var_names['flux_time_var']]
        flux = []
        for j,var in enumerate(var_names['flux_vars']):
            print FluxData[var][:].shape
            if j==0:
                flux = FluxData[var][:]
            else:
                flux += FluxData[var][:]
        
        # Loop over all UH Grids
        for j,pointDict in stateDict.iteritems():
            if options['verbose']:
                print 'starting point:',j+1, 'of',len(stateDict), '(flux file #',i+1,'of',len(fluxFiles),')'
            if i ==0:
                lons[j]=pointDict['outlet_lon']
                lats[j]=pointDict['outlet_lat']
            if j==0:
                luh = len(pointDict['unit_hydrograph'])
                td = time[-1]-time[-2]
                lt = len(time)+luh-1
                data[i]=np.zeros((len(flux)+luh-1,len(stateDict)))
                times[i] = np.arange(time[0],time[0]+lt*td,td)
                if i==0:
                    start = time[0]
                    end = time[0]+lt*td
                else:
                    start = np.minimum(start,time[0])
                    end = np.maximum(end,time[0]+lt*td)

            data[i][:,j] = convolution(flux,pointDict['unit_hydrograph'],pointDict['fraction'],pointDict['yi'],pointDict['xi'],GridData['mask'],area)
            if data[i][:,j].max()>10e10 and options['debug']:
                print '+++++++++++++++++++++++++++'
                print i,j,data[i][:,j].max()
                print pointDict['fraction'].shape
                print pointDict['fraction'], pointDict['fraction'].min(), pointDict['fraction'].max()
                print pointDict['yi'],pointDict['xi']
                print flux[:,pointDict['yi'],pointDict['xi']].max(),flux[:,pointDict['yi'],pointDict['xi']].min()
                print pointDict['file_name']
                print '==========================='
                tm.sleep(5)

    # Make a empty table to put the flows in
    print start,end, td, end-start, 'days'
    outTimes = np.arange(start,end+td,td)
    outFlows = np.zeros((len(outTimes),len(stateDict)))

    # combine the dictionary values from each flux file
    for i in xrange(len(fluxFiles)):
        s = times[i][0]; e = times[i][-1] #get the start and end values
        si = find_nearest(outTimes,s); ei=find_nearest(outTimes,e)+1 # get the start and end inds
        outFlows[si:ei,:] += data[i]
        print outFlows.shape, s,si,e,ei,data[i].shape
    
    # Adjust the flow units
    outFlows[:,:] /= options['flow_div']
    
    ## Write to NC file ##
    print 'initialized ', outFile
    f = Dataset(outFile, 'w', format='NETCDF4_CLASSIC')

    # set dimensions
    time = f.createDimension('time', None)
    points = f.createDimension('point', len(stateDict))
    
    # initialize variables
    time = f.createVariable('time','f8',('time',))
    lon = f.createVariable('longitudes','f8',('point',))
    lat = f.createVariable('latitudes','f8',('point',))
    flow = f.createVariable('Streamflow','f8',('time','point',))
    
    # write attributes for netcdf
    f.description = 'Streamflow.'
    f.history = 'Created ' + tm.ctime(tm.time())
    f.source = 'Streamflow convolution program - convolve_offline.py'

    lat.standard_name = "latitude"
    lat.long_name = "latitude of grid cell center"
    lat.units = "degrees_north"
    
    lon.standard_name = "longitude"
    lon.long_name = "longitude of grid cell center"
    lon.units = "degrees_east"

    time.units = FluxAttrs['time']['units']
    time.calendar = FluxAttrs['time']['calendar']
    time.longname = 'time'
    time.type_prefered = 'double'
    
    flow.description = 'Streamflow'
    flow.units = 'm^3/s'
    
    flow[:,:] = outFlows
    time[:] = outTimes
    lon[:] = lons
    lat[:] = lats
    
    f.close()

    return

def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type convolve.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--gridFile", type=str, help="Input grid netcdf including area variable")
    parser.add_argument("--fluxFiles", type=str, help="Input grid(s) netcdf including runoff/baseflow fluxes",nargs='+')
    parser.add_argument("--uhFiles", type=str, help="Input unit(s) hydrograph grid netcdf",nargs='+')
    parser.add_argument("--stateFiles", type=str, help="State files including index/uh/fraction variables",nargs='+')    
    parser.add_argument("--outFile", type=str, help="Input unit(s) hydrograph grid netcdf",default='out.nc')
    parser.add_argument("--area_var", type=str, help="variable name for area in the grid_nc",default='area')
    parser.add_argument("--flux_vars", type=str, help="If only one var should be used in the convolution, input here",default=['Runoff','Baseflow'],nargs='+')
    parser.add_argument("--flux_time_var", type=str, help="Time variable for ",default='time')
    parser.add_argument("--uh_var",type=str, help="variable name for unit hydrograph",default='unit_hydrograph')
    parser.add_argument("--uh_time_var", type=str, help="Time variable for ",default='time')
    parser.add_argument("--frac_var",type=str, help="variable name for unit hydrograph",default='fraction')
    parser.add_argument("--flow_div",type=float, help="multiplier value to convert to m3/s",default=86400*1000)
    parser.add_argument("--verbose",help="Make script verbose",action="store_true")
    parser.add_argument("--clean",help="Clean up temporary files", action='store_true')
    parser.add_argument("--dryrun",help="Check all files to have appropriate variables and lengths", action='store_true')
    parser.add_argument("--debug",help="Turn on Debugging", action='store_true')
    args = parser.parse_args()

    gridFile = args.gridFile
 
    try:
        stateFiles = glob.glob(args.stateFiles[0])
        uhFiles = False
        print 'found',len(stateFiles), 'state files'
    except:
        uhFiles = glob.glob(args.uhFiles[0])
        stateFiles = False
        print 'found',len(uhFiles), 'UH files'

    fluxFiles = glob.glob(args.fluxFiles[0])
    print 'found',len(fluxFiles),'Flux Files'
    outFile = args.outFile
    print 'will write out to:',outFile
    
    var_names = {}
    var_names['area_var'] = args.area_var
    var_names['flux_vars'] = args.flux_vars
    var_names['uh_var'] = args.uh_var
    var_names['frac_var'] = args.frac_var
    var_names['uh_time_var'] = args.uh_time_var
    var_names['flux_time_var'] = args.flux_time_var
    
    options = {}
    options['verbose'] = args.verbose
    options['clean'] = args.clean
    options['dryrun'] = args.dryrun
    options['debug'] = args.debug
    options['flow_div'] = args.flow_div
    
    return gridFile,fluxFiles,uhFiles,stateFiles,outFile,var_names,options

##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile,vars = [],coords = False, verbose = False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """
    if verbose:
        print 'Reading input data vars:', vars, ', from file:',ncFile
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    d={}
    a={}
    g={}
    if coords:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][coords]
            a[vars] = f.variables[vars].__dict__
        else:
            for var in vars:
                d[var] = f.variables[var][coords]
                a[var] = f.variables[var].__dict__
    else:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][:]
            a[vars] = f.variables[vars].__dict__
        else:
            for var in vars:
                d[var] = f.variables[var][:]
                a[var] = f.variables[var].__dict__
    
    for attr in f.ncattrs():
        g[attr] = getattr(f,attr)
    
    f.close()
    
    return d,a,g

def convolution(flux,uh,frac,yi,xi,mask,area):   
    """
    Do the convolution by looping over all uh grid cells and flux timesteps.
    """
    q = np.zeros((len(flux)+len(uh)-1))
    for i in xrange(len(xi)):
        y=yi[i];x=xi[i]
        if mask[y,x]==1:
            q += np.convolve(flux[:,y,x],uh[:,i])*area[y,x]*frac[i]
    return q

def find_nearest(array,value):
    """ Find the index location in (array) with value nearest to (value)"""
    idx = (np.abs(array-value)).argmin()
    return idx

def read_state_files(state_files):
    """
    Make a dictionary of points and within each dictionary store a dictionary containing
    the state/uh/fraction arrays
    """
    state_dict = {}
    fails = []
    j = 0
    print 'reading',len(state_files), 'flattened state files and storing all the data in a nice dictionary'

    for i,sf in enumerate(state_files):
        try:
            state_dict[j],a,g = read_netcdf(sf)
            state_dict[j]['outlet_lon']=g['outlet_lon']
            state_dict[j]['outlet_lat']=g['outlet_lat']
            state_dict[j]['file_name']=sf
            j += 1
        except:
            fails.append(sf)

    print 'opened all the state files'
    print 'there was a problem with these files: %i' %len(fails)
    print fails

    return state_dict

def make_state_dict(uh_files):
    """
    In the abscence of state files, we can flatten all the grids once,
    and store in a python dictionary
    (this will help with speed for long convolutions)
    """
    state_dict = {}
    print 'reading',len(uh_files), 'unit hydrograph grid files and storing all the grids in a nice dictionary'
    for i,uf in enumerate(uh_files):
        inData,inAttrs,inGlobs = read_netcdf(uf)
        
        # Find inds and flatten uh time-series
        y,x = np.nonzero(inData['fraction'])
        state_dict[i]['unit_hydrograph'] = inData['unit_hydrograph'][:,y,x]
        state_dict[i]['fraction'] = inData['fraction'][y,x] 
        state_dict[i]['time'] = inData['time']
        state_dict[i]['state_streamflow'] = np.zeros(len(time))
        state_dict[i]['xc'] = inData['xc'][y,x]
        state_dict[i]['yc'] = inData['yc'][y,x]
        state_dict[i]['xi'] = x
        state_dict[i]['yi'] = y
        state_dict[i]['file_name']=uf
        
    return state_dict
##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
