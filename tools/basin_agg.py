#!/usr/local/bin/python
"""
Aggregate high resolution Unit Hydrograph NetCDF Grids

Joe Hamman March, 2013
"""
import os
import sys
import numpy as np
import numpy.ma as ma
from collections import defaultdict
from netCDF4 import Dataset
from cdo import *
import time as tm
import argparse
cdo = Cdo()


##################################################################################

def main():

    Rvars,paths,options = process_command_line()
    Mvars = ('src_address','dst_address','src_grid_center_lat','dst_grid_center_lat','src_grid_center_lon','dst_grid_center_lon')
    Cvars = ('src_grid_center_lat','dst_grid_center_lat','src_grid_center_lon','dst_grid_center_lon')

    print paths
    #Load Inputs
    Inputs = read_netcdf(paths['mapFile'],Mvars,options['verbose'])
    
    #Convert lats and lons for src/dst to decimal degrees
    Coords = make_degrees(Cvars,Inputs,options['verbose'])
    
    #Find Unique Destination grid cells
    d = defaultdict(list)
    if options['verbose']:
        print 'Finding addresses now...'
    for i in xrange(len(Inputs['src_address'])):
        # parse links
        di = Inputs['dst_address'][i]-1
        si = Inputs['src_address'][i]-1
        key = (Coords['src_grid_center_lon'][si]-180, Coords['src_grid_center_lat'][si])
        coord =  (Coords['dst_grid_center_lon'][di]-180, Coords['dst_grid_center_lat'][di])
        d[key].append(coord)

    # Aggregate all basins in each key
    count = 0
    if options['verbose']:
        print 'Aggregating now...'
    for i,key in enumerate(d):
        flag  = 0
        Flist = []
        aggFile = filename(options['outPrefix'],key)
        for j,loc in enumerate(d[key]):
            inFile = filename(options['inPrefix'],loc)
            if os.path.isfile(paths['srcDir']+inFile):
                if options['verbose']:
                    print 'just opened', paths['srcDir']+inFile
                Flist.append(inFile)
                f = Dataset(paths['srcDir']+inFile,'r')
                if flag == 0:
                    aggData={}
                    for var in Rvars:
                        aggData[var] = f.variables[var][:]
                    if count == 0:
                        velocity = f.velocity
                        diffusion = f.diffusion
                else:
                    for var in Rvars:
                        inData[var] = f.variables[var][:] 
                    aggData = agg(inData,aggData,options['resolution'],options['fill_value'],options['verbose'])
                flag = 1
                count += 1
        # Add pad to final file
        if (len(Flist)>0 and flag==1):
             aggData = agg([],aggData,options['resolution'],options['verbose'],options['fill_value'],pad=options['pad'])
             # Write out to netCDF
             write_netcdf(paths['aggDir']+aggFile,aggData['lon'],aggData['lat'],aggData['time'],aggData['unit_hydrograph'],
                          aggData['fraction'],loc,Flist,velocity,diffusion,options['fill_value'],options['verbose'])
        if (flag == 1 and options['remap']):
            remap_file(paths['gridFile'],aggFile,paths['aggDir'],paths['remapDir'],options['verbose'])
            if options['verbose']:
                print 'Finished',key, i,'of',len(d), 'and placed', count, 'files'
            if options['clean']:
                #clean agg directory
                clean(paths['aggDir'])
    print 'done with everything'
    return

##################################################################################
def read_netcdf(nc_str,vars,verbose):
    #if verbose == True:
    #    print 'Reading input data vars:', vars, 'from file:',nc_str
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
    f.close()
    return d

##################################################################################
def make_degrees(vars,Inputs,verbose):  
    d={}
    for var in vars:
        d[var]=np.rad2deg(Inputs[var])
    if verbose == True:
        print 'Converted', vars, 'to degrees'
    return d

##################################################################################
def filename(prefix,tup):
    lon = ('%.5f' % tup[0])[:-2]
    lat = ('%.5f' % tup[1])[:-2]
    string = prefix+lon+'_'+lat+'.nc'

    return string

##################################################################################
def remap_file(gridFile,aggFile,aggDir,remapDir,verbose):
    aggFile = aggDir+aggFile
    remapFile = remapDir+aggFile
    cdo.remapcon(gridFile, input = aggFile, output = remapFile, options = '-f nc4')
    if verbose:
        print 'remapped to out file:', outFile
    return

##################################################################################
##  Write output to netCDF
##  Writes out a netCD4 data file containing the UH_S and fractions
##################################################################################
def write_netcdf(file,lons,lats,times,hydrographs,fractions,loc,Flist,velocity,diffusion,fill_value,verbose):
    """
    Write output to netCDF.  Writes out a netCDF4-64BIT data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    f = Dataset(file,'w', format='NETCDF4_CLASSIC')

    # set dimensions
    time = f.createDimension('time', None)
    lon = f.createDimension('lon', (len(lons)))
    lat = f.createDimension('lat', (len(lats)))

    # initialize variables
    time = f.createVariable('time','f8',('time',))
    lon = f.createVariable('lon','f8',('lon',))
    lat = f.createVariable('lat','f8',('lat',))
    fraction = f.createVariable('fraction','f8',('lat','lon',),fill_value=fill_value)
    UHS = f.createVariable('unit_hydrograph','f8',('time','lat','lon',),fill_value=fill_value)

    # write attributes for netcdf
    f.description = 'Aggregated UH_S and Fraction Vars'
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0] # prints the name of script used
    f.velocity = velocity
    f.diffusion = diffusion
    f.outlet_lon = loc[0]
    f.outlet_lat = loc[1]
    f.includes = ', '.join(Flist)

    lat.long_name = 'latitude coordinate'
    lat.standard_name = 'latitude'
    lat.units = 'degrees_north'

    lon.long_name = 'longitude coordinate'
    lon.standard_name = 'longitude'
    lon.units = 'degrees_east'

    time.standard_name = 'time'
    time.units = 'seconds'
    time.description = 'Seconds since initial impulse'
    time.calendar = 'proleptic_gregorian'

    UHS.units = 'unitless'
    UHS.description = 'unit hydrograph for each grid cell with respect to downstream grid location'
    
    fraction.units = 'unitless'
    fraction.description = 'fraction of grid cell contributing to guage location'

    # write data to variables initialized above
    time[:]= times
    lon[:] = lons
    lat[:] = lats
    UHS[:,:,:] = hydrographs
    fraction[:,:]= fractions
    f.close()
##################################################################################
def agg(inData,aggData,resolution,verbose,fill_value,pad=0):
    """
    Add the two data sets together and return the combined arrays.
    The two data sets must include the coordinate variables lon,lat, and time
    """
    if (inData and aggData):
        # find range of vectors
        Range = [np.minimum(inData['time'].min(),aggData['time'].min()),
                 np.maximum(inData['time'].max(),aggData['time'].max()),
                 np.minimum(inData['lat'].min(),aggData['lat'].min()),
                 np.maximum(inData['lat'].max(),aggData['lat'].max()),
                 np.minimum(inData['lon'].min(),aggData['lon'].min()),
                 np.maximum(inData['lon'].max(),aggData['lon'].max())]
        tStep = inData['time'][1] - inData['time'][0]
        try: yres = np.absolute(aggData['lat'][1] - aggData['lat'][0])
        except:
            try: yres = np.absolute(inData['lat'][1] - inData['lat'][0])
            except: yres = resolution
        try: xres = np.absolute(aggData['lon'][1] - aggData['lon'][0])
        except:
            try: xres = np.absolute(inData['lon'][1] - inData['lon'][0])
            except: xres = resolution
    elif inData:
        Range = [inData['time'].min(),inData['time'].max(),
                 inData['lat'].min(),inData['lat'].max(),
                 inData['lon'].min(),inData['lon'].max()]
        tStep = inData['time'][1] - inData['time'][0]
        try: yres = np.absolute(inData['lat'][1] - inData['lat'][0])
        except: yres = resolution
        try: xres = np.absolute(inData['lon'][1] - inData['lon'][0])
        except: xres = resolution
    elif aggData:
        Range = [aggData['time'].min(),aggData['time'].max(),
                 aggData['lat'].min(),aggData['lat'].max(),
                 aggData['lon'].min(),aggData['lon'].max()]
        tStep = aggData['time'][1] - aggData['time'][0]
        try: yres = np.absolute(aggData['lat'][1] - aggData['lat'][0])
        except: yres = resolution
        try: xres = np.absolute(aggData['lon'][1] - aggData['lon'][0])
        except: xres = resolution
    else:
        raise IOError('no inputs to agg function')
    # make output arrays for lons/lats and initialize fractions/hydrographs
    # pad output arrays so there is a space = pad around inputs
    times = np.arange(Range[0],Range[1]+tStep,tStep)
    lats = np.arange(Range[2]-yres*(pad),Range[3]+yres*(1+pad),yres)[::-1]
    lons = np.arange(Range[4]-xres*(pad),Range[5]+xres*(1+pad),xres)
    fractions = np.zeros((lats.shape[0],lons.shape[0]))
    hydrographs = np.zeros((times.shape[0],lats.shape[0],lons.shape[0]))
    
    # find target index locations of all corners for both datasets
    if inData:
        In = [find_nearest(times,np.min(inData['time'])), find_nearest(times,np.max(inData['time']))+1,
              find_nearest(lats,np.max(inData['lat'])), find_nearest(lats,np.min(inData['lat']))+1,
              find_nearest(lons,np.min(inData['lon'])), find_nearest(lons,np.max(inData['lon']))+1]
    if aggData:
        Ex = [find_nearest(times,np.min(aggData['time'])), find_nearest(times,np.max(aggData['time']))+1,
              find_nearest(lats,np.max(aggData['lat'])), find_nearest(lats,np.min(aggData['lat']))+1,
              find_nearest(lons,np.min(aggData['lon'])), find_nearest(lons,np.max(aggData['lon']))+1]

    # Make sure all values in the unit hydrograph are zero (no mask)
    if inData: inData['unit_hydrograph'][inData['unit_hydrograph']<0] = 0.0
    if aggData: aggData['unit_hydrograph'][aggData['unit_hydrograph']<0] = 0.0

    # Place data
    if inData:
        fractions[In[2]:In[3],In[4]:In[5]] += inData['fraction']
        hydrographs[In[0]:In[1],In[2]:In[3],In[4]:In[5]] += inData['unit_hydrograph']
    if aggData:
        fractions[Ex[2]:Ex[3],Ex[4]:Ex[5]] += aggData['fraction']
        hydrographs[Ex[0]:Ex[1],Ex[2]:Ex[3],Ex[4]:Ex[5]] += aggData['unit_hydrograph']
    
    # Mask the hydrographs and make sure they sum to 1 at each grid cell
    if (inData == [] or aggData == []):
        mask = np.broadcast_arrays(hydrographs, fractions<=0.0)[1]
        hydrographs = ma.masked_array(hydrographs,mask=mask)
        ma.set_fill_value(hydrographs, fill_value)
        
        # Normalize the hydrographs (each cell should sum to 1)
        hydrographs /= hydrographs.sum(axis=0)
        hydrographs = ma.filled(hydrographs, fill_value)
    # Put all the data into aggData variable and return to main
    aggData['lon']  = lons
    aggData['lat']  = lats
    aggData['fraction']  = fractions
    aggData['unit_hydrograph']  = hydrographs
    aggData['time']  = times

    return aggData
    
##################################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##################################################################################
def find_nearest(array,value):
    """ Find the index location in (array) with value nearest to (value)"""
    idx = (np.abs(array-value)).argmin()
    return idx

##################################################################################
def process_command_line():
    """
    Process command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("mapFile", type=str, help="Input mapping weights beween target grid and routing grid, example of how to make these weights: cdo gennn,Wu_routing_inputs.nc domain.lnd.wr50a_ar9v4.100920.nc out_weights.nc")
    parser.add_argument("srcDir", type=str, help="Directory containing Unit Hydrograph grids to be aggregated")
    parser.add_argument("--gridFile", type=str, help="Input netCDF target grid")
    parser.add_argument("--remapDir", type=str, help="Directory containing Output Unit Hydrograph grids")
    parser.add_argument("--aggDir", type=str, help="Directory where to store aggregated files (before remap)")
    parser.add_argument("--inPrefix", type=str, help="Input Unit Hydrograph File Prefix",default='UH_')
    parser.add_argument("--outPrefix", type=str, help="Output Unit Hydrograph File Prefix", default="Agg_UH_")
    parser.add_argument("--time", type=str, help="Input Unit Hydrograph time variable name",default='time')
    parser.add_argument("--lon", type=str, help="Input Unit Hydrograph longitude variable name",default='lon')
    parser.add_argument("--lat", type=str, help="Input Unit Hydrograph latitude variable name",default='lat')
    parser.add_argument("--fraction", type=str, help="Input Unit Hydrograph fraction variable name",default='fraction')
    parser.add_argument("--unit_hydrograph",type=str, help="Input Unit Hydrograph unit hydrograph variable name",default='unit_hydrograph')
    parser.add_argument("--cdoDebug",help="Enable CDO debuging (prings each step to screen)",action="store_true")
    parser.add_argument("--cdoForce",help="Enable CDO force output (will overwrite existing files during remap)",action="store_true")
    parser.add_argument("--verbose",help="Make script verbose",action="store_true")
    parser.add_argument("--remap",help="Remap the aggregated Unit Hydrographs to outDir and put the aggregated files in the tempDir",action='store_true')
    parser.add_argument("--fill_value",type=float,help="value to use as masked value",default = 9.96920996839e+36)
    parser.add_argument("--pad",type=int,help="Set number of empty cells to include around each aggregated basin",default=10)
    parser.add_argument("--resolution",type=float,help="Set resolution of input Unit Hydrographs",default=1/16.)
    parser.add_argument("--clean",help="Clean up aggregated Unit Hydrograph grids if remapping", action='store_true')
    args = parser.parse_args()

    options = {}
    paths = {}
    # parse the basics
    Rvars = (args.time,args.lon,args.lat,args.fraction,args.unit_hydrograph)
    paths['srcDir'] = args.srcDir
    paths['mapFile'] = args.mapFile
    paths['gridFile'] = args.gridFile

    if args.aggDir:
        paths['aggDir'] = args.aggDir
    else:
        paths['aggDir'] = paths['srcDir']+'../aggregated/'
        if not os.path.exists(paths['aggDir']):
            os.makedirs(paths['aggDir'])

    options['verbose'] = args.verbose
    options['fill_value'] = args.fill_value
    options['pad'] = args.pad
    options['resolution'] = args.resolution
    options['inPrefix'] = args.inPrefix
    options['outPrefix'] = args.outPrefix

    if args.remap:
        options['remap']=True
        options['clean']=args.clean
        cdo.debug=args.cdoDebug
        cdo.forceOutput=args.cdoForce
        if args.remapDir:
            paths['remapDir'] = args.remapDir
        else:
            paths['remapDir'] = paths['srcDir']+'../remaped/'
            if not os.path.exists(paths['remapDir']):
                os.makedirs(paths['remapDir'])
        print paths['remapDir']

    return Rvars,paths,options

def clean(aggDir):
    for file in os.listdir(aggDir):
        file_path = os.path.join(aggDir, file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception, e:
            print e
    return
##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
