#!/usr/local/bin/python

import os
import sys
import shutil
import numpy as np
import numpy.ma as ma
from collections import defaultdict
from netCDF4 import Dataset
from cdo import *
cdo = Cdo()
import time as tm

nc_str = '/usr1/jhamman/Dropbox/RASM_Joe/python_scripts/rout_rasm/rev_weights.nc'
gridFile = '/raid/jhamman/RASM_masks/domain.lnd.wr50a_ar9v4.100920.nc'
in_dir = '/raid/jhamman/temp_uh_files/run1_RASM/'
out_dir = '/raid/jhamman/temp_uh_files/run1_RASM_remap/'
weight_dir = '/raid/jhamman/remap_weights/'
temp_dir =  '/raid/jhamman/temp_uh_files/temp/'
in_prefix = 'UH_'
out_prefix = 'RASM_'
Mvars = ('src_address','dst_address','src_grid_center_lat','dst_grid_center_lat','src_grid_center_lon','dst_grid_center_lon')
Cvars = ('src_grid_center_lat','dst_grid_center_lat','src_grid_center_lon','dst_grid_center_lon')
Rvars = ('time','lon','lat','fraction','unit_hydrograph')

##################################################################################
cdo.debug = False
cdo.forceOutput = True
verbose=True
remap = True
NODATA = -9999.0
pad = 10
resolution = 1/16.

##################################################################################
def main():
    #Load Inputs
    Inputs = read_netcdf(nc_str,Mvars,verbose)
    
    #Convert lats and lons for src/dst to decimal degrees
    Coords = make_degrees(Cvars,Inputs,verbose)
    
    #Find Unique Destination grid cells
    d = defaultdict(list)
    if verbose == True:
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
    if verbose:
        print 'Remapping and Aggregating now...'
    for i,key in enumerate(d):
        flag  = 0
        Flist = []
        aggFile = filename(out_prefix,key)
        for j,loc in enumerate(d[key]):
            inFile = filename(in_prefix,loc)
            if os.path.isfile(in_dir+inFile):
                Flist.append(inFile)
                if flag == 0:
                    aggData = read_netcdf(in_dir+inFile,Rvars,verbose)
                    if count == 0:
                        velocity = Dataset(in_dir+inFile,'r').velocity
                        diffusion = Dataset(in_dir+inFile,'r').diffusion
                else:
                    inData = read_netcdf(in_dir+inFile,Rvars,verbose)
                    aggData = agg(inData,aggData,0,verbose)
                flag = 1
                count += 1
        # Add pad to final file
        if (len(Flist)>0 and flag==1):
             aggData = agg([],aggData,pad,verbose)
             # Write out to netCDF
             write_netcdf(temp_dir+aggFile,aggData['lon'],aggData['lat'],aggData['time'],aggData['unit_hydrograph'],
                          aggData['fraction'],loc,Flist,velocity,diffusion,NODATA,verbose)
        if (flag == 1 and remap):
            remap_file(aggFile,temp_dir,out_dir,verbose)
            if verbose == True:
                print 'Finished',key, i,'of',len(d), 'and placed', count, 'files'

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
    lon = ('%.5f' % tup[0])[:-3]
    lat = ('%.5f' % tup[1])[:-3]
    string = prefix+lon+'_'+lat+'.nc'
    return string
##################################################################################
def remap_file(aggFile,temp_dir,out_dir,verbose):
    inFile = temp_dir+aggFile
    outFile = out_dir+aggFile
    cdo.remapcon(gridFile, input = inFile, output = outFile, options = '-f nc4')
    #if verbose:
    #    print 'remapped to out file:', outFile
##################################################################################
##  Write output to netCDF
##  Writes out a netCD4 data file containing the UH_S and fractions
##################################################################################
def write_netcdf(file,lons,lats,times,hydrographs,fractions,loc,Flist,velocity,diffusion,NODATA,verbose):
    """
    Write output to netCDF.  Writes out a netCDF4-64BIT data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    f = Dataset(file,'w', format='NETCDF4')

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
def agg(inData,aggData,pad,verbose):
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
        print 'no inputs to agg function'
        raise 
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
        ma.set_fill_value(hydrographs, NODATA)
        
        # Normalize the hydrographs (each cell should sum to 1)
        hydrographs /= hydrographs.sum(axis=0)
        hydrographs = ma.filled(hydrographs, NODATA)
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
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
