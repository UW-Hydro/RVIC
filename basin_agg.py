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
weight_prefix = 'Weights_'
Mvars = ('src_address','dst_address','src_grid_center_lat','dst_grid_center_lat','src_grid_center_lon','dst_grid_center_lon')
Cvars = ('src_grid_center_lat','dst_grid_center_lat','src_grid_center_lon','dst_grid_center_lon')
Rvars = ('time','lon','lat','fraction','unit_hydrograph')

##################################################################################
cdo.debug = False
cdo.forceOutput = True
verbose=True
resolution = 1/16.
tStep = 86400
velocity = 1.0
diffusion = 2000.0
NODATA = -9999.0
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
    if verbose == True:
        print 'Remapping and Aggregating now...'
    for i,key in enumerate(d):
        flag  = 0
        Flist = []
        #---> initiate 
        for j,loc in enumerate(d[key]):
            aggFile = filename(out_prefix,key)
            inFile = filename(in_prefix,loc)
            if os.path.isfile(in_dir+inFile):
                Flist.append(inFile)
                if (flag == 0 and os.path.isfile(temp_dir+aggFile)==False):
                    shutil.copy2(in_dir+inFile,temp_dir+aggFile)
                else:
                    times,lons,lats,fractions,hydrographs = agg(in_dir+inFile,temp_dir+aggFile,verbose)
                    # Write out to netCDF
                    write_netcdf(temp_dir+aggFile,lons,lats,times,hydrographs,
                                 fractions,loc,Flist,velocity,diffusion,NODATA,verbose)
                flag = 1
                count += 1
        if flag == 1:
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
    string = prefix+('%.5f' % tup[0])+'_'+('%.5f' % tup[1])+'.nc'
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
def agg(inFile,aggFile,verbose):
    inData = read_netcdf(inFile,Rvars,verbose)
    exData = read_netcdf(aggFile,Rvars,verbose)
    
    # find range of vectors
    min_x = np.minimum(inData['lon'].min(),exData['lon'].min())
    max_x = np.maximum(inData['lon'].max(),exData['lon'].max())
    min_y = np.minimum(inData['lat'].min(),exData['lat'].min())
    max_y = np.maximum(inData['lat'].max(),exData['lat'].max())
    min_t = np.minimum(inData['time'].min(),exData['time'].min())
    max_t = np.maximum(inData['time'].max(),exData['time'].min())

    # make output arrays for lons/lats and initialize fractions/hydrographs
    lons = np.arange(min_x,max_x+resolution,resolution)
    lats = np.arange(min_y,max_y+resolution,resolution)
    times = np.arange(min_t,max_t+tStep,tStep)
    fractions = np.zeros((lats.shape[0],lons.shape[0]))
    hydrographs = np.zeros((times.shape[0],lats.shape[0],lons.shape[0]))
    
    # find target index locations of all corners for both datasets
    In = [find_nearest(times,np.min(inData['time'])), find_nearest(times,np.max(inData['time']))+1,
          find_nearest(lats,np.min(inData['lat'])), find_nearest(lats,np.max(inData['lat']))+1,
          find_nearest(lons,np.min(inData['lon'])), find_nearest(lons,np.max(inData['lon']))+1]
    Ex = [find_nearest(times,np.min(exData['time'])), find_nearest(times,np.max(exData['time']))+1,
          find_nearest(lats,np.min(exData['lat'])), find_nearest(lats,np.max(exData['lat']))+1,
          find_nearest(lons,np.min(exData['lon'])), find_nearest(lons,np.max(exData['lon']))+1]

    inData['unit_hydrograph'][inData['unit_hydrograph']<0] = 0.0
    exData['unit_hydrograph'][exData['unit_hydrograph']<0] = 0.0

    # Place data
    fractions[In[2]:In[3],In[4]:In[5]] += inData['fraction']
    fractions[Ex[2]:Ex[3],Ex[4]:Ex[5]] += exData['fraction']

    hydrographs[In[0]:In[1],In[2]:In[3],In[4]:In[5]] += inData['unit_hydrograph']
    hydrographs[Ex[0]:Ex[1],Ex[2]:Ex[3],Ex[4]:Ex[5]] += exData['unit_hydrograph']

    mask = np.broadcast_arrays(hydrographs, fractions<=0.0)[1]
    hydrographs = ma.masked_array(hydrographs,mask=mask)
    ma.set_fill_value(hydrographs, NODATA)
    
    # Normalize the hydrographs (each cell should sum to 1)
    hydrographs /= hydrographs.sum(axis=0)
    hydrographs = ma.filled(hydrographs, NODATA)

    return times,lons,lats,fractions,hydrographs
    
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
