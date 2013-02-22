#!/usr/local/bin/python

import os
import sys
import numpy as np
import numpy.ma as ma
from collections import defaultdict
from netCDF4 import Dataset
import time as tm
import glob

##################################################################################
in_dir = '/raid/jhamman/temp_uh_files/run1_RASM_remap/'
out_dir = '/raid/jhamman/temp_uh_files/'
outFile = 'All_RASM.nc'
verbose = False
Agvars = ['fraction','unit_hydrograph']
NODATA = -9999

##################################################################################
def main():
    # Get list of files
    trouble = []
    files = glob.glob(in_dir+'RASM*.nc')
    print 'aggregating',str(len(files)), 'into', outFile
    flag = 0
    # Loop over all files
    for i,file in enumerate(files):
        if flag == 0:
            data = read_netcdf(file,[],verbose)
            velocity = Dataset(file,'r').velocity
            diffusion = Dataset(file,'r').diffusion
            for var in Agvars:
                data[var] = ma.filled(data[var],0.0)
                data[var][data[var]<0] = 0.0
            flag = 1
        elif flag == 1:
            temp = read_netcdf(file,Agvars,verbose)
            for var in Agvars:
                try:
                    temp[var] = ma.filled(temp[var],0.0)
                    #data[var] = data[var] + temp[var]
                except:
                    print 'something went wrong with',file
                    trouble.append(file)
            data['fraction'] = data['fraction']+temp['fraction']
            data['unit_hydrograph'] = data['unit_hydrograph']+temp['fraction']*temp['unit_hydrograph']
            print 'done with', i+1, 'of', len(files)
            temp = []
            #if i > 20:
            #flag =2
    print 'done with aggregation, now normalizing hydrographs'
    print 'had trouble with:', trouble
    # Do some work on the hydrographs data
    mask = np.broadcast_arrays(data['unit_hydrograph'],data['fraction']<=0.001)[1]
    data['unit_hydrograph'] = ma.masked_array(data['unit_hydrograph'],mask=mask)
    ma.set_fill_value(data['unit_hydrograph'], NODATA)
    # Normalize the hydrographs (each cell should sum to 1)
    data['unit_hydrograph'] /= data['unit_hydrograph'].sum(axis=0)
    data['unit_hydrograph'] = ma.filled(data['unit_hydrograph'], NODATA)
    # Write the data out
    print 'writing data out now'
    write_netcdf(out_dir+outFile,data['xc'],data['xc_bnds'],data['yc'],data['yc_bnds'],data['time'],data['unit_hydrograph'],
                 data['fraction'],('RASM','RASM'),files,velocity,diffusion,NODATA,verbose)
    
##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile,vars,verbose):
    """
    Read data from input netCDF. Will read all variables if none provided.
    """
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    if verbose: print 'Reading input data vars:', vars, ',from file:',ncFile
    d={}
    for var in vars:
        try:
            d[var] = f.variables[var][:]
        except:
            d[var] = []
    f.close()
    return d

##################################################################################
##  Write output to netCDF
##  Writes out a netCD4 data file containing the UH_S and fractions
##################################################################################
def write_netcdf(file,xc,xc_bnd,yc,yc_bnd,times,hydrographs,fractions,loc,Flist,velocity,diffusion,NODATA,verbose):
    """
    Write output to netCDF.  Writes out a netCDF4 data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    f = Dataset(file,'w', format='NETCDF4')

    # set dimensions
    time = f.createDimension('time', None)
    x = f.createDimension('x',xc.shape[1])
    y = f.createDimension('y',xc.shape[0])
    nv4 = f.createDimension('nv4',4)

    # initialize variables
    time = f.createVariable('time','f8',('time',))
    xcs = f.createVariable('xc','f8',('y','x',))
    ycs = f.createVariable('yc','f8',('y','x',))
    xc_bnds = f.createVariable('xc_bnds','f8',('y','x','nv4',))
    yc_bnds = f.createVariable('yc_bnds','f8',('y','x','nv4',))
    fraction = f.createVariable('fraction','f8',('y','x',),fill_value=NODATA)
    UHS = f.createVariable('unit_hydrograph','f8',('time','y','x',),fill_value=NODATA)

    # write attributes for netcdf
    f.description = 'Aggregated UH_S and Fraction Vars for full RASM domain'
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0] # prints the name of script used
    f.velocity = velocity
    f.diffusion = diffusion
    f.outlet_lon = loc[0]
    f.outlet_lat = loc[1]
    f.includes = str(len(Flist))+' files'

    ycs.long_name = 'latitude of grid cell center'
    ycs.standard_name = 'latitude'
    ycs.units = 'degrees_north'
    ycs._CoordinateAxisType = 'Lat'
    ycs.bounds = 'yc_bnds'

    xcs.long_name = 'longitude of grid cell center'
    xcs.standard_name = 'longitude'
    xcs.units = 'degrees_east'
    xcs._CoordinateAxisType = 'Lon'
    xcs.bounds = 'xc_bnds'

    time.standard_name = 'time'
    time.units = 'seconds'
    time.description = 'Seconds since initial impulse'
    time.calendar = 'proleptic_gregorian'

    UHS.units = 'unitless'
    UHS.description = 'unit hydrograph for each grid cell with respect to basin outlet location'
    
    fraction.units = 'unitless'
    fraction.description = 'fraction of grid cell contributing to guage location'

    # write data to variables initialized above
    time[:]= times
    xcs[:,:] = xc
    ycs[:,:] = yc
    xc_bnds[:,:,:] = xc_bnd
    yc_bnds[:,:,:] = yc_bnd
    UHS[:,:,:] = hydrographs
    fraction[:,:]= fractions
    f.close()

##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
