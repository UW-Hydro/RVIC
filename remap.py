#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import os
import time
from cdo import *
cdo = Cdo()
import argparse

#infile='/usr1/jhamman/RASM/rout_rasm/YUKPS.uh_s'
fill_value=-9999
res=0.0625

## # Parse arguments
## parser = argparse.ArgumentParser()
## parser.add_argument("basin", help="5 letter basin name (ie. MACRR")
## args = parser.parse_args() 
## basin=args.basin

basins = ['PECUT','OBSAL','LENKU','KOLKO','YENIG','MACRR','YUKPS']
path = '/raid/jhamman/temp_uh_files/'
for basin in basins:
    print 'Starting basin - '+basin
    #Load the Header data for each point
    basin_in = path+basin+'.uh_s'
    indata1 = open(basin_in, "r").readlines()[1::2]
    indata1 = [item.rstrip() for item in indata1]
    indata1 = np.array([x.split() for x in indata1], dtype='d')

    #Load the UH data for each point
    uhs_in = open(basin_in, "r").readlines()[2::2]
    uhs_in = [item.rstrip() for item in uhs_in]
    uhs_in = np.array([x.split() for x in uhs_in], dtype='d')

    #reshape the uhs output
    #make individual arrays of lats, lons, cell nums and fracs
    lon_in=indata1[:,0]
    lat_in=indata1[:,1]
    frac_in=indata1[:,2]
    xin=indata1[:,3].astype('int32')-min(indata1[:,3].astype('int32'))
    yin=indata1[:,4].astype('int32')-min(indata1[:,4].astype('int32'))

    lons = np.arange(min(lon_in),max(lon_in)+res,res)
    lats = np.arange(min(lat_in),max(lat_in)+res,res)

    # initialize the target arrays (uhs&fra-> a&b) with the nodata value
    a = np.zeros((uhs_in.shape[1],len(lats),len(lons)))+fill_value
    b = np.zeros((len(lats),len(lons)), dtype=np.int)
    #fill in values values
    for n, i in enumerate(xin):
        b[yin[n],xin[n]] = frac_in[n]
        a[:,yin[n],xin[n]] = uhs_in[n,:]

    ################ WRITE TO NETCDF4 #########################
    basin_nc = path+basin+'_UH.nc'   
    # initialize the netcdf file
    f = Dataset(basin_nc, 'w', format='NETCDF4_CLASSIC')
    print(f.file_format)

    # set dimensions
    time = f.createDimension('time', None)
    lon = f.createDimension('lon', (len(lons)))
    lat = f.createDimension('lat', (len(lats)))
    print(f.dimensions)

    # initialize variables
    time = f.createVariable('time','f8',('time',))
    lon = f.createVariable('lon','f8',('lon',))
    lat = f.createVariable('lat','f8',('lat',))
    frac = f.createVariable('fraction','f8',('lat','lon',),fill_value=fill_value)
    uh = f.createVariable('unit_hydrograph','f8',('time','lat','lon',),fill_value=fill_value)

    import time
    # write attributes for netcdf
    f.description = 'Unit Hydrographs and Fractional Area for ' +basin+' basin . Created at resolution: '+str(res)
    f.history = 'Created ' + time.ctime(time.time())
    f.source = 'RASM routing program remap.py'
    f.grid_res = res

    lat.long_name = 'latitude coordinate'
    lat.standard_name = 'latitude'
    lat.units = 'degrees_north'

    lon.long_name = 'longitude coordinate'
    lon.standard_name = 'longitude'
    lon.units = 'degrees_east'

    time
    time.units = 'day'
    time.description = 'Days to runoff at downstream outlet'
    frac.units = 'Fraction of contributing area'

    uh.units = 'Unit Hydrograph.  Contributing fraction of runoff from each grid cell to the downstream point of interest.' 

    # write data to variables initialized above
    lon[:] = lons
    lat[:] = lats
    frac[:,:] = b
    uh[:,:,:] = a

    f.close()

    print '*** SUCCESS writing', basin_nc, '***'

    print 'remapping ', basin_nc, ' to RASM grid'
    # remap the uhs.nc file using cdo
    grid = '/usr1/jhamman/RASM/cdo_interp/domain.lnd.wr50a_ar9v4.100920.nc'
    ifile = basin_nc
    ofile = path+basin+'_RASM_UH.nc'
    cdo.remapcon(grid,input = ifile,output = ofile)

print 'done remapping'
