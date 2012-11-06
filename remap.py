#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import os
import time

infile='/usr1/jhamman/RASM/vic_routing/route_rasm/arctic/UH_S_files/YUKON.uh_s'
fill_value=-9999
res=0.0625

#Load the Header data for each point
indata1 = open( infile, "r").readlines()[1::2]
indata1 = [item.rstrip() for item in indata1]
indata1 = np.array([x.split() for x in indata1], dtype='float')

#Load the UH data for each point
uhs_in = open( infile, "r").readlines()[2::2]
uhs_in = [item.rstrip() for item in uhs_in]
uhs_in = np.array([x.split() for x in uhs_in], dtype='float')

#reshape the uhs output
#make individual arrays of lats, lons, cell nums and fracs
lon_in=indata1[:,0]
lat_in=indata1[:,1]
frac_in=indata1[:,2]
xcell_in=indata1[:,3]
ycell_in=indata1[:,4]

print 'lon_in = '
print lon_in
print 'lat_in = '
print lat_in

xin=xcell_in-min(xcell_in)
yin=ycell_in-min(ycell_in)
xin=xin.astype('int32')
yin=yin.astype('int32')

lons = np.arange(min(lon_in),max(lon_in)+res,res)
lats = np.arange(min(lat_in),max(lat_in)+res,res)

# initialize the target arrays (uhs&fra-> a&b) with the nodata value
a = np.zeros((uhs_in.shape[1],len(lats),len(lons)))+fill_value
b = np.zeros((len(lats),len(lons)), dtype=np.int)+fill_value
#fill in values values
for n, i in enumerate(xin):
    b[yin[n],xin[n]] = frac_in[n]
    a[:,yin[n],xin[n]] = uhs_in[n,:]

   ############# WRITE TO NETCDF4 ##################
   
# initialize the netcdf file
f = Dataset('YUKON_UH.nc', 'w', format='NETCDF4_CLASSIC')
print(f.file_format)


# set dimensions
time = f.createDimension('time', None)
lon = f.createDimension('lon', (len(lons)))
lat = f.createDimension('lat', (len(lats)))
print(f.dimensions)

# initialize variables
times = f.createVariable('time','f8',('time',))
lon = f.createVariable('lon','f4',('lon',))
lat = f.createVariable('lat','f4',('lat',))
frac = f.createVariable('fraction','f8',('lat','lon',),fill_value=fill_value)
uh = f.createVariable('unit_hydrograph','f8',('time','lat','lon',),fill_value=fill_value)

import time
# write attributes for netcdf
f.description = 'test_script'
f.history = 'Created ' + time.ctime(time.time())
f.source = 'RASM routing program remap.py'
lat.units = 'degrees_north'
lon.units = 'degrees_east'
times.units = 'Days'
times.units = 'Days to runoff at downstream outlet'

# data
lon[:] = lons
lat[:] = lats
frac[:,:] = b
uh[:,:,:] = a

f.close()

print '*** SUCCESS writing test.nc!'

# remap the uhs.nc file using cdo
#outfile = "outtest.nc"
#os.system("cdo .......")

