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
uhs = open( infile, "r").readlines()[2::2]
uhs = [item.rstrip() for item in uhs]
uhs = np.array([x.split() for x in uhs], dtype='float')

#reshape the uhs output
#make individual arrays of lats, lons, cell nums and fracs
lat_in=indata1[:,0]
lon_in=indata1[:,1]
frac_in=indata1[:,2]
xcell_in=indata1[:,3]
ycell_in=indata1[:,4]

xin=xcell_in-min(xcell_in)
yin=ycell_in-min(ycell_in)

xmin=min(xin)
xmax=max(xin)
ymin=min(yin)
ymax=max(yin)

xs = np.arange(xmin,xmax,res)
ys = np.arange(ymin,ymax,res)

# initialize the target arrays (uhs&fra-> a&b) with the nodata value
#a = np.zeros((len(xs),len(ys),uhs.shape[1]))+fill_value
b = np.zeros((1,len(xs),len(ys)))+fill_value

#fill in values values
for n, i in enumerate(xin):
   b[0,xin[i],yin[i]] = frac_in[i]
   #a[xin[i],yin[i],j] = uhs_in[:,i]


   ############# WRITE TO NETCDF4 ##################
   
# initialize the netcdf file
f = Dataset('test.nc', 'w', format='NETCDF4_CLASSIC')
print(f.file_format)


# set dimensions
time = f.createDimension('time', None)
lon = f.createDimension('lon', (len(xs)))
lat = f.createDimension('lat', (len(ys)))
#print(f.dimensions)

# initialize variables
times = f.createVariable('time','f8',('time',))
longitudes = f.createVariable('longitude','f4',('lon',))
latitudes = f.createVariable('latitude','f4',('lat',))
frac = f.createVariable('fraction','f8',('lon','lat','time'),fill_value=fill_value)
#uh = f.createVariable('unit_hydrograph','f8',('lon','lat','time'),fill_value=fill_value)

import time
# write attributes for netcdf
f.description = 'test script'
f.history = 'Created ' + time.ctime(time.time())
f.source = 'RASM routing program remap.py'
latitudes.units = 'degrees north'
longitudes.units = 'degrees east'
times.units = 'Days'

#f.close()

#print attributes
#for name in f.ncattrs():
#   print 'Global attr', name, '=', getattr(f.name)

f = Dataset('test.nc','a')
# data
lons1 = np.arange(xmin,xmax,res)
longitudes[:] = lons1
latitudes[:] = np.arange(ymin,ymax,res)

print np.shape(b)

frac[:,:,:] = b[0,:,:]
#uhs[:] = a

f.close()

print '*** SUCCESS writing tst.nc!'

# remap the uhs.nc file using cdo
#outfile = "outtest.nc"
#os.system("cdo .......")

