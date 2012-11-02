#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import os
import time

infile='/Users/jhamman/Dropbox/RASM_Joe/YUKON.uh_s'
nodata=-9999


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
lats=indata1[:,0]
lons=indata1[:,1]
frac=indata1[:,2]
xcell=indata1[:,3]
ycell=indata1[:,4]

# initialize the target arrays (uhs&fra-> a&b) with the nodata value
a = np.zeros((len(lons),len(lats),uhs.shape[1]))+nodata
b = np.zeros((lons.length,lats.length))+nodata

xin=xcell-min(xcell)+1
yin=ycell-min(ycell)+1

#fill in values values
for n, i in enumerate(xin):
   b[xin[i],yin[i]] = fracin[i]
   a[xin[i],yin[i],:] = uhs_in[:,i]


   ############# WRITE TO NETCDF4 ##################
   
# dimensions
lon = f.createDimension('lon', xdim)
lat = f.createDimension('lat', ydim)
time = f.createDimension('time',tdim)
print(f.dimensions)

# variables
longitudes = f.createVariable('longitude','f4',('lon',))
latitudes = f.createVariable('latitude','f4',('lat',))
times = f.createVariable('time','i4',('time',))
frac = f.createVariable('fraction','f8',('lon','lat','time'))
uh = f.createVariable('unit_hydrograph','f8',('lon','lat','time'))

#write attributes for netcdf
f.description = 'test script'
f.history = 'Created ' + create_t
f.source = 'RASM routing program'
latitudes.units = 'degrees north'
longitudes.units = 'degrees east'
times.units = 'Days'

f.close()

#print attributes
for name in f.ncattrs():
   print 'Global attr', name, '=', getattr(f.name)

f = Dataset('test.nc','a')
# data
lon[:] = np.arange(xmin,xmax,res)
lat[:] = np.arange(ymin,ymax,res)
fracs[:] = b
uhs[:] = a

f.close()

# remap the uhs.nc file using cdo
#outfile = "outtest.nc"
#os.system("cdo .......")

