#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import time

#load UHS file in column form
path='/raid/jhamman/UH_S_files/'
res=0.0625
cell=12098
infile='0.0625_YUKON.uh_s'

indata = np.loadtxt(path+infile,delimiter=' ', skiprows=2)
print indata

xdim=5760
ydim=2240
tdim= 107
nodata=-9999

xmin=-179.969
xmax=179.969
ymin=-55.9688
ymax=83.9688


# initialize the target arrays (uhs&fra-> a&b) with the nodata value
#a = np.zeros((xdim,ydim,tdim))+nodata
b = np.zeros((xdim,ydim))+nodata

# parse the input array
#lonin=indata[0,:]
#latin=indata[1,:]
fracin=indata[2,:]
xin=indata[3,:]-1
yin=indata[4,:]-1
uhs_in=indata[3:,:]

#fill in values values
for n, i in enumerate(xin):
   b[xin[i],yin[i]] = fracin[i]
   #a[xin[i],yin[i],:] = uhs_in[i,:]


# initialize the netcdf file
rootgrp = Dataset('test.nc', 'w', format='NETCDF4_CLASSIC')
print rootgrp.file_format

rootgrp.close()

# dimensions
lon = rootgrp.createDimension('lon', xdim)
lat = rootgrp.createDimension('lat', ydim)
time = rootgrp.createDimension('time',tdim
print rootgrp.dimensions

# variables
longitudes = rootgrp.createVariable('longitude','f4',('lon',))
latitudes = rootgrp.createVariable('latitude','f4',('lat',))
times = rootgrp.createVariable('time','i4',('time',))
frac = rootgrp.createVariable('fraction','f8',('lon','lat','time'))
uh = rootgrp.createVariable('unit_hydrograph','f8',('lon','lat','time'))

#write attributes for netcdf
rootgrp.description = 'test script'
rootgrp.history = 'Created ' + time.ctime(time.time())
rootgrp.source = 'RASM routing program'
latitudes.units = 'degrees north'
longitudes.units = 'degrees east'
times.units = 'Days'

#print attributes
for name in rootgrp.ncattrs():
   print 'Global attr', name, '=', getattr(rootgrp.name)

# data
lon[:] = np.linspace(xmin, xmax, xdim)
lat[:] = np.linspace(ymin, ymax, ydim)
fracs[:] = b
#uhs[:] = a


