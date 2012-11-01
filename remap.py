#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import os
import time

path='/raid/jhamman/UH_S_files/'
res=0.0625
cell=12098
infile='0.0625_YUKON.uh_s'
nodata=-9999

#load UHS file in column form
indata = np.loadtxt(path+infile,delimiter=' ')

#reshape the uhs output
#make individual arrays of lats, lons, cell nums and uhs
lats=
lons=
xcell=
ycell=
frac=
uhs=

# initialize the target arrays (uhs&fra-> a&b) with the nodata value
a = np.zeros((lons.length,lats.length,uhs.length))+nodata
b = np.zeros((lons.length,lats.length))+nodata

#fill in values values
for n, i in enumerate(xin):
   b[xin[i],yin[i]] = fracin[i]
   a[xin[i],yin[i],:] = uhs_in[i,:]

print 'Global attr', name, '=', getattr(f.name)

f = Dataset('test.nc','a')

# data
lon[:] = np.arange(min.lons,max.lons,res)
lat[:] = np.arange(min.lats,max.lats,res)
fracs[:] = b
uhs[:] = a

f.close()

# remap the uhs.nc file using cdo
outfile = "outtest.nc"
os.system("cdo .......")

