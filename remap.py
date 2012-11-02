#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import os
import time

infile='0.0625_YUKON.uh_s'
nodata=-9999


#Load the Header data for each point
indata1 = open( '/Users/jhamman/Dropbox/RASM_Joe/YUKON.uh_s', "r").readlines()[1::2]
indata1 = [item.rstrip() for item in indata1]
indata1 = np.array([x.split() for x in indata1], dtype='float')

#Load the UH data for each point
uhs = open( '/Users/jhamman/Dropbox/RASM_Joe/YUKON.uh_s', "r").readlines()[2::2]
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

