#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

basin_nc ='/usr1/jhamman/RASM/rout_rasm/YUKON_RASM_UH.nc'
# open the UH netcdf file

f1 = Dataset(basin_nc)
uh = f1.variables['unit_hydrograph'][:]
frac = f1.variables['fraction'][:]
t_uh = f1.variables['time'][:]

rasm_nc = '/nfs/hydro6/raid/nijssen/rasm/r33RBVIC70/lnd/hist/r33RBVIC70.vic.ha.1998.daily.nc'
f2 = Dataset(rasm_nc)
runoff = f2.variables['Runoff']
baseflow = f2.variables['Baseflow'][:]
time = f2.variables['time'][:]
flux = runoff + baseflow

q = np.zeros((len(flux)+len(t_uh)))
area = 1

for i in xrange(len(flux)):
    for j in xrange(len(t_uh)):
        tstep = i+j
        q[tstep]= q[tstep] + np.sum(flux[i,:,:] * uh[j,:,:] * frac[:,:]* area)

plt.plot(q)
plt.show()
