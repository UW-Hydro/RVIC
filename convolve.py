
#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import argparse

# Parse arguments
## parser = argparse.ArgumentParser()
## parser.add_argument("basin", help="5 letter basin name (ie. MACRR")
## args = parser.parse_args() 
## basin=args.basin

basins = ['PECUT','OBSAL','LENKU','KOLKO','YENIG','MACRR','YUKPS']
runs = ['r33RBVIC60', 'r33RBVIC70']

grid_nc = '/nfs/hydro6/raid/nijssen/rasm/masks/racm_masks_121108.nc'
f3 = Dataset(grid_nc)
area  = f3.variables['DOMA_AREA'][:]

re = 6.37122e6
area = area[:,:]*re*re
for run in runs:

    rasm_nc = '/nfs/thermal/raid/jhamman/'+run+'.vic.ha.1990s.qflux.nc'
    f2 = Dataset(rasm_nc)
    runoff = f2.variables['Runoff'][:]
    baseflow = f2.variables['Baseflow'][:]
    time1 = f2.variables['time'][:]
    flux = runoff + baseflow

    for basin in basins:
        print 'starting basin', basin
        # open the UH netcdf file
        basin_nc ='/raid/jhamman/temp_uh_files/'+basin+'_RASM_UH.nc'
        f1 = Dataset(basin_nc)
        uh = f1.variables['unit_hydrograph'][:]
        frac = f1.variables['fraction'][:]
        t_uh = f1.variables['time'][:]

        q = np.zeros((len(flux)+len(t_uh)))

        #Do the convolution
        for i in xrange(len(flux)):
            for j in xrange(len(t_uh)):
                tstep = i+j
                q[tstep]= q[tstep] + np.sum(flux[i,:,:]*uh[j,:,:]*area[:,:]*frac[:,:])

        #Change units to m3/sec
        q = q[:]/1000/86400

        #Write to NC file
        out_nc = basin+'_'+run+'_streamflow.nc'
        f = Dataset(out_nc, 'w', format='NETCDF4_CLASSIC')
        print(f.file_format)

        # set dimensions
        time = f.createDimension('time', None)
        print(f.dimensions)

        # initialize variables
        time = f.createVariable('time','f8',('time',))
        Q = f.createVariable('Q','f8',('time',))

        import time as t
        # write attributes for netcdf
        f.description = 'Streamflow for RASM'
        f.history = 'Created ' + t.ctime(t.time())
        f.source = 'RASM routing program - convolve'

        time.units = 'day'
        time.description = 'Days since 1-1-1'

        Q.description = 'Daily runoff at station'
        Q.units = 'm^3/s'

        # write data to variables initialized above
        time[:] = time1
        Q[:] = q[0:len(time1)]

        f.close()
        print '*** SUCCESS writing streamflows for ',run+'_'+basin, '***'

print 'Done'    

