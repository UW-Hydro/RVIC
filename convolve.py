
#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np

basins = ['PECUT','OBSAL','LENKU','KOLKO','YENIG','MACRR','YUKPS','DVIUP']
runs = ['run2','run3','run4','run5','run6','run7','run8']

sims = ['r33RBVIC60', 'r33RBVIC70']

grid_nc = '/nfs/hydro6/raid/nijssen/rasm/masks/racm_masks_121108.nc'
f3 = Dataset(grid_nc)
area  = f3.variables['DOMA_AREA'][:]

re = 6.37122e6
area = area[:,:]*re*re
for sim in sims:

    rasm_nc = '/nfs/thermal/raid/jhamman/RASM_results/r33/'+sim+'.vic.ha.1990s.qflux.nc'
    f2 = Dataset(rasm_nc)
    runoff = f2.variables['Runoff'][:]
    baseflow = f2.variables['Baseflow'][:]
    time1 = f2.variables['time'][:]
    flux = runoff + baseflow

    for run in runs:

        for basin in basins:
            print 'starting basin', basin, ' in run ', run
            # open the UH netcdf file
            basin_nc ='/raid/jhamman/temp_uh_files/'+run+'/'+basin+'_RASM_UH.nc'
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
            out_nc = 'Streamflows/'+basin+'_'+run+'_'+sim+'_streamflow.nc'
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
            f.description = 'Streamflow for RASM,created using UH file from '+run+' and fluxes from '+sim+'.'
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
            print '*** SUCCESS writing streamflows for ',out_nc, '***'

print 'Done'    

