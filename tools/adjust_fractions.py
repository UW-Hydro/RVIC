#!/usr/local/bin/python

import os
import sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import time as tm
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import argparse

def main():

    files, options = process_command_line()
    # read gridFile
    try:
        gridData, gridAttrs, gridGlobs = read_netcdf(files['gridFile'],verbose = options['verbose'])
    except:
        raise IOError('Need gridFile')
    
    # read in all the inputFiles, just store the flattened data (not the full / mostly empty grid)
    lf = len(files['inputFiles'])
    fracs = {}; uhs = {}; times = {}
    xcs = {}; ycs = {}
    xs = {}; ys = {}
    attrs = {}; globs = {}
    for i, fname in enumerate(files['inputFiles']):
        inputFile = os.path.join(files['inPath'],fname)
        if options['verbose']:
            print 'On file %i of %i' %(i+1,lf)

        fData, attrs[fname], globs[fname] = read_netcdf(inputFile,verbose=options['verbose'])

        if i==0:
            aggFracs = np.zeros(fData['fraction'].shape)
            aggUHs = np.zeros(fData['unit_hydrograph'].shape)

        y,x = ((fData['fraction']>0.0)*(gridData['mask']==1)).nonzero()
        
        aggFracs[y,x] +=fData['fraction'][y,x]
        aggUHs[:,y,x] +=fData['unit_hydrograph'][:,y,x]
        
        fracs[fname] = fData['fraction'][y,x]
        uhs[fname] = fData['unit_hydrograph'][:,y,x]
        times[fname] = fData['time']
        xcs[fname] = fData['xc'][y,x]
        ycs[fname] = fData['yc'][y,x]
        xs[fname] = x
        ys[fname] = y

    if options['verbose']:
        print 'Done reading input files'

    gridFracs = gridData['frac']
    diffFracs = aggFracs-gridFracs
    yi,xi = np.nonzero(aggFracs>gridFracs)
    ratioFracs = np.ones(diffFracs.shape)
    ratioFracs[yi,xi] = gridFracs[yi,xi]/aggFracs[yi,xi]

    if files['diagPath']:
        make_plot(gridFracs,"gridFracs",aggFracs,"aggFracs",diffFracs,'diffFracs=aggFracs-gridFracs',ratioFracs
                  ,'ratioFracs=gridFracs/aggFracs','Diagnostic Plot-0 for Aggregated Fractions'
                  ,os.path.join(files['diagPath'],'diag_0.png'))

        make_plot(gridData['mask'],"Land Mask",gridFracs,"gridFracs",aggUHs.sum(axis=0),'Sum of UHs'
                  ,gridData['mask']-aggUHs.sum(axis=0),'land+mask-uh_sum',
                  'Diagnostic Plot-1 for Aggregated Fractions',os.path.join(files['diagPath'],'diag_1.png'))

    if (files['outPath']):
        aggFracs2 = np.zeros(aggFracs.shape)
        
        for i, fname in enumerate(files['inputFiles']):
            if options['verbose']:
                print 'On file %i of %i' %(i,lf)
            x = xs[fname]
            y = ys[fname]
            out_fracs = fracs[fname]*ratioFracs[y,x]  #Adjust fracs based on ratioFracs 
            aggFracs2[y,x] += out_fracs
    
            # write flat uh file 
            if files['outPath']:
                outFile = os.path.join(files['outPath'],fname)
                if options['verbose']:
                    print 'Making flat uh file %s' %outFile

                write_flat_netcdf(outFile,times[fname],out_fracs,uhs[fname],
                                  x,y,xcs[fname],ycs[fname],globs[fname],attrs[fname])

        if files['diagPath']:
            make_plot(gridFracs,"gridFracs",aggFracs,"aggFracs",aggFracs2,"aggFracs2",aggFracs2-gridFracs
                      ,"aggFracs2-gridFracs",'Diagnostic Plot-2 for Aggregated Fractions'
                      ,os.path.join(files['diagPath'],'diag_2.png'))
    return

##################################################################################
## Read netCDF Inputs
## Read data from input netCDF.
##################################################################################
def read_netcdf(ncFile,vars = [],coords = False, verbose = False):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """
    if verbose:
        print 'Reading input data vars:', vars, ', from file:',ncFile
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    d={}
    a={}
    g={}
    if coords:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][coords]
            a[vars] = f.variables[vars].__dict__
        else:
            for var in vars:
                d[var] = f.variables[var][coords]
                a[var] = f.variables[var].__dict__
    else:
        if isinstance(vars,str):
            d[vars] = f.variables[vars][:]
            a[vars] = f.variables[vars].__dict__
        else:
            for var in vars:
                d[var] = f.variables[var][:]
                a[var] = f.variables[var].__dict__
    
    for attr in f.ncattrs():
        g[attr] = getattr(f,attr)
    
    f.close()
    
    return d,a,g

##################################################################################
##  Write output to netCDF
##  Writes out a netCDF4 data file containing the UH_S and fractions
##################################################################################
def write_netcdf(file,xc,xc_bnd,yc,yc_bnd,times,hydrographs,fractions,loc,Flist,velocity,diffusion,NODATA,verbose):
    """
    Write output to netCDF.  Writes out a netCDF4 data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    
    f = Dataset(file,'w', format='NETCDF4')

    # set dimensions
    time = f.createDimension('time', None)
    x = f.createDimension('x',xc.shape[1])
    y = f.createDimension('y',xc.shape[0])
    nv4 = f.createDimension('nv4',4)

    # initialize variables
    time = f.createVariable('time','f8',('time',))
    xcs = f.createVariable('xc','f8',('y','x',))
    ycs = f.createVariable('yc','f8',('y','x',))
    xc_bnds = f.createVariable('xc_bnds','f8',('y','x','nv4',))
    yc_bnds = f.createVariable('yc_bnds','f8',('y','x','nv4',))
    fraction = f.createVariable('fraction','f8',('y','x',),fill_value=NODATA)
    UHS = f.createVariable('unit_hydrograph','f8',('time','y','x',),fill_value=NODATA)

    # write attributes for netcdf
    f.description = 'Aggregated UH_S and Fraction Vars for full RASM domain'
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0] # prints the name of script used
    f.velocity = velocity
    f.diffusion = diffusion
    f.outlet_lon = loc[0]
    f.outlet_lat = loc[1]
    f.includes = str(len(Flist))+' files'

    ycs.long_name = 'latitude of grid cell center'
    ycs.standard_name = 'latitude'
    ycs.units = 'degrees_north'
    ycs._CoordinateAxisType = 'Lat'
    ycs.bounds = 'yc_bnds'

    xcs.long_name = 'longitude of grid cell center'
    xcs.standard_name = 'longitude'
    xcs.units = 'degrees_east'
    xcs._CoordinateAxisType = 'Lon'
    xcs.bounds = 'xc_bnds'

    time.standard_name = 'time'
    time.units = 'seconds'
    time.description = 'Seconds since initial impulse'
    time.calendar = 'proleptic_gregorian'

    UHS.units = 'unitless'
    UHS.description = 'unit hydrograph for each grid cell with respect to basin outlet location'
    
    fraction.units = 'unitless'
    fraction.description = 'fraction of grid cell contributing to guage location'

    # write data to variables initialized above
    time[:]= times
    xcs[:,:] = xc
    ycs[:,:] = yc
    xc_bnds[:,:,:] = xc_bnd
    yc_bnds[:,:,:] = yc_bnd
    UHS[:,:,:] = hydrographs
    fraction[:,:]= fractions
    f.close()

    return

def write_flat_netcdf(outFile,time,frac,uh,x,y,xc,yc,inGlobs,inAttrs):
    """
    Write a flattened uh file that includes the fractions, uhs, and flows
    """
    f = Dataset(outFile, 'w', format='NETCDF4')

    # set dimensions
    times = f.createDimension('time', len(time))
    npoints = f.createDimension('npoints', len(frac))
    
    # initialize variables
    times = f.createVariable('time','f8',('time',))
    fracs = f.createVariable('fraction','f8',('npoints',))
    xis = f.createVariable('xi','i4',('npoints',))
    yis = f.createVariable('yi','i4',('npoints',))
    xcs = f.createVariable('xc','f8',('npoints',))
    ycs = f.createVariable('yc','f8',('npoints',))
    uhs = f.createVariable('unit_hydrograph','f8',('time','npoints',))
    
    # deal with attributes
    f.description = 'Flattened uh/fraction grid file'
    f.history = 'Created ' + tm.ctime(tm.time())
    f.velocity = inGlobs['velocity']
    f.diffusion = inGlobs['diffusion']
    f.outlet_lon = inGlobs['outlet_lon']
    f.outlet_lat = inGlobs['outlet_lat']
    f.outlet_y = int(inGlobs['outlet_y'])
    f.outlet_x = int(inGlobs['outlet_x'])
    f.outlet_id = int(inGlobs['outlet_id'])
    try:
        f.includes = inGlobs['includes']
    except:
        pass
    
    times.standard_name = inAttrs['time']['standard_name']
    times.units = inAttrs['time']['units']
    times.calendar = inAttrs['time']['calendar']
    
    try:
        fracs.units = inAttrs['fraction']['units']
    except:
        fracs.units = '%'
    fracs.description = inAttrs['fraction']['description']
    
    uhs.units = inAttrs['unit_hydrograph']['units']
    uhs.description = inAttrs['unit_hydrograph']['description']
    
    xis.standard_name = 'x_ind'
    xis.description = 'x index location'
    
    yis.standard_name = 'y_ind'
    yis.description = 'y index location'
    
    xcs.standard_name =inAttrs['xc']['standard_name']
    xcs.long_name = inAttrs['xc']['long_name']
    xcs.units =inAttrs['xc']['units']
    
    ycs.standard_name =inAttrs['yc']['standard_name']
    ycs.long_name = inAttrs['yc']['long_name']
    ycs.units =inAttrs['yc']['units']
    
    times[:] = time
    fracs[:] = frac
    uhs[:,:] = uh
    xis[:] = x
    yis[:] = y
    xcs[:] = xc
    ycs[:] = yc

    f.close()
    
    return

def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type adjust_fractions.py -f
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--outPath", type=str, help="Output path for flat uh files with adjusted fractions")
    parser.add_argument("--gridFile", type=str, help="Grid File containing full domain fractions variable ")
    parser.add_argument("--inputFiles", help="Input netcdf grid(s) containing fraction/uh data")
    parser.add_argument("--verbose",help="Make script verbose",action="store_true")
    parser.add_argument("--diagPath",type=str,help="Path to place diagnostic outputs")

    args = parser.parse_args()

    options = {}
    options['verbose'] = args.verbose
            
    files={}
    temp = glob.glob(args.inputFiles)
    try:
        files['gridFile'] = args.gridFile
    except:
        files['gridFile'] = False
    try:
        files['diagPath'] = args.diagPath
        if not os.path.exists(files['diagPath']):
            print 'making diagnostic directory'
            os.makedirs(files['diagPath'])
    except:
        files['diagPath'] = False

    files['inputFiles'] = []
    for fi in temp:
        files['inputFiles'].append(os.path.basename(fi))
    files['inPath'] = os.path.dirname(fi)

    try:
        files['outPath'] = args.outPath
        if not os.path.exists(files['outPath']):
            print 'output directory'
            os.makedirs(files['outPath'])
    except:
        files['outPath'] = False
    
    return files,options

def make_plot(d0,t0,d1,t1,d2,t2,d3,t3,suptitle,path_out):
    """
    Make a 4 pannel plot that shows map views of the domain of interest
    No projection is defined.  
    """
    fig0 = plt.figure()
    grid = ImageGrid(fig0, 111, # similar to subplot(111)
                     nrows_ncols = (2, 2), # creates 2x2 grid of axes
                     axes_pad=0.4, # pad between axes in inch.
                     share_all=True, # share axes
                     cbar_mode='each')
    
    im = grid[0].pcolor(d0)
    grid[0].set_title(t0)
    grid.cbar_axes[0].colorbar(im)
    grid[0].axis([0,d0.shape[1],0,d0.shape[0]])
    
    im = grid[1].pcolor(d1)
    grid[1].set_title(t1)
    grid.cbar_axes[1].colorbar(im)
    grid[0].axis([0,d1.shape[1],0,d1.shape[0]])
     
    im = grid[2].pcolor(d2)
    grid[2].set_title(t2)
    grid.cbar_axes[2].colorbar(im)
    grid[0].axis([0,d2.shape[1],0,d2.shape[0]])
    
    im = grid[3].pcolor(d3)
    grid[3].set_title(t3) 
    grid.cbar_axes[3].colorbar(im)
    grid[0].axis([0,d3.shape[1],0,d3.shape[0]])

    
    fig0.suptitle(suptitle)
    
    fig0.savefig(path_out)
    fig0.clf()
    return

##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
