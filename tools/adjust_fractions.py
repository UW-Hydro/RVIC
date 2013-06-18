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
        gridData, gridAttrs, gridGlobs = read_netcdf(files['gridFile'], verbose = options['verbose'])
    except:
        raise IOError('Need gridFile')
    
    # read in all the inputFiles, just store the flattened data (not the full / mostly empty grid)
    lf = len(files['inputFiles'])
    fracs = {}; uhs = {}; times = {}
    xcs = {}; ycs = {}
    xs = {}; ys = {}
    attrs = {}; globs = {}
    for i, fname in enumerate(files['inputFiles']):
        inputFile = os.path.join(files['inPath'], fname)
        if options['verbose']:
            print 'On file %i of %i' %(i+1,lf)

        fData, attrs[fname], globs[fname] = read_netcdf(inputFile,verbose = options['verbose'])

        if i==0:
            aggFracs = np.zeros(fData['fraction'].shape)
            aggUHs = np.zeros(fData['unit_hydrograph'].shape)

        y,x = ((fData['fraction']>0.0)*(gridData['mask']==1)).nonzero()
        
        aggFracs[y,x] += fData['fraction'][y,x]
        aggUHs[:,y,x] += fData['unit_hydrograph'][:,y,x]
        
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
        make_plot(gridFracs, "gridFracs", aggFracs, "aggFracs", diffFracs, 
                  'diffFracs=aggFracs-gridFracs',ratioFracs, 
                  'ratioFracs=gridFracs/aggFracs', 
                  'Diagnostic Plot-0 for Aggregated Fractions', 
                  os.path.join(files['diagPath'],'diag_0.png'))

        make_plot(gridData['mask'], "Land Mask",gridFracs, "gridFracs", 
                  aggUHs.sum(axis=0),'Sum of UHs', gridData['mask']-aggUHs.sum(axis=0),
                  'land+mask-uh_sum', 'Diagnostic Plot-1 for Aggregated Fractions',
                  os.path.join(files['diagPath'],'diag_1.png'))

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

                # Subset if needed
                if options['subset']:
                    offset, out_uh, full_length = subset(uhs[fname], options['subset'], 
                                                         options['threshold'])
                    time = np.arange(options['subset'])
                else:
                    offset = np.zeros(len(out_fracs))
                    full_length = uhs[fname].shape[0]
                    out_uh = uhs[fname]
                    time = np.arange(full_length)

                write_flat_netcdf(outFile, time, offset, full_length, out_fracs, 
                                  out_uh, x, y, xcs[fname], ycs[fname],
                                  options['out_format'], globs[fname], attrs[fname])

        if files['diagPath']:
            make_plot(gridFracs, "gridFracs", aggFracs,"aggFracs", aggFracs2,
                      "aggFracs2",aggFracs2-gridFracs ,"aggFracs2-gridFracs",
                      'Diagnostic Plot-2 for Aggregated Fractions'
                      ,os.path.join(files['diagPath'],'diag_2.png'))
    return

def subset(uh, subset, threshold):
    full_length = uh.shape[0]
    offset = np.empty(uh.shape[1])
    out_uh = np.zeros((subset, uh.shape[1]))
    for i in xrange(uh.shape[1]):
        offset[i]  = np.nonzero(uh[:,i]>threshold)[0][0]     # find position of first index > threshold
        end = np.minimum(uh.shape[0],offset[i]+subset)        # find end point
        out_uh[:end-offset[i],i] = uh[offset[i]:end,i]/uh[offset[i]:end,i].sum() # clip and normalize
    return offset, out_uh, full_length

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

def write_flat_netcdf(outFile, time, offset, full_length, frac, uh, 
                      x, y, xc, yc, out_format, inGlobs, inAttrs):
    """
    Write a flattened uh file that includes the fractions, uhs, and flows
    """
    f = Dataset(outFile, 'w', format=out_format)

    # set dimensions
    times = f.createDimension('time', len(time))
    npoints = f.createDimension('npoints', len(frac))
    
    # initialize variables
    times = f.createVariable('time','f8',('time',))
    time_offset = f.createVariable('time_offset', 'i8', ('npoints'))
    fracs = f.createVariable('fraction','f8',('npoints',))
    xis = f.createVariable('xi','i8',('npoints',))
    yis = f.createVariable('yi','i8',('npoints',))
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
    
    times.standard_name = 'time'
    times.units = 'timesteps'
    times.full_length = full_length 

    time_offset.standard_name = 'time_offset'
    time_offset.units = 'timesteps'
    time_offset.description = 'Number of leading timesteps'
    
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
    xcs.units = inAttrs['xc']['units']
    
    ycs.standard_name =inAttrs['yc']['standard_name']
    ycs.long_name = inAttrs['yc']['long_name']
    ycs.units =inAttrs['yc']['units']
    
    times[:] = time
    time_offset[:] = offset
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
    parser.add_argument("--outPath", type = str, help = "Output path for flat uh files with adjusted fractions")
    parser.add_argument("--gridFile", type = str, help = "Grid File containing full domain fractions variable ")
    parser.add_argument("--inputFiles", help = "Input netcdf grid(s) containing fraction/uh data")
    parser.add_argument("--verbose", help = "Make script verbose", action = "store_true")
    parser.add_argument("--diagPath", type = str,help="Path to place diagnostic outputs")
    parser.add_argument("--subset", type = int, default = 0, help = "Number of timesteps to include of each UH after first point > threshold.")
    parser.add_argument("--threshold", type = float, default = 0, help = "Threshold to use to select the unit hydrograph.")
    parser.add_argument("--out_format", type = str, help = "Output netcdf format for flat uh files", default = "NETCDF4_CLASSIC")

    args = parser.parse_args()

    options = {}
    options['verbose'] = args.verbose
    options['out_format'] = args.out_format
            
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

    if (args.subset and args.threshold):
        options['subset'] = args.subset
        options['threshold'] = args.threshold
    else:
        options['subset'] = False
        options['threshold'] = False

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
    
    return files, options

def make_plot(d0, t0, d1, t1, d2, t2, d3, t3, suptitle, path_out):
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
