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

earthRadius = 6.37122e6 # meters

def main():

    files, options = process_command_line()
    # read gridFile
    try:
        gridData, gridAttrs, gridGlobs = read_netcdf(files['gridFile'], verbose = options['verbose'])
    except:
        raise IOError('Need gridFile')

    cell_id_grid = np.flipud(np.arange(gridData['yc'].size).reshape(gridData['yc'].shape))

    # read in all the inputFiles, just store the flattened data (not the full grid)
    lf = len(files['inputFiles'])
    attrs_dict = {}
    globs_dict = {}

    # point specific inputs
    fraction_dict = {}
    unit_hydrograph_dict = {}
    time_dict = {}
    lon_point_dict = {}
    lat_point_dict = {}
    x_ind_point_dict = {}
    y_ind_point_dict = {}
    cell_id_point_dict = {}
    
    # outlet specific inputs
    cell_id_outlet_dict = {}
    lon_outlet_dict = {}
    lat_outlet_dict = {}
    x_ind_outlet_dict = {}
    y_ind_outlet_dict = {}

    # Now loop over all input files, storing data in above dictionaries
    for i, fname in enumerate(files['inputFiles']):
        inputFile = os.path.join(files['inPath'], fname)
        if options['verbose']:
            print 'On file %i of %i' % (i+1, lf)

        fData, attrs_dict[fname], globs_dict[fname] = read_netcdf(inputFile, verbose = options['verbose'])
        if i == 0:
            aggFracs = np.zeros(fData['fraction'].shape)
            aggUHs = np.zeros(fData['unit_hydrograph'].shape)

        y, x = ((fData['fraction'] > 0.0)*(gridData['mask'] == 1)).nonzero()
        
        aggFracs[y, x] += fData['fraction'][y, x]
        aggUHs[:, y, x] += fData['unit_hydrograph'][:, y, x]
        
        fraction_dict[fname] = fData['fraction'][y, x]
        unit_hydrograph_dict[fname] = fData['unit_hydrograph'][:, y, x]
        time_dict[fname] = fData['time']
        lat_point_dict[fname] = fData['xc'][y, x]
        lon_point_dict[fname] = fData['yc'][y, x]
        cell_id_point_dict[fname] = cell_id_grid[y, x]
        x_ind_point_dict[fname] = x
        y_ind_point_dict[fname] = y

        cell_id_outlet_dict[fname] = globs_dict[fname]['outlet_id']
        x_ind_outlet_dict[fname] = globs_dict[fname]['outlet_x']
        y_ind_outlet_dict[fname] = globs_dict[fname]['outlet_y']
        lon_outlet_dict[fname] = globs_dict[fname]['outlet_lon']
        lat_outlet_dict[fname] = globs_dict[fname]['outlet_lat']

    if options['verbose']:
        print 'Done reading input files'

    gridFracs = gridData['frac']
    diffFracs = aggFracs-gridFracs
    yi, xi = np.nonzero(aggFracs > gridFracs)
    ratioFracs = np.ones(diffFracs.shape)
    ratioFracs[yi, xi] = gridFracs[yi, xi]/aggFracs[yi,xi]

    if files['diagPath']:
        make_plot(gridFracs, "gridFracs", aggFracs, "aggFracs", diffFracs, 
                  'diffFracs=aggFracs-gridFracs', ratioFracs, 
                  'ratioFracs=gridFracs/aggFracs', 
                  'Diagnostic Plot-0 for Aggregated Fractions', 
                  os.path.join(files['diagPath'], 'diag_0.png'))

        make_plot(gridData['mask'], "Land Mask", gridFracs, "gridFracs", 
                  aggUHs.sum(axis=0), 'Sum of UHs', gridData['mask']-aggUHs.sum(axis=0),
                  'land+mask-uh_sum', 'Diagnostic Plot-1 for Aggregated Fractions',
                  os.path.join(files['diagPath'], 'diag_1.png'))

    if (files['outPath']):
        aggFracs2 = np.zeros(aggFracs.shape)
        
        for i, fname in enumerate(files['inputFiles']):
            if options['verbose']:
                print 'On file %i of %i' %(i,lf)
            x = x_ind_point_dict[fname]
            y = y_ind_point_dict[fname]
            out_fracs = fraction_dict[fname]*ratioFracs[y, x]  #Adjust fracs based on ratioFracs 
            aggFracs2[y, x] += out_fracs
    
            # write flat uh file 
            if files['outPath']:
                
                # Subset if needed
                if options['subset']:
                    offset, out_uh, full_length = subset(unit_hydrograph_dict[fname], options['subset'], 
                                                         options['threshold'])
                    time = np.arange(options['subset'])
                else:
                    offset = np.zeros(len(out_fracs))
                    full_length = unit_hydrogaph_dict[fname].shape[0]
                    out_uh = unit_hydrogaph_dict[fname]
                    time = np.arange(full_length)

                if options['flat']:
                    outFile = os.path.join(files['outPath'],fname)
                    if options['verbose']:
                        print 'Making flat uh file %s' %outFile
                    write_flat_netcdf(outFile, time, offset, full_length, out_fracs, 
                                      out_uh, x, y, xcs[fname], ycs[fname],
                                      options['out_format'], globs[fname], attrs[fname])

                if options['rvic_params']:
                    j = i + 1
                    if i == 0:
                        uh_points = out_uh
                        frac_points = out_fracs
                        times = np.arange(options['subset'])
                        lon_points = lon_point_dict[fname]
                        lat_points = lat_point_dict[fname]
                        x_ind_points = x_ind_point_dict[fname]
                        y_ind_points = y_ind_point_dict[fname]
                        cell_id_points = cell_id_point_dict[fname]
                        t_offset = offset
                        point2outlet_inds = np.zeros(len(out_fracs))
                        
                        # outlet specific inputs
                        cell_id_outlets = np.array(cell_id_outlet_dict[fname])
                        lon_outlets = np.array(lon_outlet_dict[fname])
                        lat_outlets = np.array(lat_outlet_dict[fname])
                        x_ind_outlets = np.array(x_ind_outlet_dict[fname])
                        y_ind_outlets = np.array(y_ind_outlet_dict[fname])
                        outlet_nums = np.array(j)

                        # get a few global values
                        t_subset_length = options['subset']
                        t_full_length = full_length
                        t_timestep = 'placeholder'

                    else:
                        uh_points = np.append(uh_points, out_uh, axis=1)
                        frac_points = np.append(frac_points, out_fracs)
                        lon_points = np.append(lon_points, lon_point_dict[fname])
                        lat_points = np.append(lat_points, lat_point_dict[fname])
                        x_ind_points = np.append(x_ind_points, x_ind_point_dict[fname])
                        y_ind_points = np.append(y_ind_points, y_ind_point_dict[fname])
                        cell_id_points = np.append(cell_id_points, cell_id_point_dict[fname])
                        t_offset = np.append(t_offset, offset)
                        point2outlet_inds = np.append(point2outlet_inds, np.zeros_like(offset)+j)
                        
                        # outlet specific inputs
                        cell_id_outlets = np.append(cell_id_outlets, cell_id_outlet_dict[fname])
                        lon_outlets = np.append(lon_outlets, lon_outlet_dict[fname])
                        lat_outlets = np.append(lat_outlets, lat_outlet_dict[fname])
                        x_ind_outlets = np.append(x_ind_outlets, x_ind_outlet_dict[fname])
                        y_ind_outlets = np.append(y_ind_outlets, y_ind_outlet_dict[fname])     
                        outlet_nums = np.append(outlet_nums, j)        

        if options['rvic_params']:
            print 'writing param file...'

            outFile = os.path.join(files['outPath'], options['rvic_params'])
            
            uh_points *= frac_points * gridData['area'][y_ind_points,x_ind_points] * earthRadius * earthRadius

            write_param_file(outFile, options['out_format'], times, t_subset_length, t_full_length, 
                             t_timestep, uh_points, frac_points, lon_points, lat_points,
                             cell_id_points, x_ind_points, y_ind_points,
                             t_offset, point2outlet_inds, cell_id_outlets, lon_outlets,
                             lat_outlets, x_ind_outlets, y_ind_outlets, outlet_nums)
        if files['diagPath']:
            make_plot(gridFracs, "gridFracs", aggFracs, "aggFracs", aggFracs2,
                      "aggFracs2", aggFracs2-gridFracs , "aggFracs2-gridFracs",
                      'Diagnostic Plot-2 for Aggregated Fractions', 
                      os.path.join(files['diagPath'], 'diag_2.png'))
    return

def write_param_file(outFile, out_format, times, t_subset_length, t_full_length,
                     t_timestep, uh_points, frac_points, lon_points, lat_points,
                     cell_id_points, x_ind_points, y_ind_points,
                     time_offset_points, point2outlet_inds, cell_id_outlets, lon_outlets,
                     lat_outlets, x_ind_outlets, y_ind_outlets, outlet_nums):
    
    f = Dataset(outFile, 'w', format = out_format)

    # set dimensions
    time = f.createDimension('time', len(times))
    n_points = f.createDimension('n_points', len(frac_points))
    n_outlets = f.createDimension('n_outlets', len(outlet_nums))
    
    # Variables
    time = f.createVariable('time','i8',('time',))
    time.standard_name = 'time'    
    time.units = 'timesteps'
    time.subset_length = t_subset_length
    time.full_time_length = t_full_length
    time.timestep = t_timestep
    time[:] = times

    uh_point = f.createVariable('uh_point', 'f8', ('time', 'n_points',))
    uh_point.long_name = 'Unit Hydrographs'
    uh_point.units = 'm^2'
    uh_point.description = 'Subset and flattened unit hydrograph'
    uh_point[:,:] = uh_points

    # frac_point = f.createVariable('frac_point', 'f8', ('n_points',))
    # frac_point.long_name = 'Fraction'
    # frac_point.units = 'unitless'
    # frac_point.description = 'Fraction of grid cell contributing to outlet point'
    # frac_point[:] = fractions

    cell_id_point = f.createVariable('cell_id_point', 'i8', ('n_points',))
    cell_id_point.long_name = 'Cell ID Point'
    cell_id_point.units = 'unitless'
    cell_id_point.description = 'Land Model Grid Cell ID'
    cell_id_point[:] = cell_id_points

    y_ind_point = f.createVariable('y_ind_point', 'i8',('n_points',))
    y_ind_point.long_name = 'Y Index Location'
    y_ind_point.units = 'unitless'
    y_ind_point.description = 'Y Index Location of Origin Grid Cell'
    y_ind_point[:] = y_ind_points

    x_ind_point = f.createVariable('x_ind_point', 'i8', ('n_points',))
    x_ind_point.long_name = 'X Index Location'
    x_ind_point.units = 'unitless'
    x_ind_point.description = 'X Index Location of Origin Grid Cell'
    x_ind_point[:] = x_ind_points

    lon_point = f.createVariable('lon_point', 'f8', ('n_points',))
    lon_point.long_name = 'longitude coordinate'
    lon_point.units = 'degrees_east'
    lon_point.description = 'Longitude Coordinate of Origin Grid Cell'
    lon_point[:] = lon_points

    lat_point = f.createVariable('lat_point', 'f8', ('n_points',))
    lat_point.long_name = 'latitude coordinate'
    lat_point.units = 'degrees_north'
    lat_point.description = 'Latitude Coordinate of Origin Grid Cell'
    lat_point[:] = lat_points

    time_offset_point = f.createVariable('t_offset_point', 'i8', ('n_points',))
    time_offset_point.long_name = 'time_offset'
    time_offset_point.units = 'timesteps'
    time_offset_point.description = 'Number of ommited leading timesteps'
    time_offset_point[:] = time_offset_points

    point2outlet_index = f.createVariable('point2outlet_index', 'i8', ('n_points',))
    point2outlet_index.long_name = 'Point to outlet index mapping'
    point2outlet_index.description = '1D outlet index associated with source point'
    point2outlet_index[:] = point2outlet_inds

    cell_id_outlet = f.createVariable('cell_id_outlet', 'i8', ('n_outlets',))
    cell_id_outlet.long_name = 'Outlet ID Point'
    cell_id_outlet.units = 'unitless'
    cell_id_outlet.description = 'Outlet Grid Cell ID'
    cell_id_outlet[:] = cell_id_outlets

    outlet_num = f.createVariable('outlet_num', 'i8', ('n_outlets',))
    outlet_num.long_name = 'Outlet Index'
    outlet_num.description = 'Outlet Point Index'
    outlet_num[:] = outlet_nums

    x_ind_outlet = f.createVariable('x_ind_outlet', 'i8', ('n_outlets',))
    x_ind_outlet.long_name = 'X Index Location'
    x_ind_outlet.units = 'unitless'
    x_ind_outlet.description = 'X Index Location of Outlet Grid Cell'
    x_ind_outlet[:] = x_ind_outlets
       
    y_ind_outlet = f.createVariable('y_ind_outlet', 'i8', ('n_outlets',))
    y_ind_outlet.long_name = 'Y Index Location'
    y_ind_outlet.units = 'unitless'
    y_ind_outlet.description = 'Y Index Location of Outlet Grid Cell'
    y_ind_outlet[:] = y_ind_outlets

    lon_outlet = f.createVariable('lon_outlet', 'f8', ('n_outlets',))
    lon_outlet.long_name = 'longitude coordinate'
    lon_outlet.units = 'degrees_east'
    lon_outlet.description = 'Longitude Coordinate of Outlet Grid Cell'
    lon_outlet[:] = lon_outlets

    lat_outlet = f.createVariable('lat_outlet', 'f8', ('n_outlets',))
    lat_outlet.long_name = 'latitude coordinate'
    lat_outlet.units = 'degrees_north'
    lat_outlet.description = 'Latitude Coordinate of Outlet Grid Cell'
    lat_outlet[:] = lat_outlets

    # Globals
    f.description = 'Flattened uh/fraction RVIC parameter file'
    f.history = 'Created ' + tm.ctime(tm.time())
    f.subset_length = t_subset_length
    f.full_time_length = t_full_length
    f.timestep = t_timestep

    f.close()

def subset(uh, subset, threshold):
    full_length = uh.shape[0]
    offset = np.empty(uh.shape[1], dtype = int)
    out_uh = np.zeros((subset, uh.shape[1]))
    for i in xrange(uh.shape[1]):
        offset[i]  = np.nonzero(uh[:, i] > threshold)[0][0]     # find position of first index > threshold
        end = np.minimum(uh.shape[0],offset[i]+subset)        # find end point
        out_uh[:end-offset[i], i] = uh[offset[i]:end, i]/uh[offset[i]:end, i].sum() # clip and normalize
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
    parser.add_argument("--rvic_param", type = str, help = "RVIC paramameter file name", default = False)
    parser.add_argument("--flat", type = bool, default = False, help = "write individual flat files")

    args = parser.parse_args()

    options = {}
    options['verbose'] = args.verbose
    options['out_format'] = args.out_format
    options['rvic_params'] = args.rvic_param
    options['flat'] = args.flat
            
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
