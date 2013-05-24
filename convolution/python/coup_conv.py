#!/usr/local/bin/python
"""
This is the convolution routine developed in preparation in coupling RVIC to RASM

Written by Joe Hamman, May 2013
"""
from netCDF4 import Dataset, num2date, date2num
import numpy as np
import glob, os, shutil
import ConfigParser
import argparse
import time as tm
from collections import OrderedDict, deque

def main(re = 6.37122e6,rho_h20=1000):
    """
    The basic workflow of the main routine of the coupled model convolution is:
    1.  Read configuration file
    2.  Read Grid File
    3.  Load the unit hydrograph files (put into point_dict)
        - Load the initial state file and put in convolution ring
    4.  Loop over flux files
        - Combine the Baseflow and Runoff Variables
        - Adjust units as necessary
        - Do convolution
        - Write output files
    """
    
    # read command line and config file
    config_file = process_command_line()
    Config, uh_files,flux_files,grid_file,out_path,initial_state,outputs,options = process_config_file(config_file)

    # find the grid cell area in square meters
    try:
        f=Dataset(grid_file,'r')
        if (f.variables['area'].units in ["rad2", "radians2", "radian2","rad^2",
                                          "radians^2", "rads^2","radians squared",
                                          "square-radians"]):
            area = f.variables['area'][:]*re*re
        elif f.variables['area'].units in ["m2","m^2","meters^2","meters2",
                                           "square-meters","meters squared"]:
            area = f.variables['area'][:]
        elif f.variables['area'].units in ["km2","km^2","kilometers^2",
                                           "kilometers2","square-kilometers",
                                           "kilometers squared"]:
            area = f.variables['area'][:]/1000./1000.
        elif f.variables['area'].units in ["mi2","mi^2","miles^2","miles",
                                           "square-miles","miles squared"]:
            area = f.variables['area'][:]*2.59e+6
        elif f.variables['area'].units in ["acres", "ac","ac."]:
            area = f.variables['area'][:]*4046.86
        else:
            print "WARNING: UNKNOWN AREA UNITS (%s), ASSUMING THEY ARE IN SQUARE METERS" %f.variables['area'].units
        out_dict = {}
        if outputs['out_type']=='grid':
            out_dict['longitudes'] = f.variables['xc'][:]
            out_dict['latitudes'] = f.variables['yc'][:]
            shape = area.shape
        else:
            shape = (len(uh_files),)
        f.close()
    except:
        e = sys.exc_info()[0]
        print "Either no grid_file with area variable was provided or there was a problem loading the file"
        print "In the future, we can calculate the area of the grid cells based on their spacing (if on regular grid)"
        raise e

    if options['verbose']:
        print 'reading input files'
    point_dict, out_dict = make_point_dict(uh_files,area,out_dict,initial_state)

    if options['verbose']:
        print 'starting convolution now'
    time_dict = {}
    i = 0
    while flux_files:
        ff = flux_files.popleft()
        i += 1
        # Check to see if it's time to save a state file
        if any(date in ff for date in outputs['state']):
            print 'making statefile from %s' %ff
            return_state = True
        else:
            return_state = False

        f=Dataset(ff,'r')
        # read time step
        time_dict['time_step'] = f.variables['time'][:]
        if i == 1:
            time_dict['units'] = f.variables['time'].units
            time_dict['cal'] = f.variables['time'].calendar
            time_dict['long_name'] = f.variables['time'].long_name

            # convert to m3/s
            if f.output_frequency=='hourly' and f.output_mode=='instantaneous':
                div = 1200 # assumes vic timestep of 20min
            elif f.output_frequency=='hourly' and f.output_mode=='averaged':
                div = 3600  #averaged should really mean accumulated here
                time_dict['output_frequency'] = 'hourly'
            elif f.output_frequency=='dailyy' and f.output_mode=='averaged': # note the typo in dailyy
                div = 86400
            else:
                print 'Unexpected flux output frequency %s, assuming it is hourly accumulated' %f.output_frequency
                div = 3600
                
            if f.variables['Runoff'].units == 'mm':
                div *= 1000
            elif f.variables['Runoff'].units == 'cm':
                div *= 100
            else:
                print 'Unexpected flux units %s, assuming they are mm' %f.variables['Runoff'].units
                div *= 1000
        
        # Get the fluxes and convert to m3
        f.variables["Runoff"].set_auto_maskandscale(False)
        f.variables['Baseflow'].set_auto_maskandscale(False)
        flux=(f.variables['Runoff'][:]+f.variables['Baseflow'][:])/div

        f.close()

        # do the covolutions for this timestep
        point_dict,out_flow,out_state,time_dict = convolve(point_dict,time_dict,flux,
                                                           return_state,shape=shape)

        # write this timestep's streamflows to out_name
        if outputs['out_type'] != "false":
            out_name = os.path.join(out_path,os.path.split(ff)[1])
            write_output(out_name,out_flow,out_dict,time_dict,"streamflow",options,shape=shape)
        
        # write this timestep's state
        if return_state:
            state_name = os.path.join(out_path,'state_'+os.path.split(ff)[1])
            restart_name = os.path.join(out_path,'restart_'+os.path.split(ff)[1][:-2]+'cfg')
            write_output(state_name,out_state,out_dict,time_dict,"state",options,shape=shape)

            # make an associated restart file
            Config = write_restart(Config,state_name,restart_name,flux_files)
    return

def write_restart(Config,state_name,restart_name,flux_files):
    """
    Write a restart configuration file for startup of next simulation, use current state
    """
    # use the configuration parser to update the fields in restart_file
    Config.set('Paths', 'flux_files', ",".join(flux_files))
    Config.set('Paths', 'initial_state',state_name)
    # Writing our configuration file to 'example.cfg'
    with open(restart_name, 'wb') as configfile:
        Config.write(configfile)
    return Config

def write_output(out_name,out_flow,out_dict,time_dict,out_type,options,shape):
    """
    Write output file
    This routine is setup to handle the creation of streamflow or state files
    in grid or array formats.
    """
    if options['verbose']:
        print 'writing %s' %out_name
    f = Dataset(out_name,'w',format='NETCDF4')
    
    # Items that apply for all cases
    f.history = 'Created ' + tm.ctime(tm.time())
    f.source = 'Streamflow convolution program - coup_conv.py'
    time = f.createDimension('time', None)
    time = f.createVariable('time','f8',('time',))
    time.units = time_dict['units']
    time.calendar = time_dict['cal']
    time.longname = time_dict['long_name']
    time.type_prefered = 'double'
    if out_type=='state':
        time[:] = time_dict['out_state_time']
    else:
        time[:] = time_dict['time_step']
    
    if len(shape)==1:
        points = f.createDimension('point', shape[0])
        
        xis = f.createVariable('xi','i8',('point',))
        xis.standard_name = "x_outlet"
        xis.long_name = "X Grid Location of Outlet"
        xis.units = "grid_cells"
        xis[:] = out_dict['outlet_xs']
        
        yis = f.createVariable('yi','i8',('point',))
        yis.standard_name = "y_outlet"
        yis.long_name = "Y Grid Location of Outlet"
        yis.units = "grid_cells"
        yis[:] = out_dict['outlet_ys']
        
        lat = f.createVariable('latitudes','f8',('point',))
        lat.standard_name = "latitude"
        lat.long_name = "latitude of outlet grid cell center"
        lat.units = "degrees_north"
        lat[:] = out_dict['lats']
        
        lon = f.createVariable('longitudes','f8',('point',))
        lon.standard_name = "longitude"
        lon.long_name = "longitude of outlet grid cell center"
        lon.units = "degrees_east"
        lon[:] = out_dict['lons']
        
        flow = f.createVariable('Streamflow','f8',('time','point',))
        flow.description = 'Streamflow'
        flow.units = 'm^3/s'
        if out_type=='state':
            flow[:,:] = out_flow
        else:
            flow[0,:] = out_flow
    else:
        # Put all data into a grid
        x = f.createDimension('x', shape[1])
        y = f.createDimension('y', shape[0])
        
        lat = f.createVariable('latitudes','f8',('y','x',))
        lat.standard_name = "latitude"
        lat.long_name = "latitude of grid cell center"
        lat.units = "degrees_north"
        lat[:,:] = out_dict['latitudes']
        
        lon = f.createVariable('longitudes','f8',('y','x',))
        lon.standard_name = "longitude"
        lon.long_name = "longitude of grid cell center"
        lon.units = "degrees_east"
        lon[:,:] = out_dict['longitudes']
        
        flow = f.createVariable('Streamflow','f8',('time','y','x',))
        flow.description = 'Streamflow'
        flow.units = 'm^3/s'
        if out_type=='state':
            flow[:,:,:] = out_flow
        else:
            flow[0,:,:] = out_flow
            
    # write attributes for netcdf
    if out_type=='state':
        f.description = 'Streamflow state'
    else:
        f.description = 'Streamflow'
    f.close()

def make_point_dict(uh_files,area,out_dict,initial_state=None):
    """
    Read the initial state file if present
    Open all the unit hydrograph grids and store in dictionary
    Include location indecies
    Turn IRFs to true Unit Hydrographs 
    Setup convolution structures in save dictionary
    """
    # Create an ordered dictionary so that we can trust that the outputs will always be the same
    point_dict = OrderedDict()
    
    for i,uh_file in enumerate(uh_files):
        d = {}
        f=Dataset(uh_file,'r')
        if f.variables['time'].units =="seconds since 0-01-01 00:00:00":
            # convert to ordinal of day since...
            d['time'] = f.variables['time'][:]/86400
        else:
            # for now assume they are in days since 0-01-01 00:00:00
            d['time'] = f.variables['time'][:]
        
        # Get basin indicies
        d['xi'] = f.variables['xi'][:]
        d['yi'] = f.variables['yi'][:]

        # make unit hydrograph (no longer has volume of 1)
        frac = f.variables['fraction'][:]
        d['uh']=f.variables['unit_hydrograph'][:]*frac*area[d['yi'],d['xi']]

        # Get grid outlet locations
        d['x'] = f.outlet_x
        d['y'] = f.outlet_y
        d['lat'] = f.outlet_lat
        d['lon'] = f.outlet_lon
        f.close()
        
        d['end'] = len(d['uh'])-1
        d['now'] = 0

        # If there is an inital state, put that in the ring, if now, make a ring of zeros
        if not initial_state:
            d['ring'] = np.zeros(len(d['uh']))

        # store each individual point dictionary in the larger point_dict
        key = (d['y'],d['x'])
        point_dict[key]=d

    if initial_state:
        print "Reading Initial State File: %s" %initial_state
        f=Dataset(initial_state,'r')
        state = f.variables['Streamflow'][:]

        if state.ndim==3:
            for i,(key,d) in enumerate(point_dict.iteritems()):
                point_dict[key]['ring']=state[:,key[0],key[1]]
        else:
            for i in xrange(state.shape[1]):
                x_outlets = f.variables['xi'][:]
                y_outlets = f.variables['yi'][:]
                key = (y_outlets[i],x_outlets[i])
                point_dict[key]['ring']=state[:,i]
        f.close()

    # Now make a few numpy arrays from the point dict that will go in each output file
    out_dict['lats'] = np.zeros(len(point_dict))
    out_dict['lons'] = np.zeros(len(point_dict))
    out_dict['outlet_xs'] = np.zeros(len(point_dict))
    out_dict['outlet_ys'] = np.zeros(len(point_dict))
    for i,(key,d) in enumerate(point_dict.iteritems()):
        out_dict['lats'][i] = d['lat']
        out_dict['lons'][i] = d['lon']
        out_dict['outlet_ys'][i] = d['y']
        out_dict['outlet_xs'][i] = d['x']
        
    return point_dict, out_dict

def convolve(point_dict,time_dict,flux,return_state,shape):
    """
    This convoluition funciton works by looping over all points and doing the
    convolution one timestep at a time.  This is accomplished by creating an
    convolution ring.  Contributing flow from each timestep is added to the
    convolution ring.  The convolution ring is unwrapped when state is being saved.    
    """
    out_flow = np.zeros(shape)
    if out_flow.ndim==1:
        out_type='array'
    else:
        out_type='grid'
            
    for i,(point,d) in enumerate(point_dict.iteritems()):
        if i == 0:
            if return_state and out_type=='array':
                out_state = np.zeros((len(d['uh']),shape[0]))
                time_dict['out_state_time'] = d['time']+ time_dict['time_step']
            elif return_state and out_type=='grid':
                out_state = np.zeros((len(d['uh']),shape[0],shape[1]))
                time_dict['out_state_time'] = d['time']+ time_dict['time_step']
            else:
                out_state = None
                time_dict['out_state_time'] = None
        
        # Get the convolved hydrograph from the flux and add to convolution ring
        d['ring'] += (flux[:,d['yi'],d['xi']]*d['uh']).sum(axis=1)
        
        # Store the streamflow for this timestep
        if out_type=='array':
            out_flow[i] = d['ring'][0]
        elif out_type=='grid':
            out_flow[d['y'],d['x']]=d['ring'][0]
        
        #Set the current ring value to 0
        d['ring'][d['now']]=0

        # Shift the ring 
        d['ring'] = shift(d['ring'],1)
                
        #get the starting state for the next timestep from the ring
        if return_state and out_type=='array':
            out_state[:,i] = d['ring']
        elif return_state and out_type=='grid':
            out_state[:,point[0],point[1]] = d['ring']
        point_dict[point]=d
        
    return  point_dict, out_flow, out_state, time_dict

def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type coup_conv.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("configFile", type=str, help="Input Configuration File")
    args = parser.parse_args()
    config_file = args.configFile
    return config_file

def process_config_file(config_file):
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)

    # Read Options Section
    options = {}
    try:
        options['verbose'] = Config.getboolean("Options", "verbose")
    except:
        options['verbose'] = False
    # Read Outputs Section
    outputs={}
    try:
        outputs["out_type"] = Config.get("Outputs","out_type")
    except:
        print "WARNING:  outputs[out_type] not found in configuration file, output type will be array"
        outputs["out_type"] = "array"

    try:
        outputs["state"] = Config.get("Outputs","state").split(',')
    except:
        outputs["state"] = False
    try:
        outputs["case_name"] = Config.get("Outputs","case_name")
    except:
        pass
    # Read Paths section
    uh_files = deque(sorted(glob.glob(Config.get("Paths","uh_files"))))

    f = Config.get("Paths","flux_files").split(',')
    if len(f)>1:
        flux_files = deque(f) # this is a list of files to read
    else:
        flux_files = deque(sorted(glob.glob(f[0]))) # this is a string for glob

    try:
        grid_file = Config.get("Paths","grid_file")
    except:
        raise IOerror('REQUIRED FILE NOT PROVIDED:  grid_file')

    try:
        out_path = Config.get("Paths","out_path")
    except:
        raise IOerror('REQUIRED PATH NOT PROVIDED:  out_path')

    try:
        f = Config.get("Paths","initial_state")
        if os.path.exists(f):  initial_state = Config.get("Paths","initial_state")
        else: raise IOerror('Initial State file %s does not exist' %initial_state)
    except:
        initial_state = None
        
    return Config, uh_files, flux_files, grid_file, out_path, initial_state, outputs, options

def shift(l, offset): 
    """
    see F90  cshift with offset=-offset
    """
    offset %= len(l)
    return np.concatenate((l[offset:], l[:offset]))

##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
