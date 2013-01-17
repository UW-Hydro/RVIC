#!/usr/bin/python

# Written by Joe Hamman, University of Washington, January 2012
"""
This script launches batch streamflow routing jobs on the Hydra Cluster (SGE).
It does this by:
1.  Loading input points to (infile) route to and input grids (nc_in)
2.  Splits global netcdf to basin size (uses ncks)
3.  Checks to see how many jobs are already running, if there is room in the queue it contintues,
    otherwise it waits until there is room. 
4.  Writes small qsub shell script
5.  Executes qsub script (submits job)
6.  Inprogress...Cleanup
"""
######################################
import numpy as np
import os
import getpass
from datetime import datetime
from netCDF4 import Dataset
import py_compile
import subprocess
######################################
# Inpaths
temp_path = '/raid2/jhamman/route_rasm/inputs/run1_RASM/temp/'
out_path = '/raid2/jhamman/route_rasm/outputs/run1_RASM/'
infile = '/raid2/jhamman/route_rasm/inputs/run1_RASM/Arctic_stations.csv'
UH_BOX = '/raid2/jhamman/route_rasm/inputs/run1_RASM/UH_RASM_hourly.csv'
nc_in = '/raid2/jhamman/route_rasm/inputs/run1_RASM/Wu_routing_inputs.nc'
scriptname = '/home/jhamman/rout_rasm/rout.py'
stdout = '/raid2/jhamman/route_rasm/outputs/run1_RASM/stdout/'
qsub_template = '/home/jhamman/rout_rasm/rout.scr'

#compile the script
#scriptname = py_compile.compile(scriptname)

# Input Parameters
velocity = 1
diffusion = 2000

#Available Queues and nodes
Avail_Queues = {'forecast.q':60,'default.q':20}
#Avail_Queues = {'forecast.q':60, 'ben.q':48, 'default.q':20}

######################################
######################################
def main():
    username = getpass.getuser()

    # Load input files
    (lon,lat,basins) = load_points(infile)
    Global = read_netcdf(nc_in,('lon','lat'))

    # Run routing program for each basin
    Q_ind = 0
    for i, basin in enumerate(basins):
        basin_string = 'Basin_'+str("%06d" % basin)
        flag = 0
        while flag==0:
            QUEUE = Avail_Queues.keys()[Q_ind]
            cmd = 'qstat -u ' + username + ' -q ' + QUEUE + ' | grep ' + scriptname + ' | wc -l'
            print cmd
            used = os.system(cmd)
            print 'used = ',used
            if used > Avail_Queues[QUEUE]:
                time.sleep(20)
            else:
                flag = 1
            if Q_ind==len(Avail_Queues)-1:
                Q_ind = 0
            else:
                Q_ind = Q_ind+1
        # Write qsub job string
        qsub_cmd = qsub_string(lat[i],lon[i],basins[i],QUEUE)
        print qsub_cmd
        # Run qsub script
        os.system(qsub_cmd)
        print 'just started '+basin_string+' ('+str(i+1)+' of '+str(len(basins))+') at time: '+ str(datetime.now())
    print 'Finished submitting jobs'
##     print 'Will start cleaning up when all jobs under userid: ', username, ' are finished'
##     # Remove all files from temporary folder (subset netcdfs and qsub shell scripts)
##     running_flag = 0
##     while running_flag == 0:
##         cmd = 'qstat -u ' + username + ' | grep ' + scriptname + ' | wc -l'
##         running = os.system(cmd)
##         if running > 0:
##             print 'waiting on ' + str(running) + ' jobs'
##             time.sleep(600)
##         else:
##             running_flag = 1
##     print 'Cleaning up'
##     cleanup()
    print 'Done with everything at time: ' + str(datetime.now())
    
#################################################################################
## Read Pour Points
##   format: (lon,lat,basins,area,lon_min,lat_min, lon_max,lat_max)
#################################################################################
def load_points(infile):
    print 'Reading Pour Points from file: ', infile
    (lon,lat,basins) = np.genfromtxt(infile,delimiter=',',unpack=True)
    #Sort by area (start with larger basins)
    #ii = np.argsort(area)[::-1]
    #lon = lon[ii]
    #lat = lat[ii]
    #basins = basins[ii].astype(int)
    return (lon,lat,basins)

##################################################################################
##  Read netCDF Inputs [just read lons and lats]
##  Read data from input netCDF.  Input netCDF should include the following vars:
##  ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat')
##  Velocity and Diffusion are assumed to be constant in this version
##################################################################################
def read_netcdf(nc_str,vars):
    print 'Reading input data vars:', vars
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
    f.close()
    return d

##################################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##################################################################################
## def find_nearest(array,value):
##     idx = (np.abs(array-value)).argmin()
##     return idx

##################################################################################
##  Extract a subset of a netcdf using the ncks command
##  Requires in indecies of the subset (uses find_nearest function)
##################################################################################
## def split_netCDF(y_max_ind, y_min_ind, x_min_ind, x_max_ind,basin_string):
##     cmd = 'ncks -d lat,' + str(y_max_ind) +','+str(y_min_ind)+' -d lon,'+str(x_min_ind)+','+str(x_max_ind)+' '+nc_in+ ' '+ temp_path + basin_string+'.nc'
##     os.system(cmd)

##################################################################################
##  Writes the qsub string
##################################################################################
## def write_qsub(lat,lon,basin,basin_string,QUEUE):
##     arg1 = temp_path+basin_string
##     arg2 = str(lat)
##     arg3 = str(lon)
##     arg4 = str(basin)
##     arg5 = str(velocity)
##     arg6 = str(diffusion)
##     arg7 = 'path/to/UH_BOX'
##     args = arg1+' '+arg2+' '+arg3+' '+arg4+' '+arg5+' '+arg6+' '+arg7
##     QSUB_SCRIPT = QSUB_SCRIPT.replace("_scriptname_", scriptname)
##     QSUB_SCRIPT = QSUB_SCRIPT.replace("_WorkingDir_", out_path)
##     QSUB_SCRIPT = QSUB_SCRIPT.replace("_args_", args)
##     QSUB_SCRIPT = QSUB_SCRIPT.replace("_queue_", QUEUE)
##     qsfname = temp_path+'rn.'+basin_string+'.sh'
##     qsub = open(qsfname,'w')
##     qsub.write(QSUB_SCRIPT)
##     qsub.close()
##     QSUB_SCRIPT = QSUB_TEMPLATE # Resets QSUB_TEMPLATE
##     return qsfname

def qsub_string(lat,lon,basin_id,QUEUE):
    arg1 = nc_in
    arg2 = UH_BOX
    arg3 = str(lon)
    arg4 = str(lat)
    arg5 = str(velocity)
    arg6 = str(diffusion)
    args = arg1+' '+arg2+' '+arg3+' '+arg4+' '+arg5+' '+arg6
    string = 'qsub -wd '+out_path +' -o ' +stdout +' -q '+QUEUE+' -N rt.'+str("%06d" % basin_id) + ' ' + qsub_template+' '+scriptname+' '+args
    return string

##################################################################################
##  Removes all files in temporary directory
##################################################################################
def cleanup ():
    cmd = "rm "+temp_path+"*"
    #os.system(cmd)
    print cmd

main()
