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
import commands
import time
######################################
# Inpaths
temp_path = '/raid2/jhamman/route_rasm/inputs/run2_RASM/temp/'
out_path = '/raid2/jhamman/route_rasm/outputs/run2_RASM/'
infile = '/raid2/jhamman/route_rasm/inputs/run2_RASM/POP_rout2points2.csv'
UH_BOX = '/raid2/jhamman/route_rasm/inputs/run2_RASM/UH_RASM_hourly.csv'
nc_in = '/raid2/jhamman/route_rasm/inputs/run2_RASM/Wu_routing_inputs.nc'
scriptname = '/home/jhamman/Streamflow_Routing/routing/rout.py'
configfile = '/home/jhamman/Streamflow_Routing/routing/rout.cfg'
stdout = '/raid2/jhamman/route_rasm/outputs/run2_RASM/stdout/'
qsub_template = '/home/jhamman/Streamflow_Routing/tools/rout.scr'

#compile the script
#scriptname = py_compile.compile(scriptname)

## # Input Parameters
## velocity = 1
## diffusion = 2000

#Available Queues and nodes
Avail_Queues = {'forecast.q':60,'default.q':20,'ben.q':20}
#Avail_Queues = {'forecast.q':60, 'ben.q':48, 'default.q':20}

######################################
######################################
def main():
    username = getpass.getuser()

    # Load input files
    (lon,lat,basins) = load_points(infile)

    # Run routing program for each basin
    Q_ind = 0
    for i, basin in enumerate(basins):
        basin_string = 'Basin_'+str("%06d" % basin)
        flag = 0
        while flag==0:
            QUEUE = Avail_Queues.keys()[Q_ind]
            cmd = 'qstat -u ' + username + ' -q ' + QUEUE + ' | grep ' + username + ' | wc -l'
            #print cmd
            used = int(commands.getoutput(cmd))
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
        #print qsub_cmd
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
##   format: (lon,lat,basins)
#################################################################################
def load_points(infile):
    print 'Reading Pour Points from file: ', infile
    (lon,lat,basins) = np.genfromtxt(infile,delimiter=',',unpack=True)
    return (lon,lat,basins)

def qsub_string(lat,lon,basin_id,QUEUE):
    arg1 = str(lon)
    arg2 = str(lat)
    args = '-C '+configfile+' -lon '+str(lon)+' -lat '+str(lat)
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
