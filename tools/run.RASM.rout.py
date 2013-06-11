#!/usr/bin/python

# Written by Joe Hamman, University of Washington, January 2012
"""
This script launches batch streamflow routing jobs on the Hydra Cluster (SGE).
It does this by:
1.  Loading input points to (infile) route to
2.  Checks to see how many jobs are already running, if there is room in the queue it contintues,otherwise it waits until there is room. 
5.  Executes qsub script (submits job)
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

out_path = '/raid2/jhamman/route_rasm/outputs/run4_RASM/uh_files/'
infile = '/raid2/jhamman/route_rasm/inputs/run4_RASM/POP_rout2points2.csv'
scriptname = '/home/jhamman/Streamflow_Routing/routing/rout.py'
configfile = '/home/jhamman/Streamflow_Routing/routing/rout.cfg'
stdout = '/raid2/jhamman/route_rasm/outputs/run4_RASM/stdout/'
qsub_template = '/home/jhamman/Streamflow_Routing/tools/rout.scr'

#Available Queues and nodes
Avail_Queues = {'default.q':100}

######################################
######################################
def main():
    username = getpass.getuser()

    # Load input files
    (lons,lats) = load_points(infile)

    # Run routing program for each basin
    Q_ind = 0
    for (lon,lat), basin in enumerate(zip(lons,lats)):
        basin_string = 'Basin_'+str("%06d" % basin)
        flag = 0
        while flag==0:
            QUEUE = Avail_Queues.keys()[Q_ind]
            cmd = 'qstat -u ' + username + ' -q ' + QUEUE + ' | grep ' + username + ' | wc -l'
            #print cmd
            used = int(commands.getoutput(cmd))
            print 'used = ',used
            if used > Avail_Queues[QUEUE]:
                time.sleep(5)
            else:
                flag = 1
            if Q_ind==len(Avail_Queues)-1:
                Q_ind = 0
            else:
                Q_ind = Q_ind+1
        # Write qsub job string
        qsub_cmd = qsub_string(lat,lon,basin,QUEUE)
        #print qsub_cmd
        # Run qsub script
        os.system(qsub_cmd)
        print 'just started '+basin_string+' ('+str(basin+1)+' of '+str(len(basins))+') at time: '+ str(datetime.now())
    print 'Finished submitting jobs'
    print 'Done with everything at time: ' + str(datetime.now())
    
#################################################################################
## Read Pour Points
##   format: (lon,lat,basins)
#################################################################################
def load_points(infile):
    print 'Reading Pour Points from file: ', infile
    (lon,lat) = np.genfromtxt(infile,delimiter=',',unpack=True)
    return (lon,lat)

def qsub_string(lat,lon,basin_id,QUEUE):
    args = '-C='+configfile+' -lon='+str(lon)+' -lat='+str(lat)
    string = 'qsub -wd '+out_path +' -o ' +stdout +' -q '+QUEUE+' -N rt.'+str("%06d" % basin_id) + ' ' + qsub_template+' '+scriptname+' '+args
    return string

main()
