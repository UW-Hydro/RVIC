#!/usr/local/bin/python
"""
Note that this isn't really a unit test, rather a system test.

This is designed to test the functionality of the full system as it runs
"""
import argparse
import coup_conv  
import cProfile, pstats
import os, glob
from collections import deque
from subprocess import Popen, PIPE

def main():
    config_file, ProfileFullRun, TestBinaryRestart = process_command_line()
    if TestBinaryRestart:
        Test_Binary_Resart(config_file)
    if ProfileFullRun:
        Profile_Full_Run(config_file)
    return

def Test_Binary_Resart(config_file, N = 24):
    """
    Run 2 parallel convolutions
    1.  Do a convolution without stoping (continuous)
    2.  Do a convolution, stopping every N timesteps
    """
    
    # Do first run 
    Config,uh_files,flux_files,grid_file,out_path,initial_state,outputs,options = coup_conv.process_config_file(config_file)
    out_path1 = os.path.join(out_path,'continuous')
    if not os.path.exists(out_path1):
        os.makedirs(out_path1)
    point_dict,out_dict,area,shape,counts = coup_conv.init(uh_files,flux_files,grid_file,
                                                           initial_state,outputs,options)
    out_name,state_name,restart_name,counts = coup_conv.run(Config,flux_files,out_path1,outputs,
                                                            options,point_dict,out_dict,area,shape,counts)
    coup_conv.final(counts,outputs,out_path)

    print "Done with run 1, starting run 2 now..."
    print "-----------------------------------------------------------"

    # Do second run, this time control the IO more closely
    # and limit each run to a certain length
    Config,uh_files,flux_files,grid_file,out_path,initial_state,outputs,options = coup_conv.process_config_file(config_file)
    out_path2 = os.path.join(out_path,'restart')
    if not os.path.exists(out_path2):
        os.makedirs(out_path2)
    while flux_files:
        new_list = deque()
        for n in xrange(N):
            # pop the next N elements
            try:
                new_list.append(flux_files.popleft())
            except: pass
        Config,uh_files,junk,grid_file,out_path,initial_state,outputs,options = coup_conv.process_config_file(config_file)
        point_dict,out_dict,area,shape,counts = coup_conv.init(uh_files,new_list,grid_file,
                                                               initial_state,outputs,options)
        out_name,state_name,restart_name,counts = coup_conv.run(Config,new_list,out_path2,outputs,
                                                                options,point_dict,out_dict,area,shape,counts)
        config_file = restart_name
    coup_conv.final(counts,outputs,out_path2)

    print "Done with run2"
    
    # now compare the outputs using nccmp
    print "-----------------------------------------------------------"
    print 'Comparing the Outputs of runs 1 and 2'
    print 'Run 1 Directory:  %s' %out_path1
    print 'Run 2 Directory:  %s' %out_path2
    errors = 0
    files1 = glob.glob(os.path.join(out_path1,"*.nc"))
    for f1 in files1:
        d1, f = os.path.split(f1)
        stdout,stdout = Popen(["nccmp", "-d", "-f",os.path.join(out_path1,f), os.path.join(out_path2,f)], stdout=PIPE, stderr=PIPE).communicate()
        if stdout:
            errors += 1
            print "Diff of %s : %s" %f %stdout

    print "Finished Comparing Run1 and Run2...Results..."
    print "Total Comparisons that Returned an Error: %i" %errors
    print "Percentage of Files that Differ: ", round((float(errors)/float(len(files1))*100),1),"%"
    print "-----------------------------------------------------------"
    return

def Profile_Full_Run(config_file):
    """
    Run coup_conv and use cProfile to do profiling
    """
    cProfile.run('coup_conv.main("%s")' %config_file, 'cstats')
    p = pstats.Stats('cstats')
    p.strip_dirs().sort_stats('cumulative').print_stats(25)
    return

def Check_Mass_Ballance(config_file):
    pass
    
def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type ConvolutionSystemTests -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-TBR","--TestBinaryRestart", help="Test Binary Restart",action='store_true')
    parser.add_argument("-PFR","--ProfileFullRun", help="Run the profiler over the full run",action='store_true')
    parser.add_argument("configFile", type=str, help="Input Configuration File")
    args = parser.parse_args()
    TestBinaryRestart = args.TestBinaryRestart
    ProfileFullRun = args.ProfileFullRun
    config_file = args.configFile
    
    return config_file, ProfileFullRun, TestBinaryRestart

if __name__ == '__main__':
    main()
