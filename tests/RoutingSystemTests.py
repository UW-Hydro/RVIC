#!/usr/local/bin/python
"""
Note that this isn't really a unit test, rather a system test.

This is designed to test the functionality of the full system as it runs
"""
import argparse
import rout  
import cProfile, pstats

def main():
    config_file, ProfileFullRun, TestBinaryRestart = process_command_line()
    if TestBinaryRestart:
        Test_Binary_Resart(config_file)
    if ProfileFullRun:
        Profile_Full_Run(config_file)
    return

def Profile_Full_Run(config_file):
    """
    Run coup_conv and use cProfile to do profiling
    """
    cProfile.run('rout.main("%s")' %config_file, 'cstats')
    p = pstats.Stats('cstats')
    p.strip_dirs().sort_stats('cumulative').print_stats(25)
    return

    
def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type ConvolutionSystemTests -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-PFR","--ProfileFullRun", help="Run the profiler over the full run",action='store_true')
    parser.add_argument("configFile", type=str, help="Input Configuration File")
    args = parser.parse_args()
    ProfileFullRun = args.ProfileFullRun
    config_file = args.configFile
    
    return config_file, ProfileFullRun, TestBinaryRestart

if __name__ == '__main__':
    main()
