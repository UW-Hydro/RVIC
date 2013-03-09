
#!/usr/bin/python

from netCDF4 import Dataset
import numpy as np
import time as tm
import argparse


re = 6.37122e6
area = area[:,:]**re

def main():
    grid_nc,area_var,flux_nc,flux_vars,uh_nc,uh_var = process_command_line()
    q = convolve(grid_nc,area_var,flux_nc,flux_vars,uh_nc,uh_var)
    write_netcdf()

def convolve(grid_nc,flux_nc,uh_nc):
    Grid = read_netcdf(grid_nc)
    Flux = read_netcdf(flux_nc)
    UH = read_netcdf(uh_nc)

    convolution():

    return q
    
    
def process_command_line():
    """
    Parse arguments and assign flags for further loading of variables, for
    information on input arguments, type rout.py -h
    """
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("grid_nc", type=str, help="Input grid netcdf including area variable")
    parser.add_argument("flux_nc", type=str, help="Input grid netcdf including runoff/baseflow fluxes")
    parser.add_argument("uh_nc", type=str, help="Input unit hydrograph grid netcdf")
    parser.add_argument("-av","--area_var", type=str, help="variable name for area",default='area')
    parser.add_argument("-fv","--flux_vars", type=str, help="If only one var should be used in the convolution, input here",default=['runoff','baseflow'],nargs='+')
    parser.add_argument("-uv", "--uh_var",type=str, help="variable name for unit hydrograph")
    args = parser.parse_args()

def write_netcdf(out_nc,q,time):
    #Write to NC file
    f = Dataset(out_nc, 'w', format='NETCDF4_CLASSIC')
    
    # set dimensions
    time = f.createDimension('time', None)
    
    # initialize variables
    time = f.createVariable('time','f8',('time',))
    Q = f.createVariable('Q','f8',('time',))
    
    # write attributes for netcdf
    f.description = 'Streamflow.'
    f.history = 'Created ' + tm.ctime(tm.time())
    f.source = 'RASM routing program - convolve'
    
    time.units = 'day'
    time.description = 'Days since 1-1-1'
    
    Q.description = 'Daily runoff at station'
    Q.units = 'm^3/s'
    
    # write data to variables initialized above
    time[:] = time1
    Q[:] = q[0:len(time1)]
    
    f.close()
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
    if verbose: print 'Reading input data vars:', vars, ', from file:',ncFile
    f = Dataset(ncFile,'r')
    if vars==[]: vars = f.variables.keys()
    d={}
    a={}
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
    f.close()
    return d,a

def convolution(flux,uh,area,frac,full = False):   
    """
    Do the convolution by looping over all uh and flux timesteps.
    Return a 1d array that matches the time length of the flux variable.
    If full = True, the length will be the combined length of the flux and uh.  
    """
    q = np.zeros((len(flux)+len(t_uh)))

    #Do the convolution
    for i in xrange(len(flux)):
        for j in xrange(len(t_uh)):
            tstep = i+j
            q[tstep]= q[tstep] + np.sum(flux[i,:,:]*uh[j,:,:]*area[:,:]*frac[:,:])
            
    #Change units to m3/sec
    q /= 1000/86400

    # Truncate if not true
    if not Full:
        q = q[:len(flux)]

    return q

##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
