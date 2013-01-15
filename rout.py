#!/usr/bin/python

# PROGRAM rout, Python-Version, written by Joe Hamman winter 2012/2013
# Routing algorithm developed by D. Lohmann.

## INPUTS
## arg1 = input grids (netcdf)
## arg2 = Pour Point Latitude
## arg3 = Pour Point Longitude
## arg4 = Basin name (string)
## arg5 = velocity
## arg6 = diffusion
## arg7 = UH_BOX (txt file)

## OUTPUTS
## netCDF with gridded UH_S and fraction data

####################################################################################
import numpy as np
import argparse
from netCDF4 import Dataset
import time as Time
import os 
import sys
from datetime import datetime

(mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(sys.argv[0])
print "last modified: %s" % Time.ctime(mtime)
####################################################################################
NODATA = -9999
Cell_flowtime = 4     # Max time for flow to pass through single cell (days)
UH_len =  60          # Max time to outlet (days)
PREC = 1e-30
out_TS = 86400

####################################################################################
# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("basin_x", type=float, help="input longitude of point to route to")
parser.add_argument("basin_y", type=float, help="input latitude of point to route to")
args = parser.parse_args()

#Replace these values with input args
nc_in = '/raid/jhamman/temp_input_files/Wu_routing_inputs.nc'
basin_x = args.basin_x
basin_y = args.basin_y
#basin_x = 127.21875000
#basin_y = 26.66546249
velocity = 1.
diffusion = 2000.
UH_BOX_file = '/raid/jhamman/temp_input_files/UH_RASM_hourly.csv'

##################################################################################
###############################  MAIN PROGRAM ####################################
##################################################################################
def main():
        print 'Start time: ', datetime.now()

        # Load input arrays, store in python list.  (Format - Global['var'])
	Global = read_netcdf(nc_in,('Basin_ID','lon','lat'))

	# Find boudnds of basin and basin_id
	global basin_id
	(basin_id,x_min,x_max,y_min,y_max) = read_global_inputs(Global['lon'],Global['lat'],Global['Basin_ID'],basin_x,basin_y)
	print 'basin_id: ', basin_id
	# Load input arrays, store in python list.  (Format - Basin['var'])
        Basin = clip_netcdf(nc_in,('Basin_ID','Flow_Direction','Flow_Distance','lon','lat'),x_min,x_max,y_min,y_max)
        
        # Load UH_BOX input
        (uh_t,UH_Box) = load_uh(UH_BOX_file)
        
        # Find timestep (timestep is determined from UH_BOX input file)
	global DELTA_T
	global T_Cell
	global T_UH
        DELTA_T = find_TS(uh_t)
        T_Cell = Cell_flowtime*86400/DELTA_T
        T_UH = UH_len*86400/DELTA_T
        
        # Read direction grid and find to_col (to_x) and to_row (to_y)
        (to_x,to_y) = ReadDirection(Basin['Flow_Direction'],Basin['lon'],Basin['lat'],Basin['Basin_ID'])
        
        # Find row/column indicies of lat/lon inputs
        x_ind = find_nearest(Basin['lon'],basin_x)
        y_ind = find_nearest(Basin['lat'],basin_y)
        
        # Find all grid cells upstream of pour point
        (CATCH, fractions) = SearchCatchment(to_x,to_y,x_ind,y_ind,Basin['lon'],Basin['lat'], Basin['Basin_ID'],Basin['Flow_Distance'])
        
        # Make UH for each grid cell upstream of basin pour point (linear routing model - Saint-Venant equation)
        UH = MakeUH(CATCH['x_inds'],CATCH['y_inds'],Basin['Flow_Distance'],Basin['lon'],Basin['lat'])

        # Make UH_RIVER by incrementally moving upstream comining UH functions
        UH_RIVER = MakeGridUHRIVER(UH,to_x,to_y,x_ind,y_ind,CATCH['x_inds'], CATCH['y_inds'], CATCH['count_ds'])
        
        # Make UH_S for each grid cell upstream of basin pour point (combine IRFs for all grid cells in flow path) 
        UH_S = MakeGridUH(UH_RIVER,UH_Box,to_x,to_y,CATCH['x_inds'], CATCH['y_inds'], CATCH['count_ds'])

	# Agregate to output timestep
	UH_out = aggregate(UH_S,DELTA_T,out_TS,CATCH['x_inds'], CATCH['y_inds'])

        #Write to output netcdf
        write_netcdf(Basin['lon'],Basin['lat'], np.linspace(0,out_TS*UH_out.shape[0],UH_out.shape[0],endpoint=False), UH_out, fractions)
       
        print 'routing program finished.'
        print 'End time: ', datetime.now()

##################################################################################
############### Routines #########################################################
##################################################################################

##################################################################################
##  Read netCDF Inputs
##  Read data from input netCDF.  Input netCDF should include the following vars:
##  ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat')
##  Velocity and Diffusion are assumed to be constant in this version
##################################################################################
def read_netcdf(nc_str,vars):
    print 'Reading input data vars:', vars, 'from file:',nc_str
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
    f.close()
    return d

##################################################################################
## Find Basin Dims and ID
##################################################################################
def read_global_inputs(lons,lats,basins,basin_x,basin_y):
	print 'reading global inputs'
	basin_id = basins[find_nearest(lats,basin_y),find_nearest(lons,basin_x)]
	inds = np.nonzero(basins==basin_id)
	x,y = np.meshgrid(np.arange(len(lons)), np.arange(len(lats)))
	x_inds = x[inds]
	y_inds = y[inds]
	x_min = min(x_inds)
	x_max = max(x_inds)
	y_min = min(y_inds)
	y_max = max(y_inds)	
	return (basin_id,x_min,x_max,y_min,y_max)

##################################################################################
##  Clip netCDF Inputs
##  Read data from input netCDF.  Input netCDF should include the following vars:
##  ('Basin_ID','Flow_Direction','Flow_Distance','Land_Mask','lon','lat')
##  Velocity and Diffusion are assumed to be constant in this version
##################################################################################
def clip_netcdf(nc_str,vars,x_min,x_max,y_min,y_max):
    print 'Reading input data vars:', vars
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        temp = f.variables[var][:]
	if np.rank(temp) > 1:
		d[var] = temp[y_min:y_max+1,x_min:x_max+1]
	elif var == 'lon':
		d['lon'] = temp[x_min:x_max+1]
	elif var == 'lat':
		d['lat'] = temp[y_min:y_max+1]
	else:
		print 'error clipping ', var
    f.close()
    return d

##################################################################################
##  Read UH_BOX
##  Read the UH_BOX timeseries.  Save both the UH timeseries and the timestamps
##################################################################################
def load_uh(infile):
    print 'Reading UH_Box from file: ', infile
    #(uh_t,uh) = np.loadtxt(infile, skiprows=1, unpack=True)
    (uh_t,uh) = np.genfromtxt(infile, delimiter=',', skip_header=1, unpack=True)
    uh_t = uh_t.astype(int)
    return (uh_t, uh)

##################################################################################
##  Read Timestep from UH_BOX
##  In this version of the code, the timestep of the UH_BOX must match DELTA_T
##################################################################################
def find_TS(uh_t):
    DELTA_T = uh_t[1]-uh_t[0]
    print 'Timestep = '+ str(DELTA_T) + ' seconds'
    return DELTA_T

##################################################################################
##  Read Direction
##  Determines downstream row/col numbers from flow direction input
##      1 = north      [-1][0]
##      2 = northeast  [-1][+1]
##      3 = east       [0][+1]
##      4 = southeast  [+1][+1]
##      5 = south      [+1][0]
##      6 = southwest  [+1][-1]
##      7 = west       [0][-1]
##      8 = northwest  [-1][-1]
##################################################################################
def ReadDirection(fdr,lons,lats,basin_ids):
    print 'Reading Direction input and finding target row/columns'
    to_x = np.zeros((fdr.shape[0],fdr.shape[1]),dtype=int)
    to_y = np.zeros((fdr.shape[0],fdr.shape[1]),dtype=int)
    for y in xrange(len(lats)):
        for x in xrange(len(lons)):
            if basin_ids[y][x]==basin_id:
                if fdr[y][x] == 1:
                    to_x[y][x] = x
                    to_y[y][x] = y-1
                elif fdr[y][x] == 2:
                    to_x[y][x] = x+1
                    to_y[y][x] = y-1
                elif fdr[y][x] == 3:
                    to_x[y][x] = x+1
                    to_y[y][x] = y
                elif fdr[y][x] == 4:
                    to_x[y][x] = x+1
                    to_y[y][x] = y+1
                elif fdr[y][x] == 5:
                    to_x[y][x] = x
                    to_y[y][x] = y+1
                elif fdr[y][x] == 6:
                    to_x[y][x] = x-1
                    to_y[y][x] = y+1
                elif fdr[y][x] == 7:
                    to_x[y][x] = x-1
                    to_y[y][x] = y
                elif fdr[y][x] == 8:
                    to_x[y][x] = x-1
                    to_y[y][x] = y-1
                else:
                    to_x[y][x] = 0
                    to_y[y][x] = 0
            x = 0
        y = 0
    return (to_x,to_y)

##################################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##################################################################################
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

##################################################################################
## Seach Catchment
## Find all cells upstream of pour point.  Retrun a dictionary with x_inds, yinds,
## and #of cell to downstream pour point.  All are sorted the by the latter.  
##################################################################################
def SearchCatchment(to_x,to_y,x_ind,y_ind,lons,lats,basin_ids,xmask):
    print 'Searching Catchment'
    COUNT = 0
    x_inds = np.zeros(to_x.shape[0]*to_x.shape[1],dtype=int)
    y_inds = np.zeros(to_x.shape[0]*to_x.shape[1],dtype=int)
    count_ds = np.zeros(to_x.shape[0]*to_x.shape[1],dtype=int)
    fractions = np.zeros((to_x.shape[0],to_x.shape[1]))
    for y in xrange(len(lats)):
        for x in xrange(len(lons)):
            if basin_ids[y][x]==basin_id:
                yy = y
                xx = x
                op=0
                cells = 0
                while op == 0:
                    if (yy==y_ind and xx==x_ind):
                        op = 1
                        x_inds[COUNT] = x
                        y_inds[COUNT] = y
                        count_ds[COUNT] = cells
                        COUNT = COUNT+1
                        fractions[y][x] = 1.
                    elif (to_x[yy][xx] > -1 and to_y[yy][xx] > -1):
                        yyy = to_y[yy][xx]
                        xxx = to_x[yy][xx]
                        yy=yyy
                        xx=xxx
                        cells = cells + 1
                        if (xx>len(lons)-1 or xx<0 or yy>len(lats)-1 or yy<0):
                            op =-1
                    else:
                        op = -1
            x = 0
        y = 0
    print "Upstream grid cells from present station: ", COUNT
    # Save to sorted dictionary
    CATCH = {}
    ii = np.argsort(count_ds[:COUNT])
    x_inds = x_inds[:COUNT]
    y_inds = y_inds[:COUNT]
    CATCH['x_inds'] = x_inds[ii]
    CATCH['y_inds'] = y_inds[ii]
    CATCH['count_ds'] = count_ds[ii]
    return (CATCH, fractions)

##################################################################################
##  MakeUH
##  Calculate impulse response function for grid cells using equation (15) from
##  Lohmann, et al. (1996) Tellus article.  Return 3d UH grid.  
###################################################################################
def MakeUH(x_inds,y_inds,xmask,lons,lats):
	print 'Making UH for each cell'
	UH = np.zeros((T_Cell, xmask.shape[0],xmask.shape[1]))
	for j in xrange(len(x_inds)):
		#print 'Making UH: ', j, 'of ', len(x_inds)
		x = x_inds[j]
		y = y_inds[j]
		time = 0.
		green = np.zeros(T_Cell)
		for t in xrange(T_Cell):
			time = time+DELTA_T
			exponent = -1*np.power(velocity*time-xmask[y][x],2)/(4*diffusion*time)
			if exponent > -69.:
				green[t] = xmask[y][x]/(2*time*np.sqrt(np.pi*time*diffusion))*np.exp(exponent)
			else:
				green[t] = 0.
		sum = np.sum(green)
		if sum>0.:
			UH[:,y,x] = green[:]/sum
	#print UH[:,y_inds[0],x_inds[0]]
	return UH

##################################################################################
##  Make RIVER UH
##  Calculate impulse response function for river routing
##  Steps upstream combining unit hydrographs
###################################################################################        
def MakeGridUHRIVER(UH,to_x,to_y,x_ind,y_ind,x_inds,y_inds,count_ds):
	print "Making UH_RIVER grid.... It takes a while...\n"
	UH_RIVER = np.zeros((T_UH,UH.shape[1],UH.shape[2]))
	for j in xrange(len(x_inds)):
		d = count_ds[j]
		if d>0:
			x = x_inds[j]
			y = y_inds[j]
			yy = to_y[y,x]
			xx = to_x[y,x]
			#print j, 'of', len(x_inds)
			#print d, x, y, xx, yy
			IRF_temp = np.zeros(T_UH)
			for t in np.nonzero(UH_RIVER[:,yy,xx]>PREC)[0]: # not sure this is ok but its really fast
			#for t in np.nonzero(UH_RIVER[:,yy,xx])[0]:
				for l in xrange(T_Cell):
					IRF_temp[t+l] = IRF_temp[t+l] + UH[l,y,x]*UH_RIVER[t,yy,xx]
			sum = np.sum(IRF_temp[:T_UH])
			if sum>0:
				UH_RIVER[:,y,x] = IRF_temp[:]/sum
		elif d==0:
			UH_RIVER[:T_Cell,y_ind,x_ind] = UH[:,y_ind,x_ind]
	return UH_RIVER

##################################################################################
## Make Grid UH
## Combines the UH_BOX with downstream cell UH_River IRF.
## Cell [0] is given the UH_Box without river routing
##################################################################################
def MakeGridUH(UH_RIVER,UH_BOX,to_x,to_y,x_inds, y_inds,count_ds):
	print "Making UH_S grid" 
	UH_S = np.zeros((T_UH,UH_RIVER.shape[1],UH_RIVER.shape[2]))+NODATA
	for j in xrange(len(x_inds)):
		#print 'Making UH_S ', j, 'of ', len(x_inds)
		x = x_inds[j]
		y = y_inds[j]
		d = count_ds[j]
		IRF_temp = np.zeros(T_UH)
		if d > 0:
			yy = to_y[y,x]
			xx = to_x[y,x]
			for t in np.nonzero(UH_RIVER[:,yy,xx]>PREC)[0]: 
				for l in xrange(len(UH_BOX)):
					IRF_temp[t+l] = IRF_temp[t+l] + UH_BOX[l]*UH_RIVER[t,yy,xx]
			sum = np.sum(IRF_temp)
			if sum>0:
				UH_S[:,y,x] = IRF_temp[:]/sum
		elif d==0:		
			UH_S[:,y,x] = IRF_temp
			UH_S[:len(UH_BOX),y,x] = UH_BOX[:]
	return UH_S

##################################################################################
## Aggregate to larger timestep 
##################################################################################
def aggregate(UH_S,in_TS,out_TS,x_inds,y_inds):
	if out_TS<in_TS:
		print 'error in aggregation step'
	elif out_TS==in_TS:
		print 'no need to aggregate, skipping this step (UH_out = UH_S'
		UH_out = UH_S
	elif np.remainder(out_TS,in_TS)==0:
		print 'aggregating to ', out_TS, 'from ', in_TS, 'seconds'
		fac = int(out_TS/in_TS)
		T_UH_out = int(T_UH/fac)
		UH_out = np.zeros((T_UH_out,UH_S.shape[1],UH_S.shape[2]))+NODATA
		for j in xrange(len(x_inds)):
			x = x_inds[j]
			y = y_inds[j]
			for t in xrange(T_UH_out):
				UH_out[t,y,x] = np.sum(UH_S[t*fac:t*fac+fac,y,x])
	else:
		print 'error in aggregation step'
	return UH_out 

##################################################################################
##  Write output to netCDF
##  Writes out a netCDF3-64BIT data file containing the UH_S and fractions
##################################################################################
def write_netcdf(lons, lats, times, UH_S,fractions):
	string = 'BasinUH_'+('%.8f' % basin_x)+'_'+('%.8f' % basin_y)+'.nc'
	f = Dataset(string,'w', format='NETCDF3_64BIT')
	
	# set dimensions
	time = f.createDimension('time', None)
	lon = f.createDimension('lon', (len(lons)))
	lat = f.createDimension('lat', (len(lats)))
	
	# initialize variables
	time = f.createVariable('time','f8',('time',))
	lon = f.createVariable('lon','f8',('lon',))
	lat = f.createVariable('lat','f8',('lat',))
	fraction = f.createVariable('fraction','f8',('lat','lon',),fill_value=NODATA)
	UHS = f.createVariable('unit_hydrograph','f8',('time','lat','lon',),fill_value=NODATA)
	
	# write attributes for netcdf
	f.description = 'UH_S for point '+('%.8f' % basin_x)+', '+('%.8f' % basin_y)
	f.history = 'Created ' + Time.ctime(Time.time())
	f.source = sys.argv[0] # prints the name of script used
	f.velocity = velocity
	f.diffusion = diffusion
	f.outlet_lon = basin_x
	f.outlet_lat = basin_y
	f.global_basin_id = basin_id
	
	lat.long_name = 'latitude coordinate'
	lat.standard_name = 'latitude'
	lat.units = 'degrees_north'
	
	lon.long_name = 'longitude coordinate'
	lon.standard_name = 'longitude'
	lon.units = 'degrees_east'

	time.units = 'seconds'
	time.description = 'Seconds since initial impulse'
	
	UHS.units = 'unitless'
	UHS.description = 'unit hydrograph for each grid cell with respect to downstream grid location - ' +('%.8f' % basin_x)+','+('%.8f' % basin_y)

	fraction.description = 'fraction of grid cell contributing to guage location'
	
	# write data to variables initialized above
	time[:]= times
	lon[:] = lons
	lat[:] = lats
	UHS[:,:,:] = UH_S
	fraction[:,:]= fractions
	f.close()
	print 'wrote UH_S netCDF for ' + ('%.8f' % basin_x)+'_'+('%.8f' % basin_y)

##################################################################################
# Run Program
##################################################################################
main()
