#!/usr/local/bin/python
"""
Aggregate high resolution Unit Hydrograph NetCDF Grids

Joe Hamman March, 2013
"""
import os
import sys
import numpy as np
import numpy.ma as ma
from collections import defaultdict
from netCDF4 import Dataset
from cdo import *
import time as tm
import argparse
cdo = Cdo()
from scipy.spatial import cKDTree

##################################################################################

def main():

    Rvars,paths,options = process_command_line()

    Cvars = ('yc','xc')
    Inputs=read_netcdf(paths['gridFile'],Cvars,options['verbose'])

    #Find Destination grid cells 
    if options['verbose']:
        print 'Finding addresses now...'

    #1.  Make a list of source points (regular gid lat/lons)
    file_list = os.listdir(paths['srcDir'])
    lons, lats = [], []
    for file in file_list:
        f,suffix = os.path.splitext(file)
        prefix,lon,lat = f.split('_')
        lons.append(float(lon))
        lats.append(float(lat))

    if (min(lons)<0 and Inputs['xc'].min()>=0):
        posinds = np.nonzero(Inputs['xc']>180)
        Inputs['xc'][posinds] -= 360
        if options['verbose']:
            print 'adjusted xc lon minimum'

    if options['testAgg']:
        # limit the inputs arrays to a single point
        # so that all points are mapped to just one location
        combined = np.dstack(([Inputs['xc'][0,0].ravel(),Inputs['yc'][0,0].ravel()]))[0]
    else:
        combined = np.dstack(([Inputs['xc'].ravel(),Inputs['yc'].ravel()]))[0]

    points=list(np.vstack((np.array(lons),np.array(lats))).transpose())

    mytree = cKDTree(combined)
    dist,indexes =mytree.query(points,k=1)
    indexes = np.array(indexes)

    d = defaultdict(list)
    for i,ind in enumerate(indexes):
        key = (combined[ind][0],combined[ind][1])
        d[key].append(file_list[i])

    if options['dryrun']:
        print 'Starting dry run now'
        lcount = 0
        kcount = 0
        file_list = os.listdir(paths['srcDir'])
        num = len(file_list)
        for i,key in enumerate(d):
            kcount += 1
            for j,fname in enumerate(d[key]):
                if os.path.isfile(os.path.join(paths['srcDir'],fname)):
                    lcount += 1
                    file_list.remove(fname)
        print 'NUMBER OF FILES IN SOURCE DIRECTORY:',num
        print 'NUMBER OF POINTS TO AGGREGATE TO:', kcount
        print 'NUMBER OF SOURCE FILES FOUND:',lcount
        print 'NUMBER OF FILES STILL IN FILE LIST:',len(file_list)
        print 'EFFECIENCY OF:', 100.*lcount/num, '%'
        print 'FILES LEFT IN FILE LIST:\n',file_list
        return
        

    # Aggregate all basins in each key
    count = 0
    if options['verbose']:
        print 'Aggregating now...'
        if options['testAgg']:
            print 'Cool, youre doing the aggTest...'
    for i,key in enumerate(d):
        if options['verbose']:
            print '  Started',key,i+1,'of',len(d)
        flag  = 0
        Flist = []
        for j,inFile in enumerate(d[key]):
            if options['verbose']:
                sys.stdout.write('\r')
                sys.stdout.write('Last opened file '+str(j+1)+" of "+str(len(d[key]))+" - "+inFile)
                sys.stdout.flush()
            Flist.append(inFile)
            f = Dataset(os.path.join(paths['srcDir'],inFile),'r')
            if flag == 0:
                aggData={}
                aggFile = filename(options['outPrefix'],key)
                for var in Rvars:
                    aggData[var] = f.variables[var][:]
                if count == 0:
                    velocity = f.velocity
                    diffusion = f.diffusion
            else:
                inData={}
                for var in Rvars:
                    inData[var] = f.variables[var][:] 
                aggData = agg(inData,aggData,options['resolution'],
                              options['fill_value'],options['verbose'])
            f.close()
            flag = 1
            count += 1
        # Add pad to final file
        if (len(Flist)>0 and flag==1):
             aggData = agg([],aggData,options['resolution'],options['verbose'],
                           options['fill_value'],pad=options['pad'])
             # Write out to netCDF
             write_netcdf(os.path.join(paths['aggDir'],aggFile),aggData['lon'],
                          aggData['lat'],aggData['time'],aggData['unit_hydrograph'],
                          aggData['fraction'],key,Flist,velocity,diffusion,
                          options['fill_value'],options['verbose'])
        if (flag == 1 and options['remap']):
            remap_file(paths['gridFile'],aggFile,paths['aggDir'],paths['remapDir'],
                       options['verbose'])
            if options['verbose']:
                print 'Finished',key, i+1,'of',len(d), 'and placed', count, 'files\n'
            if options['clean']:
                #clean agg directory
                clean(paths['aggDir'])
    print 'done with everything'
    return

##################################################################################
def read_netcdf(nc_str,vars,verbose):
    #if verbose == True:
    #    print 'Reading input data vars:', vars, 'from file:',nc_str
    f = Dataset(nc_str,'r')
    d={}
    for var in vars:
        d[var] = f.variables[var][:]
    f.close()
    return d

##################################################################################
def make_degrees(vars,Inputs,verbose):  
    d={}
    for var in vars:
        d[var]=np.rad2deg(Inputs[var])
    if verbose == True:
        print 'Converted', vars, 'to degrees'
    return d

##################################################################################
def filename(prefix,tup):
    lon = str(round(tup[0],2))
    lat = str(round(tup[1],2))
    string = prefix+lon+'_'+lat+'.nc'
    return string

##################################################################################
def remap_file(gridFile,aggFile,aggDir,remapDir,verbose):
    remapFile = os.path.join(remapDir,aggFile)
    aggFile = os.path.join(aggDir,aggFile)
    cdo.remapcon(gridFile, input = aggFile, output = remapFile, options = '-f nc4 -z zip')
    # cdo.remapcon(gridFile, input = aggFile, output = remapFile, options = '-f nc4')
    if verbose:
        print '\nremapped to out file:', remapFile
    return

##################################################################################
##  Write output to netCDF
##  Writes out a netCD4 data file containing the UH_S and fractions
##################################################################################
def write_netcdf(file,lons,lats,times,hydrographs,fractions,loc,Flist,velocity,diffusion,fill_value,verbose):
    """
    Write output to netCDF.  Writes out a netCDF4 data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    f = Dataset(file,'w', format='NETCDF4')

    # set dimensions
    time = f.createDimension('time', None)
    lon = f.createDimension('lon', (len(lons)))
    lat = f.createDimension('lat', (len(lats)))

    # initialize variables
    time = f.createVariable('time','f8',('time',))
    lon = f.createVariable('lon','f8',('lon',))
    lat = f.createVariable('lat','f8',('lat',))
    fraction = f.createVariable('fraction','f8',('lat','lon',),fill_value=fill_value)
    UHS = f.createVariable('unit_hydrograph','f8',('time','lat','lon',),fill_value=fill_value)

    # write attributes for netcdf
    f.description = 'Aggregated UH_S and Fraction Vars'
    f.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    f.history += ' '.join(sys.argv) + '\n'
    f.source = sys.argv[0] # prints the name of script used
    f.velocity = velocity
    f.diffusion = diffusion
    f.outlet_lon = loc[0]
    f.outlet_lat = loc[1]
    f.includes = ', '.join(Flist)

    lat.long_name = 'latitude coordinate'
    lat.standard_name = 'latitude'
    lat.units = 'degrees_north'

    lon.long_name = 'longitude coordinate'
    lon.standard_name = 'longitude'
    lon.units = 'degrees_east'

    time.units = 'seconds since 0000-1-1 0:0:0'
    time.calendar = 'noleap'
    time.longname = 'time'
    time.type_prefered = 'float'
    time.description = 'Seconds since initial impulse'

    UHS.units = 'unitless'
    UHS.description = 'unit hydrograph for each grid cell with respect to downstream grid location'
    
    fraction.units = 'unitless'
    fraction.description = 'fraction of grid cell contributing to guage location'

    # write data to variables initialized above
    time[:]= times
    lon[:] = lons
    lat[:] = lats
    UHS[:,:,:] = hydrographs
    fraction[:,:]= fractions
    f.close()
##################################################################################
def agg(inData,aggData,resolution,verbose,fill_value,pad=0):
    """
    Add the two data sets together and return the combined arrays.
    The two data sets must include the coordinate variables lon,lat, and time
    """
    if (inData and aggData):
        # find range of vectors
        Range = [np.minimum(inData['time'].min(),aggData['time'].min()),
                 np.maximum(inData['time'].max(),aggData['time'].max()),
                 np.minimum(inData['lat'].min(),aggData['lat'].min()),
                 np.maximum(inData['lat'].max(),aggData['lat'].max()),
                 np.minimum(inData['lon'].min(),aggData['lon'].min()),
                 np.maximum(inData['lon'].max(),aggData['lon'].max())]
        tStep = inData['time'][1] - inData['time'][0]
        try: yres = np.absolute(aggData['lat'][1] - aggData['lat'][0])
        except:
            try: yres = np.absolute(inData['lat'][1] - inData['lat'][0])
            except: yres = resolution
        try: xres = np.absolute(aggData['lon'][1] - aggData['lon'][0])
        except:
            try: xres = np.absolute(inData['lon'][1] - inData['lon'][0])
            except: xres = resolution
    elif inData:
        Range = [inData['time'].min(),inData['time'].max(),
                 inData['lat'].min(),inData['lat'].max(),
                 inData['lon'].min(),inData['lon'].max()]
        tStep = inData['time'][1] - inData['time'][0]
        try: yres = np.absolute(inData['lat'][1] - inData['lat'][0])
        except: yres = resolution
        try: xres = np.absolute(inData['lon'][1] - inData['lon'][0])
        except: xres = resolution
    elif aggData:
        Range = [aggData['time'].min(),aggData['time'].max(),
                 aggData['lat'].min(),aggData['lat'].max(),
                 aggData['lon'].min(),aggData['lon'].max()]
        tStep = aggData['time'][1] - aggData['time'][0]
        try: yres = np.absolute(aggData['lat'][1] - aggData['lat'][0])
        except: yres = resolution
        try: xres = np.absolute(aggData['lon'][1] - aggData['lon'][0])
        except: xres = resolution
    else:
        raise IOError('no inputs to agg function')
    # make output arrays for lons/lats and initialize fractions/hydrographs
    # pad output arrays so there is a space =pad around inputs
    times = np.arange(Range[0],Range[1]+tStep,tStep)
    lats = np.arange(Range[2]-yres*(pad),Range[3]+yres*(1+pad),yres)[::-1]
    lons = np.arange(Range[4]-xres*(pad),Range[5]+xres*(1+pad),xres)
    fractions = np.zeros((lats.shape[0],lons.shape[0]))
    hydrographs = np.zeros((times.shape[0],lats.shape[0],lons.shape[0]))
    
    # find target index locations of all corners for both datasets
    if inData:
        In = [find_nearest(times,np.min(inData['time'])), find_nearest(times,np.max(inData['time']))+1,
              find_nearest(lats,np.max(inData['lat'])), find_nearest(lats,np.min(inData['lat']))+1,
              find_nearest(lons,np.min(inData['lon'])), find_nearest(lons,np.max(inData['lon']))+1]
    if aggData:
        Ex = [find_nearest(times,np.min(aggData['time'])), find_nearest(times,np.max(aggData['time']))+1,
              find_nearest(lats,np.max(aggData['lat'])), find_nearest(lats,np.min(aggData['lat']))+1,
              find_nearest(lons,np.min(aggData['lon'])), find_nearest(lons,np.max(aggData['lon']))+1]

    # Make sure all values in the unit hydrograph are zero (no mask)
    if inData:
        inData['unit_hydrograph'][inData['unit_hydrograph']<0] = 0.0
        try:
            inData['unit_hydrograph'] = inData['unit_hydrograph'].filled(fill_value=0)
        except:
            pass
    if aggData:
        aggData['unit_hydrograph'][aggData['unit_hydrograph']<0] = 0.0
        try:
            aggData['unit_hydrograph']=aggData['unit_hydrograph'].filled(fill_value=0)
        except:
            pass
        
    # Place data
    # First the fractions
    if inData:
        fractions[In[2]:In[3],In[4]:In[5]] += inData['fraction']
    if aggData:
        fractions[Ex[2]:Ex[3],Ex[4]:Ex[5]] += aggData['fraction']

    # If there is a chance that there is overlap between basins, this method will need to be used.
    # Otherwise, the simplier method below should work fine  
    # # Then the hydrographs 
    # if inData:
    #     pvals = np.nonzero(fractions[In[2]:In[3],In[4]:In[5]]>0)
    #     hydrographs[In[0]:In[1],In[2]:In[3],In[4]:In[5]][:,pvals[0],pvals[1]] += inData['unit_hydrograph'][:,pvals[0],pvals[1]]*(inData['fraction'][pvals]/fractions[In[2]:In[3],In[4]:In[5]][pvals])
    # if aggData:
    #     pvals = np.nonzero(fractions[Ex[2]:Ex[3],Ex[4]:Ex[5]]>0)
    #     hydrographs[Ex[0]:Ex[1],Ex[2]:Ex[3],Ex[4]:Ex[5]][:,pvals[0],pvals[1]] += aggData['unit_hydrograph'][:,pvals[0],pvals[1]]*(aggData['fraction'][pvals]/fractions[Ex[2]:Ex[3],Ex[4]:Ex[5]][pvals])
    if inData:
        hydrographs[In[0]:In[1],In[2]:In[3],In[4]:In[5]] += inData['unit_hydrograph']
    if aggData:
        hydrographs[Ex[0]:Ex[1],Ex[2]:Ex[3],Ex[4]:Ex[5]] += aggData['unit_hydrograph']
    
    # Mask the hydrographs and make sure they sum to 1 at each grid cell
    if (inData == [] or aggData == []):
        ym,xm = np.nonzero((fractions<=0)*(hydrographs.sum(axis=0)<=0))
        fractions[ym,xm] = 0
        hydrographs[:,ym,xm] = fill_value
        
        # Normalize the hydrographs (each cell should sum to 1)
        yv,xv = np.nonzero(fractions>0)
        # print '\n'
        # print hydrographs[:,yv,xv].sum(axis=0)
        hydrographs[:,yv,xv] /= hydrographs[:,yv,xv].sum(axis=0)
        # print 'just normalized the uh grid'
        # print hydrographs[:,yv,xv].sum(axis=0)

    # Put all the data into aggData variable and return to main
    
    aggData['lon'] = lons
    aggData['lat'] = lats
    aggData['fraction'] = fractions
    aggData['unit_hydrograph'] = hydrographs
    aggData['time'] = times

    return aggData
    
##################################################################################
##  Find Indicies
##  Given an input lat or lon, the function returns the nearest index location
##################################################################################
def find_nearest(array,value):
    """ Find the index location in (array) with value nearest to (value)"""
    idx = (np.abs(array-value)).argmin()
    return idx

##################################################################################
def process_command_line():
    """
    Process command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("srcDir", type=str, help="Directory containing Unit Hydrograph grids to be aggregated")
    parser.add_argument("--gridFile", type=str, help="Input netCDF target grid")
    parser.add_argument("--remapDir", type=str, help="Directory containing Output Unit Hydrograph grids")
    parser.add_argument("--aggDir", type=str, help="Directory where to store aggregated files (before remap)")
    parser.add_argument("--inPrefix", type=str, help="Input Unit Hydrograph File Prefix",default='UH_')
    parser.add_argument("--outPrefix", type=str, help="Output Unit Hydrograph File Prefix", default="Agg_UH_")
    parser.add_argument("--time", type=str, help="Input Unit Hydrograph time variable name",default='time')
    parser.add_argument("--lon", type=str, help="Input Unit Hydrograph longitude variable name",default='lon')
    parser.add_argument("--lat", type=str, help="Input Unit Hydrograph latitude variable name",default='lat')
    parser.add_argument("--fraction", type=str, help="Input Unit Hydrograph fraction variable name",default='fraction')
    parser.add_argument("--unit_hydrograph",type=str, help="Input Unit Hydrograph unit hydrograph variable name",default='unit_hydrograph')
    parser.add_argument("--testAgg",help="Do a test aggregation, where all inpoint points are aggregated into one file, remapping can be done afterwards using the --remap flag",action="store_true")
    parser.add_argument("--cdoDebug",help="Enable CDO debuging (prings each step to screen)",action="store_true")
    parser.add_argument("--cdoForce",help="Enable CDO force output (will overwrite existing files during remap)",action="store_true")
    parser.add_argument("--verbose",help="Make script verbose",action="store_true")
    parser.add_argument("--remap",help="Remap the aggregated Unit Hydrographs to outDir and put the aggregated files in the tempDir",action='store_true')
    parser.add_argument("--fill_value",type=float,help="value to use as masked value",default = 9.96920996839e+36)
    parser.add_argument("--pad",type=int,help="Set number of empty cells to include around each aggregated basin",default=10)
    parser.add_argument("--resolution",type=float,help="Set resolution of input Unit Hydrographs",default=1/16.)
    parser.add_argument("--clean",help="Clean up aggregated Unit Hydrograph grids if remapping", action='store_true')
    parser.add_argument("--dryrun",help="Do the mapping between the source and target grid based on the files in the input directory, return the performance stats for the run", action='store_true')
    args = parser.parse_args()

    options = {}
    paths = {}
    # parse the basics
    Rvars = (args.time,args.lon,args.lat,args.fraction,args.unit_hydrograph)
    paths['srcDir'] = args.srcDir
    paths['gridFile'] = args.gridFile

    if args.aggDir:
        paths['aggDir'] = args.aggDir
    else:
        paths['aggDir'] = os.path.join(paths['srcDir'],'../aggregated/')
        if not os.path.exists(paths['aggDir']):
            os.makedirs(paths['aggDir'])

    options['verbose'] = args.verbose
    options['fill_value'] = args.fill_value
    options['pad'] = args.pad
    options['resolution'] = args.resolution
    options['inPrefix'] = args.inPrefix
    options['outPrefix'] = args.outPrefix
    options['dryrun'] = args.dryrun
    options['testAgg'] = args.testAgg
    options['clean']=args.clean
    options['remap']=args.remap
    
    if options['remap']:
        cdo.debug=args.cdoDebug
        cdo.forceOutput=args.cdoForce
        if args.remapDir:
            paths['remapDir'] = args.remapDir
        else:
            paths['remapDir'] = os.path.join(paths['srcDir'],'../remaped/')
            if not os.path.exists(paths['remapDir']):
                os.makedirs(paths['remapDir'])
        print paths['remapDir']        

    return Rvars,paths,options

def clean(aggDir):
    for file in os.listdir(aggDir):
        file_path = os.path.join(aggDir, file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception, e:
            print e
    return
##################################################################################
# Run Program
##################################################################################
if __name__ == "__main__":
    main()
