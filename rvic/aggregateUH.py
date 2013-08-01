#!/usr/local/bin/python
# encoding: utf-8
"""
aggregate.py
"""

import numpy as np
import RvicShare import fillValue_f
import logging

def agg(inData = None, aggData = None, resolution = 0, pad = 0):
    """
    Add the two data sets together and return the combined arrays.
    Expand the horizontal dimensions as necessary to fit inData with aggData.
    The two data sets must include the coordinate variables lon,lat, and time.
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
            inData['unit_hydrograph'] = inData['unit_hydrograph'].filled(fill_value=0.0)
        except:
            pass
    if aggData:
        aggData['unit_hydrograph'][aggData['unit_hydrograph']<0] = 0.0
        try:
            aggData['unit_hydrograph']=aggData['unit_hydrograph'].filled(fill_value=0.0)
        except:
            pass
        
    # Place data
    # First the fractions
    if inData:
        fractions[In[2]:In[3],In[4]:In[5]] += inData['fraction']
    if aggData:
        fractions[Ex[2]:Ex[3],Ex[4]:Ex[5]] += aggData['fraction']

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
        hydrographs[:,yv,xv] /= hydrographs[:,yv,xv].sum(axis=0)

    # Put all the data into aggData variable and return to main
    
    aggData['lon'] = lons
    aggData['lat'] = lats
    aggData['fraction'] = fractions
    aggData['unit_hydrograph'] = hydrographs
    aggData['time'] = times

    return aggData