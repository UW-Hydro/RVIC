#!/usr/local/bin/python
"""
aggregate.py
"""

import numpy as np
from scipy.spatial import cKDTree
from rvic.share import fillValue_f
from rvic.log import log_name
from rvic.utilities import find_nearest
from rvic.vars import Point
import logging

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(log_name)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find target cells for pour points
def MakeAggPairs(lons, lats, DomLon, DomLat, DomIds, agg_type='agg'):
    """
    Group pour points by domain grid outlet cell
    """

    # ---------------------------------------------------------------- #
    #Find Destination grid cells
    log.info('Finding addresses now...')

    if (min(lons) < 0 and DomLon.min() >= 0):
        posinds = np.nonzero(DomLon > 180)
        DomLon[posinds] -= 360
        log.info('adjusted domain lon minimum')

    if agg_type != 'test':
        combined = np.dstack(([DomLat.ravel(), DomLon.ravel()]))[0]
    else:
        # limit the inputs arrays to a single point
        # so that all points are mapped to just one location
        combined = np.dstack(([DomLat[0, 0].ravel(), DomLon[0, 0].ravel()]))[0]
    points = list(np.vstack((np.array(lats), np.array(lons))).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points, k=1)
    indexes = np.array(indexes)
    yinds, xinds = np.unravel_index(np.array(indexes), DomLat.shape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Do the aggregation
    outlets = {}

    for i, ind in enumerate(indexes):
        cell_id = DomIds[yinds[i], xinds[i]]
        if cell_id in outlets:
            outlets[cell_id].PourPoints.append(Point(lat=points[i][0], lon=points[i][1]))
        else:
            outlets[cell_id] = Point(y=yinds[i], x=xinds[i],
                                     lat=combined[ind][0], lon=combined[ind][1])
            outlets[cell_id].PourPoints = [Point(lat=points[i][0], lon=points[i][1])]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Count the pairs
    pp_count = 0
    key_count = 0
    num = len(lons)
    for i, key in enumerate(outlets):
        key_count += 1
        pp_count += len(outlets[key].PourPoints)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Print Debug Results
    log.debug('\n------------------ SUMMARY OF MakeAggPairs ------------------')
    log.debug('NUMBER OF POUR POINTS IN INPUT LIST: %i' % num)
    log.debug('NUMBER OF POINTS TO AGGREGATE TO: %i' % key_count)
    log.debug('NUMBER OF POUR POINTS AGGREGATED: %i' % pp_count)
    log.debug('EFFECIENCY OF: %.2f %%' % (100. * pp_count / num))
    log.debug('UNASSIGNED POUR POINTS: %i \n' % (num - pp_count))
    # ---------------------------------------------------------------- #
    return outlets
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Aggregate the UH grids
def agg(inData, aggData, res=0, pad=0, maskandnorm=False):
    """
    Add the two data sets together and return the combined arrays.
    Expand the horizontal dimensions as necessary to fit inData with aggData.
    The two data sets must include the coordinate variables lon,lat, and time.
    """

    # ---------------------------------------------------------------- #
    # find range of coordinates
    if aggData:
        lat_min = np.minimum(inData['lat'].min(), aggData['lat'].min())
        lat_max = np.maximum(inData['lat'].max(), aggData['lat'].max())
        lon_min = np.minimum(inData['lon'].min(), aggData['lon'].min())
        lon_max = np.maximum(inData['lon'].max(), aggData['lon'].max())
        tshape = inData['UHgrid'].shape[0]
    else:
        lat_min = inData['lat'].min()
        lat_max = inData['lat'].max()
        lon_min = inData['lon'].min()
        lon_max = inData['lon'].max()
        tshape = inData['UHgrid'].shape[0]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # make output arrays for lons/lats and initialize fractions/hydrographs
    # pad output arrays so there is a space =pad around inputs
    lats = np.arange((lat_min - res * (pad)), (lat_max + res * (1 + pad)), res)[::-1]
    lons = np.arange((lon_min - res * (pad)), (lon_max + res * (1 + pad)), res)

    fractions = np.zeros((lats.shape[0], lons.shape[0]))
    hydrographs = np.zeros((tshape, lats.shape[0], lons.shape[0]))
    log.debug('fractions shape %s' % str(fractions.shape))
    log.debug('hydrographs shape %s' % str(hydrographs.shape))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # find target index locations of all corners for both datasets
    # Not that the lat inds are inverted
    ilat_min_ind = find_nearest(lats, np.max(inData['lat']))
    ilat_max_ind = find_nearest(lats, np.min(inData['lat'])) + 1
    ilon_min_ind = find_nearest(lons, np.min(inData['lon']))
    ilon_max_ind = find_nearest(lons, np.max(inData['lon'])) + 1

    log.info('___inData fractions shape: %s' % str(inData['fractions'].shape))

    if aggData:
        alat_min_ind = find_nearest(lats, np.max(aggData['lat']))
        alat_max_ind = find_nearest(lats, np.min(aggData['lat'])) + 1
        alon_min_ind = find_nearest(lons, np.min(aggData['lon']))
        alon_max_ind = find_nearest(lons, np.max(aggData['lon'])) + 1

        log.info('___aggData fractions shape: %s' % str(aggData['fractions'].shape))

        print 'alat_max_ind: %i' % alat_max_ind
        print 'alat_min_ind: %i' % alat_min_ind

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Place data
    fractions[ilat_min_ind:ilat_max_ind, ilon_min_ind:ilon_max_ind] += inData['fractions']
    hydrographs[:, ilat_min_ind:ilat_max_ind, ilon_min_ind:ilon_max_ind] += inData['UHgrid']

    if aggData:
        fractions[alat_min_ind:alat_max_ind, alon_min_ind:alon_max_ind] += aggData['fractions']
        hydrographs[:, alat_min_ind:alat_max_ind, alon_min_ind:alon_max_ind] += aggData['UHgrid']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Mask and Normalize the unit hydrographs
    if (maskandnorm):
        # Normalize the hydrographs (each cell should sum to 1)
        yv, xv = np.nonzero(fractions > 0)
        hydrographs[:, yv, xv] /= hydrographs[:, yv, xv].sum(axis=0)

        # Mask the hydrographs and make sure they sum to 1 at each grid cell
        yv, xv = np.nonzero(((fractions <= 0) * (hydrographs.sum(axis=0) <= 0)))
        hydrographs[:, yv, xv] = fillValue_f
        #hydrographs = np.ma.masked_where(hydrographs == fillValue_f, hydrographs, copy=False)

        log.info('Done with Aggregation')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Put all the data into aggData variable and return to main

    aggData['timesteps'] = inData['timesteps']
    aggData['unit_hydrogaph_dt'] = inData['unit_hydrogaph_dt']
    aggData['lon'] = lons
    aggData['lat'] = lats
    aggData['fractions'] = fractions
    aggData['UHgrid'] = hydrographs
    # ---------------------------------------------------------------- #
    return aggData
# -------------------------------------------------------------------- #
