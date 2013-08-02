#!/usr/local/bin/python
"""
aggregate.py
"""

import numpy as np
from scipy.spatial import cKDTree
from share import FILLVALUE_F
from utilities import find_nearest
from variables import Point
from logging import getLogger
from log import LOG_NAME

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find target cells for pour points
def make_agg_pairs(lons, lats, dom_lon, dom_lat, dom_ids, agg_type='agg'):
    """
    Group pour points by domain grid outlet cell
    """

    # ---------------------------------------------------------------- #
    #Find Destination grid cells
    log.info('Finding addresses now...')

    if (min(lons) < 0 and dom_lon.min() >= 0):
        posinds = np.nonzero(dom_lon > 180)
        dom_lon[posinds] -= 360
        log.info('adjusted domain lon minimum')

    if agg_type != 'test':
        combined = np.dstack(([dom_lat.ravel(), dom_lon.ravel()]))[0]
    else:
        # limit the inputs arrays to a single point
        # so that all points are mapped to just one location
        combined = np.dstack(([dom_lat[0, 0].ravel(), dom_lon[0, 0].ravel()]))[0]
    points = list(np.vstack((np.array(lats), np.array(lons))).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points, k=1)
    indexes = np.array(indexes)
    yinds, xinds = np.unravel_index(np.array(indexes), dom_lat.shape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Do the aggregation
    outlets = {}

    for i, ind in enumerate(indexes):
        cell_id = dom_ids[yinds[i], xinds[i]]
        if cell_id in outlets:
            outlets[cell_id].pour_points.append(Point(lat=points[i][0], lon=points[i][1]))
        else:
            outlets[cell_id] = Point(y=yinds[i], x=xinds[i],
                                     lat=combined[ind][0], lon=combined[ind][1])
            outlets[cell_id].pour_points = [Point(lat=points[i][0], lon=points[i][1])]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Count the pairs
    pp_count = 0
    key_count = 0
    num = len(lons)
    for i, key in enumerate(outlets):
        key_count += 1
        pp_count += len(outlets[key].pour_points)
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
def aggregate(in_data, agg_data, res=0, pad=0, maskandnorm=False):
    """
    Add the two data sets together and return the combined arrays.
    Expand the horizontal dimensions as necessary to fit in_data with agg_data.
    The two data sets must include the coordinate variables lon,lat, and time.
    """

    # ---------------------------------------------------------------- #
    # find range of coordinates
    if agg_data:
        lat_min = np.minimum(in_data['lat'].min(), agg_data['lat'].min())
        lat_max = np.maximum(in_data['lat'].max(), agg_data['lat'].max())
        lon_min = np.minimum(in_data['lon'].min(), agg_data['lon'].min())
        lon_max = np.maximum(in_data['lon'].max(), agg_data['lon'].max())
        tshape = in_data['UHgrid'].shape[0]
    else:
        lat_min = in_data['lat'].min()
        lat_max = in_data['lat'].max()
        lon_min = in_data['lon'].min()
        lon_max = in_data['lon'].max()
        tshape = in_data['UHgrid'].shape[0]
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
    ilat_min_ind = find_nearest(lats, np.max(in_data['lat']))
    ilat_max_ind = find_nearest(lats, np.min(in_data['lat'])) + 1
    ilon_min_ind = find_nearest(lons, np.min(in_data['lon']))
    ilon_max_ind = find_nearest(lons, np.max(in_data['lon'])) + 1

    log.info('___in_data fractions shape: %s' % str(in_data['fractions'].shape))

    if agg_data:
        alat_min_ind = find_nearest(lats, np.max(agg_data['lat']))
        alat_max_ind = find_nearest(lats, np.min(agg_data['lat'])) + 1
        alon_min_ind = find_nearest(lons, np.min(agg_data['lon']))
        alon_max_ind = find_nearest(lons, np.max(agg_data['lon'])) + 1

        log.info('___agg_data fractions shape: %s' % str(agg_data['fractions'].shape))

        print 'alat_max_ind: %i' % alat_max_ind
        print 'alat_min_ind: %i' % alat_min_ind

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Place data
    fractions[ilat_min_ind:ilat_max_ind, ilon_min_ind:ilon_max_ind] += in_data['fractions']
    hydrographs[:, ilat_min_ind:ilat_max_ind, ilon_min_ind:ilon_max_ind] += in_data['UHgrid']

    if agg_data:
        fractions[alat_min_ind:alat_max_ind, alon_min_ind:alon_max_ind] += agg_data['fractions']
        hydrographs[:, alat_min_ind:alat_max_ind, alon_min_ind:alon_max_ind] += agg_data['UHgrid']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Mask and Normalize the unit hydrographs
    if (maskandnorm):
        # Normalize the hydrographs (each cell should sum to 1)
        yv, xv = np.nonzero(fractions > 0)
        hydrographs[:, yv, xv] /= hydrographs[:, yv, xv].sum(axis=0)

        # Mask the hydrographs and make sure they sum to 1 at each grid cell
        yv, xv = np.nonzero(((fractions <= 0) * (hydrographs.sum(axis=0) <= 0)))
        hydrographs[:, yv, xv] = FILLVALUE_F
        #hydrographs = np.ma.masked_where(hydrographs == FILLVALUE_F, hydrographs, copy=False)

        log.info('Done with Aggregation')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Put all the data into agg_data variable and return to main

    agg_data['timesteps'] = in_data['timesteps']
    agg_data['unit_hydrogaph_dt'] = in_data['unit_hydrogaph_dt']
    agg_data['lon'] = lons
    agg_data['lat'] = lats
    agg_data['fractions'] = fractions
    agg_data['UHgrid'] = hydrographs
    # ---------------------------------------------------------------- #
    return agg_data
# -------------------------------------------------------------------- #
