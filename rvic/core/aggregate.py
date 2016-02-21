# -*- coding: utf-8 -*-
'''
aggregate.py
'''

import numpy as np
from .share import FILLVALUE_F
from .utilities import find_nearest, latlon2yx
from .variables import Point
from logging import getLogger
from .log import LOG_NAME
from .pycompat import OrderedDict, iteritems, pyzip

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find target cells for pour points
def make_agg_pairs(pour_points, dom_data, fdr_data, config_dict):
    '''
    Group pour points by domain grid outlet cell
    '''
    lons = pour_points['lons']
    lats = pour_points['lats']
    dom_lon = dom_data[config_dict['DOMAIN']['LONGITUDE_VAR']]
    dom_lat = dom_data[config_dict['DOMAIN']['LATITUDE_VAR']]
    dom_ids = dom_data['cell_ids']
    fdr_lons = fdr_data[config_dict['ROUTING']['LONGITUDE_VAR']]
    fdr_lats = fdr_data[config_dict['ROUTING']['LATITUDE_VAR']]
    fdr_srcarea = fdr_data[config_dict['ROUTING']['SOURCE_AREA_VAR']]

    # ---------------------------------------------------------------- #
    # Find Destination grid cells
    log.info('Finding addresses now...')
    routys, routxs = latlon2yx(plats=lats,
                               plons=lons,
                               glats=fdr_lats,
                               glons=fdr_lons)

    domys, domxs = latlon2yx(plats=lats,
                             plons=lons,
                             glats=dom_lat,
                             glons=dom_lon)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Do the aggregation
    outlets = OrderedDict()

    for i, (lat, lon) in enumerate(pyzip(lats, lons)):
        # Define pour point object
        pour_point = Point(lat=lat,
                           lon=lon,
                           domx=domxs[i],
                           domy=domys[i],
                           routx=routxs[i],
                           routy=routys[i],
                           name=None,
                           cell_id=dom_ids[domys[i], domxs[i]])
        pour_point.source_area = fdr_srcarea[pour_point.routy,
                                             pour_point.routx]

        cell_id = dom_ids[domys[i], domxs[i]]

        if cell_id in outlets:
            outlets[cell_id].pour_points.append(pour_point)
            outlets[cell_id].upstream_area += pour_point.source_area
        else:
            # define outlet grid cell (on domain grid)
            outlets[cell_id] = Point(domy=domys[i],
                                     domx=domxs[i],
                                     lat=dom_lat[domys[i], domxs[i]],
                                     lon=dom_lon[domys[i], domxs[i]])
            outlets[cell_id].pour_points = [pour_point]
            outlets[cell_id].cell_id = cell_id
            outlets[cell_id].upstream_area = pour_point.source_area
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Sort based on outlet total source area pour_point.source_area
    outlets = OrderedDict(sorted(list(iteritems(outlets)),
                          key=lambda t: t[1].upstream_area,
                          reverse=True))
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
    log.info('\n------------------ SUMMARY OF MakeAggPairs ------------------')
    log.info('NUMBER OF POUR POINTS IN INPUT LIST: %i', num)
    log.info('NUMBER OF POINTS TO AGGREGATE TO: %i', key_count)
    log.info('NUMBER OF POUR POINTS AGGREGATED: %i', pp_count)
    log.info('EFFECIENCY OF: %.2f %%', (100. * pp_count / num))
    log.info('UNASSIGNED POUR POINTS: %i \n', (num - pp_count))
    log.info('-------------------------------------------------------------\n')

    # ---------------------------------------------------------------- #
    return outlets
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Aggregate the UH grids
def aggregate(in_data, agg_data, res=0, pad=0, maskandnorm=False):
    '''
    Add the two data sets together and return the combined arrays.
    Expand the horizontal dimensions as necessary to fit in_data with agg_data.
    The two data sets must include the coordinate variables lon,lat, and time.
    '''

    # ---------------------------------------------------------------- #
    # find range of coordinates
    if agg_data:
        lat_min = np.minimum(in_data['lat'].min(), agg_data['lat'].min())
        lat_max = np.maximum(in_data['lat'].max(), agg_data['lat'].max())
        lon_min = np.minimum(in_data['lon'].min(), agg_data['lon'].min())
        lon_max = np.maximum(in_data['lon'].max(), agg_data['lon'].max())
        tshape = in_data['unit_hydrograph'].shape[0]
    else:
        lat_min = in_data['lat'].min()
        lat_max = in_data['lat'].max()
        lon_min = in_data['lon'].min()
        lon_max = in_data['lon'].max()
        tshape = in_data['unit_hydrograph'].shape[0]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # make output arrays for lons/lats and initialize fraction/unit_hydrograph
    # pad output arrays so there is a space =pad around inputs
    lats = np.arange((lat_min - res * pad),
                     (lat_max + res * (pad + 1)), res)[::-1]
    lons = np.arange((lon_min - res * pad),
                     (lon_max + res * (pad + 1)), res)

    fraction = np.zeros((len(lats), len(lons)), dtype=np.float64)
    unit_hydrograph = np.zeros((tshape, len(lats), len(lons)),
                               dtype=np.float64)
    log.debug('fraction shape %s', fraction.shape)
    log.debug('unit_hydrograph shape %s', unit_hydrograph.shape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # find target index locations of all corners for both datasets
    # Not that the lat inds are inverted
    ilat_min_ind = find_nearest(lats, np.max(in_data['lat']))
    ilat_max_ind = find_nearest(lats, np.min(in_data['lat'])) + 1
    ilon_min_ind = find_nearest(lons, np.min(in_data['lon']))
    ilon_max_ind = find_nearest(lons, np.max(in_data['lon'])) + 1

    log.debug('in_data fraction shape: %s', in_data['fraction'].shape)

    if agg_data:
        alat_min_ind = find_nearest(lats, np.max(agg_data['lat']))
        alat_max_ind = find_nearest(lats, np.min(agg_data['lat'])) + 1
        alon_min_ind = find_nearest(lons, np.min(agg_data['lon']))
        alon_max_ind = find_nearest(lons, np.max(agg_data['lon'])) + 1

        log.debug('agg_data fraction shape: %s', agg_data['fraction'].shape)

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Place data
    fraction[ilat_min_ind:ilat_max_ind,
             ilon_min_ind:ilon_max_ind] += in_data['fraction']
    unit_hydrograph[:,
                    ilat_min_ind:ilat_max_ind,
                    ilon_min_ind:ilon_max_ind] += in_data['unit_hydrograph']

    if agg_data:
        fraction[alat_min_ind:alat_max_ind,
                 alon_min_ind:alon_max_ind] += agg_data['fraction']
        unit_hydrograph[:, alat_min_ind:alat_max_ind,
                        alon_min_ind:alon_max_ind] += \
            agg_data['unit_hydrograph']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Mask and Normalize the unit unit_hydrograph
    if (maskandnorm):
        # Normalize the unit_hydrograph (each cell should sum to 1)

        yv, xv = np.nonzero(fraction > 0.0)
        unit_hydrograph[:, yv, xv] /= unit_hydrograph[:, yv, xv].sum(axis=0)

        # Mask the unit_hydrograph and make sure they sum to 1 at each grid
        # cell
        ym, xm = np.nonzero(fraction <= 0.0)
        unit_hydrograph[:, ym, xm] = FILLVALUE_F

        log.info('Completed final aggregation step')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Put all the data into agg_data variable and return to main
    agg_data['unit_hydrograph_dt'] = in_data['unit_hydrograph_dt']
    agg_data['lon'] = lons
    agg_data['lat'] = lats
    agg_data['fraction'] = fraction
    agg_data['unit_hydrograph'] = unit_hydrograph
    # ---------------------------------------------------------------- #
    return agg_data
# -------------------------------------------------------------------- #
