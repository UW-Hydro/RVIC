# -*- coding: utf-8 -*-
'''
make_uh.py

PROGRAM rout, Python-Version, written by Joe Hamman winter 2012/2013
Routing algorithm developed by D. Lohmann.

Terminology:
- basin refers to the basin mask given in the RVIC parameter file
- catachment refers to the set of points upstream from the outlet
__________________________________
REVISION HISTORY
--------
July 2013, Joe Hamman
Removed read and write functions to make more modular.
Now called from make_parameters.py
'''

import numpy as np
import logging
from scipy.interpolate import interp1d
from .utilities import latlon2yx
from .share import SECSPERDAY
from .log import LOG_NAME
from .pycompat import pyzip, pyrange

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def rout(pour_point, uh_box, fdr_data, fdr_atts, rout_dict):
    '''
    Make the Unit Hydrograph Grid
    '''
    log.info('Starting routing program for point: %s', pour_point)
    # ---------------------------------------------------------------- #
    # Unpack a few structures
    uh_t = uh_box['time']
    uh_box = uh_box['func']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find Basin Dims and ID
    basin_id = fdr_data[rout_dict['BASIN_ID_VAR']][pour_point.routy,
                                                   pour_point.routx]

    log.info('Input Latitude: %f', pour_point.lat)
    log.info('Input Longitude: %f', pour_point.lon)
    log.info('Global Basid ID: %i', basin_id)

    y_inds, x_inds = \
        np.nonzero(fdr_data[rout_dict['BASIN_ID_VAR']] == basin_id)
    y = np.arange(len(fdr_data[rout_dict['LATITUDE_VAR']]))
    x = np.arange(len(fdr_data[rout_dict['LONGITUDE_VAR']]))

    x_min = min(x[x_inds])
    x_max = max(x[x_inds]) + 1
    y_min = min(y[y_inds])
    y_max = max(y[y_inds]) + 1
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create the Basin Dictionary, a subset of the fdr_data
    basin = {}
    basin['lat'] = fdr_data[rout_dict['LATITUDE_VAR']][y_min:y_max]
    basin['lon'] = fdr_data[rout_dict['LONGITUDE_VAR']][x_min:x_max]
    basin['basin_id'] = fdr_data[rout_dict['BASIN_ID_VAR']][y_min:y_max,
                                                            x_min:x_max]
    basin['flow_direction'] = \
        fdr_data[rout_dict['FLOW_DIRECTION_VAR']][y_min:y_max, x_min:x_max]
    basin['flow_distance'] = \
        fdr_data[rout_dict['FLOW_DISTANCE_VAR']][y_min:y_max, x_min:x_max]
    basin['velocity'] = fdr_data['velocity'][y_min:y_max, x_min:x_max]
    basin['diffusion'] = fdr_data['diffusion'][y_min:y_max, x_min:x_max]

    log.debug('Grid cells in subset: %i', basin['velocity'].size)

    pour_point.basiny, pour_point.basinx = latlon2yx(plats=pour_point.lat,
                                                     plons=pour_point.lon,
                                                     glats=basin['lat'],
                                                     glons=basin['lon'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create the rout_data Dictionary
    rout_data = {'lat': basin['lat'], 'lon': basin['lon'],
                 'basiny': pour_point.basiny[0],
                 'basinx': pour_point.basinx[0]}
    log.debug('Clipping indicies:')
    log.debug('rout_x_min: %s', x_min)
    log.debug('rout_x_max: %s', x_max)
    log.debug('rout_y_min: %s', y_min)
    log.debug('rout_y_max: %s', y_max)
    log.debug('pour_point.basiny: %s', pour_point.basiny)
    log.debug('pour_point.basinx: %s', pour_point.basinx)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine/Set flow direction syntax
    # Flow directions {north, northeast, east, southeast,
    # south, southwest, west, northwest}
    if 'VIC' in fdr_atts[rout_dict['FLOW_DIRECTION_VAR']] or \
            fdr_data[rout_dict['FLOW_DIRECTION_VAR']].max() < 10:
        # VIC Directions: http://www.hydro.washington.edu/Lettenmaier/Models/\
        # VIC/Documentation/Routing/FlowDirection.shtml
        dy = {1: -1, 2: -1, 3: 0, 4: 1,
              5: 1, 6: 1, 7: 0, 8: -1}
        dx = {1: 0, 2: 1, 3: 1, 4: 1,
              5: 0, 6: -1, 7: -1, 8: - 1}
        log.debug('Using VIC flow directions (1-8).')
    else:
        # ARCMAP Directions: http://webhelp.esri.com/arcgisdesktop/9.2/\
        # index.cfm?TopicName=flow_direction
        dy = {64: -1, 128: -1, 1: 0, 2: 1,
              4: 1, 8: 1, 16: 0, 32: -1}
        dx = {64: 0, 128: 1, 1: 1, 2: 1,
              4: 0, 8: -1, 16: -1, 32: -1}
        log.debug('Using ARCMAP flow directions (1-128).')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find timestep (timestep is determined from uh_BOX input file)
    input_interval = find_ts(uh_t)
    rout_data['unit_hydrograph_dt'] = input_interval
    t_cell = int(rout_dict['CELL_FLOWDAYS'] * SECSPERDAY / input_interval)
    t_uh = int(rout_dict['BASIN_FLOWDAYS'] * SECSPERDAY / input_interval)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read direction grid and find to_col (to_x) and to_row (to_y)
    to_y, to_x = read_direction(basin['flow_direction'], dy, dx)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find all grid cells upstream of pour point
    catchment, rout_data['fraction'] = search_catchment(to_y, to_x, pour_point,
                                                        basin['basin_id'],
                                                        basin_id)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make uh for each grid cell upstream of basin pour point
    # (linear routing model - Saint-Venant equation)
    uh = make_uh(input_interval, t_cell, catchment['y_inds'],
                 catchment['x_inds'], basin['velocity'], basin['diffusion'],
                 basin['flow_distance'])
    # ---------------------------------------------------------------- #
    # Make uh_river by incrementally moving upstream comining uh functions
    uh_river = make_grid_uh_river(t_uh, t_cell, uh, to_y, to_x, pour_point,
                                  catchment['y_inds'], catchment['x_inds'],
                                  catchment['count_ds'])
    # ---------------------------------------------------------------- #
    # Make uh_s for each grid cell upstream of basin pour point
    # (combine IRFs for all grid cells in flow path)
    uh_s = make_grid_uh(t_uh, t_cell, uh_river, uh_box, to_y, to_x,
                        catchment['y_inds'], catchment['x_inds'],
                        catchment['count_ds'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Agregate to output timestep
    rout_data['unit_hydrograph'] = adjust_uh_timestep(
        uh_s, t_uh, input_interval, rout_dict['OUTPUT_INTERVAL'],
        catchment['x_inds'], catchment['y_inds'])
    # ---------------------------------------------------------------- #
    return rout_data
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the Timestep from the uh box
def find_ts(uh_t):
    '''
    Determines the (input_interval) based on the timestep given in uhfile
    '''
    input_interval = uh_t[1] - uh_t[0]
    log.debug('Input Timestep = %i seconds', input_interval)
    return input_interval
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the flow direction file
def read_direction(fdr, dy, dx):
    '''
    Reads the direction file and makes two grids (to_x) and (to_y).
    The input grids follow the 1-8 or 1-128 grid directions as shown below.
    val = direction  [to_y][to_x]
    '''
    log.debug('Reading direction input and finding target row/columns')

    to_y = np.zeros_like(fdr, dtype=np.int16)
    to_x = np.zeros_like(fdr, dtype=np.int16)

    valid_dirs = list(dy.keys())

    for (y, x), d in np.ndenumerate(fdr):
        if d in valid_dirs:
            to_y[y, x] = y + dy[d]
            to_x[y, x] = x + dx[d]
        else:
            to_y[y, x] = -9999
            to_x[y, x] = -9999

    return to_y, to_x
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Search the catchment
def search_catchment(to_y, to_x, pour_point, basin_ids, basin_id):
    '''
    Find all cells upstream of pour point.  Retrun a dictionary with x_inds,
    yinds, and #of cell to downstream pour point.  All are sorted the by the
    latter. For each x,y pair, the flow path is followed until either the
    catchment outlet is encountered
    (if (yy==pour_point.basiny and xx==pour_point.basinx):)
    or the flowpath leads outside of grid.
    *** Does not handle wrapped coordinates. ***
    '''
    log.debug('Searching catchment')

    (len_y, len_x) = to_x.shape

    byinds, bxinds = np.nonzero(basin_ids == basin_id)
    bsize = len(byinds)

    # out fractions
    catch_fracs = np.zeros_like(to_x, dtype=np.float64)

    # temporary 2d arrays
    count_ds = np.zeros_like(to_x, dtype=np.uint32)

    # In catch array (starts as -1)
    # -1 - not in catchment
    # 0 - unknown
    # 1 - in catchment
    in_catch = np.zeros_like(to_x, dtype=np.int16) - 1
    in_catch[byinds, bxinds] = 0  # set basin inds as 0

    # temporary variables for tracking flow path
    path_count = np.zeros(bsize, dtype=np.uint32)
    pathy = np.zeros(bsize, dtype=np.uint32)
    pathx = np.zeros(bsize, dtype=np.uint32)

    cells = 0

    for yy, xx in pyzip(byinds, bxinds):
        if in_catch[yy, xx] >= 0:
            # set the old path to zero
            pathy[:cells + 1] = 0
            pathx[:cells + 1] = 0

            # reset the cells counter
            cells = 0

            while True:
                pathy[cells] = yy
                pathx[cells] = xx
                path_count[cells] = cells

                # stop if you reach the pour point or if you reach a point you
                # know is in the catchment
                if ((yy == pour_point.basiny) and (xx == pour_point.basinx)) \
                        or in_catch[yy, xx] == 1:

                    # get the path inds
                    py = pathy[:cells + 1]
                    px = pathx[:cells + 1]

                    # Add path to in_catch and count_ds
                    in_catch[py, px] = 1

                    # reverse the path count
                    temp_path = path_count[:cells + 1][::-1]
                    count_ds[py, px] = (count_ds[yy, xx] + temp_path)

                    break
                else:
                    # Move downstream
                    yy, xx = to_y[yy, xx], to_x[yy, xx]
                    cells += 1
                    if (xx == len_x) or (xx < 0) or (yy == len_y) or (yy < 0) \
                            or (in_catch[yy, xx] == -1):
                        # That path is not in the catchment
                        py = pathy[:cells]
                        px = pathx[:cells]
                        in_catch[py, px] = -1  # set to -1
                        break

    catchment = {}
    cyinds, cxinds = np.nonzero(in_catch == 1)
    catchment['count_ds'] = count_ds[cyinds, cxinds]
    catch_fracs[cyinds, cxinds] = 1.0

    count = len(cyinds)

    log.debug('Found %i upstream grid cells from present station', count)
    log.debug('Expected at most %i upstream grid cells from present station',
              bsize)

    if count > bsize:
        log.exception('Error, too many points found.')
        raise

    # ---------------------------------------------------------------- #
    # sort catchment
    ii = np.argsort(catchment['count_ds'])
    catchment['count_ds'] = catchment['count_ds'][ii]
    catchment['x_inds'] = cxinds[ii]
    catchment['y_inds'] = cyinds[ii]
    # ---------------------------------------------------------------- #
    return catchment, catch_fracs
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make uh Grid
def make_uh(dt, t_cell, y_inds, x_inds, velocity, diffusion, xmask):
    '''
    Calculate the impulse response function for grid cells using equation 15
    from Lohmann, et al. (1996) Tellus article.  Return 3d uh grid.
    '''
    log.debug('Making uh for each cell')

    uh = np.zeros((t_cell, xmask.shape[0], xmask.shape[1]), dtype=np.float64)
    time = np.arange(dt, t_cell * dt + dt, dt, dtype=np.float64)

    for y, x in pyzip(y_inds, x_inds):
        xm = xmask[y, x]
        v = velocity[y, x]
        d = diffusion[y, x]

        exponent = -1 * np.power(v * time - xm, 2) / (4 * d * time)
        green = xm / (2 * time * np.sqrt(np.pi * time * d)) * np.exp(exponent)

        # Normalize
        uh[:, y, x] = green / green.sum()
    return uh
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make uh river
def make_grid_uh_river(t_uh, t_cell, uh, to_y, to_x, pour_point, y_inds,
                       x_inds, count_ds):
    '''
    Calculate impulse response function for river routing.  Starts at
    downstream point incrementally moves upstream.
    '''
    log.debug('Making uh_river grid....')
    y_ind = pour_point.basiny
    x_ind = pour_point.basinx

    uh_river = np.zeros((t_uh, uh.shape[1], uh.shape[2]), dtype=np.float64)

    for (y, x, d) in pyzip(y_inds, x_inds, count_ds):
        if d > 0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            irf_temp = np.convolve(uh_river[:, yy, xx], uh[:, y, x])

            # Normalize
            uh_river[:, y, x] = irf_temp[:t_uh] / irf_temp[:t_uh].sum()
        elif d == 0:
            # Just use the UH calculated previously
            uh_river[:t_cell, y_ind, x_ind] = uh[:, y_ind, x_ind]
        else:
            raise ValueError('Got negative value ({0}) for count_ds at y={1}'
                             'x={2}'.format(d, y, x))

    return uh_river
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make grid uh
def make_grid_uh(t_uh, t_cell, uh_river, uh_box, to_y, to_x, y_inds, x_inds,
                 count_ds):
    '''
    Combines the uh_box with downstream cell uh_river.  Cell [0] is given the
    uh_box without river routing
    '''
    log.debug('Making unit_hydrograph grid')

    unit_hydrograph = np.zeros((t_uh, uh_river.shape[1], uh_river.shape[2]))
    irf_temp = np.zeros(t_uh + t_cell, dtype=np.float64)

    for (y, x, d) in pyzip(y_inds, x_inds, count_ds):
        irf_temp[:] = 0.0
        if d > 0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            irf_temp = np.convolve(uh_box, uh_river[:, yy, xx])
            unit_hydrograph[:, y, x] = \
                irf_temp[:t_uh] / np.sum(irf_temp[:t_uh])
        else:
            irf_temp[:len(uh_box)] = uh_box[:]
            unit_hydrograph[:, y, x] = irf_temp[:t_uh]

    return unit_hydrograph
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Adjust the timestep
def adjust_uh_timestep(unit_hydrograph, t_uh, input_interval, output_interval,
                       x_inds, y_inds):
    '''
    Aggregates to timestep (output_interval).  output_interval must be a
    multiple of input_interval.  This function is not setup to disaggregate
    the Unit Hydrographs to a output_interval<input_interval.
    '''
    if output_interval == input_interval:
        log.debug('No need to aggregate in time (output_interval = '
                  'input_interval) Skipping the adjust_uh_timestep step')
        uh_out = unit_hydrograph
    elif np.remainder(output_interval, input_interval) == 0:
        log.debug('Aggregating to %i from %i seconds', output_interval,
                  input_interval)
        fac = int(output_interval / input_interval)
        t_uh_out = int(t_uh / fac)
        uh_out = np.zeros((t_uh_out, unit_hydrograph.shape[1],
                          unit_hydrograph.shape[2]), dtype=np.float64)
        for (y, x) in pyzip(y_inds, x_inds):
            for t in pyrange(t_uh_out):
                uh_out[t, y, x] = unit_hydrograph[t * fac:t * fac + fac,
                                                  y, x].sum()
    else:
        log.debug('Interpolating unit hydrograph from input_interval: %i to '
                  'output_interval: %i', input_interval, output_interval)
        fac = int(input_interval / output_interval)
        t_uh_out = int(t_uh * fac)
        uh_out = np.zeros((t_uh_out,
                          unit_hydrograph.shape[1],
                          unit_hydrograph.shape[2]), dtype=np.float64)
        ts_orig = np.linspace(0, t_uh, t_uh)
        ts_new = np.linspace(0, t_uh, t_uh_out)
        for (y, x) in pyzip(y_inds, x_inds):
            f = interp1d(ts_orig, unit_hydrograph[:, y, x])
            uh_out[:, y, x] = f(ts_new)

        # normalize
        uh_out /= uh_out.sum(axis=0)
    return uh_out
# -------------------------------------------------------------------- #
