"""
make_uh.py

PROGRAM rout, Python-Version, written by Joe Hamman winter 2012/2013
Routing algorithm developed by D. Lohmann.
__________________________________
REVISION HISTORY
--------
July 2013, Joe Hamman
Removed read and write functions to make more modular.
Now called from make_parameters.py
"""

import numpy as np
import logging
from scipy.interpolate import interp1d
from utilities import find_nearest
from variables import Point
from share import SECSPERDAY, PRECISION
from log import LOG_NAME

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def rout(pour_point, uh_box, fdr_data, fdr_atts, rout_dict):
    """
    Make the Unit Hydrograph Grid
    """
    log.info("Starting routing program...")
    # ---------------------------------------------------------------- #
    # Unpack a few structures
    uh_t = uh_box['time']
    uh_box = uh_box['func']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find Basin Dims and ID
    pour_point.x = find_nearest(fdr_data[rout_dict['LONGITUDE_VAR']], pour_point.lon)
    pour_point.y = find_nearest(fdr_data[rout_dict['LATITUDE_VAR']], pour_point.lat)
    basin_id = fdr_data[rout_dict['BASIN_ID_VAR']][pour_point.y, pour_point.x]

    log.info('Input Latitude: %f' % pour_point.lat)
    log.info('Input Longitude: %f' % pour_point.lon)
    log.info('Global Basid ID: %i' % basin_id)

    y_inds, x_inds = np.nonzero(fdr_data[rout_dict['BASIN_ID_VAR']] == basin_id)
    y = np.arange(len(fdr_data[rout_dict['LATITUDE_VAR']]))
    x = np.arange(len(fdr_data[rout_dict['LONGITUDE_VAR']]))

    x_min = min(x[x_inds])
    x_max = max(x[x_inds])+1
    y_min = min(y[y_inds])
    y_max = max(y[y_inds])+1
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create the Basin Dictionary, a subset of the fdr_data
    basin = {}
    basin['lat'] = fdr_data[rout_dict['LATITUDE_VAR']][y_min:y_max]
    basin['lon'] = fdr_data[rout_dict['LONGITUDE_VAR']][x_min:x_max]
    basin['basin_id'] = fdr_data[rout_dict['BASIN_ID_VAR']][y_min:y_max, x_min:x_max]
    basin['flow_direction'] = fdr_data[rout_dict['FLOW_DIRECTION_VAR']][y_min:y_max, x_min:x_max]
    basin['flow_distance'] = fdr_data[rout_dict['FLOW_DISTANCE_VAR']][y_min:y_max, x_min:x_max]
    basin['velocity'] = fdr_data['velocity'][y_min:y_max, x_min:x_max]
    basin['diffusion'] = fdr_data['diffusion'][y_min:y_max, x_min:x_max]

    log.debug('Grid cells in subset: %i' % basin['velocity'].size)

    basin_point = Point(x=find_nearest(basin['lon'], pour_point.lon),
                        y=find_nearest(basin['lat'], pour_point.lat),
                        lon=pour_point.lon, lat=pour_point.lat)

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create the rout_data Dictionary
    rout_data = {'lat': basin['lat'], 'lon': basin['lon']}

    # ---------------------------------------------------------------- #
    # Determine low direction syntax
    if 'VIC' in fdr_atts[rout_dict['FLOW_DIRECTION_VAR']]:
        # VIC Directions: http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/Routing/FlowDirection.shtml
        dy = {1: -1, 2: -1, 3: 0, 4: 1, 5: 1, 6: 1, 7: 0, 8: -1}
        dx = {1: 0, 2: 1, 3: 1, 4: 1, 5: 0, 6: -1, 7: -1, 8: - 1}
        log.debug('Using VIC flow directions (1-8).')
    else:
        # ARCMAP Directions: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=flow_direction
        dy = {1: 0, 2: 1, 4: 1, 8: 1, 16: 0, 32: -1, 64: -1, 128: -1}
        dx = {1: 1, 2: 1, 4: 0, 8: -1, 16: -1, 32: -1, 64: 0, 128: 1}
        log.debug('Using ARCMAP flow directions (1-128).')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find timestep (timestep is determined from uh_BOX input file)
    input_interval = find_ts(uh_t)
    rout_data['unit_hydrograph_dt'] = input_interval
    t_cell = int(rout_dict['CELL_FLOWDAYS']*SECSPERDAY/input_interval)
    t_uh = int(rout_dict['BASIN_FLOWDAYS']*SECSPERDAY/input_interval)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read direction grid and find to_col (to_x) and to_row (to_y)
    to_y, to_x = read_direction(basin['flow_direction'], dy, dx)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find all grid cells upstream of pour point
    catchment, rout_data['fraction'] = search_catchment(to_y, to_x, basin_point,
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
    uh_river = make_grid_uh_river(t_uh, t_cell, uh, to_y, to_x, basin_point,
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
    rout_data['unit_hydrograph'], rout_data['timesteps'] = adjust_uh_timestep(uh_s, t_uh,
                                                                     input_interval,
                                                                     rout_dict['OUTPUT_INTERVAL'],
                                                                     catchment['x_inds'],
                                                                     catchment['y_inds'])
    # ---------------------------------------------------------------- #
    return rout_data
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the Timestep from the uh box
def find_ts(uh_t):
    """
    Determines the (input_interval) based on the timestep given in uhfile
    """
    input_interval = uh_t[1]-uh_t[0]
    log.debug('Input Timestep = %i seconds' % input_interval)
    return input_interval
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the flow direction file
def read_direction(fdr, dy, dx):
    """
    Reads the direction file and makes two grids (to_x) and (to_y).
    The input grids follow the 1-8 or 1-128 grid directions as shown below.
    val = direction  [to_y][to_x]
    """
    log.debug('Reading direction input and finding target row/columns')

    to_y = np.zeros(fdr.shape, dtype=int)
    to_x = np.zeros(fdr.shape, dtype=int)

    for (y, x), d in np.ndenumerate(fdr):
        try:
            to_y[y, x] = y+dy[d]
            to_x[y, x] = x+dx[d]
        except:
            to_y[y, x] = -9999
            to_x[y, x] = -9999

    return to_y, to_x
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Search the catchment
def search_catchment(to_y, to_x, basin_point, basin_ids, basin_id):
    """
    Find all cells upstream of pour point.  Retrun a dictionary with x_inds,
    yinds, and #of cell to downstream pour point.  All are sorted the by the
    latter. For each x,y pair, the flow path is followed until either the
    catchment outlet is encountered (if (yy==basin_point.y and xx==basin_point.x):)
    or the flowpath leads outside of grid.
    *** Does not handle wrapped coordinates. ***
    """
    log.debug('Searching catchment')

    count = 0
    (len_y, len_x) = to_x.shape
    catchment = {}

    yinds, xinds = np.nonzero(basin_ids == basin_id)

    fractions = np.zeros((len_y, len_x))
    fractions[yinds, xinds] = 1.0
    catchment['count_ds'] = np.empty(len(yinds))

    for i, (y, x) in enumerate(zip(yinds, xinds)):
        yy, xx = y, x
        cells = 0
        while True:
            if (yy == basin_point.y and xx == basin_point.x):
                catchment['count_ds'][i] = cells
                count += 1
                break
            else:
                yy, xx = to_y[yy, xx], to_x[yy, xx]
                cells += 1
                if ((xx > (len_x - 1)) or (xx < 0) or (yy > (len_y - 1)) or (yy < 0)):
                    break

    log.debug("Found %i upstream grid cells from present station" % count)
    log.debug("Expected at most %i upstream grid cells from present station" % len(yinds))
    if count>len(yinds):
        log.exception('Error, too many points found.')
        raise

    # ---------------------------------------------------------------- #
    # sort catchment
    ii = np.argsort(catchment['count_ds'])
    catchment['count_ds'] = catchment['count_ds'][ii]
    catchment['x_inds'] = xinds[ii]
    catchment['y_inds'] = yinds[ii]
    # ---------------------------------------------------------------- #
    return catchment, fractions
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make uh Grid
def make_uh(dt, t_cell, y_inds, x_inds, velocity, diffusion, xmask):
    """
    Calculate the impulse response function for grid cells using equation 15
    from Lohmann, et al. (1996) Tellus article.  Return 3d uh grid.
    """
    log.debug('Making uh for each cell')

    uh = np.zeros((t_cell, xmask.shape[0], xmask.shape[1]))
    for (y, x) in zip(y_inds, x_inds):
        time = dt
        flag = 0
        t = 0
        green = np.zeros(t_cell)
        while (t < t_cell and flag == 0):
            exponent = -1*np.power(velocity[y, x]*time-xmask[y, x], 2)/(4*diffusion[y, x]*time)
            if exponent > np.log(PRECISION):
                green[t] = xmask[y, x]/(2*time*np.sqrt(np.pi*time*diffusion[y, x]))*np.exp(exponent)
                t += 1
                time = time+dt
            else:
                flag = 1
        tot = np.sum(green)
        if tot > 0.:
            uh[:, y, x] = green[:]/tot
    return uh
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Make uh river
def make_grid_uh_river(t_uh, t_cell, uh, to_y, to_x, basin_point, y_inds,
                       x_inds, count_ds):
    """
    Calculate impulse response function for river routing.  Starts at
    downstream point incrementally moves upstream.
    """
    log.debug("Making uh_river grid.... It takes a while...")
    y_ind = basin_point.y
    x_ind = basin_point.x

    uh_river = np.zeros((t_uh, uh.shape[1], uh.shape[2]))
    for (y, x, d) in zip(y_inds, x_inds, count_ds):
        if d > 0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            irf_temp = np.zeros(t_uh+t_cell)
            active_timesteps = np.nonzero(uh_river[:, yy, xx] > PRECISION)[0]
            for t in active_timesteps:
                for l in xrange(t_cell):
                    irf_temp[t+l] = irf_temp[t + l] + uh[l, y, x] * uh_river[t, yy, xx]
            tot = np.sum(irf_temp[:t_uh])
            if tot > 0:
                uh_river[:, y, x] = irf_temp[:t_uh] / tot
        elif d == 0:
            uh_river[:t_cell, y_ind, x_ind] = uh[:, y_ind, x_ind]

    return uh_river
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make grid uh
def make_grid_uh(t_uh, t_cell, uh_river, uh_BOX, to_y, to_x, y_inds, x_inds,
                 count_ds):
    """
    Combines the uh_BOX with downstream cell uh_river.  Cell [0] is given the
    uh_box without river routing
    """
    log.debug("Making unit_hydrograph grid")

    unit_hydrograph = np.zeros((t_uh, uh_river.shape[1], uh_river.shape[2]))
    for (y, x, d) in zip(y_inds, x_inds, count_ds):
        irf_temp = np.zeros(t_uh+t_cell)
        if d > 0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            active_timesteps = np.nonzero(uh_river[:, yy, xx] > PRECISION)[0]
            for t in active_timesteps:
                for l in xrange(len(uh_BOX)):
                    irf_temp[t + l] = irf_temp[t + l] + uh_BOX[l] * uh_river[t, yy, xx]
            tot = np.sum(irf_temp[:t_uh])
            if tot > 0:
                unit_hydrograph[:, y, x] = irf_temp[:t_uh] / tot
        else:
            irf_temp[:len(uh_BOX)] = uh_BOX[:]
            unit_hydrograph[:, y, x] = irf_temp[:t_uh]
    return unit_hydrograph
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Adjust the timestep
def adjust_uh_timestep(unit_hydrograph, t_uh, input_interval, output_interval, x_inds, y_inds):
    """
    Aggregates to timestep (output_interval).  output_interval must be a
    multiple of input_interval.  This function is not setup to disaggregate
    the Unit Hydrographs to a output_interval<input_interval.
    """
    if output_interval == input_interval:
        log.debug('No need to aggregate (output_interval = input_interval) \
                      Skipping the adjust_uh_timestep step')
        uh_out = unit_hydrograph
        ts_new = np.arange(t_uh)
    elif np.remainder(output_interval, input_interval) == 0:
        log.debug('Aggregating to %i from %i seconds' % (output_interval, input_interval))
        fac = int(output_interval/input_interval)
        t_uh_out = int(t_uh/fac)
        ts_new = np.arange(t_uh_out)
        uh_out = np.zeros((t_uh_out, unit_hydrograph.shape[1], unit_hydrograph.shape[2]))
        for (y, x) in zip(y_inds, x_inds):
            for t in xrange(t_uh_out):
                uh_out[t, y, x] = np.sum(unit_hydrograph[t*fac:t*fac+fac, y, x])
    elif np.remainder(input_interval, output_interval):
        log.debug('Interpolating unit hydrograph from input_interval: %i to output_interval: %i' % (input_interval, output_interval))
        fac = int(input_interval / output_interval)
        t_uh_out = int(t_uh * fac)
        uh_out = np.zeros((t_uh_out, unit_hydrograph.shape[1], unit_hydrograph.shape[2]))
        ts_orig = np.linspace(0, t_uh, t_uh)
        ts_new = np.linspace(0, t_uh, t_uh_out)
        for (y, x) in zip(y_inds, x_inds):
            f = interp1d(ts_orig, unit_hydrograph[:, y, x])
            uh_out[:, y, x] = f(ts_new)
    return uh_out, ts_new
# -------------------------------------------------------------------- #

