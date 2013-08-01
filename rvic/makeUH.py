"""
PROGRAM rout, Python-Version, written by Joe Hamman winter 2012/2013
Routing algorithm developed by D. Lohmann.
"""

import numpy as np
import logging
from scipy.interpolate import interp1d
from rvic.share import secsPerDay, precision
from rvic.log import log_name
from rvic.utilities import find_nearest
from rvic.vars import Point

# create logger
log = logging.getLogger(log_name)


def rout(PourPoint, UHbox, FdrData, FdrAtts, RoutDict):
    """
    Make the Unit Hydrograph Grid
    """
    log.info("Starting routing program...")
    # ---------------------------------------------------------------- #
    # Unpack a few structures
    UH_t = UHbox['time']
    UH_Box = UHbox['func']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find Basin Dims and ID
    PourPoint.x = find_nearest(FdrData[RoutDict['longitude_var']], PourPoint.lon)
    PourPoint.y = find_nearest(FdrData[RoutDict['latitude_var']], PourPoint.lat)
    basin_id = FdrData[RoutDict['basin_id_var']][PourPoint.y, PourPoint.x]

    log.info('Input Latitude: %f' % PourPoint.lat)
    log.info('Input Longitude: %f' % PourPoint.lon)
    log.info('Global Basid ID: %i' % basin_id)

    y_inds, x_inds = np.nonzero(FdrData[RoutDict['basin_id_var']] == basin_id)
    y = np.arange(len(FdrData[RoutDict['latitude_var']]))
    x = np.arange(len(FdrData[RoutDict['longitude_var']]))

    x_min = min(x[x_inds])
    x_max = max(x[x_inds])+1
    y_min = min(y[y_inds])
    y_max = max(y[y_inds])+1
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create the Basin Dictionary, a subset of the FdrData
    Basin = {}
    Basin['lat'] = FdrData[RoutDict['latitude_var']][y_min:y_max]
    Basin['lon'] = FdrData[RoutDict['longitude_var']][x_min:x_max]
    Basin['Basin_ID'] = FdrData[RoutDict['basin_id_var']][y_min:y_max, x_min:x_max]
    Basin['Flow_Direction'] = FdrData[RoutDict['flow_direction_var']][y_min:y_max, x_min:x_max]
    Basin['Flow_Distance'] = FdrData[RoutDict['flow_distance_var']][y_min:y_max, x_min:x_max]
    Basin['Velocity'] = FdrData['velocity'][y_min:y_max, x_min:x_max]
    Basin['Diffusion'] = FdrData['diffusion'][y_min:y_max, x_min:x_max]

    log.info('Grid cells in subset: %i' % Basin['Velocity'].size)

    BasinPoint = Point(x=find_nearest(Basin['lon'], PourPoint.lon),
                       y=find_nearest(Basin['lat'], PourPoint.lat),
                       lon=PourPoint.lon, lat=PourPoint.lat)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Create the OutData Dictionary
    OutData = {'lat': Basin['lat'], 'lon': Basin['lon']}

    # ---------------------------------------------------------------- #
    # Determine low direction syntax
    if 'VIC' in FdrAtts[RoutDict['flow_direction_var']]:
        # VIC Directions: http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/Routing/FlowDirection.shtml
        dy = {1: -1, 2: -1, 3: 0, 4: 1, 5: 1, 6: 1, 7: 0, 8: -1}
        dx = {1: 0, 2: 1, 3: 1, 4: 1, 5: 0, 6: -1, 7: -1, 8: - 1}
        log.info('Using VIC flow directions (1-8).')
    else:
        # ARCMAP Directions: http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=flow_direction
        dy = {1: 0, 2: 1, 4: 1, 8: 1, 16: 0, 32: -1, 64: -1, 128: -1}
        dx = {1: 1, 2: 1, 4: 0, 8: -1, 16: -1, 32: -1, 64: 0, 128: 1}
        log.info('Using ARCMAP flow directions (1-128).')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Check latitude order, flip if necessary.
    if Basin['lat'][-1] > Basin['lat'][0]:
        log.info('Inputs came in upside down, flipping everything now.')
        Basin_vars = Basin.keys()
        Basin_vars.remove('lon')
        for var in Basin_vars:
            Basin[var] = np.flipud(Basin[var])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find timestep (timestep is determined from UH_BOX input file)
    InputInterval = find_TS(UH_t)
    OutData['unit_hydrogaph_dt'] = InputInterval
    Tcell = int(RoutDict['cell_flowdays']*secsPerDay/InputInterval)
    Tuh = int(RoutDict['basin_flowdays']*secsPerDay/InputInterval)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read direction grid and find to_col (to_x) and to_row (to_y)
    to_y, to_x = read_direction(Basin['Flow_Direction'], Basin['Basin_ID'],
                                dy, dx, basin_id)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Find all grid cells upstream of pour point
    Catchment, OutData['fractions'] = search_catchment(to_y, to_x, BasinPoint,
                                                       Basin['Basin_ID'],
                                                       basin_id)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make UH for each grid cell upstream of basin pour point
    # (linear routing model - Saint-Venant equation)
    UH = make_UH(InputInterval, Tcell, Catchment['y_inds'],
                 Catchment['x_inds'], Basin['Velocity'], Basin['Diffusion'],
                 Basin['Flow_Distance'])

    # ---------------------------------------------------------------- #
    # Make UH_RIVER by incrementally moving upstream comining UH functions
    UH_RIVER = make_grid_UH_river(Tuh, Tcell, UH, to_y, to_x, BasinPoint,
                                  Catchment['y_inds'], Catchment['x_inds'],
                                  Catchment['count_ds'])

    # ---------------------------------------------------------------- #
    # Make UH_S for each grid cell upstream of basin pour point
    # (combine IRFs for all grid cells in flow path)
    UH_S = make_grid_UH(Tuh, Tcell, UH_RIVER, UH_Box, to_y, to_x,
                        Catchment['y_inds'], Catchment['x_inds'],
                        Catchment['count_ds'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Agregate to output timestep
    OutData['UHgrid'], OutData['timesteps'] = adjust_uh_timestep(UH_S, Tuh,
                                                                 InputInterval,
                                                                 RoutDict['output_interval'],
                                                                 Catchment['x_inds'],
                                                                 Catchment['y_inds'])
    # ---------------------------------------------------------------- #

    return OutData

# -------------------------------------------------------------------- #
# Find the Timestep from the uh box
def find_TS(uh_t):
    """
    Determines the (InputInterval) based on the timestep given in UHfile
    """
    InputInterval = uh_t[1]-uh_t[0]
    log.info('Input Timestep = %i seconds' % InputInterval)
    return InputInterval
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Read the flow direction file
def read_direction(fdr, basin_ids, dy, dx, basin_id):
    """
    Reads the direction file and makes two grids (to_x) and (to_y).
    The input grids follow the 1-8 or 1-128 grid directions as shown below.
    val = direction  [to_y][to_x]
    """
    log.info('Reading direction input and finding target row/columns')

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
def search_catchment(to_y, to_x, BasinPoint, basin_ids, basin_id):
    """
    Find all cells upstream of pour point.  Retrun a dictionary with x_inds,
    yinds, and #of cell to downstream pour point.  All are sorted the by the
    latter. For each x,y pair, the flow path is followed until either the
    catchment outlet is encountered (if (yy==BasinPoint.y and xx==BasinPoint.x):)
    or the flowpath leads outside of grid.
    *** Does not handle wrapped coordinates. ***
    """
    log.info('Searching Catchment')

    COUNT = 0
    (len_y, len_x) = to_x.shape
    CATCH = {}
    CATCH['x_inds'] = np.array([], dtype=int)
    CATCH['y_inds'] = np.array([], dtype=int)
    CATCH['count_ds'] = np.array([], dtype=int)
    fractions = np.zeros((len_y, len_x))
    for (y, x), basin_num in np.ndenumerate(basin_ids):
        if basin_num == basin_id:
            yy, xx = y, x
            op = 0
            cells = 0
            while op == 0:
                if (yy == BasinPoint.y and xx == BasinPoint.x):
                    op = 1
                    CATCH['x_inds'] = np.append(CATCH['x_inds'], x)
                    CATCH['y_inds'] = np.append(CATCH['y_inds'], y)
                    CATCH['count_ds'] = np.append(CATCH['count_ds'], cells)
                    COUNT += 1
                    fractions[y, x] = 1.
                else:
                    yy, xx = to_y[yy, xx], to_x[yy, xx]
                    cells += 1
                    if ((xx > (len_x - 1)) or (xx < 0) or (yy > (len_y - 1)) or (yy < 0)):
                        op = -1
    log.info("Upstream grid cells from present station: %i" % COUNT)

    # ---------------------------------------------------------------- #
    # sort CATCH
    ii = np.argsort(CATCH['count_ds'])
    CATCH['count_ds'] = CATCH['count_ds'][ii]
    CATCH['x_inds'] = CATCH['x_inds'][ii]
    CATCH['y_inds'] = CATCH['y_inds'][ii]
    # ---------------------------------------------------------------- #
    return (CATCH, fractions)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make UH Grid
def make_UH(DELTA_T, T_Cell, y_inds, x_inds, velocity, diffusion, xmask):
    """
    Calculate the impulse response function for grid cells using equation 15
    from Lohmann, et al. (1996) Tellus article.  Return 3d UH grid.
    """
    log.info('Making UH for each cell')

    UH = np.zeros((T_Cell, xmask.shape[0], xmask.shape[1]))
    for (y, x) in zip(y_inds, x_inds):
        time = DELTA_T
        flag = 0
        t = 0
        green = np.zeros(T_Cell)
        while (t < T_Cell and flag == 0):
            exponent = -1*np.power(velocity[y, x]*time-xmask[y, x], 2)/(4*diffusion[y, x]*time)
            if exponent > np.log(precision):
                green[t] = xmask[y, x]/(2*time*np.sqrt(np.pi*time*diffusion[y, x]))*np.exp(exponent)
                t += 1
                time = time+DELTA_T
            else:
                flag = 1
        tot = np.sum(green)
        if tot > 0.:
            UH[:, y, x] = green[:] / tot
    return UH
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Make UH river
def make_grid_UH_river(T_UH, T_Cell, UH, to_y, to_x, BasinPoint, y_inds,
                       x_inds, count_ds):
    """
    Calculate impulse response function for river routing.  Starts at
    downstream point incrementally moves upstream.
    """
    log.info("Making UH_RIVER grid.... It takes a while...")
    y_ind = BasinPoint.y
    x_ind = BasinPoint.x

    UH_RIVER = np.zeros((T_UH, UH.shape[1], UH.shape[2]))
    for (y, x, d) in zip(y_inds, x_inds, count_ds):
        if d > 0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            IRF_temp = np.zeros(T_UH+T_Cell)
            active_timesteps = np.nonzero(UH_RIVER[:, yy, xx] > precision)[0]
            for t in active_timesteps:
                for l in xrange(T_Cell):
                    IRF_temp[t+l] = IRF_temp[t + l] + UH[l, y, x] * UH_RIVER[t, yy, xx]
            tot = np.sum(IRF_temp[:T_UH])
            if tot > 0:
                UH_RIVER[:, y, x] = IRF_temp[:T_UH] / tot
        elif d == 0:
            UH_RIVER[:T_Cell, y_ind, x_ind] = UH[:, y_ind, x_ind]

    return UH_RIVER
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make grid UH
def make_grid_UH(T_UH, T_Cell, UH_RIVER, UH_BOX, to_y, to_x, y_inds, x_inds,
                 count_ds):
    """
    Combines the UH_BOX with downstream cell UH_RIVER.  Cell [0] is given the
    UH_Box without river routing
    """
    log.info("Making UH_S grid")

    UH_S = np.zeros((T_UH, UH_RIVER.shape[1], UH_RIVER.shape[2]))
    for (y, x, d) in zip(y_inds, x_inds, count_ds):
        IRF_temp = np.zeros(T_UH+T_Cell)
        if d > 0:
            yy = to_y[y, x]
            xx = to_x[y, x]
            active_timesteps = np.nonzero(UH_RIVER[:, yy, xx] > precision)[0]
            for t in active_timesteps:
                for l in xrange(len(UH_BOX)):
                    IRF_temp[t + l] = IRF_temp[t + l] + UH_BOX[l] * UH_RIVER[t, yy, xx]
            tot = np.sum(IRF_temp[:T_UH])
            if tot > 0:
                UH_S[:, y, x] = IRF_temp[:T_UH] / tot
        else:
            IRF_temp[:len(UH_BOX)] = UH_BOX[:]
            UH_S[:, y, x] = IRF_temp[:T_UH]
    return UH_S
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Adjust the timestep
def adjust_uh_timestep(UH_S, T_UH, INPUT_INTERVAL, OUTPUT_INTERVAL, x_inds, y_inds):
    """
    Aggregates to timestep (OUTPUT_INTERVAL).  OUTPUT_INTERVAL must be a
    multiple of INPUT_INTERVAL.  This function is not setup to disaggregate
    the Unit Hydrographs to a OUTPUT_INTERVAL<INPUT_INTERVAL.
    """
    if OUTPUT_INTERVAL == INPUT_INTERVAL:
        log.info('No need to aggregate (OUTPUT_INTERVAL = INPUT_INTERVAL) \
                      Skipping the adjust_uh_timestep step')
        UH_out = UH_S
        ts_new = np.arange(T_UH)
    elif np.remainder(OUTPUT_INTERVAL, INPUT_INTERVAL) == 0:
        log.info('Aggregating to %i from %i seconds' % (OUTPUT_INTERVAL, INPUT_INTERVAL))
        fac = int(OUTPUT_INTERVAL/INPUT_INTERVAL)
        T_UH_out = int(T_UH/fac)
        ts_new = np.arange(T_UH_out)
        UH_out = np.zeros((T_UH_out, UH_S.shape[1], UH_S.shape[2]))
        for (y, x) in zip(y_inds, x_inds):
            for t in xrange(T_UH_out):
                UH_out[t, y, x] = np.sum(UH_S[t*fac:t*fac+fac, y, x])
    elif np.remainder(INPUT_INTERVAL, OUTPUT_INTERVAL):
        log.info('Interpolating unit hydrograph from input_interval: %i to output_interval: %i' % (INPUT_INTERVAL, OUTPUT_INTERVAL))
        fac = int(INPUT_INTERVAL / OUTPUT_INTERVAL)
        T_UH_out = int(T_UH * fac)
        UH_out = np.zeros((T_UH_out, UH_S.shape[1], UH_S.shape[2]))
        ts_orig = np.linspace(0, T_UH, T_UH)
        ts_new = np.linspace(0, T_UH, T_UH_out)
        for (y, x) in zip(y_inds, x_inds):
            f = interp1d(ts_orig, UH_S[:, y, x])
            UH_out[:, y, x] = f(ts_new)
    return UH_out, ts_new
# -------------------------------------------------------------------- #

