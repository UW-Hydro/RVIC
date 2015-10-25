# -*- coding: utf-8 -*-
'''
RVIC parameter file development driver
'''
import os
import numpy as np
import pandas as pd
from logging import getLogger
from .core.log import init_logger, close_logger, LOG_NAME
from .core.multi_proc import error
from .core.utilities import make_directories, copy_inputs, strip_invalid_char
from .core.utilities import read_netcdf, tar_inputs, latlon2yx
from .core.utilities import check_ncvars, clean_file, read_domain
from .core.utilities import search_for_channel
from .core.aggregate import make_agg_pairs, aggregate
from .core.make_uh import rout
from .core.share import NcGlobals
from .core.write import write_agg_netcdf
from .core.variables import Point
from .core.param_file import finish_params
from .core.config import read_config
from .core.pycompat import OrderedDict, iteritems, pyrange, basestring

try:
    from .core.remap import remap
    remap_available = True
except ImportError:
    remap_available = False

# global multiprocessing results container
results = {}


# -------------------------------------------------------------------- #
# Top level driver
def parameters(config_file, numofproc=1):
    '''
    Top level function for RVIC parameter generation function.

    Parameters
    ----------
    config_file : str
        Path to RVIC parameters configuration file.
    numofproc : int
        Number of processors to use when developing RVIC parameters.
    '''

    # ---------------------------------------------------------------- #
    # Initilize
    uh_box, fdr_data, fdr_vatts, dom_data, \
        outlets, config_dict, directories = gen_uh_init(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get main logger
    log = getLogger(LOG_NAME)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run
    if numofproc > 1:
        from multiprocessing import Pool
        pool = Pool(processes=numofproc)
        status = []

        for i, (key, outlet) in enumerate(iteritems(outlets)):
            log.info('On Outlet #%s of %s', i + 1, len(outlets))
            stat = pool.apply_async(gen_uh_run,
                                    (uh_box, fdr_data, fdr_vatts,
                                     dom_data, outlet, config_dict,
                                     directories),
                                    callback=store_result,
                                    error_callback=error)
            # Store the result
            status.append(stat)

        # Close the pool
        pool.close()

        # Check that everything worked
        [stat.get() for stat in status]

        pool.join()

        outlets = OrderedDict(sorted(results.items(), reverse=True))
    else:
        for i, (key, outlet) in enumerate(iteritems(outlets)):
            outlet = gen_uh_run(uh_box, fdr_data, fdr_vatts, dom_data, outlet,
                                config_dict, directories)

    if not outlets:
        raise ValueError('outlets in parameters are empty')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Finally, make the parameter file
    gen_uh_final(outlets, dom_data, config_dict, directories)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


def gen_uh_init(config_file):
    '''Initialize RVIC parameters script.

    This function:
        - Reads the configuration file
        - Sets up the RVIC case directories
        - Copies all input files to the case directory
        - Initializes the logging
        - Reads the pour-points, uh-box, FDR, and domain files
        - Aggregates pour points into outlet grid cells

    Parameters
    ----------
    config_file : str
        Path to RVIC parameters configuration file.

    Returns
    ----------
    uh_box : numpy.ndarray
        UH-box array.
    fdr_data : dict
        Dictionary of arrays of flow direction, velocity, diffusion, etc.
        This dictionary includes all the variables from the FDR netCDF file.
    fdr_vatts : dict
        Dictionary of attributes from the FDR netCDF file.
    dom_data : dict
        Dictionary of arrays of mask, fraction, lats, lons, etc.
        This dictionary includes all the variables from the domain netCDF file.
    outlets : dict
        Dictionary of outlet `Point` objects.
    config_dict : dict
        Dictionary of values from the configuration file.
    directories : dict
        Dictionary of directories created by this function.
    '''

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Import optional modules
    if config_dict['OPTIONS']['REMAP'] and not remap_available:
        raise ValueError('Problem importing remap module '
                         'check to make sure cdo.py is available)')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structures
    directories = make_directories(config_dict['OPTIONS']['CASE_DIR'],
                                   ['plots', 'logs', 'params', 'inputs'])
    directories.update(make_directories(config_dict['OPTIONS']['TEMP_DIR'],
                                        ['aggregated', 'remapped']))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # copy inputs to $case_dir/inputs and update configuration
    config_dict = copy_inputs(config_file, directories['inputs'])
    options = config_dict['OPTIONS']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Start Logging
    log = init_logger(directories['logs'], options['LOG_LEVEL'],
                      options['VERBOSE'])

    for direc in directories:
        log.info('%s directory is %s', direc, directories[direc])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Pour Points files
    try:
        pour_points = pd.read_csv(config_dict['POUR_POINTS']['FILE_NAME'],
                                  comment='#')
        log.info('Opened Pour Points File: %s',
                 config_dict['POUR_POINTS']['FILE_NAME'])
        if not (all(x in list(pour_points.keys()) for x in ['lons', 'lats']) or
                all(x in list(pour_points.keys()) for x in ['x', 'y'])):
            raise ValueError('Pour Points File must include '
                             'variables (lons, lats) or (x, y)')
        if 'names' in pour_points:
            pour_points.fillna(inplace=True, value='unknown')
            for i, name in enumerate(pour_points.names):
                pour_points.ix[i, 'names'] = strip_invalid_char(name)

        pour_points.drop_duplicates(inplace=True)
        pour_points.dropna()
    except Exception as e:
        log.error('Error opening pour points file: %s',
                  config_dict['POUR_POINTS']['FILE_NAME'])
        log.exception(e)
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read uh box file
    uh_file = config_dict['UH_BOX']['FILE_NAME']
    uh_header = int(config_dict['UH_BOX']['HEADER_LINES'])
    uh_box = {}
    try:
        uh_box['time'], uh_box['func'] = np.genfromtxt(uh_file,
                                                       skip_header=uh_header,
                                                       delimiter=',',
                                                       unpack=True)
        log.info('Opened UHbox File: %s',
                 config_dict['UH_BOX']['FILE_NAME'])
    except:
        log.exception('Error opening uh_box file: %s',
                      config_dict['POUR_POINTS']['FILE_NAME'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read FDR file
    fdr_file = config_dict['ROUTING']['FILE_NAME']
    fdr_var = config_dict['ROUTING']['FLOW_DIRECTION_VAR']
    fdr_lat = config_dict['ROUTING']['LATITUDE_VAR']
    fdr_lon = config_dict['ROUTING']['LONGITUDE_VAR']
    fdr_vel = config_dict['ROUTING']['VELOCITY']
    fdr_dif = config_dict['ROUTING']['DIFFUSION']
    try:
        fdr_data, fdr_vatts, _ = read_netcdf(fdr_file)
        fdr_shape = fdr_data[fdr_var].shape

        # ---------------------------------------------------------------- #
        # Check latitude order, flip if necessary.
        if fdr_data[fdr_lat][-1] > fdr_data[fdr_lat][0]:
            log.debug('Flow Direction inputs came in upside down, flipping '
                      'everything now.')

            remove_vars = []

            for var, data in iteritems(fdr_data):
                log.debug('flipping %s', var)
                if data.ndim >= 1 and var != fdr_lon:
                    fdr_data[var] = np.flipud(data)
                elif data.ndim == 0:
                    remove_vars.append(var)

            if remove_vars:
                for var in remove_vars:
                    del fdr_data[var]
        # ---------------------------------------------------------------- #

        # ---------------------------------------------------------------- #
        # Add velocity and/or diffusion grids if not present yet
        if not isinstance(fdr_vel, basestring):
            fdr_data['velocity'] = \
                np.zeros(fdr_shape, dtype=np.float64) + fdr_vel
            config_dict['ROUTING']['VELOCITY'] = 'velocity'
            log.info('Added velocity grid to fdr_data')
        if not isinstance(fdr_dif, basestring):
            fdr_data['diffusion'] = \
                np.zeros(fdr_shape, dtype=np.float64) + fdr_dif
            config_dict['ROUTING']['DIFFUSION'] = 'diffusion'
            log.info('Added diffusion grid to fdr_data')
        if ('SOURCE_AREA_VAR' not in config_dict['ROUTING'] or
                config_dict['ROUTING']['SOURCE_AREA_VAR'] not in fdr_data):
            log.warning('Upstream `SOURCE_AREA` was not provided, output '
                        'source area will be zero.')
            config_dict['ROUTING']['SOURCE_AREA_VAR'] = 'src_area'
            fdr_data[config_dict['ROUTING']['SOURCE_AREA_VAR']] = \
                np.zeros(fdr_shape, dtype=np.float64)
        # ---------------------------------------------------------------- #

        # ---------------------------------------------------------------- #
        fdr_data['resolution'] = np.abs(fdr_data[fdr_lon][1] -
                                        fdr_data[fdr_lon][0])
        check_ncvars(config_dict['ROUTING'], list(fdr_data.keys()))
        # ---------------------------------------------------------------- #

        log.info('Opened FDR File: %s', fdr_file)
    except:
        log.exception('Error opening FDR file')
        raise
        # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read domain file
    domain = config_dict['DOMAIN']
    dom_data = read_domain(domain)[0]
    log.info('Opened Domain File: %s', domain['FILE_NAME'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # If remap is False, domain coordinates needs to be in the fdr coordinates
    # We can move the unit hydrographs to the domain grid later
    if options['AGGREGATE'] and not options['REMAP']:
        log.error('RVIC parameter generation requires REMAP option to be True'
                  ' if AGGREGATE is True')
        raise ValueError('Invalid option')

    # If remap is False, then the resolution needs to match the routing data
    if not options['REMAP']:
        domain_res = np.abs(dom_data[domain['LONGITUDE_VAR']][0, 1] -
                            dom_data[domain['LONGITUDE_VAR']][0, 0])
        if not np.isclose(fdr_data['resolution'], domain_res):
            log.error('routing grid resolution: %s', fdr_data['resolution'])
            log.error('domain grid resolution: %s', domain_res)
            raise ValueError('If remap is false, domain and routing grid '
                             'resolutions must match.')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Group pour points (if aggregate)
    if options['AGGREGATE']:
        outlets = make_agg_pairs(pour_points, dom_data, fdr_data, config_dict)

        log.info('Finished making agg pairs of '
                 'pour points and outlet grid cells')

    else:
        outlets = OrderedDict()
        if all(x in list(pour_points.keys()) for x in ['x', 'y',
                                                       'lons', 'lats']):
            lats = pour_points['lats'].values
            lons = pour_points['lons'].values
            routys = pour_points['y'].values
            routxs = pour_points['x'].values
        elif all(x in list(pour_points.keys()) for x in ['x', 'y']):
            # use x and y (assume from routing inputs grid)
            # find lons and lats from xs and ys
            routys = pour_points['y'].values
            routxs = pour_points['x'].values
            lats = fdr_data[fdr_lat][routys]
            lons = fdr_data[fdr_lon][routxs]
        else:
            # use and lats to find xs and ys
            lats = pour_points['lats'].values
            lons = pour_points['lons'].values

            # find x and y on routing grid
            routys, routxs = latlon2yx(plats=lats,
                                       plons=lons,
                                       glats=fdr_data[fdr_lat],
                                       glons=fdr_data[fdr_lon])

        if options['SEARCH_FOR_CHANNEL']:
            routys, routxs = search_for_channel(
                fdr_data[config_dict['ROUTING']['SOURCE_AREA_VAR']],
                routys, routxs, tol=10, search=5)

            # update lats and lons
            lats = fdr_data[fdr_lat][routys]
            lons = fdr_data[fdr_lon][routxs]

        # Find location on domain grid
        domys, domxs = latlon2yx(plats=lats,
                                 plons=lons,
                                 glats=dom_data[domain['LATITUDE_VAR']],
                                 glons=dom_data[domain['LONGITUDE_VAR']])

        for i in pyrange(len(lats)):
            if 'names' in list(pour_points.keys()):
                name = pour_points['names'].values[i]
                name = name.replace("'", '').replace(' ', '_')
            else:
                # fill name filed with p-outlet_num
                name = 'p-{0}'.format(i)

            outlets[i] = Point(lat=lats[i],
                               lon=lons[i],
                               domx=domxs[i],
                               domy=domys[i],
                               routx=routxs[i],
                               routy=routys[i],
                               name=name,
                               cell_id=dom_data['cell_ids'][domys[i], domxs[i]])

            outlets[i].pour_points = [outlets[i]]
    # ---------------------------------------------------------------- #

    log.debug(outlets)
    log.info('Finished with gen_uh_init')
    log.info('-------------------------------------------------------------\n')

    return (uh_box, fdr_data, fdr_vatts, dom_data, outlets,
            config_dict, directories)
# -------------------------------------------------------------------- #


def gen_uh_run(uh_box, fdr_data, fdr_vatts, dom_data, outlet, config_dict,
               directories):
    '''
    Develop unit hydrographs for one outlet `Point`.

    Parameters
    ----------
    uh_box : numpy.ndarray
        UH-box array.
    fdr_data : dict
        Dictionary of arrays of flow direction, velocity, diffusion, etc.
        This dictionary includes all the variables from the FDR netCDF file.
    fdr_vatts : dict
        Dictionary of attributes from the FDR netCDF file.
    dom_data : dict
        Dictionary of arrays of mask, fraction, lats, lons, etc.
        This dictionary includes all the variables from the domain netCDF file.
    outlet : Point
        Outlet Point
    config_dict : dict
        Dictionary of values from the configuration file.
    directories : dict
        Dictionary of directories created by gen_uh_init.

    Returns
    ----------
    outlet: Point
        Point object with unit hydrographs.
    '''
    log = getLogger(LOG_NAME)

    log.info('Running outlet cell id %s', outlet.cell_id)

    agg_data = {}
    domain = config_dict['DOMAIN']
    dom_lat = domain['LATITUDE_VAR']
    dom_lon = domain['LONGITUDE_VAR']
    dom_mask = domain['LAND_MASK_VAR']
    options = config_dict['OPTIONS']

    # ------------------------------------------------------------ #
    # netCDF variable options
    ncvaropts = {}
    if 'NETCDF_ZLIB' in options:
        ncvaropts['zlib'] = options['NETCDF_ZLIB']
    if 'NETCDF_COMPLEVEL' in options:
        ncvaropts['complevel'] = options['NETCDF_COMPLEVEL']
    if 'NETCDF_SIGFIGS' in options:
        ncvaropts['least_significant_digit'] = options['NETCDF_SIGFIGS']
    # ------------------------------------------------------------ #

    # ------------------------------------------------------------ #
    # Loop over pour points
    n_pour_points = len(outlet.pour_points)
    for j, pour_point in enumerate(outlet.pour_points):

        log.info('On pour_point #%s of %s', j + 1, n_pour_points)

        # -------------------------------------------------------- #
        # Make the Unit Hydrograph Grid
        rout_data = rout(pour_point, uh_box, fdr_data, fdr_vatts,
                         config_dict['ROUTING'])

        log.debug('Done routing to pour_point')
        log.debug('rout_data: %s, %s', rout_data['unit_hydrograph'].min(),
                  rout_data['unit_hydrograph'].max())
        log.debug('rout_data sum: %s, %s', rout_data['unit_hydrograph'].sum(),
                  rout_data['fraction'].sum())

        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # aggregate
        if options['AGGREGATE']:
            if j != len(outlet.pour_points) - 1:
                agg_data = aggregate(rout_data, agg_data,
                                     res=fdr_data['resolution'])
            else:
                agg_data = aggregate(rout_data, agg_data,
                                     res=fdr_data['resolution'],
                                     pad=options['AGG_PAD'],
                                     maskandnorm=True)

                log.debug('agg_data: %s, %s',
                          agg_data['unit_hydrograph'].min(),
                          agg_data['unit_hydrograph'].max())
        else:
            agg_data = rout_data
        # -------------------------------------------------------- #

    # ------------------------------------------------------------ #
    # write temporary file #1
    if options['REMAP']:
        glob_atts = NcGlobals(
            title='RVIC Unit Hydrograph Grid File',
            RvicPourPointsFile=os.path.split(
                config_dict['POUR_POINTS']['FILE_NAME'])[1],
            RvicUHFile=os.path.split(config_dict['UH_BOX']['FILE_NAME'])[1],
            RvicFdrFile=os.path.split(config_dict['ROUTING']['FILE_NAME'])[1],
            RvicDomainFile=os.path.split(domain['FILE_NAME'])[1])

        temp_file_1 = os.path.join(
            directories['aggregated'],
            'aggUH_{0}.nc'.format(outlet.name.replace(' ', '_')))

        write_agg_netcdf(temp_file_1, agg_data, glob_atts,
                         options['NETCDF_FORMAT'], **ncvaropts)

        # -------------------------------------------------------- #
        # Remap temporary file #1 to temporary file #2
        temp_file_2 = os.path.join(
            directories['remapped'],
            'remapUH_{0}.nc'.format(outlet.name.replace(' ', '_')))

        remap(domain['FILE_NAME'], temp_file_1, temp_file_2)

        # -------------------------------------------------------- #
        # Read temporary file #2
        final_data = read_netcdf(
            temp_file_2, variables=['unit_hydrograph', 'fraction', dom_lat])[0]

        # -------------------------------------------------------- #
        # Check latitude order, flip if necessary.
        if final_data[dom_lat].ndim == 1:
            if final_data[dom_lat][-1] > final_data[dom_lat][0]:
                var_list = list(final_data.keys())

                log.debug('Remapped inputs came in upside down, flipping %s',
                          ', '.join(var_list))
                # flip lattiutude and fraction along y axis (axis 0)
                final_data[dom_lat] = final_data[dom_lat][::-1]
                final_data['fraction'] = final_data['fraction'][::-1, :]
                # flip unit hydrograph along y axis (axis 1)
                final_data['unit_hydrograph'] = \
                    final_data['unit_hydrograph'][:, ::-1, :]
            assert dom_data['cord_lats'][0] == final_data[dom_lat][0]
        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # Clean temporary file #2 (if applicable)
        if config_dict['OPTIONS']['CLEAN']:
            clean_file(temp_file_1)
            clean_file(temp_file_2)
        # -------------------------------------------------------- #

    else:
        # -------------------------------------------------------- #
        # Put the agg data back onto the original grid
        uh_shape = (agg_data['unit_hydrograph'].shape[0], ) + \
            dom_data[dom_mask].shape
        final_data = {}
        final_data['unit_hydrograph'] = np.zeros(uh_shape, dtype=np.float64)
        final_data['fraction'] = np.zeros(dom_data[dom_mask].shape,
                                          dtype=np.float64)

        bys, bxs = np.nonzero(agg_data['fraction'])

        ys = bys + outlet.domy - agg_data['basiny']
        xs = bxs + outlet.domx - agg_data['basinx']

        if (ys < 0).any() or (xs < 0).any():
            raise ValueError('Negative indicies found when mapping '
                             '`non-remapped` rout_data to domain grid.')

        final_data['unit_hydrograph'][:, ys, xs] = \
            agg_data['unit_hydrograph'][:, bys, bxs]
        final_data['fraction'][ys, xs] = agg_data['fraction'][bys, bxs]
        # -------------------------------------------------------- #
    # ------------------------------------------------------------ #

    # ------------------------------------------------------------ #
    # Add to 'adjust fractions structure'
    y, x = np.nonzero((final_data['fraction'] > 0.) *
                      (dom_data[dom_mask] > np.finfo(np.float).resolution))

    # From final data
    outlet.time = np.arange(final_data['unit_hydrograph'].shape[0])
    outlet.fractions = final_data['fraction'][y, x]
    outlet.unit_hydrograph = final_data['unit_hydrograph'][:, y, x]

    # From domain data
    outlet.lon_source = dom_data[dom_lon][y, x]
    outlet.lat_source = dom_data[dom_lat][y, x]
    outlet.cell_id_source = dom_data['cell_ids'][y, x]
    outlet.x_source = x
    outlet.y_source = y
    # ---------------------------------------------------------------- #
    return outlet
# -------------------------------------------------------------------- #


def gen_uh_final(outlets, dom_data, config_dict, directories):
    '''
    Make the RVIC Parameter File

    Parameters
    ----------
    outlets : dict
        Dictionary of outlet `Point` objects.
    dom_data : dict
        Dictionary of arrays of mask, fraction, lats, lons, etc.
        This dictionary includes all the variables from the domain netCDF file.
    config_dict : dict
        Dictionary of values from the configuration file.
    directories : dict
        Dictionary of directories created by gen_uh_init.
    '''
    log = getLogger(LOG_NAME)

    log.info('In gen_uh_final')

    if not len(outlets) > 0:
        raise ValueError('outlets in gen_uh_final are empty')

    # ---------------------------------------------------------------- #
    # Write the parameter file
    param_file, today = finish_params(outlets, dom_data, config_dict,
                                      directories)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    inputs_tar = tar_inputs(directories['inputs'], suffix=today)
    log_tar = tar_inputs(log.filename)

    log.info('Done with RvicGenParam.')
    log.info('Location of Inputs: %s', inputs_tar)
    log.info('Location of Log: %s', log_tar)
    log.info('Location of Parmeter File %s', param_file)

    close_logger()
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# store_result helper function
def store_result(result):
    '''
    Store values returned by a multiprocessing.pool member.

    This is called whenever foo_pool(i) returns a result.
    result_list is modified only by the main process, not the pool workers.

    Parameters
    ----------
    result : object
        Result to append to the global `results` list

    Globals
    ----------
    results : dict
        Global results container for multiprocessing results to be appended to.
    '''
    results[result.cell_id] = result
# -------------------------------------------------------------------- #
