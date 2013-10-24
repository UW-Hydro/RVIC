#!/usr/bin/env python2.7
"""
RVIC parameter file development driver
"""
import os
import numpy as np
import argparse
from collections import OrderedDict
from logging import getLogger
from rvic.log import init_logger, LOG_NAME
from rvic.mpi import LoggingPool
from rvic.utilities import read_config, make_directories, copy_inputs, read_netcdf, tar_inputs
from rvic.utilities import check_ncvars, remap, clean_file, read_domain
from rvic.aggregate import make_agg_pairs, aggregate
from rvic.make_uh import rout
from rvic.share import NcGlobals
from rvic.write import write_agg_netcdf
from rvic.variables import Point
from rvic.param_file import finish_params


# -------------------------------------------------------------------- #
# Top level driver
def main():

    # ---------------------------------------------------------------- #
    # Read command Line
    config_file, numofproc = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initilize
    uh_box, fdr_data, fdr_vatts, dom_data, outlets, config_dict, directories = gen_uh_init(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get main logger
    log = getLogger(LOG_NAME)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run
    if numofproc >1:
        pool = LoggingPool(processes=numofproc)

        for i, (cell_id, outlet) in enumerate(outlets.iteritems()):
            log.info('On Outlet #%i of %i' %(i+1, len(outlets)))
            pool.apply_async(gen_uh_run,
                             args = (uh_box, fdr_data, fdr_vatts, dom_data, outlet, config_dict, directories),
                             callback = store_result)
        pool.close()
        pool.join()

        outlets = OrderedDict(sorted(results.items(), key=lambda t: t[0]))
    else:
        for i, (cell_id, outlet) in enumerate(outlets.iteritems()):
            outlets[cell_id] = gen_uh_run(uh_box, fdr_data, fdr_vatts, dom_data, outlet, config_dict, directories)

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Finally, make the parameter file
    gen_uh_final(outlets, dom_data, config_dict, directories)
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Initialize the Genuh Program
def gen_uh_init(config_file):
    """Initialize RVIC parameter """

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
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
    log = init_logger(directories['logs'], options['LOG_LEVEL'], options['VERBOSE'])

    for direc in directories:
        log.info('%s directory is %s' % (direc, directories[direc]))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Pour Points files
    pour_points = {}
    try:
        pour_points['lons'], pour_points['lats'] = np.genfromtxt(config_dict['POUR_POINTS']['FILE_NAME'],
                                                                 skip_header=int(config_dict['POUR_POINTS']['HEADER_LINES']),
                                                                 delimiter=',', unpack=True)
        log.info('Opened Pour Points File: %s' % config_dict['POUR_POINTS']['FILE_NAME'])
    except:
        log.exception('Error opening pour points file: %s' % config_dict['POUR_POINTS']['FILE_NAME'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read uh box file
    uh_box = {}
    try:
        uh_box['time'], uh_box['func'] = np.genfromtxt(config_dict['UH_BOX']['FILE_NAME'],
                                                       skip_header=int(config_dict['UH_BOX']['HEADER_LINES']),
                                                       delimiter=',', unpack=True)
        log.info('Opened UHbox File: %s' % config_dict['UH_BOX']['FILE_NAME'])
    except:
        log.exception('Error opening uh_box file: %s' % config_dict['POUR_POINTS']['FILE_NAME'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read FDR file
    try:
        fdr_data, fdr_vatts, fdr_gatts = read_netcdf(config_dict['ROUTING']['FILE_NAME'])
        fdr_shape = fdr_data[config_dict['ROUTING']['FLOW_DIRECTION_VAR']].shape

        # ---------------------------------------------------------------- #
        # Check latitude order, flip if necessary.
        if fdr_data[config_dict['ROUTING']['LATITUDE_VAR']][-1] > fdr_data[config_dict['ROUTING']['LATITUDE_VAR']][0]:
            log.debug('Inputs came in upside down, flipping everything now.')
            var_list = fdr_data.keys()
            var_list.remove(config_dict['ROUTING']['LONGITUDE_VAR'])
            for var in var_list:
                fdr_data[var] = np.flipud(fdr_data[var])
        # ---------------------------------------------------------------- #

        # Add velocity and/or diffusion grids if not present yet
        if not type(config_dict['ROUTING']['VELOCITY']) == str:
            fdr_data['velocity'] = np.zeros(fdr_shape) + config_dict['ROUTING']['VELOCITY']
            config_dict['ROUTING']['VELOCITY'] = 'velocity'
            log.info('Added velocity grid to fdr_data')
        if not type(config_dict['ROUTING']['DIFFUSION']) == str:
            fdr_data['diffusion'] = np.zeros(fdr_shape) + config_dict['ROUTING']['DIFFUSION']
            config_dict['ROUTING']['DIFFUSION'] = 'diffusion'
            log.info('Added diffusion grid to fdr_data')

        fdr_data['resolution'] = np.abs(fdr_data[config_dict['ROUTING']['LONGITUDE_VAR']][1] -
                                       fdr_data[config_dict['ROUTING']['LONGITUDE_VAR']][0])

        check_ncvars(config_dict['ROUTING'], fdr_data.keys())

        log.info('Opened FDR File: %s' % config_dict['ROUTING']['FILE_NAME'])
    except:
        log.exception('Error opening FDR file')
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read domain file (if applicable)
    if options['REMAP']:
        dom_data, DomVats, DomGats = read_domain(config_dict['DOMAIN'])
        log.info('Opened Domain File: %s' % config_dict['DOMAIN']['FILE_NAME'])
    else:
        dom_data = {}
    # ---------------------------------------------------------------- #



    # ---------------------------------------------------------------- #
    # Group pour points (if aggregate)
    if options['AGGREGATE']:
        outlets = make_agg_pairs(pour_points['lons'], pour_points['lats'],
                                 dom_data[config_dict['DOMAIN']['LONGITUDE_VAR']],
                                 dom_data[config_dict['DOMAIN']['LATITUDE_VAR']],
                                 dom_data['cell_ids'],
                                 fdr_data[config_dict['ROUTING']['LONGITUDE_VAR']],
                                 fdr_data[config_dict['ROUTING']['LATITUDE_VAR']],
                                 fdr_data[config_dict['ROUTING']['SOURCE_AREA_VAR']],
                                 agg_type='agg')

        log.info('Finished making agg pairs of pour points and outlet grid cells')

    else:
        outlets = {}
        for i, (lon, lat) in enumerate(zip(pour_points['lons'], pour_points['lats'])):
            outlets[i] = Point(lat=lat, lon=lon)
            outlets[i].pour_points = [Point(lat=lat, lon=lon)]
    # ---------------------------------------------------------------- #

    log.info('Finished with gen_uh_init')
    log.info('--------------------------------------------------------------------\n')

    return uh_box, fdr_data, fdr_vatts, dom_data, outlets, config_dict, directories
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
#
def gen_uh_run(uh_box, fdr_data, fdr_vatts, dom_data, outlet, config_dict, directories):
    """
    Run Genuh_run
    """
    log = getLogger(LOG_NAME)

    log.info('Running outlet cell id  %i', outlet.cell_id)

    agg_data = {}
    # ------------------------------------------------------------ #
    # Loop over pour points
    for j, pour_point in enumerate(outlet.pour_points):

        log.info('On pour_point #%i of %i' %(j+1, len(outlet.pour_points)))

        # -------------------------------------------------------- #
        # Make the Unit Hydrograph Grid
        rout_data = rout(pour_point, uh_box, fdr_data, fdr_vatts,
                        config_dict['ROUTING'])

        log.debug('Done routing to pour_point')
        log.debug('rout_data: %f, %f' % (rout_data['unit_hydrograph'].min(), rout_data['unit_hydrograph'].max()))

        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # aggregate
        if config_dict['OPTIONS']['AGGREGATE']:
            if j != len(outlet.pour_points)-1:
                agg_data = aggregate(rout_data, agg_data, res=fdr_data['resolution'])
            else:
                agg_data = aggregate(rout_data, agg_data, res=fdr_data['resolution'],
                                     pad=config_dict['OPTIONS']['AGG_PAD'], maskandnorm=True)

                log.debug('agg_data: %f, %f' % (agg_data['unit_hydrograph'].min(), agg_data['unit_hydrograph'].max()))
        elif config_dict['OPTIONS']['REMAP']:
            agg_data = rout_data
        else:
            remap_data = rout_data
        # -------------------------------------------------------- #

    # ------------------------------------------------------------ #
    # write temporary file #1
    if config_dict['OPTIONS']['REMAP']:
        glob_atts = NcGlobals(title='RVIC Unit Hydrograph Grid File',
                              RvicPourPointsFile=os.path.split(config_dict['POUR_POINTS']['FILE_NAME'])[1],
                              RvicUHFile=os.path.split(config_dict['UH_BOX']['FILE_NAME'])[1],
                              RvicFdrFile=os.path.split(config_dict['ROUTING']['FILE_NAME'])[1],
                              RvicDomainFile=os.path.split(config_dict['DOMAIN']['FILE_NAME'])[1])

        temp_file_1 = os.path.join(directories['aggregated'], 'aggUH_%i.nc' % outlet.cell_id)

        write_agg_netcdf(temp_file_1, agg_data, glob_atts,
                         config_dict['OPTIONS']['NETCDF_FORMAT'])
    elif config_dict['OPTIONS']['AGGREGATE']:
        remap_data = agg_data
    else:
        pass
    # ------------------------------------------------------------ #

    # ------------------------------------------------------------ #
    # Remap temporary file #1 to temporary file #2
    if config_dict['OPTIONS']['REMAP']:

        temp_file_2 = os.path.join(directories['remapped'], 'remapUH_%i.nc' % outlet.cell_id)

        remap(config_dict['DOMAIN']['FILE_NAME'], temp_file_1, temp_file_2)

        # -------------------------------------------------------- #
        # Clean temporary file #1 (if applicable)
        if config_dict['OPTIONS']['CLEAN']:
            clean_file(temp_file_1)
        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # Read temporary file #2
        remap_data, remap_vatts, remap_gatts = read_netcdf(temp_file_2,
                                                           variables=['unit_hydrograph', 'fraction'])
        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # Clean temporary file #2 (if applicable)
        if config_dict['OPTIONS']['CLEAN']:
            clean_file(temp_file_2)
        # -------------------------------------------------------- #

    else:
        remap_data = agg_data
    # ------------------------------------------------------------ #

    # ------------------------------------------------------------ #
    # Add to adjust fractions Structure
    if config_dict['OPTIONS']['REMAP']:
        y, x = np.nonzero((remap_data['fraction'] > 0.0) * (dom_data[config_dict['DOMAIN']['LAND_MASK_VAR']] == 1))

        outlet.fractions = remap_data['fraction'][y, x]
        outlet.unit_hydrograph = remap_data['unit_hydrograph'][:, y, x]
        outlet.time = np.arange(remap_data['unit_hydrograph'].shape[0])
        outlet.lon_source = dom_data[config_dict['DOMAIN']['LONGITUDE_VAR']][y, x]
        outlet.lat_source = dom_data[config_dict['DOMAIN']['LATITUDE_VAR']][y, x]
        outlet.cell_id_source = dom_data['cell_ids'][y, x]
        outlet.x_source = x
        outlet.y_source = y
    # ---------------------------------------------------------------- #
    return outlet
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def gen_uh_final(outlets, dom_data, config_dict, directories):
    """
    Make the RVIC Parameter File
    """
    log = getLogger(LOG_NAME)

    log.info('In gen_uh_final')

    # ---------------------------------------------------------------- #
    # Write the parameter file
    param_file, today = finish_params(outlets, dom_data, config_dict, directories)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    inputs_tar = tar_inputs(directories['inputs'], suffix=today)
    log_tar = tar_inputs(log.filename)

    log.info('Done with RvicGenParam.')
    log.info('Location of Inputs: %s' % inputs_tar)
    log.info('Location of Log: %s' % log_tar)
    log.info('Location of Parmeter File %s' % param_file)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def process_command_line():
    """
    Get the path to the config_file
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate RVIC parameter files.')
    parser.add_argument("config_file", type=str,
                        help="Input configuration file")
    parser.add_argument("-np", "--numofproc", type=int,
                        help="Number of processors used to run job", default=1)

    args = parser.parse_args()

    return args.config_file, args.numofproc
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# store_result helper function
results = {}
def store_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    results[result.cell_id] = result
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
