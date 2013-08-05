#!/opt/local/bin/python
"""
RVIC parameter file development driver
"""
import os
import numpy as np
import argparse
from datetime import date
from logging import getLogger
from rvic.log import init_logger, LOG_NAME
from rvic.utilities import read_config, make_directories, copy_inputs, read_netcdf, tar_inputs
from rvic.utilities import check_ncvars, remap, clean_file, subset, read_domain
from rvic.aggregate import make_agg_pairs, aggregate
from rvic.make_uh import rout
from rvic.share import NcGlobals
from rvic.write import write_agg_netcdf, write_param_file
from rvic.variables import Point
from multiprocessing import Pool


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
    # Setup the pool of processors
    pool = Pool(processes=numofproc)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get main logger
    log = getLogger(LOG_NAME)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run
    for i, (cell_id, outlet) in enumerate(outlets.iteritems()):
        log.info('On Outlet #%i of %i' %(i+1, len(outlets)))
        pool.apply_async(gen_uh_run,
                         args = (uh_box, fdr_data, fdr_vatts, dom_data, outlet, config_dict, directories),
                         callback = store_result)
    pool.close()
    pool.join()

    outlets = results
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
    directories = make_directories(config_dict['options']['case_dir'],
                                   ['plots', 'logs', 'params', 'inputs'])
    directories.update(make_directories(config_dict['options']['temp_dir'],
                                        ['aggregated', 'remapped']))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # copy inputs to $case_dir/inputs and update configuration
    config_dict = copy_inputs(config_file, directories['inputs'])
    options = config_dict['options']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Start Logging
    log = init_logger(directories['logs'], options['log_level'], options['verbose'])

    for direc in directories:
        log.info('%s directory is %s' % (direc, directories[direc]))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Pour Points files
    pour_points = {}
    try:
        pour_points['lons'], pour_points['lats'] = np.genfromtxt(config_dict['pour_points']['file_name'],
                                                                 skip_header=int(config_dict['pour_points']['header_lines']),
                                                                 delimiter=',', unpack=True)
        log.info('Opened Pour Points File: %s' % config_dict['pour_points']['file_name'])
    except:
        log.exception('Error opening pour points file: %s' % config_dict['pour_points']['file_name'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read uh box file
    uh_box = {}
    try:
        uh_box['time'], uh_box['func'] = np.genfromtxt(config_dict['uh_box']['file_name'],
                                                       skip_header=int(config_dict['uh_box']['header_lines']),
                                                       delimiter=',', unpack=True)
        log.info('Opened UHbox File: %s' % config_dict['uh_box']['file_name'])
    except:
        log.exception('Error opening uh_box file: %s' % config_dict['pour_points']['file_name'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read FDR file
    try:
        fdr_data, fdr_vatts, fdr_gatts = read_netcdf(config_dict['routing']['file_name'])
        fdr_shape = fdr_data[config_dict['routing']['flow_direction_var']].shape

        # ---------------------------------------------------------------- #
        # Check latitude order, flip if necessary.
        if fdr_data[config_dict['routing']['latitude_var']][-1] > fdr_data[config_dict['routing']['latitude_var']][0]:
            log.debug('Inputs came in upside down, flipping everything now.')
            var_list = fdr_data.keys()
            var_list.remove(config_dict['routing']['longitude_var'])
            for var in var_list:
                fdr_data[var] = np.flipud(fdr_data[var])
        # ---------------------------------------------------------------- #

        # Add velocity and/or diffusion grids if not present yet
        if not type(config_dict['routing']['velocity']) == str:
            fdr_data['velocity'] = np.zeros(fdr_shape) + config_dict['routing']['velocity']
            config_dict['routing']['velocity'] = 'velocity'
            log.info('Added velocity grid to fdr_data')
        if not type(config_dict['routing']['diffusion']) == str:
            fdr_data['diffusion'] = np.zeros(fdr_shape) + config_dict['routing']['diffusion']
            config_dict['routing']['diffusion'] = 'diffusion'
            log.info('Added diffusion grid to fdr_data')

        fdr_data['resolution'] = np.abs(fdr_data[config_dict['routing']['longitude_var']][1] -
                                       fdr_data[config_dict['routing']['longitude_var']][0])

        check_ncvars(config_dict['routing'], fdr_data.keys())

        log.info('Opened FDR File: %s' % config_dict['routing']['file_name'])
    except:
        log.exception('Error opening FDR file')
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read domain file (if applicable)
    if options['remap']:
        dom_data, DomVats, DomGats = read_domain(config_dict['domain'])
        log.info('Opened Domain File: %s' % config_dict['domain']['file_name'])
    else:
        dom_data = {}
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Group pour points (if aggregate)
    if options['aggregate']:
        outlets = make_agg_pairs(pour_points['lons'], pour_points['lats'],
                                 dom_data[config_dict['domain']['longitude_var']],
                                 dom_data[config_dict['domain']['latitude_var']],
                                 dom_data['cell_ids'], agg_type='agg')

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
    Run Genuh_run on a single processor (slow)
    """
    log = getLogger(LOG_NAME)

    log.info('Running outlet cell id  %i' % outlet.cell_id)

    agg_data = {}
    # ------------------------------------------------------------ #
    # Loop over pour points
    for j, pour_point in enumerate(outlet.pour_points):

        log.info('On pour_point #%i of %i' %(j+1, len(outlet.pour_points)))

        # -------------------------------------------------------- #
        # Make the Unit Hydrograph Grid
        rout_data = rout(pour_point, uh_box, fdr_data, fdr_vatts,
                        config_dict['routing'])

        log.debug('Done routing to pour_point')
        log.debug('rout_data: %f, %f' % (rout_data['unit_hydrograph'].min(), rout_data['unit_hydrograph'].max()))

        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # aggregate
        if config_dict['options']['aggregate']:
            if j != len(outlet.pour_points)-1:
                agg_data = aggregate(rout_data, agg_data, res=fdr_data['resolution'])
            else:
                agg_data = aggregate(rout_data, agg_data, res=fdr_data['resolution'],
                                     pad=config_dict['options']['agg_pad'], maskandnorm=True)

                log.debug('agg_data: %f, %f' % (agg_data['unit_hydrograph'].min(), agg_data['unit_hydrograph'].max()))
        elif config_dict['options']['remap']:
            agg_data = rout_data
        else:
            remap_data = rout_data
        # -------------------------------------------------------- #

    # ------------------------------------------------------------ #
    # write temporary file #1
    if  config_dict['options']['remap']:
        glob_atts = NcGlobals(title='RVIC Unit Hydrograph Grid File',
                              RvicPourPointsFile=os.path.split(config_dict['pour_points']['file_name'])[1],
                              RvicUHFile=os.path.split(config_dict['uh_box']['file_name'])[1],
                              RvicFdrFile=os.path.split(config_dict['routing']['file_name'])[1],
                              RvicDomainFile=os.path.split(config_dict['domain']['file_name'])[1])

        temp_file_1 = os.path.join(directories['aggregated'], 'aggUH_%i.nc' % outlet.cell_id)

        write_agg_netcdf(temp_file_1, agg_data, glob_atts,
                         config_dict['options']['netcdf_format'])
    elif config_dict['options']['aggregate']:
        remap_data = agg_data
    else:
        pass
    # ------------------------------------------------------------ #

    # ------------------------------------------------------------ #
    # Remap temporary file #1 to temporary file #2
    if config_dict['options']['remap']:

        temp_file_2 = os.path.join(directories['remapped'], 'remapUH_%i.nc' % outlet.cell_id)

        remap(config_dict['domain']['file_name'], temp_file_1, temp_file_2)

        # -------------------------------------------------------- #
        # Clean temporary file #1 (if applicable)
        if config_dict['options']['clean']:
            clean_file(temp_file_1)
        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # Read temporary file #2
        remap_data, remap_vatts, remap_gatts = read_netcdf(temp_file_2,
                                                           variables=['unit_hydrograph', 'fraction'])
        # -------------------------------------------------------- #

        # -------------------------------------------------------- #
        # Clean temporary file #2 (if applicable)
        if config_dict['options']['clean']:
            clean_file(temp_file_2)
        # -------------------------------------------------------- #

    else:
        remap_data = agg_data
    # ------------------------------------------------------------ #

    # ------------------------------------------------------------ #
    # Add to adjust fractions Structure
    if config_dict['options']['remap']:
        y, x = np.nonzero((remap_data['fraction'] > 0.0) * (dom_data[config_dict['domain']['land_mask_var']] == 1))

        outlet.fractions = remap_data['fraction'][y, x]
        outlet.unit_hydrograph = remap_data['unit_hydrograph'][:, y, x]
        outlet.time = np.arange(remap_data['unit_hydrograph'].shape[0])
        outlet.lon_source = dom_data[config_dict['domain']['longitude_var']][y, x]
        outlet.lat_source = dom_data[config_dict['domain']['latitude_var']][y, x]
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

    log.info('Starting gen_uh_final now.')

    options = config_dict['options']


    # ---------------------------------------------------------------- #
    # Aggregate the fractions
    fractions = np.zeros(dom_data[config_dict['domain']['fraction_var']].shape)
    for i, cell_id in enumerate(outlets):
        y = outlets[cell_id].y_source
        x = outlets[cell_id].x_source

        fractions[y, x] += outlets[cell_id].fractions
        # unit_hydrograph[:, y, x] += outlets[cell_id].unit_hydrograph
    # ---------------------------------------------------------------- #


    # ---------------------------------------------------------------- #
    # Determin how to adjust the fractions
    grid_fractions = dom_data[config_dict['domain']['fraction_var']]
    diff_fractions = fractions - grid_fractions
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Only adjust fractions where the fractions are gt the domain fractions
    yi, xi = np.nonzero(fractions > grid_fractions)
    ratio_fraction = np.ones(diff_fractions.shape)
    adjust_fractions = np.zeros(diff_fractions.shape)
    ratio_fraction[yi, xi] = grid_fractions[yi, xi]/fractions[yi, xi]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Subset
    for i, cell_id in enumerate(outlets):

        y = outlets[cell_id].y_source
        x = outlets[cell_id].x_source

        offset, out_uh, full_length = subset(outlets[cell_id].unit_hydrograph,
                                             options['subset_length'],
                                             options['subset_threshold'])
        outlets[cell_id].fractions *= ratio_fraction[y, x]  # Adjust fracs based on ratio_fraction
        adjust_fractions[y, x] += outlets[cell_id].fractions

        if i == 0:
            # -------------------------------------------------------- #
            # Source specific values
            unit_hydrograph = out_uh
            frac_sources = outlets[cell_id].fractions
            source_lon = outlets[cell_id].lon_source
            source_lat = outlets[cell_id].lat_source
            source_x_ind = x
            source_y_ind = y
            source_decomp_ind = outlets[cell_id].cell_id_source
            source_time_offset = offset
            source2outlet_ind = np.zeros(len(outlets[cell_id].fractions))
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_decomp_ind = np.array(cell_id)
            outlet_lon = np.array(outlets[cell_id].lon)
            outlet_lat = np.array(outlets[cell_id].lat)
            outlet_x_ind = np.array(outlets[cell_id].x)
            outlet_y_ind = np.array(outlets[cell_id].y)
            outlet_number = np.array(i)

            # -------------------------------------------------------- #
            # get a few global values
            subset_length = options['subset_length']
            full_time_length = full_length
            unit_hydrograph_dt = config_dict['routing']['output_interval']
            # -------------------------------------------------------- #

        else:
            # -------------------------------------------------------- #
            # Point specific values
            unit_hydrograph = np.append(unit_hydrograph, out_uh, axis=1)
            frac_sources = np.append(frac_sources, outlets[cell_id].fractions)
            source_lon = np.append(source_lon, outlets[cell_id].lon_source)
            source_lat = np.append(source_lat, outlets[cell_id].lat_source)
            source_x_ind = np.append(source_x_ind, x)
            source_y_ind = np.append(source_y_ind, y)
            source_decomp_ind = np.append(source_decomp_ind, outlets[cell_id].cell_id_source)
            source_time_offset = np.append(source_time_offset, offset)
            source2outlet_ind = np.append(source2outlet_ind, np.zeros_like(offset) + i)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_decomp_ind = np.append(outlet_decomp_ind, cell_id)
            outlet_lon = np.append(outlet_lon, outlets[cell_id].lon)
            outlet_lat = np.append(outlet_lat, outlets[cell_id].lat)
            outlet_x_ind = np.append(outlet_x_ind, outlets[cell_id].x)
            outlet_y_ind = np.append(outlet_y_ind, outlets[cell_id].y)
            outlet_number = np.append(outlet_number, i)
            # -------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust Unit Hydrographs for differences in source/outlet areas
    unit_hydrograph *= frac_sources
    unit_hydrograph *= dom_data[config_dict['domain']['area_var']][source_y_ind, source_x_ind]

    for p, ind in enumerate(source2outlet_ind):
        unit_hydrograph[:, p] /= dom_data[config_dict['domain']['area_var']][outlet_y_ind[ind], outlet_x_ind[ind]]
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    # Write parameter file
    today = date.today().strftime('%Y%m%d')
    param_file = os.path.join(directories['params'],
                             '%s.rvic.prm.%s.%s.nc' % (options['caseid'],
                                                       options['gridid'], today))

    write_param_file(param_file,
                         nc_format = options['netcdf_format'],
                         glob_atts = NcGlobals(title='RVIC parameter file',
                                               RvicPourPointsFile=os.path.split(config_dict['pour_points']['file_name'])[1],
                                               RvicUHFile=os.path.split(config_dict['uh_box']['file_name'])[1],
                                               RvicFdrFile=os.path.split(config_dict['routing']['file_name'])[1],
                                               RvicDomainFile=os.path.split(config_dict['domain']['file_name'])[1]),
                         full_time_length = full_time_length,
                         subset_length = subset_length,
                         unit_hydrograph_dt = unit_hydrograph_dt,
                         outlet_lon = outlet_lon,
                         outlet_lat = outlet_lat,
                         outlet_x_ind = outlet_x_ind,
                         outlet_y_ind = outlet_y_ind,
                         outlet_decomp_ind = outlet_decomp_ind,
                         outlet_number = outlet_number,
                         source_lon = source_lon,
                         source_lat = source_lat,
                         source_x_ind = source_x_ind,
                         source_y_ind = source_y_ind,
                         source_decomp_ind = source_decomp_ind,
                         source_time_offset = source_time_offset,
                         source2outlet_ind = source2outlet_ind,
                         unit_hydrograph = unit_hydrograph)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write a summary of what was done to the log file.
    log.info('Parameter file includes %i outlets' % (len(outlets)))
    log.info('Parameter file includes %i Source Points' % (len(outlets)))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    InputsTar = tar_inputs(directories['inputs'], suffix=today)
    LogTar = tar_inputs(log.filename)

    log.info('Done with RvicGenParam.')
    log.info('Location of Inputs: %s' % InputsTar)
    log.info('Location of Log: %s' % LogTar)
    log.info('Location of Parmeter File %s' % param_file)
    # ---------------------------------------------------------------- #
    return
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
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
