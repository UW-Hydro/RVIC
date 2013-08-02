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
# try:
#     from IPython.parallel import Client
# except:
#     pass
# try:
#     import SGEpy
# except:
#     pass

# -------------------------------------------------------------------- #
# Top level driver
def main():

    uh_box, fdr_data, fdr_vatts, dom_data, outlets, config_dict, directories = gen_uh_init()

    outlets, fractions, uh = gen_uh_run(uh_box, fdr_data, fdr_vatts, dom_data, outlets, config_dict, directories)

    gen_uh_final(outlets, fractions, uh, dom_data, config_dict, directories)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Initialize the Genuh Program
def gen_uh_init(configFile=None):
    """Initialize RVIC parameter """

    # ---------------------------------------------------------------- #
    # Read command Line
    if not configFile:
        configFile = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(configFile)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structure
    directories = make_directories(config_dict['options']['case_dir'],
                        ['aggregated', 'remapped', 'plots',
                        'logs', 'params', 'inputs'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # copy inputs to $case_dir/inputs and update configuration
    config_dict = copy_inputs(configFile, directories['inputs'])
    options = config_dict['options']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Start Logging
    log = init_logger(directories['logs'], options['log_level'], options['verbose'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Pour Points files
    pour_points = {}
    try:
        pour_points['lons'], pour_points['lats'] = np.genfromtxt(config_dict['pour_points']['file_name'],
                                                                 skip_header=int(config_dict['pour_points']['header_lines']),
                                                                 delimiter=',', unpack=True)
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
    except:
        log.exception('Error opening uh_box file: %s' % config_dict['pour_points']['file_name'])
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read FDR file
    try:
        fdr_data, fdr_vatts, fdr_gatts = read_netcdf(config_dict['routing']['file_name'])
        fdr_shape = fdr_data[config_dict['routing']['flow_direction_var']].shape
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
    except:
        log.exception('Error opening FDR file')
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read domain file (if applicable)
    if options['remap']:
        dom_data, DomVats, DomGats = read_domain(config_dict['domain'])
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
    else:
        outlets = {}
        for i, (lon, lat) in enumerate(zip(pour_points['lons'], pour_points['lats'])):
            outlets[i] = Point(lat=lat, lon=lon)
            outlets[i].pour_points = [Point(lat=lat, lon=lon)]

    # ---------------------------------------------------------------- #
    return uh_box, fdr_data, fdr_vatts, dom_data, outlets, config_dict, directories
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# The standard single processor driver for gerating parameter files
def gen_uh_run(uh_box, fdr_data, fdr_vatts, dom_data, outlets, config_dict, directories):
    """
    Run Genuh_run on a single processor (slow)
    """
    log = getLogger(LOG_NAME)

    log.info('Starting gen_uh_run now.')

    # ---------------------------------------------------------------- #
    # Loop over agg points
    for i, cell_id in enumerate(outlets):

        agg_data = {}
        # ------------------------------------------------------------ #
        # Loop over pour points
        for j, pour_point in enumerate(outlets[cell_id].pour_points):

            # -------------------------------------------------------- #
            # Make the Unit Hydrograph Grid
            routData = rout(pour_point, uh_box, fdr_data, fdr_vatts,
                            config_dict['routing'])

            log.debug('routData: %f, %f' % (routData['uhgrid'].min(), routData['uhgrid'].max()))

            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # aggregate
            if config_dict['options']['aggregate']:
                if j != len(outlets[cell_id].pour_points)-1:
                    agg_data = aggregate(routData, agg_data, res=fdr_data['resolution'])
                else:
                    agg_data = aggregate(routData, agg_data, res=fdr_data['resolution'],
                                         pad=config_dict['options']['agg_pad'], maskandnorm=True)

                    log.debug('agg_data: %f, %f' % (agg_data['uhgrid'].min(), agg_data['uhgrid'].max()))
            elif config_dict['options']['remap']:
                agg_data = routData
            else:
                remap_data = routData
            # -------------------------------------------------------- #

        # ------------------------------------------------------------ #
        # write temporary file #1
        if  config_dict['options']['remap']:
            glob_atts = NcGlobals(title='RVIC Unit Hydrograph Grid File',
                                RvicPourPointsFile=os.path.split(config_dict['pour_points']['file_name'])[1],
                                RvicUhFile=os.path.split(config_dict['uh_box']['file_name'])[1],
                                RvicFdrFile=os.path.split(config_dict['routing']['file_name'])[1],
                                RvicDomainFile=os.path.split(config_dict['domain']['file_name'])[1])

            temp_file_1 = os.path.join(directories['aggregated'], 'aggUH_%i.nc' % cell_id)

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

            temp_file_2 = os.path.join(directories['remapped'], 'remapUH_%i.nc' % cell_id)

            remap(config_dict['domain']['file_name'], temp_file_1, temp_file_2)

            # -------------------------------------------------------- #
            # Clean temporary file #1 (if applicable)
            if config_dict['options']['clean']:
                clean_file(temp_file_1)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # Read temporary file #2
            remap_data, remap_vatts, remap_gatts = read_netcdf(temp_file_2)
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
            if i == 0:
                fractions = np.zeros(remap_data['fractions'].shape)
                unit_hydrographs = np.zeros(remap_data['unit_hydrographs'].shape)

            y, x = np.nonzero((remap_data['fractions'] > 0.0) * (dom_data[config_dict['domain']['land_mask_var']] == 1))

            outlets[cell_id].fractions = remap_data['fractions'][y, x]
            outlets[cell_id].unit_hydrographs = remap_data['unit_hydrographs'][:, y, x]
            outlets[cell_id].time = np.arange(remap_data['unit_hydrographs'].shape[0])
            outlets[cell_id].lon_source = remap_data[config_dict['domain']['longitude_var']][y, x]
            outlets[cell_id].lat_source = remap_data[config_dict['domain']['latitude_var']][y, x]
            outlets[cell_id].cell_id_source = dom_data['cell_ids'][y, x]
            outlets[cell_id].x_source = x
            outlets[cell_id].y_source = y

            fractions[y, x] += remap_data['fractions'][y, x]
            unit_hydrographs[:, y, x] += remap_data['unit_hydrographs'][:, y, x]
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
    return outlets, fractions, unit_hydrographs
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# The Poor Mans Parallel Driver for generating parmeter files
def gen_uh_run_ppp():
    """
    Run Genuh_run in pour mans parallel (SGE only)
    """
    # ---------------------------------------------------------------- #
    # Loop over agg points
    # ---------------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # Loop over pour points
        # ------------------------------------------------------------ #
            # -------------------------------------------------------- #
            # rout(PourPoint, uh, fdr_data, RoutDict)
            # -------------------------------------------------------- #
            # -------------------------------------------------------- #
            # aggregate
            # -------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # write temporary file #1
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Remap temporary file #1 to temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #1 (if applicable)
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Read temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Add to adjust fractions Structure
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #2 (if applicable)
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
    # Adjust fractions (if applicable)
    # ---------------------------------------------------------------- #
    # Write parameter file
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def gen_uh_run_ipp():
    """
    Run GenUH_run using ipython parallel
    """
    # ---------------------------------------------------------------- #
    # Loop over agg points
    # ---------------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # Loop over pour points
        # ------------------------------------------------------------ #
            # -------------------------------------------------------- #
            # RvicMakeUH()
            # -------------------------------------------------------- #
            # -------------------------------------------------------- #
            # aggregate
            # -------------------------------------------------------- #
        # ------------------------------------------------------------ #
        # write temporary file #1
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Remap temporary file #1 to temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #1 (if applicable)
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Read temporary file #2
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Add to adjust fractions Structure
        # ------------------------------------------------------------ #
        # ------------------------------------------------------------ #
        # Clean temporary file #2 (if applicable)
        # ------------------------------------------------------------ #
    # ---------------------------------------------------------------- #
    # Adjust fractions (if applicable)
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    # Write parameter file
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def gen_uh_final(outlets, fractions, unit_hydrographs, dom_data, config_dict, directories):
    """
    Make the RVIC Parameter File
    """

    log = getLogger(LOG_NAME)

    log.info('Starting gen_uh_final now.')

    options = config_dict['options']
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

        offset, out_uh, full_length = subset(outlets[cell_id].unit_hydrographs,
                                             options['subset_length'],
                                             options['subset_threshold'])
        outlets[cell_id].fractions *= ratio_fraction[y, x]  # Adjust fracs based on ratio_fraction
        adjust_fractions[y, x] += outlets[cell_id].fractions

        if i == 0:
            # -------------------------------------------------------- #
            # Source specific values
            unit_hydrographs = out_uh
            frac_sources = outlets[cell_id].fractions
            lon_sources = outlets[cell_id].lon_source
            lat_sources = outlets[cell_id].lat_source
            x_ind_source = x
            y_ind_source = y
            source_decomp_id = outlets[cell_id].cell_id_source
            time_offset_source = offset
            source2outlet_index = np.zeros(len(outlets[cell_id].fractions))
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_decomp_id = np.array(cell_id)
            lon_outlet = np.array(outlets[cell_id].lon)
            lat_outlet = np.array(outlets[cell_id].lat)
            x_ind_outlet = np.array(outlets[cell_id].x)
            y_ind_outlet = np.array(outlets[cell_id].y)
            outlet_nums = np.array(i)

            # -------------------------------------------------------- #
            # get a few global values
            subset_length = options['subset_length']
            full_time_length = full_length
            unit_hydrogaph_dt = config_dict['routing']['output_interval']
            # -------------------------------------------------------- #

        else:
            # -------------------------------------------------------- #
            # Point specific values
            unit_hydrographs = np.append(unit_hydrographs, out_uh, axis=1)
            frac_sources = np.append(frac_sources, outlets[cell_id].fractions)
            lon_sources = np.append(lon_sources, outlets[cell_id].lon_source)
            lat_sources = np.append(lat_sources, outlets[cell_id].lat_source)
            x_ind_source = np.append(x_ind_source, x)
            y_ind_source = np.append(y_ind_source, y)
            source_decomp_id = np.append(source_decomp_id, outlets[cell_id].cell_id_source)
            time_offset_source = np.append(time_offset_source, offset)
            source2outlet_index = np.append(source2outlet_index, np.zeros_like(offset) + i)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_decomp_id = np.append(outlet_decomp_id, cell_id)
            lon_outlet = np.append(lon_outlet, outlets[cell_id].lon)
            lat_outlet = np.append(lat_outlet, outlets[cell_id].lat)
            x_ind_outlet = np.append(x_ind_outlet, outlets[cell_id].x)
            y_ind_outlet = np.append(y_ind_outlet, outlets[cell_id].y)
            outlet_nums = np.append(outlet_nums, i)
            # -------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust Unit Hydrographs for differences in source/outlet areas
    unit_hydrographs *= frac_sources
    unit_hydrographs *= dom_data[config_dict['domain']['area_var']][y_ind_source, x_ind_source]

    for p, ind in enumerate(source2outlet_index):
        unit_hydrographs[:, p] /= dom_data[config_dict['domain']['area_var']][y_ind_outlet[ind], x_ind_outlet[ind]]
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #
    # Write parameter file
    today = date.today().strftime('%Y%m%d')
    param_file = os.path.join(directories['params'],
                             '%s.rvic.prm.%s.%s.nc' % (options['caseid'],
                                                       options['gridid'], today))

    write_param_file(param_file, options['netcdf_format'], subset_length,
                     full_time_length, unit_hydrogaph_dt, unit_hydrographs,
                     lon_sources, lat_sources, source_decomp_id, x_ind_source,
                     y_ind_source, outlet_decomp_id, time_offset_source,
                     source2outlet_index, lon_outlet, lat_outlet, x_ind_outlet,
                     y_ind_outlet, outlet_nums,
                     NcGlobals(title='RVIC parameter file',
                               RvicPourPointsFile=os.path.split(config_dict['pour_points']['file_name'])[1],
                               RvicuhFile=os.path.split(config_dict['uh_box']['file_name'])[1],
                               RvicFdrFile=os.path.split(config_dict['routing']['file_name'])[1],
                               RvicDomainFile=os.path.split(config_dict['domain']['file_name'])[1]))
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
def process_command_line():
    """
    Get the path to the configFile
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generate RVIC parameter files.')
    parser.add_argument("configFile", type=str, help="Input configuration file")

    args = parser.parse_args()
    configFile = args.configFile

    return configFile


# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
