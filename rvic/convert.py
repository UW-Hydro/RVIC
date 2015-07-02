# -*- coding: utf-8 -*-
'''
Read a set of uhs files and write an RVIC parameter file
'''

from logging import getLogger
from .core.log import init_logger, close_logger, LOG_NAME
from .core.utilities import make_directories, copy_inputs, read_domain
from .core.utilities import tar_inputs
from .core.convert import read_station_file, read_uhs_files, move_domain
from .core.param_file import finish_params
from .core.config import read_config


# -------------------------------------------------------------------- #
# Top level driver
def convert(config_file):

    # ---------------------------------------------------------------- #
    # Initilize
    dom_data, new_dom_data, outlets, config_dict, \
        directories = uhs2param_init(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get main logger
    log = getLogger(LOG_NAME)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run
    log.info('getting outlets now')
    outlets = uhs2param_run(dom_data, outlets, config_dict)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Finally, make the parameter file
    uhs2param_final(outlets, dom_data, new_dom_data, config_dict, directories)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Init
def uhs2param_init(config_file):

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structures
    directories = make_directories(config_dict['OPTIONS']['CASE_DIR'],
                                   ['plots', 'logs', 'params', 'inputs'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # copy inputs to $case_dir/inputs and update configuration
    config_dict = copy_inputs(config_file, directories['inputs'])
    options = config_dict['OPTIONS']
    config_dict['POUR_POINTS'] = {
        'FILE_NAME': config_dict['UHS_FILES']['STATION_FILE']}
    config_dict['ROUTING']['FILE_NAME'] = 'unknown'
    config_dict['UH_BOX'] = {'FILE_NAME': 'unknown'}
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Start Logging
    log = init_logger(directories['logs'], options['LOG_LEVEL'],
                      options['VERBOSE'])

    for direc in directories:
        log.info('%s directory is %s', direc, directories[direc])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read domain file (if applicable)
    dom_data = read_domain(config_dict['DOMAIN'])[0]
    log.info('Opened Domain File: %s', config_dict['DOMAIN']['FILE_NAME'])

    if 'NEW_DOMAIN' in config_dict:
        new_dom_data = read_domain(config_dict['NEW_DOMAIN'])[0]
        log.info('Opened New Domain File: %s',
                 config_dict['NEW_DOMAIN']['FILE_NAME'])
    else:
        new_dom_data = None
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read station file
    outlets = read_station_file(config_dict['UHS_FILES']['STATION_FILE'],
                                dom_data, config_dict)
    # ---------------------------------------------------------------- #

    return dom_data, new_dom_data, outlets, config_dict, directories
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# run
def uhs2param_run(dom_data, outlets, config_dict):

    # ---------------------------------------------------------------- #
    # Read uhs files
    outlets = read_uhs_files(outlets, dom_data, config_dict)
    # ---------------------------------------------------------------- #

    return outlets
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
#
def uhs2param_final(outlets, dom_data, new_dom_data, config_dict, directories):
    '''
    Make the RVIC Parameter File
    '''

    log = getLogger(LOG_NAME)

    log.info('In gen_uh_final')

    # ---------------------------------------------------------------- #
    # Move to smaller domain
    if new_dom_data:
        outlets = move_domain(dom_data, new_dom_data, outlets)
        dom_data = new_dom_data
    # ---------------------------------------------------------------- #

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
