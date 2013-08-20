#!/opt/local/bin/python
"""
This is the convolution routine developed in preparation in coupling RVIC to
CESM.  Eventually, this will be the offline RVIC model.

Written by Joe Hamman, May 2013

__________________________________
REVISION HISTORY
--------
July 2013, Joe Hamman
Changed input file type to standard RVIC parameter file
Made necessary changes to run routines to accept the new parameter file structure.
Major updates to the...
"""
import os
from argparse import ArgumentParser
from logging import getLogger
from rvic.log import init_logger, LOG_NAME
from rvic.utilities import read_config, make_directories, copy_inputs, read_domain
from rvic.utilities import write_rpointer, tar_inputs
from rvic.variables import Rvar
from rvic.time_utility import Dtime
from rvic.read_forcing import DataModel
from rvic.history import Tape
from rvic.share import NcGlobals

# -------------------------------------------------------------------- #
# Top Level Driver
def main(config_file=None):
    """
    Top level driver for RVIC model.
    """

    hist_tapes, data_model, rout_var, dom_data, time_handle, directories, config_dict = rvic_mod_init()

    time_handle, hist_tapes = rvic_mod_run(hist_tapes, data_model, rout_var, dom_data, time_handle, directories, config_dict)

    rvic_mod_final(time_handle, hist_tapes)

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Initialize RVIC
def rvic_mod_init(config_file=None):
    """
    - Read Grid File
    - Load the unit hydrograph files (put into point_dict)
    - Load the initial state file and put it in convolution rings
    """

    # ---------------------------------------------------------------- #
    # Read command Line
    if not config_file:
        config_file = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structure
    directories = make_directories(config_dict['options']['case_dir'],
                                   ['hist', 'logs', 'params', 'restarts'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Copy Inputs to $case_dir/inputs and update configuration
    config_dict = copy_inputs(config_file, directories['inputs'])
    options = config_dict['options']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Settup Logging
    log = init_logger(directories['logs'], options['log_level'], options['verbose'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Domain File
    Domain = config_dict['domain']
    dom_data, dom_vatts, dom_gatts = read_domain(Domain['file_name'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read the Parameter File
    log.info('reading parameter and state file')

    rout_var = Rvar(config_dict['PARAM_FILE']['file_name'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read the State File
    if options['RUN_TYPE'] == 'restart':
        restart = read_config(os.path.join(directories['restarts'], 'rpointer'))
        timestr = restart['restart']['timestamp']
        rout_var.initial_state(restart['restart']['file_name'])
    elif options['RUN_TYPE'] == 'startup':
        timestr = options['RUN_STARTDATE']
        rout_var.initial_state(config_dict['INITIAL_STATE']['file_name'])
    elif options['RUN_TYPE'] == 'drystart':
        timestr = options['RUN_STARTDATE']
    else:
        log.error('RUN_TYPE Option is none of these: (restart, startup, drystart)')
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup time_handle
    time_handle = Dtime(timestr, options['STOP_OPTION'], options['STOP_N'],
                        options['STOP_DATE'], options['REST_OPTION'],
                        options['REST_N'], options['REST_DATE'],
                        options['CALENDAR'], rout_var.uh_timestep)

    # ---------------------------------------------------------------- #
    # Initialize the data model
    forcings = config_dict['INPUT_FORCINGS']
    data_model = DataModel(forcings['datl_path'],
                           forcings['datl_file'],
                           forcings['time_var'],
                           forcings['datl_liq_flds'],
                           forcings['start'],
                           forcings['end'],
                           time_handle.timestamp)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup history Tape(s) and Write Initial Outputs
    history = config_dict['HISTORY']
    hist_tapes = {}

    for j in xrange(config_dict['HISTORY']['rvichist_ntapes']):
        tapename = 'Tape.%i' % j
        hist_tapes[tapename] = Tape(time_handle.timestamp,
                                    time_handle.time_ord,
                                    options['CASEID'],
                                    rout_var,
                                    fincl=['streamflow'],
                                    mfilt=history['rvichist_mfilt'][j],
                                    ndens=history['rvichist_ndens'][j],
                                    nhtfrq=history['rvichist_nhtfrq'][j],
                                    avgflag=history['rvichist_avgflag'][j],
                                    file_format=history['rvichist_ncform'][j],
                                    outtype=history['rvichist_outtype'][j],
                                    dom_lons=dom_data[Domain['longitude_var']],
                                    dom_lats=dom_data[Domain['latitude_var']],
                                    out_dir=directories['hist'],
                                    calendar=time_handle.calendar,
                                    GlobAts=NcGlobals(title='RVIC history file',
                                                      casename=options['CASEID'],
                                                      casestr=options['CASESTR'],
                                                      RvicPourPointsFile=os.path.split(rout_var.RvicPourPointsFile)[1],
                                                      RvicUHFile=os.path.split(rout_var.RvicUHFile)[1],
                                                      RvicFdrFile=os.path.split(rout_var.RvicFdrFile)[1],
                                                      RvicDomainFile=os.path.split(Domain['file_name'])[1]))
        hist_tapes[tapename].write_initial()
    # ---------------------------------------------------------------- #

    return hist_tapes, data_model, rout_var, dom_data, time_handle, directories, config_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def rvic_mod_run(hist_tapes, data_model, rout_var, dom_data, time_handle, directories, config_dict):
    """
    - Loop over flux files
    - Combine the Baseflow and Runoff Variables
    - Adjust units as necessary
    - Do convolution
    - Write output files
    """

    data2tape = {}

    # ---------------------------------------------------------------- #
    # Start log
    log = getLogger(LOG_NAME)
    log.info('Starting rvic_mod_run')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Iterate Through time_handlesteps
    while not time_handle.stop_flag:
        timestamp = time_handle.advance_timestep()
        time_ord = time_handle.time_ord

        # ------------------------------------------------------------ #
        # Get This time_handlesteps Forcing
        aggrunin = data_model(timestamp)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Do the Convolution
        rout_var.convolve(aggrunin, timestamp)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Extract the Current Variables from rout_var
        data2tape['streamflow'] = rout_var.get_rof()
        data2tape['storage'] = rout_var.get_storage()

        # Update the history Tape(s)
        for tapename, tape in hist_tapes.iteritems():
            log.debug('Updating Tape:%s' % tapename)
            tape.update(time_ord, data2tape)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write State
        if time_handle.rest_flag or time_handle.stop_flag:
            restart_file = rout_var.write_restart()

            write_rpointer(directories['restarts'], restart_file, timestamp)
        # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #
    return time_handle, hist_tapes
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Final
def rvic_mod_final(time_handle, hist_tapes):
    """ Finalize RVIC """
    # ---------------------------------------------------------------- #
    # Start log
    log = getLogger(LOG_NAME)
    log.info('Finalizing RVIC')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write final log info
    log.info("-----------------------------------------------------------")
    log.info('Done with streamflow convolution')
    log.info('Processed %i timesteps' % time_handle.timesteps)
    for name, tape in hist_tapes.iteritems():
        log.info('Wrote %i history files from %s' % (tape.files_count, name))
    log.info('Routed to %i points' % time_handle.points)
    log.info("-----------------------------------------------------------")
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    log_tar = tar_inputs(log.filename)

    log.info('Done with RvicGenParam.')
    log.info('Location of Log: %s' % log_tar)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Process Command Line
def process_command_line():
    """
    Get the path to the config_file
    """
    # Parse arguments
    parser = ArgumentParser(description='RVIC is based on the original model of Lohmann, et al., 1996, Tellus, 48(A), 708-721')
    parser.add_argument("config_file", type=str, help="Input configuration file")

    args = parser.parse_args()
    config_file = args.config_file

    return config_file
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Run Program
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
