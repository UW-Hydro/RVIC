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
from rvic.mpi import LoggingPool
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

    # ---------------------------------------------------------------- #
    # Read command Line
    config_file, numofproc = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initilize
    hist_tapes, data_model, rout_var, dom_data,\
    time_handle, directories, config_dict = rvic_mod_init(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup the pool of processors
    # if numofproc > 1:
    #     pool = LoggingPool(processes=numofproc)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run
    time_handle, hist_tapes = rvic_mod_run(hist_tapes, data_model, rout_var,
                                           dom_data, time_handle, directories,
                                           config_dict)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Finalize
    rvic_mod_final(time_handle, hist_tapes)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Initialize RVIC
def rvic_mod_init(config_file):
    """
    - Read Grid File
    - Load the unit hydrograph files (put into point_dict)
    - Load the initial state file and put it in convolution rings
    """

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structure
    directories = make_directories(config_dict['OPTIONS']['CASE_DIR'],
                                   ['hist', 'logs', 'restarts']) #'params',
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Copy Inputs to $case_dir/inputs and update configuration
    #config_dict = copy_inputs(config_file, directories['inputs'])
    options = config_dict['OPTIONS']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Settup Logging
    log = init_logger(directories['logs'], options['LOG_LEVEL'], options['VERBOSE'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Domain File
    Domain = config_dict['DOMAIN']
    dom_data, dom_vatts, dom_gatts = read_domain(Domain)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read the Parameter File
    log.info('reading parameter file %s' % config_dict['PARAM_FILE']['FILE_NAME'])

    rout_var = Rvar(config_dict['PARAM_FILE']['FILE_NAME'], options['CASEID'],
                    options['CALENDAR'], directories['restarts'], options['REST_NCFORM'])
    rout_var.check_grid_file(Domain['FILE_NAME'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine the restart options
    restart_file = None

    if options['RUN_TYPE'] == 'restart':
        restart = read_config(os.path.join(directories['restarts'], 'rpointer'))
        timestr = restart['RESTART']['TIMESTAMP']
        restart_file = restart['RESTART']['FILE_NAME']
    elif options['RUN_TYPE'] == 'startup':
        timestr = options['RUN_STARTDATE']
        restart_file = config_dict['INITIAL_STATE']['FILE_NAME']
    elif options['RUN_TYPE'] == 'drystart':
        timestr = options['RUN_STARTDATE']
    else:
        raise ValueError('RUN_TYPE option is none of these: (restart, startup, drystart)')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup time_handle
    time_handle = Dtime(timestr, options['STOP_OPTION'], options['STOP_N'],
                        options['STOP_DATE'], options['REST_OPTION'],
                        options['REST_N'], options['REST_DATE'],
                        options['CALENDAR'], rout_var.unit_hydrograph_dt)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read initial state
    rout_var.init_state(restart_file, options['RUN_TYPE'], time_handle.timestamp)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initialize the data model
    forcings = config_dict['INPUT_FORCINGS']
    data_model = DataModel(forcings['DATL_PATH'],
                           forcings['DATL_FILE'],
                           forcings['TIME_VAR'],
                           forcings['DATL_LIQ_FLDS'],
                           forcings['START'],
                           forcings['END'],
                           time_handle.timestamp)

    time_handle.end = data_model.end
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup history Tape(s) and Write Initial Outputs
    history = config_dict['HISTORY']
    hist_tapes = {}
    for j in xrange(int(config_dict['HISTORY']['RVICHIST_NTAPES'])):
        tapename = 'Tape.%i' % j
        log.info('setting up History %s' %tapename)
        hist_tapes[tapename] = Tape(time_handle.time_ord,
                                    options['CASEID'],
                                    rout_var,
                                    fincl=['streamflow'],
                                    mfilt=int(history['RVICHIST_MFILT'][j]),
                                    ndens=int(history['RVICHIST_NDENS'][j]),
                                    nhtfrq=int(history['RVICHIST_NHTFRQ'][j]),
                                    avgflag=history['RVICHIST_AVGFLAG'][j],
                                    file_format=history['RVICHIST_NCFORM'][j],
                                    outtype=history['RVICHIST_OUTTYPE'][j],
                                    grid_lons=dom_data['cord_lons'],
                                    grid_lats=dom_data['cord_lats'],
                                    out_dir=directories['hist'],
                                    calendar=time_handle.calendar,
                                    glob_ats=NcGlobals(title='RVIC history file',
                                                      casename=options['CASEID'],
                                                      casestr=options['CASESTR'],
                                                      RvicPourPointsFile=os.path.split(rout_var.RvicPourPointsFile)[1],
                                                      RvicUHFile=os.path.split(rout_var.RvicUHFile)[1],
                                                      RvicFdrFile=os.path.split(rout_var.RvicFdrFile)[1],
                                                      RvicDomainFile=os.path.split(Domain['FILE_NAME'])[1]))

    # loop over again and print summary
    for tapename, tape in hist_tapes.iteritems():
        log.info('==========%s==========' %tapename)
        log.info(tape)
        tape.write_initial()
    # ---------------------------------------------------------------- #

    return hist_tapes, data_model, rout_var, dom_data, time_handle, directories, config_dict
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def rvic_mod_run(hist_tapes, data_model, rout_var, dom_data, time_handle,
                 directories, config_dict):
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

    # ------------------------------------------------------------ #
    # Get initial time info
    time_ord = time_handle.time_ord
    timestamp = time_handle.timestamp
    stop_flag = False
    rest_flag = False
    # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #
    # Iterate Through time_handlesteps
    while True:
        # ------------------------------------------------------------ #
        # Get This time_handlesteps Forcing
        aggrunin = data_model.read(timestamp)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Do the Convolution
        # (end_timestamp is the timestamp at the end of the convolution period)
        end_timestamp = rout_var.convolve(aggrunin, time_ord)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Extract the Current Variables from rout_var
        data2tape['streamflow'] = rout_var.get_rof()
        data2tape['storage'] = rout_var.get_storage()

        # Update the history Tape(s)
        for tapename, tape in hist_tapes.iteritems():
            log.debug('Updating Tape:%s' % tapename)
            tape.update(data2tape, time_ord)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write State
        if rest_flag:
            restart_file = rout_var.write_restart()
            write_rpointer(directories['restarts'], restart_file, end_timestamp)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # advance a single timestep
        if not stop_flag:
            timestamp, time_ord, stop_flag, rest_flag = time_handle.advance_timestep()
            #check that we're still inline with convolution
            if end_timestamp != timestamp:
                raise ValueError('timestamps do not match after convolution')
        else:
            break
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
    # log.info('Routed to %i points' % time_handle.points)
    log.info("-----------------------------------------------------------")
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    log_tar = tar_inputs(log.filename)

    log.info('Done with rvic_model.')
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
    parser.add_argument("-np", "--numofproc", type=int,
                        help="Number of processors used to run job", default=1)

    args = parser.parse_args()

    return args.config_file, args.numofproc
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Run Program
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
