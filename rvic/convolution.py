# -*- coding: utf-8 -*-
'''
This is the convolution routine developed in preparation in coupling RVIC to
CESM.  Eventually, this will be the offline RVIC model.

Written by Joe Hamman, May 2013

__________________________________
REVISION HISTORY
--------
July 2013, Joe Hamman
Changed input file type to standard RVIC parameter file
Made necessary changes to run routines to accept the new parameter file
structure.
Major updates to the...
'''
import os
from collections import OrderedDict
from logging import getLogger
from .core.log import init_logger, close_logger, LOG_NAME
from .core.utilities import make_directories, read_domain
from .core.utilities import write_rpointer, tar_inputs
from .core.variables import Rvar
from .core.time_utility import Dtime
from .core.read_forcing import DataModel
from .core.history import Tape
from .core.share import NcGlobals, RVIC_TRACERS
from .core.config import read_config
from .core.pycompat import iteritems


# -------------------------------------------------------------------- #
# Top Level Driver
def convolution(config_file):
    '''
    Top level driver for RVIC convolution model.

    Parameters
    ----------
    config_file : str
        Path to RVIC convolution configuration file.
    '''

    # ---------------------------------------------------------------- #
    # Initilize
    hist_tapes, data_model, rout_var, \
        time_handle, directories = convolution_init(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup the pool of processors
    # if numofproc > 1:
    #     pool = LoggingPool(processes=numofproc)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run
    time_handle, hist_tapes = convolution_run(hist_tapes, data_model, rout_var,
                                              time_handle, directories)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Finalize
    convolution_final(time_handle, hist_tapes)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Initialize RVIC
def convolution_init(config_file):
    '''
    Initialize the RVIC convolution routine

    This function performs these main tasks:
        - Reads the configuration file
        - Sets up the RVIC case directories
        - Copies all input files to the case directory
        - Initializes the logging
        - Read domain file
        - Load the RVIC parameter file
        - Load the initial state file and put it in convolution rings
        - Setup time and history file objects

    Parameters
    ----------
    config_file : str
        Path to RVIC convolution configuration file.

    Returns
    ----------
    hist_tapes : OrderedDict
        Ordered dictionary of History objects
    data_model : DataModel
        DataModel instance containing the forcings
    rout_var : Rvar
        Rvar instance containing the RVIC parameters / unit hydrographs
    dom_data : dict
        Dictionary of arrays of mask, fraction, lats, lons, etc.
        This dictionary includes all the variables from the domain netCDF file.
    time_handle : Dtime
        Dtime instance containing information about run length, time
    directories : dict
        Dictionary of directories created by this function.
    config_dict : dict
        Dictionary of values from the configuration file.
    '''

    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = read_config(config_file)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structure
    directories = make_directories(config_dict['OPTIONS']['CASE_DIR'],
                                   ['hist', 'logs', 'restarts'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # unpack options
    options = config_dict['OPTIONS']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Settup Logging
    log = init_logger(directories['logs'], options['LOG_LEVEL'],
                      options['VERBOSE'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initialize the data model
    forcings = config_dict['INPUT_FORCINGS']
    data_model = DataModel(forcings['DATL_PATH'],
                           forcings['DATL_FILE'],
                           forcings['TIME_VAR'],
                           forcings['LATITUDE_VAR'],
                           forcings['DATL_LIQ_FLDS'],
                           forcings['START'],
                           forcings['END'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Domain File
    domain = config_dict['DOMAIN']
    dom_data = read_domain(domain, lat0_is_min=data_model.lat0_is_min)[0]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read the Parameter File
    log.info('reading parameter file %s',
             config_dict['PARAM_FILE']['FILE_NAME'])

    rout_var = Rvar(config_dict['PARAM_FILE']['FILE_NAME'], options['CASEID'],
                    options['CALENDAR'], directories['restarts'],
                    options['REST_NCFORM'])

    rout_var.set_domain(dom_data, domain, data_model.lat0_is_min)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine the restart options
    restart_file = None

    if options['RUN_TYPE'] == 'restart':
        restart = read_config(os.path.join(directories['restarts'],
                                           'rpointer'))
        timestr = restart['RESTART']['TIMESTAMP']
        restart_file = restart['RESTART']['FILE_NAME']
    elif options['RUN_TYPE'] == 'startup':
        timestr = options['RUN_STARTDATE']
        restart_file = config_dict['INITIAL_STATE']['FILE_NAME']
    elif options['RUN_TYPE'] == 'drystart':
        timestr = options['RUN_STARTDATE']
    else:
        raise ValueError('RUN_TYPE option {0} is none of these: (restart, '
                         'startup, drystart)'.format(options['RUN_TYPE']))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup time_handle
    time_handle = Dtime(timestr, options['STOP_OPTION'], options['STOP_N'],
                        options['STOP_DATE'], options['REST_OPTION'],
                        options['REST_N'], options['REST_DATE'],
                        options['CALENDAR'], data_model.secs_per_step)
    time_handle.end = data_model.end

    data_model.start(time_handle.timestamp, rout_var)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read initial state
    rout_var.init_state(restart_file, options['RUN_TYPE'],
                        time_handle.timestamp)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determine the number of aggregation timesteps
    rout_var.get_time_mode(data_model.secs_per_step)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup history Tape(s) and Write Initial Outputs
    history = config_dict['HISTORY']
    numtapes = int(history['RVICHIST_NTAPES'])
    hist_tapes = OrderedDict()

    # make sure history file fields are all in list form
    if numtapes == 1:
        for var, value in iteritems(history):
            if not isinstance(value, list):
                history[var] = list([value])

    global_atts = NcGlobals(
        title='RVIC history file',
        casename=options['CASEID'],
        casestr=options['CASESTR'],
        RvicPourPointsFile=os.path.split(rout_var.RvicPourPointsFile)[1],
        RvicUHFile=os.path.split(rout_var.RvicUHFile)[1],
        RvicFdrFile=os.path.split(rout_var.RvicFdrFile)[1],
        RvicDomainFile=os.path.split(domain['FILE_NAME'])[1])

    for j in range(numtapes):
        tapename = 'Tape.{0}'.format(j)
        log.info('setting up History %s', tapename)
        hist_tapes[tapename] = Tape(time_handle.time_ord,
                                    options['CASEID'],
                                    rout_var,
                                    tape_num=j,
                                    fincl=['streamflow'],
                                    mfilt=history['RVICHIST_MFILT'][j],
                                    ndens=int(history['RVICHIST_NDENS'][j]),
                                    nhtfrq=int(history['RVICHIST_NHTFRQ'][j]),
                                    avgflag=history['RVICHIST_AVGFLAG'][j],
                                    units=history['RVICHIST_UNITS'][j],
                                    file_format=history['RVICHIST_NCFORM'][j],
                                    outtype=history['RVICHIST_OUTTYPE'][j],
                                    grid_area=dom_data[domain['AREA_VAR']],
                                    grid_lons=dom_data['cord_lons'],
                                    grid_lats=dom_data['cord_lats'],
                                    out_dir=directories['hist'],
                                    calendar=time_handle.calendar,
                                    glob_ats=global_atts)

    # loop over again and print summary
    for tapename, tape in iteritems(hist_tapes):
        log.info('==========%s==========', tapename)
        log.info(tape)
        tape.write_initial()
    # ---------------------------------------------------------------- #

    return hist_tapes, data_model, rout_var, time_handle, directories
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def convolution_run(hist_tapes, data_model, rout_var, time_handle,
                    directories):
    '''
    Main run loop for RVIC model.

    Parameters
    ----------
    hist_tapes : OrderedDict
        Ordered dictionary of History objects
    data_model : DataModel
        DataModel instance containing the forcings
    rout_var : Rvar
        Rvar instance containing the RVIC parameters / unit hydrographs
    time_handle : Dtime
        Dtime instance containing information about run length, time
    directories : dict
        Dictionary of directories created by this function.

    Returns
    ----------
    time_handle : Dtime
        Dtime instance containing information about run length, time
    hist_tapes : OrderedDict
        Ordered dictionary of History objects
    '''

    data2tape = {}
    aggrunin = {}
    for tracer in RVIC_TRACERS:
        aggrunin[tracer] = 0.0
    aggcounter = 0

    # ---------------------------------------------------------------- #
    # Start log
    log = getLogger(LOG_NAME)
    log.info('Starting convolution_run')
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
        # Get this time_handlesteps forcing
        runin = data_model.read(timestamp)

        for tracer in RVIC_TRACERS:
            aggrunin[tracer] += runin[tracer]
        aggcounter += 1

        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Do the Convolution
        # (end_timestamp is the timestamp at the end of the convolution period)
        if aggcounter == rout_var.agg_tsteps:
            end_timestamp = rout_var.convolve(aggrunin, time_ord)
            for tracer in RVIC_TRACERS:
                aggrunin[tracer][:] = 0.0
                aggcounter = 0
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Extract the Current Variables from rout_var
        data2tape['streamflow'] = rout_var.get_rof()
        data2tape['storage'] = rout_var.get_storage()

        # Update the history Tape(s)
        for tapename, tape in iteritems(hist_tapes):
            log.debug('Updating Tape: %s', tapename)
            tape.update(data2tape, time_ord)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write State
        if rest_flag:
            # History files
            history_files = []
            history_restart_files = []
            for tapename, tape in iteritems(hist_tapes):
                log.debug('Writing Restart File for Tape: %s', tapename)
                # hist_fname, rest_fname = tape.write_restart()
                history_files.append(tape.filename)
                history_restart_files.append(tape.rest_filename)

            restart_file = rout_var.write_restart(history_files,
                                                  history_restart_files)
            write_rpointer(directories['restarts'], restart_file,
                           end_timestamp)

        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # advance a single timestep
        if not stop_flag:
            timestamp, time_ord, stop_flag, \
                rest_flag = time_handle.advance_timestep()
            # check that we're still inline with convolution
            if end_timestamp != timestamp:
                raise ValueError('timestamps do not match after convolution')
        else:
            break
        # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #
    # Make sure we write out the last history file
    for tapename, tape in iteritems(hist_tapes):
        log.debug('Closing Tape: %s', tapename)
        tape.finish()
    # ---------------------------------------------------------------- #
    return time_handle, hist_tapes
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Final
def convolution_final(time_handle, hist_tapes):
    '''Finalize RVIC Convolution

    Parameters
    ----------
    time_handle : Dtime
        Dtime instance containing information about run length, time
    hist_tapes : OrderedDict
        Ordered dictionary of History objects
    '''
    # ---------------------------------------------------------------- #
    # Start log
    log = getLogger(LOG_NAME)
    log.info('Finalizing RVIC')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write final log info
    log.info('-----------------------------------------------------------')
    log.info('Done with streamflow convolution')
    log.info('Processed %i timesteps', time_handle.timesteps)
    for name, tape in iteritems(hist_tapes):
        log.info('Wrote %i history files from %s', tape.files_count, name)
    log.info('-----------------------------------------------------------')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    log_tar = tar_inputs(log.filename)

    log.info('Done with rvic convolution.')
    log.info('Location of Log: %s', log_tar)

    close_logger()
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #
