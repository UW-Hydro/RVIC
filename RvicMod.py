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
import argparse
import os
import logging
from rvic.log import init_logger, log_name
from rvic.utilities import ReadConfig, MakeDirs, CopyInputs, Rvar, ReadDomain
from rvic.utilities import write_rpointer, TarInputs
from rvic.time_util import Dtime
from rvic.read_forcing import DataModel
from rvic.history import Tape
from rvic.share import ncGlobals


# -------------------------------------------------------------------- #
# Top Level Driver
def main(config_file=None):
    """
    Top level driver for RVIC model.
    """

    HistTapes, dataMod, routVar, DomData, Time, dirPaths, ConfigDict = rvic_init()

    Time, HistTapes = rvic_run(HistTapes, dataMod, routVar, DomData, Time, dirPaths, ConfigDict)

    rvic_final(Time, HistTapes)

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Initialize RVIC
def rvic_init(configFile=None):
    """
    - Read Grid File
    - Load the unit hydrograph files (put into point_dict)
    - Load the initial state file and put it in convolution rings
    """

    # ---------------------------------------------------------------- #
    # Read command Line
    if not configFile:
        configFile = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Configuration files
    ConfigDict = ReadConfig(configFile)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structure
    dirPaths = MakeDirs(ConfigDict['options']['case_dir'],
                        ['hist', 'logs', 'params', 'restarts'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Copy Inputs to $case_dir/inputs and update configuration
    ConfigDict = CopyInputs(configFile, dirPaths['inputs'])
    options = ConfigDict['options']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Settup Logging
    log = init_logger(dirPaths['logs'], options['log_level'], options['verbose'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read Domain File
    Domain = ConfigDict['domain']
    DomData, DomVats, DomGats = ReadDomain(Domain['file_name'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read the Parameter File
    log.info('reading parameter and state file')

    routVar = Rvar(ConfigDict['PARAM_FILE']['file_name'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Read the State File
    if options['RUN_TYPE'] == 'restart':
        restart = ReadConfig(os.path.join(dirPaths['restarts'], 'rpointer'))
        timestr = restart['restart']['timestamp']
        routVar.initial_state(restart['restart']['file_name'])
    elif options['RUN_TYPE'] == 'startup':
        timestr = options['RUN_STARTDATE']
        routVar.initial_state(ConfigDict['INITIAL_STATE']['file_name'])
    elif options['RUN_TYPE'] == 'drystart':
        timestr = options['RUN_STARTDATE']
    else:
        log.error('RUN_TYPE Option is none of these: (restart, startup, drystart)')
        raise
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Time
    Time = Dtime(timestr, options['STOP_OPTION'], options['STOP_N'],
                 options['STOP_DATE'], options['REST_OPTION'],
                 options['REST_N'], options['REST_DATE'],
                 options['CALENDAR'], routVar.uh_timestep)

    # ---------------------------------------------------------------- #
    # Initialize the data model
    Forcings = ConfigDict['INPUT_FORCINGS']
    dataMod = DataModel(Forcings['datl_path'],
                        Forcings['datl_file'],
                        Forcings['time_var'],
                        Forcings['datl_liq_flds'],
                        Forcings['start'],
                        Forcings['end'],
                        Time.timestamp)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup History Tape(s) and Write Initial Outputs
    History = ConfigDict['HISTORY']
    HistTapes = {}

    for j in xrange(ConfigDict['HISTORY']['rvichist_ntapes']):
        tapename = 'Tape.%i' % j
        HistTapes[tapename] = Tape(Time.timestamp,
                                   Time.time_ord,
                                   options['CASEID'],
                                   routVar,
                                   fincl=['streamflow'],
                                   mfilt=History['rvichist_mfilt'][j],
                                   ndens=History['rvichist_ndens'][j],
                                   nhtfrq=History['rvichist_nhtfrq'][j],
                                   avgflag=History['rvichist_avgflag'][j],
                                   file_format=History['rvichist_ncform'][j],
                                   outtype=History['rvichist_outtype'][j],
                                   dom_lons=DomData[Domain['longitude_var']],
                                   dom_lats=DomData[Domain['latitude_var']],
                                   out_dir=dirPaths['hist'],
                                   calendar=Time.calendar,
                                   GlobAts=ncGlobals(title='RVIC history file',
                                                     casename=options['CASEID'],
                                                     casestr=options['CASESTR'],
                                                     RvicPourPointsFile=os.path.split(routVar.RvicPourPointsFile)[1],
                                                     RvicUHFile=os.path.split(routVar.RvicUHFile)[1],
                                                     RvicFdrFile=os.path.split(routVar.RvicFdrFile)[1],
                                                     RvicDomainFile=os.path.split(Domain['file_name'])[1]))
        HistTapes[tapename].write_initial()
    # ---------------------------------------------------------------- #

    return HistTapes, dataMod, routVar, DomData, Time, dirPaths, ConfigDict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def rvic_run(HistTapes, dataMod, routVar, DomData, Time, dirPaths, ConfigDict):
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
    log = logging.getLogger(log_name)
    log.info('Starting rvic_run')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Iterate Through Timesteps
    while not Time.stop_flag:
        timestamp = Time.advance_timestep()
        time_ord = Time.time_ord

        # ------------------------------------------------------------ #
        # Get This Timesteps Forcing
        aggrunin = dataMod(timestamp)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Do the Convolution
        routVar.convolve(aggrunin, timestamp)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Extract the Current Variables from routVar
        data2tape['streamflow'] = routVar.get_rof()
        data2tape['storage'] = routVar.get_storage()

        # Update the History Tape(s)
        for tapename, tape in HistTapes.iteritems():
            log.debug('Updating Tape:%s' % tapename)
            tape.update(time_ord, data2tape)
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------ #
        # Write State
        if Time.rest_flag or Time.stop_flag:
            restart_file = routVar.write_restart()

            write_rpointer(dirPaths['restarts'], restart_file, timestamp)
        # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #
    return Time, HistTapes
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Final
def rvic_final(Time, HistTapes):
    """ Finalize RVIC """
    # ---------------------------------------------------------------- #
    # Start log
    log = logging.getLogger(log_name)
    log.info('Finalizing RVIC')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write final log info
    log.info("-----------------------------------------------------------")
    log.info('Done with streamflow convolution')
    log.info('Processed %i timesteps' % Time.timesteps)
    for name, tape in HistTapes.iteritems():
        log.info('Wrote %i history files from %s' % (tape.files_count, name))
    log.info('Routed to %i points' % Time.points)
    log.info("-----------------------------------------------------------")
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # tar the inputs directory / log file
    LogTar = TarInputs(log.filename)

    log.info('Done with RvicGenParam.')
    log.info('Location of Log: %s' % LogTar)
    # ---------------------------------------------------------------- #
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Process Command Line
def process_command_line():
    """
    Get the path to the configFile
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description='Do the RVIC Convoluition')
    parser.add_argument("configFile", type=str, help="Input configuration file")

    args = parser.parse_args()
    configFile = args.configFile

    return configFile
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Run Program
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
