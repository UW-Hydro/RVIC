# -*- coding: utf-8 -*-
'''
log.py
'''
import os
import sys
from time import gmtime, strftime
import logging

# -------------------------------------------------------------------- #
LOG_NAME = 'rvic'
FORMATTER = logging.Formatter('%(levelname)s:%(funcName)s>> %(message)s')
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Fake stream file to handle stdin/stdout
class StreamToFile(object):
    '''
    Fake file-like stream object that redirects writes to a logger instance.
    http://www.electricmonk.nl/log/2011/08/14/redirect-stdout-and-stderr-to-a-logger-in-python/
    '''
    def __init__(self, logger_name=LOG_NAME, log_level=logging.INFO):
        self.logger = logging.getLogger(logger_name)
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass
# -------------------------------------------------------------------- #


def init_logger(log_dir='./', log_level='DEBUG', verbose=False):
    ''' Setup the logger '''

    logger = logging.getLogger(LOG_NAME)
    logger.setLevel(log_level)
    logger.propograte = True

    # ---------------------------------------------------------------- #
    # create log file handler
    if log_dir:
        log_file = os.path.join(log_dir, 'RVIC-' + strftime('%Y%m%d-%H%M%S',
                                gmtime()) + '.log')
        fh = logging.FileHandler(log_file)
        fh.setLevel(log_level)
        fh.setFormatter(FORMATTER)
        logger.addHandler(fh)
        logger.filename = log_file
    else:
        log_file = None
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # If verbose, logging will also be sent to console
    if verbose:
        # print to console
        ch = logging.StreamHandler()
        ch.setLevel(log_level)
        ch.setFormatter(FORMATTER)
        logger.addHandler(ch)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Redirect stdout and stderr to logger
    sys.stdout = StreamToFile()
    sys.stderr = StreamToFile(log_level=logging.ERROR)
    # ---------------------------------------------------------------- #

    logger.info('-------------------- INITIALIZED RVIC LOG ------------------')
    logger.info('LOG LEVEL: %s', log_level)
    logger.info('Logging To Console: %s', verbose)
    logger.info('LOG FILE: %s', log_file)
    logger.info('----------------------------------------------------------\n')

    return logger
# -------------------------------------------------------------------- #


def close_logger():
    '''Close the handlers of the logger'''
    log = logging.getLogger(LOG_NAME)
    x = list(log.handlers)
    for i in x:
        log.removeHandler(i)
        i.flush()
        i.close()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
# -------------------------------------------------------------------- #
