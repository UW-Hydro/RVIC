# -*- coding: utf-8 -*-
'''
multi_proc.py
'''

from .log import LOG_NAME
import multiprocessing
from multiprocessing.pool import Pool
import traceback


def error(*args):
    """Error function"""
    return multiprocessing.get_logger(LOG_NAME).error(*args)
# -------------------------------------------------------------------- #


class LogExceptions(object):
    def __init__(self, func):
        self.__callable = func
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise

        # It was fine, give a normal answer
        return result
# -------------------------------------------------------------------- #


class LoggingPool(Pool):
    __doc__ = 'RVIC Logging Subclass of pool\n' + Pool.__doc__

    @property
    def has_logging(self):
        # just for testing
        return True

    def apply_async(self, func, callback=None,  error_callback=None,
                    *args, **kwds):
        """Overloaded Pool.apply_async to support Logging"""
        return Pool.apply_async(self, func, args=args, kwds=kwds,
                                callback=callback,
                                error_callback=error_callback)

# -------------------------------------------------------------------- #
