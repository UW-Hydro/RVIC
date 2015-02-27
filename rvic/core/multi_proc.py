"""
multi_proc.py
"""

from .log import LOG_NAME
import multiprocessing
from multiprocessing.pool import Pool
import traceback


def error(msg, *args):
    """ Error function"""
    return multiprocessing.get_logger(LOG_NAME).error(msg, *args)
# -------------------------------------------------------------------- #


class LogExceptions(object):
    def __init__(self, callable):
        self.__callable = callable
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
    pass
# -------------------------------------------------------------------- #


class LoggingPool(Pool):
    """Subclass of pool"""
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds,
                                callback)
# -------------------------------------------------------------------- #
