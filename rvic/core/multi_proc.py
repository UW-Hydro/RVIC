"""
multi_proc.py
"""

from .log import LOG_NAME
import multiprocessing
from multiprocessing.pool import Pool
import traceback


def error(*args):
    """ Error function"""
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
    """Subclass of pool"""
    def apply_async(self, func, callback=None, *args, **kwargs):
        return Pool.apply_async(self, LogExceptions(func), args, kwargs,
                                callback)
# -------------------------------------------------------------------- #
