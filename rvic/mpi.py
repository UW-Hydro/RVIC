"""
mpi.py
"""

from log import LogExceptions
from multiprocessing.pool import Pool

class LoggingPool(Pool):
    def apply_async(self, func, args=(), kwds={}, callback=None):
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)
