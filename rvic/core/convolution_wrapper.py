"""
convolution_wrapper.py

ctypes wrapper for rvic_convolution.c

gcc -shared -o rvic_convolution.so rvic_convolution.c
"""

import os
import numpy as np
import ctypes

SHAREDOBJECT = 'rvic_convolution.so'
LIBPATH = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))

try:
    _convolution = np.ctypeslib.load_library(SHAREDOBJECT, LIBPATH)
except ImportError as ie:
    print('looking for shared object {0} in {1}'.format(SHAREDOBJECT, LIBPATH))
    raise ImportError(ie)
except OSError as oe:
    print('looking for shared object {0} in {1}'.format(SHAREDOBJECT, LIBPATH))
    raise ImportError(oe)

_args = [ctypes.c_int,
         ctypes.c_int,
         ctypes.c_int,
         ctypes.c_int,
         np.ctypeslib.ndpointer(np.int32),
         np.ctypeslib.ndpointer(np.int32),
         np.ctypeslib.ndpointer(np.int32),
         np.ctypeslib.ndpointer(np.int32),
         np.ctypeslib.ndpointer(np.float64),
         np.ctypeslib.ndpointer(np.float64),
         np.ctypeslib.ndpointer(np.float64)]
_convolution.convolve.argtypes = _args
_convolution.convolve.restype = None


def rvic_convolve(*args):
    """args:

       nsources,
       noutlets,
       subset_length,
       xsize,
       source2outlet_ind,
       source_y_ind,
       source_x_ind,
       source_time_offset,
       unit_hydrograph,
       aggrunin,
       ring
    """
    _convolution.convolve(*args)

    return
