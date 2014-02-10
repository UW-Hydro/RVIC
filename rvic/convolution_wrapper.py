"""
convolution_wrapper.py

ctypes wrapper for convolve.c

gcc -shared -o c_convolve.so c_convolve.c
"""

import numpy as np
import ctypes

import os
# os.system('gcc -shared -o c_convolve.so c_convolve.c')

try:
    path_to_file = os.path.split(__file__)[0]

    _convolution = np.ctypeslib.load_library('c_convolve.so', path_to_file)
    args = [ctypes.c_int,
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
    _convolution.c_convolve.argtypes = args
    _convolution.c_convolve.restype = None

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
        _convolution.c_convolve(*args)

        return
except:
    print('Using brodcasting convolution method because there was a problem '
          'loading c_convolve.c')

    def rvic_convolve(nsources, noutlets, subset_length, xsize,
                      source2outlet_ind, source_y_ind, source_x_ind,
                      source_time_offset, unit_hydrograph, aggrunin, ring):
        # numpy brodcasting method
        for s, outlet in enumerate(source2outlet_ind):
            y = source_y_ind[s]
            x = source_x_ind[s]
            i = source_time_offset[s]
            j = i + subset_length
            ring[i:j, outlet] += np.squeeze(unit_hydrograph[:, s] * aggrunin[y, x])

        # # pure python convolution
        # for s, outlet in enumerate(source2outlet_ind):
        #     y = source_y_ind[s]
        #     x = source_x_ind[s]
        #     for i in xrange(subset_length):
        #         j = i + source_time_offset[s]
        #         ring[j, outlet] += (unit_hydrograph[i, s] * aggrunin[y, x])
        return


def test():
    nsources = 20
    subset_length = 10
    full_time_length = 15
    noutlets = 4
    source2outlet_ind = np.linspace(0, noutlets, nsources,
                                    endpoint=False).astype(np.int32)

    ysize = 12
    xsize = 15

    source_y_ind = np.random.randint(0, ysize-1, nsources).astype(np.int32)
    source_x_ind = np.random.randint(0, xsize-1, nsources).astype(np.int32)

    source_time_offset = np.random.randint(0, 4, nsources).astype(np.int32)

    aggrunin = np.random.uniform(size=(ysize, xsize))
    unit_hydrograph = np.zeros((subset_length, nsources), dtype=np.float64)
    unit_hydrograph[0:4, :] = 0.25
    ring = np.zeros((full_time_length, noutlets), dtype=np.float64)

    for i in xrange(10):
        aggrunin = np.random.uniform(size=(ysize, xsize))
        rvic_convolve(nsources, noutlets, subset_length, xsize,
                      source2outlet_ind, source_y_ind, source_x_ind,
                      source_time_offset, unit_hydrograph, aggrunin, ring)
