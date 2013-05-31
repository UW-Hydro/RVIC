#!/usr/local/bin/python

import numpy as np
from coup_conv import *
import unittest

class TestConvolutionFuctions(unittest.TestCase):

    def test_ring_shift(self):
        # Make sure the j+1 element ends up in the j location
        ring1 = np.random.random(10)
        ring2=shift(ring1,1)
        self.assertEqual(ring1[1],ring2[0])

    def test_ring_pop(self):
        # Make sure first element is replaced by a zero and moved to last element
        ring1 = np.arange(10)
        ring1[0] = 0
        ring2=shift(ring1,1)
        self.assertEqual(ring2[-1],0)
        
    def test_grid_convolve_length(self):
        # Make sure that doing a single timestep convolution
        # returns a hydrograph with length equal to the unit hydrograph
        points = np.random.randint(10,50)
        tlen = 500
        xs = np.random.randint(0,10,points)
        ys = np.random.randint(0,10,points)
        UH = np.random.randint(500,1000,(tlen,points))
        flux = np.random.random(size=(10,10))
        self.flow = flux[ys,xs]*UH
        self.assertEqual(len(self.flow),tlen)

    def test_grid_convolve_dims(self):
        # Make sure that doing a single timestep convolution
        # returns a hydrograph with 1 dimension
        points = np.random.randint(10,50)
        tlen = 500
        xs = np.random.randint(0,10,points)
        ys = np.random.randint(0,10,points)
        UH = np.random.randint(500,1000,(tlen,points))
        flux = np.random.random(size=(10,10))
        self.flow = (flux[ys,xs]*UH).sum(axis=1)
        self.assertEqual(self.flow.ndim,1)
        
suite = unittest.TestLoader().loadTestsFromTestCase(TestConvolutionFuctions)
unittest.TextTestRunner(verbosity=2).run(suite)
