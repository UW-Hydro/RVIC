#!/usr/local/bin/python

import numpy as np
from rout import *
import unittest

class TestConvolutionFuctions(unittest.TestCase):

    def test_blank():
    	pass
    
suite = unittest.TestLoader().loadTestsFromTestCase(TestConvolutionFuctions)
unittest.TextTestRunner(verbosity=2).run(suite)
