import pytest
import os
import numpy as np
from rvic.core.variables import Rvar
from rvic.core.history import Tape


@pytest.mark.skipif(os.getenv('TRAVIS', 'false').lower() == 'true',
                    reason='on travis')
@pytest.fixture()
def rvar(scope='function'):
    dirname = os.path.dirname(__file__)
    infile = os.path.join(dirname, 'unit_test_data',
                          'gunnison_parameters_01.rvic.prm.BLMSA.20150226.nc')
    rv = Rvar(infile, 'test_case', 'noleap', dirname, 'NETCDF4')
    return rv


@pytest.mark.skipif(os.getenv('TRAVIS', 'false').lower() == 'true',
                    reason='on travis')
def test_create_tape_instance(rvar):
    history_tape = Tape(1.25, 'test', rvar, grid_area=np.zeros((10, 11)),
                        outtype='array')
