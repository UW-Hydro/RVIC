import pytest
import numpy as np
from rvic.core.variables import Rvar
from rvic.core.history import Tape


@pytest.fixture()
def rvar(scope='function'):
    rv = Rvar(
        'unit_test_data/gunnison_parameters_01.rvic.prm.BLMSA.20150226.nc',
        'test_case', 'noleap', './', 'NETCDF4')
    return rv


def test_create_tape_instance(rvar):
    history_tape = Tape(1.25, 'test', rvar, grid_area=np.zeros((10, 11)),
                        outtype='array')
