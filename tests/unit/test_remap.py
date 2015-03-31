import pytest
from cdo import CDOException
from rvic.core.remap import remap


def test_cdo_raises_exception():
    with pytest.raises(CDOException):
        remap('grid_file', 'in_file', 'out_file')
