import pytest
try:
    from rvic.core.remap import remap
    from cdo import CDOException
    cdo_unavailable = False
except ImportError:
    cdo_unavailable = True


@pytest.mark.skipif(cdo_unavailable, reason='cdo not installed')
def test_cdo_raises_exception():
    with pytest.raises(CDOException):
        remap('grid_file', 'in_file', 'out_file')
