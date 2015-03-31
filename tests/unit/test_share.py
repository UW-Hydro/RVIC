from rvic.core.share import NcGlobals, NcVar


def test_nc_globals():
    ncg = NcGlobals()
    ncg.update()


def test_nc_var():
    ncv = NcVar(test='test')
    assert ncv.test == 'test'
