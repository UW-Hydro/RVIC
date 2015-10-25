from rvic.core.config import config_type, isfloat, isint


def test_config_type_int():
    assert config_type('1') == 1


def test_config_type_float():
    assert config_type('1.75') == 1.75


def test_config_type_bool():
    assert config_type('True')


def test_isfloat():
    assert isfloat(4.3)
    assert isfloat('4.3')


def test_isint():
    assert isint(4)
    assert isint('4')
    assert not isint(4.3)
    assert not isint('4.3')
