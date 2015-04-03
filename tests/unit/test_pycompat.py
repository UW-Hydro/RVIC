# -*- coding: utf-8 -*-
import pytest
from rvic.core.pycompat import (pyrange, pyzip, iteritems, itervalues,
                                OrderedDict, SafeConfigParser)


def test_pyrange():
    x = pyrange(10)
    assert hasattr(x, '__iter__')
    assert len(x) == 10
    assert x[9] == 9


def test_pyzip():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']
    x = pyzip(a, b)
    assert hasattr(x, '__iter__')
    for i, (aa, bb) in enumerate(x):
        assert aa == a[i]
        assert bb == b[i]
    with pytest.raises(TypeError):
        len(x)
        x[1]


def test_iteritems():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']
    d = dict(pyzip(a, b))
    for k, v in iteritems(d):
        assert k in a
        assert v in b
        assert d[k] == v


def test_iteritems_ordered_dict():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']
    d = OrderedDict(pyzip(a, b))
    for i, (k, v) in enumerate(iteritems(d)):
        assert k in a
        assert v in b
        assert d[k] == v
        assert b[i] == v


def test_itervalues():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']
    d = dict(pyzip(a, b))
    for v in itervalues(d):
        assert v in b


def test_itervalues_ordered_dict():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']
    d = OrderedDict(pyzip(a, b))
    for i, v in enumerate(itervalues(d)):
        assert v in b
        assert b[i] == v


def test_ordered_dict_order():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']
    d = OrderedDict(pyzip(a, b))
    assert list(d.keys()) == a
    assert list(d.values()) == b


def test_ordered_dict_append():
    a = [1, 2, 3]
    b = ['a', 'b', 'c']
    d = OrderedDict(pyzip(a, b))
    d[4] = 'd'
    assert len(d) == 4
    assert list(d.keys())[-1] == 4


def test_init_safe_config_parser():
    SafeConfigParser()
