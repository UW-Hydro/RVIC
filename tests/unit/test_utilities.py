# -*- coding: utf-8 -*-
import numpy as np
import os
import datetime
from rvic.core.utilities import (latlon2yx, search_for_channel, find_nearest,
                                 write_rpointer, strip_non_ascii,
                                 strip_invalid_char)
from rvic.core.share import RPOINTER


def test_latlon2yx():
    glats, glons = np.meshgrid(np.arange(0, 20, 2), np.arange(10, 30, 2))
    plats = np.array([0, 2, 17])
    plons = np.array([28, 17, 12])
    y, x = latlon2yx(plats, plons, glats, glons)
    assert len(y) == len(x) == len(plats)
    assert y[0] == 9
    assert x[0] == 0


def test_search_for_channel():
    source_area = np.array([[0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00],
                            [0.00, 0.00, 0.00, 2.00, 0.00, 0.00, 0.00, 0.00],
                            [0.00, 0.00, 0.00, 3.00, 0.00, 0.00, 8.00, 9.00],
                            [51.0, 52.0, 53.0, 0.00, 4.00, 0.00, 7.00, 10.0],
                            [50.0, 1.00, 54.0, 0.00, 0.00, 5.00, 6.00, 11.0],
                            [1.00, 77.0, 1.00, 0.00, 0.00, 0.00, 12.0, 0.00],
                            [78.0, 20.0, 19.0, 0.00, 0.00, 13.0, 0.00, 0.00],
                            [79.0, 21.0, 18.0, 0.00, 0.00, 14.0, 0.00, 0.00],
                            [80.0, 22.0, 0.00, 17.0, 16.0, 15.0, 0.00, 0.00]])

    routys = np.array([4, 8, 0, 1])
    routxs = np.array([1, 0, 3, 7])

    new_ys, new_xs = search_for_channel(source_area, routys, routxs, search=1)

    print(new_ys, new_xs)

    np.testing.assert_equal(new_ys, [5, 8, 0, 2])
    np.testing.assert_equal(new_xs, [1, 0, 3, 7])
    np.testing.assert_equal(source_area[new_ys, new_xs], [77., 80., 1., 9.])


def test_write_rpointer():
    write_rpointer('./', 'test_restart_file.nc', datetime.datetime.now())
    testfile = os.path.join('./', RPOINTER)
    assert os.path.isfile(testfile)
    os.remove(testfile)


def test_find_nearest():
    assert find_nearest(np.array([8, 19, 39, 100, 399]), 20) == 1


def test_find_nearest_max():
    assert find_nearest(np.array([8, 19, 39, 100, 399]), 500) == 4


def test_find_nearest_min():
    assert find_nearest(np.array([8, 19, 39, 100, 399]), 1) == 0


def test_strip_non_ascii():
    test = u'éáé123456tgreáé@€'
    new = strip_non_ascii(test)
    assert new == '123456tgre@'
    test = 'standard'
    new = strip_non_ascii(test)
    assert new == test


def test_strip_invalid_char():
    test = '123$%^789'
    new = strip_invalid_char(test)
    assert new == '123789'
