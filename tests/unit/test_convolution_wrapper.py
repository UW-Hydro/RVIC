import pytest
import numpy as np
from rvic.core.convolution_wrapper import rvic_convolve


def test_convolution_wraper():
    n_sources = 4
    n_outlets = 2
    time_length = 6
    ysize = 10
    xsize = 10
    source2outlet_ind = np.array([0, 0, 1, 1], dtype=np.int32)
    source_y_ind = np.array([2, 4, 6, 8], dtype=np.int32)
    source_x_ind = np.array([2, 6, 3, 5], dtype=np.int32)
    source_time_offset = np.zeros(n_sources, dtype=np.int32)
    unit_hydrograph = np.array([[0., 1., 0.5, 0.25],
                                [1., 0., 0.5, 0.25],
                                [0., 0., 0., 0.25],
                                [0., 0., 0., 0.25],
                                [0., 0., 0., 0.],
                                [0., 0., 0., 0.]])
    aggrunin = np.ones((ysize, xsize), dtype=np.float64)
    ring = np.zeros((time_length, n_outlets))

    rvic_convolve(n_sources,
                  n_outlets,
                  time_length,
                  xsize,
                  source2outlet_ind,
                  source_y_ind,
                  source_x_ind,
                  source_time_offset,
                  unit_hydrograph,
                  aggrunin,
                  ring)
    assert ring.sum() == 4
    np.testing.assert_almost_equal(ring.sum(axis=1),
                                   unit_hydrograph.sum(axis=1))
    np.testing.assert_almost_equal(ring.sum(axis=0),
                                   [2, 2])


def test_convolution_wraper_1_outlet():
    n_sources = 4
    n_outlets = 1
    time_length = 6
    ysize = 10
    xsize = 10
    source2outlet_ind = np.array([0, 0, 0, 0], dtype=np.int32)
    source_y_ind = np.array([2, 4, 6, 8], dtype=np.int32)
    source_x_ind = np.array([2, 6, 3, 5], dtype=np.int32)
    source_time_offset = np.zeros(n_sources, dtype=np.int32)
    unit_hydrograph = np.array([[0., 1., 0.5, 0.25],
                                [1., 0., 0.5, 0.25],
                                [0., 0., 0., 0.25],
                                [0., 0., 0., 0.25],
                                [0., 0., 0., 0.],
                                [0., 0., 0., 0.]])
    aggrunin = np.ones((ysize, xsize), dtype=np.float64)
    ring = np.zeros((time_length, n_outlets))

    rvic_convolve(n_sources,
                  n_outlets,
                  time_length,
                  xsize,
                  source2outlet_ind,
                  source_y_ind,
                  source_x_ind,
                  source_time_offset,
                  unit_hydrograph,
                  aggrunin,
                  ring)
    assert ring.sum() == 4
    np.testing.assert_almost_equal(ring.sum(axis=1),
                                   unit_hydrograph.sum(axis=1))
    np.testing.assert_almost_equal(ring.sum(axis=0),
                                   [4])


def test_convolution_wraper_1_source_1_outlet():
    n_sources = 1
    n_outlets = 1
    time_length = 6
    ysize = 10
    xsize = 10
    source2outlet_ind = np.array([0, 0, 1, 1], dtype=np.int32)
    source_y_ind = np.array([2], dtype=np.int32)
    source_x_ind = np.array([5], dtype=np.int32)
    source_time_offset = np.zeros(n_sources, dtype=np.int32)
    unit_hydrograph = np.array([[0., 0.25, 0.5, 0.25, 0., 0.]]).T
    aggrunin = np.ones((ysize, xsize), dtype=np.float64)
    ring = np.zeros((time_length, n_outlets))

    rvic_convolve(n_sources,
                  n_outlets,
                  time_length,
                  xsize,
                  source2outlet_ind,
                  source_y_ind,
                  source_x_ind,
                  source_time_offset,
                  unit_hydrograph,
                  aggrunin,
                  ring)
    assert ring.sum() == 1
    np.testing.assert_almost_equal(ring.sum(axis=1),
                                   unit_hydrograph.sum(axis=1))
    np.testing.assert_almost_equal(ring.sum(axis=0),
                                   [1])
