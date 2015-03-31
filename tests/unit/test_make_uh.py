import pytest
import numpy as np
from rvic.core.make_uh import (find_ts, read_direction, search_catchment,
                               make_uh, make_grid_uh_river, make_grid_uh,
                               adjust_uh_timestep)
from rvic.core.variables import Point


@pytest.fixture()
def fdr_vic(scope='function'):
    a = [[0, 0, 0, 5, 5, 5, 6, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 4, 0, 4, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 5, 3, 5, 5, 6, 7, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 5, 4, 0, 5, 5, 5, 5, 6, 7, 5, 5, 5, 7, 5, 0, 0, 0, 0, 0, 0],
         [0, 0, 5, 4, 0, 5, 4, 5, 5, 7, 7, 5, 6, 7, 5, 6, 0, 0, 0, 0, 0, 0],
         [4, 0, 5, 7, 4, 0, 3, 5, 5, 7, 5, 5, 7, 4, 5, 6, 0, 0, 0, 0, 0, 0],
         [4, 0, 3, 4, 0, 4, 0, 4, 3, 5, 6, 7, 7, 4, 5, 6, 0, 0, 0, 0, 0, 0],
         [0, 4, 0, 4, 0, 3, 1, 2, 4, 0, 4, 5, 7, 5, 5, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 3, 4, 3, 1, 8, 4, 5, 4, 0, 5, 5, 3, 5, 6, 5, 0, 0, 0, 0, 0],
         [0, 0, 0, 3, 4, 5, 5, 4, 5, 4, 3, 5, 6, 7, 7, 5, 5, 7, 0, 0, 0, 0],
         [0, 0, 0, 4, 3, 5, 5, 3, 4, 0, 4, 3, 4, 3, 5, 5, 7, 7, 0, 0, 0, 0],
         [0, 0, 0, 4, 0, 5, 4, 0, 4, 3, 4, 0, 4, 5, 3, 5, 7, 7, 0, 0, 0, 0],
         [0, 0, 0, 0, 4, 0, 4, 0, 4, 0, 4, 5, 3, 6, 7, 7, 5, 5, 6, 0, 5, 7],
         [0, 0, 0, 0, 0, 4, 0, 3, 4, 3, 5, 3, 6, 7, 6, 5, 7, 7, 7, 7, 7, 7],
         [0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 5, 7, 5, 6, 7, 1, 7, 1, 7, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 3, 4, 0, 5, 5, 6, 7, 7, 7, 6, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 5, 7, 1, 7, 7, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 4, 5, 7, 7, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
    return np.array(a, dtype=np.int)


@pytest.fixture()
def fdr_vic_small(scope='function'):
    a = [[0, 0, 0, 5, 5],
         [0, 0, 0, 4, 0],
         [0, 0, 5, 3, 5],
         [0, 0, 5, 4, 0],
         [0, 0, 5, 4, 0]]
    return np.array(a, dtype=np.int)


@pytest.fixture()
def dy_vic(scope='module'):
    dy = {1: -1, 2: -1, 3: 0, 4: 1, 5: 1, 6: 1, 7: 0, 8: -1}
    return dy


@pytest.fixture()
def dx_vic(scope='module'):
    dx = {1: 0, 2: 1, 3: 1, 4: 1, 5: 0, 6: -1, 7: -1, 8: - 1}
    return dx


@pytest.fixture()
def pour_point(scope='module'):
    p = Point()
    p.basinx = 4
    p.basiny = 2
    return p


def test_find_ts():
    assert find_ts(np.array([0, 86400, 172800])) == 86400


def test_find_ts_raises_when_scalar():
    with pytest.raises(TypeError):
        find_ts(4)


def test_read_direction(fdr_vic, dy_vic, dx_vic):
    to_y, to_x = read_direction(fdr_vic, dy_vic, dx_vic)
    np.testing.assert_equal(to_y.shape, to_x.shape)
    assert to_y.max() <= fdr_vic.shape[0] + 1
    assert to_x.max() <= fdr_vic.shape[1] + 1


def test_search_catchment(fdr_vic_small, dy_vic, dx_vic, pour_point):
    basin_ids = np.ones_like(fdr_vic_small, dtype=np.int)
    basin_id = 1
    to_y, to_x = read_direction(fdr_vic_small, dy_vic, dx_vic)
    catchment, catch_fracs = search_catchment(to_y, to_x, pour_point,
                                              basin_ids, basin_id)
    assert catch_fracs.min() <= 0.
    assert catch_fracs.max() == 1.
    assert type(catchment) == dict
    assert all([k in catchment for k in ['count_ds', 'x_inds', 'y_inds']])
    assert len(catchment['count_ds']) > 0
    assert len(catchment['count_ds']) == len(catchment['x_inds'])
    assert len(catchment['count_ds']) == len(catchment['y_inds'])


def test_make_uh(fdr_vic_small):
    ndays = 4
    y_inds, x_inds = np.nonzero(fdr_vic_small)
    velocity = np.zeros(fdr_vic_small.shape, dtype=np.float) + 2.
    diffusion = np.zeros(fdr_vic_small.shape, dtype=np.float) + 3000
    xmask = np.ones(fdr_vic_small.shape, dtype=np.float)
    uh = make_uh(86400, ndays, y_inds, x_inds, velocity, diffusion, xmask)
    assert uh.shape[0] == ndays
    assert uh.shape[1:] == fdr_vic_small.shape
    assert uh.min() >= 0.
    assert uh.max() <= 1.
    np.testing.assert_almost_equal(uh.sum(axis=0)[y_inds, x_inds], 1)


def test_make_grid_uh_river(fdr_vic_small, dy_vic, dx_vic, pour_point):
    ndays = 4
    t_uh = 40
    basin_ids = np.ones_like(fdr_vic_small, dtype=np.int)
    basin_id = 1
    to_y, to_x = read_direction(fdr_vic_small, dy_vic, dx_vic)
    uh = np.zeros((ndays, ) + fdr_vic_small.shape)
    uh[0, :, :] = 1.
    catchment, _ = search_catchment(to_y, to_x, pour_point,
                                    basin_ids, basin_id)
    uh_river = make_grid_uh_river(t_uh, ndays, uh, to_y, to_x, pour_point,
                                  catchment['y_inds'], catchment['x_inds'],
                                  catchment['count_ds'])
    assert uh_river.shape[0] == t_uh
    assert uh_river.max() <= 1.
    assert uh_river.min() >= 0.
    np.testing.assert_almost_equal(uh_river.sum(axis=0)[catchment['y_inds'],
                                                        catchment['x_inds']],
                                   1)


def test_make_grid_uh(fdr_vic_small, dy_vic, dx_vic, pour_point):
    ndays = 4
    t_uh = 40
    basin_ids = np.ones_like(fdr_vic_small, dtype=np.int)
    basin_id = 1
    to_y, to_x = read_direction(fdr_vic_small, dy_vic, dx_vic)
    catchment, _ = search_catchment(to_y, to_x, pour_point,
                                    basin_ids, basin_id)
    uh_river = np.zeros((t_uh, ) + fdr_vic_small.shape)
    uh_river[0] = 1.
    uh_box = np.array([1., 0, 0, 0])
    unit_hydrograph = make_grid_uh(t_uh, ndays, uh_river, uh_box, to_y, to_x,
                                   catchment['y_inds'], catchment['x_inds'],
                                   catchment['count_ds'])
    assert unit_hydrograph.shape[0] == t_uh
    assert unit_hydrograph.max() <= 1.
    assert unit_hydrograph.min() >= 0.
    np.testing.assert_almost_equal(
        unit_hydrograph.sum(axis=0)[catchment['y_inds'], catchment['x_inds']],
        1)


def test_adjust_uh_timestep_nochange(fdr_vic_small):
    t_uh = 40
    y_inds, x_inds = np.nonzero(fdr_vic_small)
    unit_hydrograph = np.zeros((t_uh, ) + fdr_vic_small.shape, dtype=np.float)
    unit_hydrograph[0] = 1.

    uh = adjust_uh_timestep(unit_hydrograph, t_uh, 3600, 3600, x_inds, y_inds)
    assert uh.shape[0] == t_uh
    np.testing.assert_almost_equal(
        uh.sum(axis=0)[y_inds, x_inds], 1)


def test_adjust_uh_timestep_downscale(fdr_vic_small):
    t_uh = 40
    y_inds, x_inds = np.nonzero(fdr_vic_small)
    unit_hydrograph = np.zeros((t_uh, ) + fdr_vic_small.shape, dtype=np.float)
    unit_hydrograph[0] = 1.

    uh = adjust_uh_timestep(unit_hydrograph, t_uh, 3600, 1800, x_inds, y_inds)
    np.testing.assert_almost_equal(
        uh.sum(axis=0)[y_inds, x_inds], 1)


def test_adjust_uh_timestep_upscale(fdr_vic_small):
    t_uh = 40
    y_inds, x_inds = np.nonzero(fdr_vic_small)
    unit_hydrograph = np.zeros((t_uh, ) + fdr_vic_small.shape, dtype=np.float)
    unit_hydrograph[0] = 1.

    uh = adjust_uh_timestep(unit_hydrograph, t_uh, 1800, 3600, x_inds, y_inds)
    np.testing.assert_almost_equal(
        uh.sum(axis=0)[y_inds, x_inds], 1)
