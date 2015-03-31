import pytest
import numpy as np
import pandas as pd
from rvic.core.aggregate import make_agg_pairs, aggregate


config_dict = {}
config_dict['DOMAIN'] = {'LONGITUDE_VAR': 'LONGITUDE',
                         'LATITUDE_VAR': 'LATITUDE'}
config_dict['ROUTING'] = {'LONGITUDE_VAR': 'LONGITUDE',
                          'LATITUDE_VAR': 'LATITUDE',
                          'SOURCE_AREA_VAR': 'SOURCE_AREA'}


@pytest.fixture()
def dom_data(scope='function'):
    dom = {}
    lons, lats = np.meshgrid(np.linspace(-10, 0., 20), np.linspace(40, 50, 20))
    dom[config_dict['DOMAIN']['LONGITUDE_VAR']] = lons
    dom[config_dict['DOMAIN']['LATITUDE_VAR']] = lats
    dom['cell_ids'] = np.arange(lons.size).reshape(lons.shape)
    return dom


@pytest.fixture()
def fdr_data(scope='function'):
    fdr = {}
    fdr['resolution'] = 1.0
    lons, lats = np.meshgrid(np.arange(-15, 15, fdr['resolution']),
                             np.arange(20, 60, fdr['resolution']))
    fdr[config_dict['ROUTING']['LONGITUDE_VAR']] = lons
    fdr[config_dict['ROUTING']['LATITUDE_VAR']] = lats
    fdr[config_dict['ROUTING']['SOURCE_AREA_VAR']] = np.ones_like(lons)
    return fdr


@pytest.fixture()
def in_data(fdr_data, scope='function'):
    rout_data = {}
    rout_data['lat'] = fdr_data[config_dict['ROUTING']['LATITUDE_VAR']]
    rout_data['lon'] = fdr_data[config_dict['ROUTING']['LONGITUDE_VAR']]
    shape = fdr_data[config_dict['ROUTING']['LATITUDE_VAR']].shape
    rout_data['unit_hydrograph'] = np.random.random((10, ) + shape)
    rout_data['fraction'] = np.ones(shape)
    rout_data['unit_hydrograph_dt'] = 3600.
    return rout_data


@pytest.fixture()
def agg_data(fdr_data, scope='function'):
    rout_data = {}
    rout_data['lat'] = fdr_data[config_dict['ROUTING']['LATITUDE_VAR']][:4]
    rout_data['lon'] = fdr_data[config_dict['ROUTING']['LONGITUDE_VAR']][:4]
    shape = rout_data['lat'].shape
    rout_data['unit_hydrograph'] = np.random.random((10, ) + shape)
    rout_data['fraction'] = np.ones(shape)
    rout_data['unit_hydrograph_dt'] = 3600.
    return rout_data


def test_make_agg_pairs_all_unique(dom_data, fdr_data):
    pour_points = pd.DataFrame({'lons': [-2.3, 0.3, 10.4, 17.8],
                                'lats': [42.1, 40.5, 49.0, 45.2]})
    outlets = make_agg_pairs(pour_points, dom_data, fdr_data, config_dict)
    assert len(outlets) == len(pour_points)


def test_make_agg_pairs_with_overlap(dom_data, fdr_data):
    pour_points = pd.DataFrame({'lons': [-2.3, 0.3, 10.4, 17.8, -2.31],
                                'lats': [42.1, 40.5, 49.0, 45.2, 42.11]})
    outlets = make_agg_pairs(pour_points, dom_data, fdr_data, config_dict)
    assert len(outlets) == len(pour_points) - 1


def test_aggregate_no_agg_data(in_data, fdr_data):
    agg_data = aggregate(in_data, {}, res=fdr_data['resolution'],
                         pad=0, maskandnorm=False)
    assert all([k in agg_data for k in in_data])


def test_aggregate_with_aggdata(in_data, agg_data, fdr_data):
    agg_data = aggregate(in_data, agg_data, res=fdr_data['resolution'],
                         pad=0, maskandnorm=False)
    assert all([k in agg_data for k in in_data])
