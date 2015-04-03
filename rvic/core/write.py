# -*- coding: utf-8 -*-
"""
write.py
"""

import numpy as np
from netCDF4 import Dataset, stringtochar
from logging import getLogger
from .log import LOG_NAME
from .share import NC_DOUBLE, NC_INT, NC_CHAR, FILLVALUE_F, NcGlobals
from .pycompat import iteritems
from . import share
# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write the agg netcdf
def write_agg_netcdf(file_name, agg_data, glob_atts, nc_format, zlib=True,
                     complevel=4, least_significant_digit=None):
    """
    Write output to netCDF.  Writes out a netCDF4 data file containing
    the UH_S and fractions and a full set of history and description
    attributes.
    """
    # ---------------------------------------------------------------- #
    # netCDF variable options
    ncvaropts = {'zlib': zlib,
                 'complevel': complevel,
                 'least_significant_digit': least_significant_digit}
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Open file
    f = Dataset(file_name, 'w', format=nc_format)
    log.info('writing aggregated netcdf: %s', file_name)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # set dimensions
    f.createDimension('timesteps', None)
    f.createDimension('lon', (len(agg_data['lon'])))
    f.createDimension('lat', (len(agg_data['lat'])))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # initialize variables
    unit_hydrograph_dt = f.createVariable('unit_hydrograph_dt', NC_INT, (),
                                          **ncvaropts)
    unit_hydrograph_dt[:] = agg_data['unit_hydrograph_dt']
    for key, val in iteritems(share.unit_hydrograph_dt.__dict__):
        if val:
            setattr(unit_hydrograph_dt, key, val)

    timesteps = f.createVariable('timesteps', NC_INT, ('timesteps',),
                                 **ncvaropts)
    timesteps[:] = np.arange(agg_data['unit_hydrograph'].shape[0])
    for key, val in iteritems(share.timesteps.__dict__):
        if val:
            setattr(timesteps, key, val)

    lon = f.createVariable('lon', NC_DOUBLE, ('lon',), **ncvaropts)
    for key, val in iteritems(share.lon.__dict__):
        if val:
            setattr(lon, key, val)

    lat = f.createVariable('lat', NC_DOUBLE, ('lat',), **ncvaropts)
    for key, val in iteritems(share.lat.__dict__):
        if val:
            setattr(lat, key, val)

    fraction = f.createVariable('fraction', NC_DOUBLE, ('lat', 'lon',),
                                **ncvaropts)
    for key, val in iteritems(share.fraction.__dict__):
        if val:
            setattr(fraction, key, val)

    unit_hydrographs = f.createVariable('unit_hydrograph', NC_DOUBLE,
                                        ('timesteps', 'lat', 'lon',),
                                        fill_value=FILLVALUE_F,
                                        **ncvaropts)
    for key, val in iteritems(share.unit_hydrograph.__dict__):
        if val:
            setattr(unit_hydrographs, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write data to variables initialized above
    lon[:] = agg_data['lon']
    lat[:] = agg_data['lat']
    unit_hydrographs[:, :, :] = agg_data['unit_hydrograph']
    fraction[:, :] = agg_data['fraction']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write globals
    for key, val in iteritems(glob_atts.atts):
        if val:
            setattr(f, key, val)
    # ---------------------------------------------------------------- #

    f.close()
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write the RVIC parameter file
def write_param_file(file_name,
                     nc_format='NETCDF3_CLASSIC',
                     glob_atts=NcGlobals(),
                     full_time_length=None,
                     subset_length=None,
                     unit_hydrograph_dt=None,
                     outlet_lon=None,
                     outlet_lat=None,
                     outlet_x_ind=None,
                     outlet_y_ind=None,
                     outlet_decomp_ind=None,
                     outlet_number=None,
                     outlet_mask=None,
                     outlet_name=None,
                     outlet_upstream_gridcells=None,
                     outlet_upstream_area=None,
                     source_lon=None,
                     source_lat=None,
                     source_x_ind=None,
                     source_y_ind=None,
                     source_decomp_ind=None,
                     source_time_offset=None,
                     source2outlet_ind=None,
                     unit_hydrograph=None,
                     zlib=True,
                     complevel=4,
                     least_significant_digit=None):

    """ Write a standard RVIC Parameter file """

    # ---------------------------------------------------------------- #
    # netCDF variable options
    ncvaropts = {'zlib': zlib,
                 'complevel': complevel,
                 'least_significant_digit': least_significant_digit}
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Open file
    f = Dataset(file_name, 'w', format=nc_format)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Time Variables

    # Timesteps
    timesteps = f.createDimension('timesteps', subset_length)
    timesteps = f.createVariable('timesteps', NC_DOUBLE, ('timesteps',),
                                 **ncvaropts)
    timesteps[:] = np.arange(subset_length)
    for key, val in iteritems(share.timesteps.__dict__):
        if val:
            setattr(timesteps, key, val)
    timesteps.timestep_length = 'unit_hydrograph_dt'

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write global attributes
    glob_atts.update()
    for key, val in iteritems(glob_atts.atts):
        if val:
            setattr(f, key, val)
    f.featureType = "timeSeries"
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 0-D variables

    # Full time length (size of ring)
    ftl = f.createVariable('full_time_length', NC_INT, (), **ncvaropts)
    ftl[:] = full_time_length
    for key, val in iteritems(share.full_time_length.__dict__):
        if val:
            setattr(ftl, key, val)

    # Subset Length
    sl = f.createVariable('subset_length', NC_INT, (), **ncvaropts)
    sl[:] = subset_length
    for key, val in iteritems(share.subset_length.__dict__):
        if val:
            setattr(sl, key, val)

    # UH timestep
    uh_dt = f.createVariable('unit_hydrograph_dt', NC_DOUBLE, (), **ncvaropts)
    uh_dt[:] = unit_hydrograph_dt
    for key, val in iteritems(share.unit_hydrograph_dt.__dict__):
        if val:
            setattr(uh_dt, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Outlet Dimensions
    if outlet_y_ind.ndim == 0:
        numoutlets = 1
        outlet_name = np.array([outlet_name])
    else:
        numoutlets = len(outlet_lon)
    ocoords = ('outlets',)
    f.createDimension(ocoords[0], numoutlets)

    nocoords = ocoords + ('nc_chars',)
    char_names = stringtochar(outlet_name)
    f.createDimension(nocoords[1], char_names.shape[1])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 1-D Outlet Variables

    # Outlet Cell Longitudes
    olon = f.createVariable('outlet_lon', NC_DOUBLE, ocoords, **ncvaropts)
    olon[:] = outlet_lon
    for key, val in iteritems(share.outlet_lon.__dict__):
        if val:
            setattr(olon, key, val)

    # Outlet Cell Latitudes
    olat = f.createVariable('outlet_lat', NC_DOUBLE, ocoords, **ncvaropts)
    olat[:] = outlet_lat
    for key, val in iteritems(share.outlet_lat.__dict__):
        if val:
            setattr(olat, key, val)

    # Outlet Cell X Indicies
    oxi = f.createVariable('outlet_x_ind', NC_INT, ocoords, **ncvaropts)
    oxi[:] = outlet_x_ind
    for key, val in iteritems(share.outlet_x_ind.__dict__):
        if val:
            setattr(oxi, key, val)

    # Outlet Cell Y Indicies
    oyi = f.createVariable('outlet_y_ind', NC_INT, ocoords, **ncvaropts)
    oyi[:] = outlet_y_ind
    for key, val in iteritems(share.outlet_y_ind.__dict__):
        if val:
            setattr(oyi, key, val)

    # Outlet Cell Decomp IDs
    odi = f.createVariable('outlet_decomp_ind', NC_INT, ocoords, **ncvaropts)
    odi[:] = outlet_decomp_ind
    for key, val in iteritems(share.outlet_decomp_ind.__dict__):
        if val:
            setattr(odi, key, val)

    # Outlet Cell Number
    on = f.createVariable('outlet_number', NC_INT, ocoords, **ncvaropts)
    on[:] = outlet_number
    for key, val in iteritems(share.outlet_number.__dict__):
        if val:
            setattr(on, key, val)

    # Outlet Mask
    om = f.createVariable('outlet_mask', NC_INT, ocoords, **ncvaropts)
    om[:] = outlet_mask
    for key, val in iteritems(share.outlet_mask.__dict__):
        if val:
            setattr(om, key, val)

    # Outlet Upstream area
    oua = f.createVariable('outlet_upstream_area', NC_DOUBLE, ocoords,
                           **ncvaropts)
    oua[:] = outlet_upstream_area
    for key, val in iteritems(share.outlet_upstream_area.__dict__):
        if val:
            setattr(oua, key, val)

    # Outlet Upstream grid cells
    oug = f.createVariable('outlet_upstream_gridcells', NC_INT, ocoords,
                           **ncvaropts)
    oug[:] = outlet_upstream_gridcells
    for key, val in iteritems(share.outlet_upstream_gridcells.__dict__):
        if val:
            setattr(oug, key, val)

    # Outlet Names
    onm = f.createVariable('outlet_name', NC_CHAR, nocoords)
    onm[:, :] = char_names
    for key, val in iteritems(share.outlet_name.__dict__):
        if val:
            setattr(onm, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Source Dimension
    scoords = ('sources',)
    f.createDimension(scoords[0], len(source_lon))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 1D Source Variables

    # Source Cell Longitudes
    slon = f.createVariable('source_lon', NC_DOUBLE, scoords, **ncvaropts)
    slon[:] = source_lon
    for key, val in iteritems(share.source_lon.__dict__):
        if val:
            setattr(slon, key, val)

    # Source Cell Latitudes
    slat = f.createVariable('source_lat', NC_DOUBLE, scoords, **ncvaropts)
    slat[:] = source_lat
    for key, val in iteritems(share.source_lat.__dict__):
        if val:
            setattr(slat, key, val)

    # Source Cell X Indicies
    sxi = f.createVariable('source_x_ind', NC_INT, scoords, **ncvaropts)
    sxi[:] = source_x_ind
    for key, val in iteritems(share.source_x_ind.__dict__):
        if val:
            setattr(sxi, key, val)

    # Source Cell Y Indicies
    syi = f.createVariable('source_y_ind', NC_INT, scoords, **ncvaropts)
    syi[:] = source_y_ind
    for key, val in iteritems(share.source_y_ind.__dict__):
        if val:
            setattr(syi, key, val)

    # Source Cell Decomp IDs
    sdi = f.createVariable('source_decomp_ind', NC_INT, scoords, **ncvaropts)
    sdi[:] = source_decomp_ind
    for key, val in iteritems(share.source_decomp_ind.__dict__):
        if val:
            setattr(sdi, key, val)

    # Source Cell Time Offset
    sto = f.createVariable('source_time_offset', NC_INT, scoords, **ncvaropts)
    sto[:] = source_time_offset
    for key, val in iteritems(share.source_time_offset.__dict__):
        if val:
            setattr(sto, key, val)

    # Source to Outlet Index Mapping
    s2o = f.createVariable('source2outlet_ind', NC_INT, scoords, **ncvaropts)
    s2o[:] = source2outlet_ind
    for key, val in iteritems(share.source2outlet_ind.__dict__):
        if val:
            setattr(s2o, key, val)

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 3-D Source Variables
    uhcords = ('timesteps',) + scoords + ('tracers',)
    f.createDimension(uhcords[2], 1)

    # Unit Hydrographs
    uhs = f.createVariable('unit_hydrograph', NC_DOUBLE, uhcords, **ncvaropts)
    uhs[:, :] = unit_hydrograph
    for key, val in iteritems(share.unit_hydrograph.__dict__):
        if val:
            setattr(uhs, key, val)
    # ---------------------------------------------------------------- #

    f.close()

    log.info('Finished writing %s', file_name)
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
