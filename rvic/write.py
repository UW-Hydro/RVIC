"""
write.py
"""

import logging
import numpy as np
from netCDF4 import Dataset
from log import log_name
from share import nc_double, nc_int, fillValue_f
import share

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(log_name)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write the agg netcdf
def write_agg_netcdf(FileName, aggData, GlobAts, format):
    """
    Write output to netCDF.  Writes out a netCDF4 data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    # ---------------------------------------------------------------- #
    # Open file
    f = Dataset(FileName, 'w', format=format)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # set dimensions
    timesteps = f.createDimension('timesteps', None)
    lon = f.createDimension('lon', (len(aggData['lon'])))
    lat = f.createDimension('lat', (len(aggData['lat'])))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # initialize variables
    unit_hydrogaph_dt = f.createVariable('unit_hydrogaph_dt', nc_int, ())
    unit_hydrogaph_dt[:] = aggData['unit_hydrogaph_dt']
    for key, val in share.unit_hydrogaph_dt.__dict__.iteritems():
        if val:
            setattr(unit_hydrogaph_dt, key, val)

    timesteps = f.createVariable('timesteps', nc_int, ('timesteps',))
    timesteps[:] = aggData['timesteps']
    for key, val in share.timesteps.__dict__.iteritems():
        if val:
            setattr(timesteps, key, val)

    lon = f.createVariable('lon', nc_double, ('lon',))
    for key, val in share.lon.__dict__.iteritems():
        if val:
            setattr(lon, key, val)

    lat = f.createVariable('lat', nc_double, ('lat',))
    for key, val in share.lat.__dict__.iteritems():
        if val:
            setattr(lat, key, val)

    fractions = f.createVariable('fractions', nc_double, ('lat', 'lon',))
    for key, val in share.fractions.__dict__.iteritems():
        if val:
            setattr(fractions, key, val)

    unit_hydrographs = f.createVariable('unit_hydrographs', nc_double,
                                        ('timesteps', 'lat', 'lon',), fill_value=fillValue_f)
    for key, val in share.unit_hydrographs.__dict__.iteritems():
        if val:
            setattr(unit_hydrographs, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write data to variables initialized above
    lon[:] = aggData['lon']
    lat[:] = aggData['lat']
    unit_hydrographs[:, :, :] = aggData['UHgrid']
    fractions[:, :] = aggData['fractions']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write globals
    for key, val in GlobAts.__dict__.iteritems():
        if val:
            setattr(f, key, val)
    # ---------------------------------------------------------------- #

    f.close()
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write the RVIC parameter file
def write_param_file(FileName, out_format, sub_length, full_t_length,
                     uh_dt, uhs, lon_s, lat_s, s_decomp_id, x_ind_s, y_ind_s,
                     o_decomp_id, t_offset_s, s2outlet_index, lon_o, lat_o,
                     x_ind_o, y_ind_o, o_num, GlobAts):
    """ Write a standard RVIC Parameter file """

    # ---------------------------------------------------------------- #
    # Open file
    f = Dataset(FileName, 'w', format=out_format)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Time Variables

    # Timesteps
    timesteps = f.createDimension('timesteps', sub_length)
    timesteps = f.createVariable('timesteps', nc_double, ('timesteps',))
    timesteps[:] = np.arange(sub_length)

    for key, val in share.timesteps.__dict__.iteritems():
        if val:
            setattr(timesteps, key, val)
    timesteps.timestep_length = 'unit_hydrogaph_dt'

    # UH timestep
    unit_hydrogaph_dt = f.createVariable('unit_hydrogaph_dt', nc_double, ())
    unit_hydrogaph_dt[:] = uh_dt
    for key, val in share.unit_hydrogaph_dt.__dict__.iteritems():
        if val:
            setattr(unit_hydrogaph_dt, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Point Dimensions
    ocoords = ('outlets',)
    scoords = ('sources',)
    tscoords = ('timesteps',) + scoords

    outlets = f.createDimension('outlets', len(lon_o))
    sources = f.createDimension('sources', len(lon_s))

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 0-D variables
    subset_length = f.createVariable('subset_length', nc_int, ())
    subset_length[:] = sub_length
    for key, val in share.subset_length.__dict__.iteritems():
        if val:
            setattr(subset_length, key, val)

    full_time_length = f.createVariable('full_time_length', nc_int, ())
    full_time_length[:] = full_t_length
    for key, val in share.full_time_length.__dict__.iteritems():
        if val:
            setattr(full_time_length, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Source Variables

    unit_hydrographs = f.createVariable('unit_hydrographs', nc_double, tscoords)
    unit_hydrographs[:, :] = uhs
    for key, val in share.unit_hydrographs.__dict__.iteritems():
        if val:
            setattr(unit_hydrographs, key, val)

    source_decomp_id = f.createVariable('source_decomp_id', nc_int, scoords)
    source_decomp_id[:] = s_decomp_id
    for key, val in share.source_decomp_id.__dict__.iteritems():
        if val:
            setattr(source_decomp_id, key, val)

    y_ind_source = f.createVariable('y_ind_source', nc_int, scoords)
    y_ind_source[:] = y_ind_s
    for key, val in share.y_ind_source.__dict__.iteritems():
        if val:
            setattr(y_ind_source, key, val)

    x_ind_source = f.createVariable('x_ind_source', nc_int, scoords)
    x_ind_source[:] = x_ind_s
    for key, val in share.x_ind_source.__dict__.iteritems():
        if val:
            setattr(x_ind_source, key, val)

    lat_source = f.createVariable('lat_source', nc_double, scoords)
    lat_source[:] = lat_s
    for key, val in share.lat_source.__dict__.iteritems():
        if val:
            setattr(lat_source, key, val)

    lon_source = f.createVariable('lon_source', nc_double, scoords)
    lon_source[:] = lon_s
    for key, val in share.lon_source.__dict__.iteritems():
        if val:
            setattr(lon_source, key, val)

    time_offset_source = f.createVariable('time_offset_source', nc_int, scoords)
    time_offset_source[:] = t_offset_s
    for key, val in share.time_offset_source.__dict__.iteritems():
        if val:
            setattr(time_offset_source, key, val)

    source2outlet_index = f.createVariable('source2outlet_index', nc_int, scoords)
    source2outlet_index[:] = s2outlet_index
    for key, val in share.source2outlet_index.__dict__.iteritems():
        if val:
            setattr(source2outlet_index, key, val)

    # ---------------------------------------------------------------- #
    # Outlet Variables
    outlet_decomp_id = f.createVariable('outlet_decomp_id', nc_int, ocoords)
    outlet_decomp_id[:] = o_decomp_id
    for key, val in share.outlet_decomp_id.__dict__.iteritems():
        if val:
            setattr(outlet_decomp_id, key, val)

    outlet_num = f.createVariable('outlet_num', nc_int, ocoords)
    outlet_num[:] = o_num
    for key, val in share.outlet_num.__dict__.iteritems():
        if val:
            setattr(outlet_num, key, val)

    x_ind_outlet = f.createVariable('x_ind_outlet', nc_int, ocoords)
    x_ind_outlet[:] = x_ind_o
    for key, val in share.x_ind_outlet.__dict__.iteritems():
        if val:
            setattr(x_ind_outlet, key, val)

    y_ind_outlet = f.createVariable('y_ind_outlet', nc_int, ocoords)
    y_ind_outlet[:] = y_ind_o
    for key, val in share.y_ind_outlet.__dict__.iteritems():
        if val:
            setattr(y_ind_outlet, key, val)

    lon_outlet = f.createVariable('lon_outlet', nc_double, ocoords)
    lon_outlet[:] = lon_o
    for key, val in share.lon_outlet.__dict__.iteritems():
        if val:
            setattr(lon_outlet, key, val)

    lat_outlet = f.createVariable('lat_outlet', nc_double, ocoords)
    lat_outlet[:] = lat_o
    for key, val in share.lat_outlet.__dict__.iteritems():
        if val:
            setattr(lat_outlet, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write global attributes
    GlobAts.update()
    for key, val in GlobAts.__dict__.iteritems():
        if val:
            setattr(f, key, val)
    # ---------------------------------------------------------------- #

    f.close()

    log.info('Finished writing %s' % FileName)
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
