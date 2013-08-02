"""
write.py
"""

import numpy as np
from netCDF4 import Dataset
from logging import getLogger
from log import LOG_NAME
from share import NC_DOUBLE, NC_INT, FILLVALUE_F
import share

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write the agg netcdf
def write_agg_netcdf(file_name, agg_data, glob_atts, format):
    """
    Write output to netCDF.  Writes out a netCDF4 data file containing
    the UH_S and fractions and a full set of history and description attributes.
    """
    # ---------------------------------------------------------------- #
    # Open file
    f = Dataset(file_name, 'w', format=format)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # set dimensions
    timesteps = f.createDimension('timesteps', None)
    lon = f.createDimension('lon', (len(agg_data['lon'])))
    lat = f.createDimension('lat', (len(agg_data['lat'])))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # initialize variables
    unit_hydrogaph_dt = f.createVariable('unit_hydrogaph_dt', NC_INT, ())
    unit_hydrogaph_dt[:] = agg_data['unit_hydrogaph_dt']
    for key, val in share.unit_hydrogaph_dt.__dict__.iteritems():
        if val:
            setattr(unit_hydrogaph_dt, key, val)

    timesteps = f.createVariable('timesteps', NC_INT, ('timesteps',))
    timesteps[:] = agg_data['timesteps']
    for key, val in share.timesteps.__dict__.iteritems():
        if val:
            setattr(timesteps, key, val)

    lon = f.createVariable('lon', NC_DOUBLE, ('lon',))
    for key, val in share.lon.__dict__.iteritems():
        if val:
            setattr(lon, key, val)

    lat = f.createVariable('lat', NC_DOUBLE, ('lat',))
    for key, val in share.lat.__dict__.iteritems():
        if val:
            setattr(lat, key, val)

    fractions = f.createVariable('fractions', NC_DOUBLE, ('lat', 'lon',))
    for key, val in share.fractions.__dict__.iteritems():
        if val:
            setattr(fractions, key, val)

    unit_hydrographs = f.createVariable('unit_hydrographs', NC_DOUBLE,
                                        ('timesteps', 'lat', 'lon',), fill_value=FILLVALUE_F)
    for key, val in share.unit_hydrographs.__dict__.iteritems():
        if val:
            setattr(unit_hydrographs, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write data to variables initialized above
    lon[:] = agg_data['lon']
    lat[:] = agg_data['lat']
    unit_hydrographs[:, :, :] = agg_data['UHgrid']
    fractions[:, :] = agg_data['fractions']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write globals
    for key, val in glob_atts.__dict__.iteritems():
        if val:
            setattr(f, key, val)
    # ---------------------------------------------------------------- #

    f.close()
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write the RVIC parameter file
def write_param_file(file_name, out_format, sub_length, full_t_length,
                     uh_dt, uhs, lon_s, lat_s, s_decomp_id, x_ind_s, y_ind_s,
                     o_decomp_id, t_offset_s, s2outlet_index, lon_o, lat_o,
                     x_ind_o, y_ind_o, o_num, glob_atts):
    """ Write a standard RVIC Parameter file """

    # ---------------------------------------------------------------- #
    # Open file
    f = Dataset(file_name, 'w', format=out_format)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Time Variables

    # Timesteps
    timesteps = f.createDimension('timesteps', sub_length)
    timesteps = f.createVariable('timesteps', NC_DOUBLE, ('timesteps',))
    timesteps[:] = np.arange(sub_length)

    for key, val in share.timesteps.__dict__.iteritems():
        if val:
            setattr(timesteps, key, val)
    timesteps.timestep_length = 'unit_hydrogaph_dt'

    # UH timestep
    unit_hydrogaph_dt = f.createVariable('unit_hydrogaph_dt', NC_DOUBLE, ())
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
    subset_length = f.createVariable('subset_length', NC_INT, ())
    subset_length[:] = sub_length
    for key, val in share.subset_length.__dict__.iteritems():
        if val:
            setattr(subset_length, key, val)

    full_time_length = f.createVariable('full_time_length', NC_INT, ())
    full_time_length[:] = full_t_length
    for key, val in share.full_time_length.__dict__.iteritems():
        if val:
            setattr(full_time_length, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Source Variables

    unit_hydrographs = f.createVariable('unit_hydrographs', NC_DOUBLE, tscoords)
    unit_hydrographs[:, :] = uhs
    for key, val in share.unit_hydrographs.__dict__.iteritems():
        if val:
            setattr(unit_hydrographs, key, val)

    source_decomp_id = f.createVariable('source_decomp_id', NC_INT, scoords)
    source_decomp_id[:] = s_decomp_id
    for key, val in share.source_decomp_id.__dict__.iteritems():
        if val:
            setattr(source_decomp_id, key, val)

    y_ind_source = f.createVariable('y_ind_source', NC_INT, scoords)
    y_ind_source[:] = y_ind_s
    for key, val in share.y_ind_source.__dict__.iteritems():
        if val:
            setattr(y_ind_source, key, val)

    x_ind_source = f.createVariable('x_ind_source', NC_INT, scoords)
    x_ind_source[:] = x_ind_s
    for key, val in share.x_ind_source.__dict__.iteritems():
        if val:
            setattr(x_ind_source, key, val)

    lat_source = f.createVariable('lat_source', NC_DOUBLE, scoords)
    lat_source[:] = lat_s
    for key, val in share.lat_source.__dict__.iteritems():
        if val:
            setattr(lat_source, key, val)

    lon_source = f.createVariable('lon_source', NC_DOUBLE, scoords)
    lon_source[:] = lon_s
    for key, val in share.lon_source.__dict__.iteritems():
        if val:
            setattr(lon_source, key, val)

    time_offset_source = f.createVariable('time_offset_source', NC_INT, scoords)
    time_offset_source[:] = t_offset_s
    for key, val in share.time_offset_source.__dict__.iteritems():
        if val:
            setattr(time_offset_source, key, val)

    source2outlet_index = f.createVariable('source2outlet_index', NC_INT, scoords)
    source2outlet_index[:] = s2outlet_index
    for key, val in share.source2outlet_index.__dict__.iteritems():
        if val:
            setattr(source2outlet_index, key, val)

    # ---------------------------------------------------------------- #
    # Outlet Variables
    outlet_decomp_id = f.createVariable('outlet_decomp_id', NC_INT, ocoords)
    outlet_decomp_id[:] = o_decomp_id
    for key, val in share.outlet_decomp_id.__dict__.iteritems():
        if val:
            setattr(outlet_decomp_id, key, val)

    outlet_num = f.createVariable('outlet_num', NC_INT, ocoords)
    outlet_num[:] = o_num
    for key, val in share.outlet_num.__dict__.iteritems():
        if val:
            setattr(outlet_num, key, val)

    x_ind_outlet = f.createVariable('x_ind_outlet', NC_INT, ocoords)
    x_ind_outlet[:] = x_ind_o
    for key, val in share.x_ind_outlet.__dict__.iteritems():
        if val:
            setattr(x_ind_outlet, key, val)

    y_ind_outlet = f.createVariable('y_ind_outlet', NC_INT, ocoords)
    y_ind_outlet[:] = y_ind_o
    for key, val in share.y_ind_outlet.__dict__.iteritems():
        if val:
            setattr(y_ind_outlet, key, val)

    lon_outlet = f.createVariable('lon_outlet', NC_DOUBLE, ocoords)
    lon_outlet[:] = lon_o
    for key, val in share.lon_outlet.__dict__.iteritems():
        if val:
            setattr(lon_outlet, key, val)

    lat_outlet = f.createVariable('lat_outlet', NC_DOUBLE, ocoords)
    lat_outlet[:] = lat_o
    for key, val in share.lat_outlet.__dict__.iteritems():
        if val:
            setattr(lat_outlet, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write global attributes
    glob_atts.update()
    for key, val in glob_atts.__dict__.iteritems():
        if val:
            setattr(f, key, val)
    # ---------------------------------------------------------------- #

    f.close()

    log.info('Finished writing %s' % file_name)
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
