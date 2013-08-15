"""
write.py
"""

import numpy as np
from netCDF4 import Dataset
from logging import getLogger
from log import LOG_NAME
from share import NC_DOUBLE, NC_INT, FILLVALUE_F, NcGlobals
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
    log.info('writing aggregated netcdf: %s' % file_name)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # set dimensions
    timesteps = f.createDimension('timesteps', None)
    lon = f.createDimension('lon', (len(agg_data['lon'])))
    lat = f.createDimension('lat', (len(agg_data['lat'])))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # initialize variables
    unit_hydrograph_dt = f.createVariable('unit_hydrograph_dt', NC_INT, ())
    unit_hydrograph_dt[:] = agg_data['unit_hydrograph_dt']
    for key, val in share.unit_hydrograph_dt.__dict__.iteritems():
        if val:
            setattr(unit_hydrograph_dt, key, val)

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

    fraction = f.createVariable('fraction', NC_DOUBLE, ('lat', 'lon',))
    for key, val in share.fraction.__dict__.iteritems():
        if val:
            setattr(fraction, key, val)

    unit_hydrographs = f.createVariable('unit_hydrograph', NC_DOUBLE,
                                        ('timesteps', 'lat', 'lon',), fill_value=FILLVALUE_F)
    for key, val in share.unit_hydrograph.__dict__.iteritems():
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
    for key, val in glob_atts.__dict__.iteritems():
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
                     source_lon=None,
                     source_lat=None,
                     source_x_ind=None,
                     source_y_ind=None,
                     source_decomp_ind=None,
                     source_time_offset=None,
                     source2outlet_ind=None,
                     source_tracer=None,
                     unit_hydrograph=None):

    """ Write a standard RVIC Parameter file """

    # ---------------------------------------------------------------- #
    # Open file
    f = Dataset(file_name, 'w', format=nc_format)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Time Variables

    # Timesteps
    timesteps = f.createDimension('timesteps', subset_length)
    timesteps = f.createVariable('timesteps', NC_DOUBLE, ('timesteps',))
    timesteps[:] = np.arange(subset_length)
    for key, val in share.timesteps.__dict__.iteritems():
        if val:
            setattr(timesteps, key, val)
    timesteps.timestep_length = 'unit_hydrograph_dt'

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write global attributes
    glob_atts.update()
    for key, val in glob_atts.__dict__.iteritems():
        if val:
            setattr(f, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 0-D variables

    # Full time length (size of ring)
    ftl = f.createVariable('full_time_length', NC_INT, ())
    ftl[:] = full_time_length
    for key, val in share.full_time_length.__dict__.iteritems():
        if val:
            setattr(ftl, key, val)

    # Subset Length
    sl = f.createVariable('subset_length', NC_INT, ())
    sl[:] = subset_length
    for key, val in share.subset_length.__dict__.iteritems():
        if val:
            setattr(sl, key, val)

    # UH timestep
    uh_dt = f.createVariable('unit_hydrograph_dt', NC_DOUBLE, ())
    uh_dt[:] = unit_hydrograph_dt
    for key, val in share.unit_hydrograph_dt.__dict__.iteritems():
        if val:
            setattr(uh_dt, key, val)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Outlet Dimension
    ocoords = ('outlets',)
    outlets = f.createDimension(ocoords[0], len(outlet_lon))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 1-D Outlet Variables

    # Outlet Cell Longitudes
    olon = f.createVariable('outlet_lon', NC_DOUBLE, ocoords)
    olon[:] = outlet_lon
    for key, val in share.outlet_lon.__dict__.iteritems():
        if val:
            setattr(olon, key, val)

    # Outlet Cell Latitudes
    olat = f.createVariable('outlet_lat', NC_DOUBLE, ocoords)
    olat[:] = outlet_lat
    for key, val in share.outlet_lat.__dict__.iteritems():
        if val:
            setattr(olat, key, val)

    # Outlet Cell X Indicies
    oxi = f.createVariable('outlet_x_ind', NC_INT, ocoords)
    oxi[:] = outlet_x_ind
    for key, val in share.outlet_x_ind.__dict__.iteritems():
        if val:
            setattr(oxi, key, val)

    # Outlet Cell Y Indicies
    oyi = f.createVariable('outlet_y_ind', NC_INT, ocoords)
    oyi[:] = outlet_y_ind
    for key, val in share.outlet_y_ind.__dict__.iteritems():
        if val:
            setattr(oyi, key, val)

    # Outlet Cell Decomp IDs
    odi = f.createVariable('outlet_decomp_ind', NC_INT, ocoords)
    odi[:] = outlet_decomp_ind
    for key, val in share.outlet_decomp_ind.__dict__.iteritems():
        if val:
            setattr(odi, key, val)

    # Outlet Cell Number
    on = f.createVariable('outlet_number', NC_INT, ocoords)
    on[:] = outlet_number
    for key, val in share.outlet_number.__dict__.iteritems():
        if val:
            setattr(on, key, val)

    # Outlet Mask
    om = f.createVariable('outlet_mask', NC_INT, ocoords)
    om[:] = outlet_mask
    for key, val in share.outlet_mask.__dict__.iteritems():
        if val:
            setattr(om, key, val)

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Source Dimension
    scoords = ('sources',)
    sources = f.createDimension(scoords[0], len(source_lon))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 1D Source Variables

    # Source Cell Longitudes
    slon = f.createVariable('source_lon', NC_DOUBLE, scoords)
    slon[:] = source_lon
    for key, val in share.source_lon.__dict__.iteritems():
        if val:
            setattr(slon, key, val)

    # Source Cell Latitudes
    slat = f.createVariable('source_lat', NC_DOUBLE, scoords)
    slat[:] = source_lat
    for key, val in share.source_lat.__dict__.iteritems():
        if val:
            setattr(slat, key, val)

    # Source Cell X Indicies
    sxi = f.createVariable('source_x_ind', NC_INT, scoords)
    sxi[:] = source_x_ind
    for key, val in share.source_x_ind.__dict__.iteritems():
        if val:
            setattr(sxi, key, val)

    # Source Cell Y Indicies
    syi = f.createVariable('source_y_ind', NC_INT, scoords)
    syi[:] = source_y_ind
    for key, val in share.source_y_ind.__dict__.iteritems():
        if val:
            setattr(syi, key, val)

    # Source Cell Decomp IDs
    sdi = f.createVariable('source_decomp_ind', NC_INT, scoords)
    sdi[:] = source_decomp_ind
    for key, val in share.source_decomp_ind.__dict__.iteritems():
        if val:
            setattr(sdi, key, val)

    # Source Cell Time Offset
    sto = f.createVariable('source_time_offset', NC_INT, scoords)
    sto[:] = source_time_offset
    for key, val in share.source_time_offset.__dict__.iteritems():
        if val:
            setattr(sto, key, val)

    # Source to Outlet Index Mapping
    s2o = f.createVariable('source2outlet_ind', NC_INT, scoords)
    s2o[:] = source2outlet_ind
    for key, val in share.source2outlet_ind.__dict__.iteritems():
        if val:
            setattr(s2o, key, val)

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # 3-D Source Variables
    uhcords = ('timesteps',) + scoords +('tracers',)
    tracers = f.createDimension(uhcords[2], 1)

    # Unit Hydrographs
    uhs = f.createVariable('unit_hydrograph', NC_DOUBLE, uhcords)
    uhs[:, :] = unit_hydrograph
    for key, val in share.unit_hydrograph.__dict__.iteritems():
        if val:
            setattr(uhs, key, val)
    # ---------------------------------------------------------------- #

    f.close()

    log.info('Finished writing %s' % file_name)
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
