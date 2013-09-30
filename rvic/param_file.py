"""
params.py
"""
import numpy as np
import logging
from log import LOG_NAME
from rvic.write import write_param_file
from rvic.share import NcGlobals, PRECISION
import os
from datetime import date

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
#
def finish_params(outlets, dom_data, config_dict, directories):
    options = config_dict['OPTIONS']

    if config_dict['NEW_DOMAIN']:
        dom_file_name = config_dict['NEW_DOMAIN']['FILE_NAME']
    else:
        dom_file_name = config_dict['DOMAIN']['FILE_NAME']

    # ---------------------------------------------------------------- #
    # adjust fractions
    if options['CONSTRAIN_FRACTIONS']:
        outlets = adjust_fractions(outlets, dom_data[config_dict['DOMAIN']['FRACTION_VAR']])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # subset
    outlets, full_time_length = subset(outlets,
                                       subset_length=options['SUBSET_LENGTH'],
                                       threshold=options['SUBSET_THRESHOLD'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Group
    unit_hydrograph, frac_sources, source_lon, source_lat, source_x_ind, \
    source_y_ind, source_decomp_ind, source_time_offset, source2outlet_ind, \
    outlet_lon, outlet_lat, outlet_x_ind, outlet_y_ind, outlet_decomp_ind, \
    outlet_number = group(outlets)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust Unit Hydrographs for differences in source/outlet areas
    unit_hydrograph *= frac_sources
    unit_hydrograph *= dom_data[config_dict['DOMAIN']['AREA_VAR']][source_y_ind, source_x_ind]

    for p, ind in enumerate(source2outlet_ind):
        unit_hydrograph[:, p] /= dom_data[config_dict['DOMAIN']['AREA_VAR']][outlet_y_ind[ind], outlet_x_ind[ind]]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # fill in some misc arrays
    outlet_mask = np.zeros(len(outlet_lon))
    newshape = unit_hydrograph.shape + (1,)
    unit_hydrograph = unit_hydrograph.reshape(newshape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write parameter file
    today = date.today().strftime('%Y%m%d')
    param_file = os.path.join(directories['params'],
                             '%s.rvic.prm.%s.%s.nc' % (options['CASEID'],
                                                       options['GRIDID'], today))

    write_param_file(param_file,
                     nc_format=options['NETCDF_FORMAT'],
                     glob_atts=NcGlobals(title='RVIC parameter file',
                                         RvicPourPointsFile=os.path.split(config_dict['POUR_POINTS']['FILE_NAME'])[1],
                                         RvicUHFile=os.path.split(config_dict['UH_BOX']['FILE_NAME'])[1],
                                         RvicFdrFile=os.path.split(config_dict['ROUTING']['FILE_NAME'])[1],
                                         RvicDomainFile=os.path.split(dom_file_name)[1]),
                     full_time_length=full_time_length,
                     subset_length=options['SUBSET_LENGTH'],
                     unit_hydrograph_dt=config_dict['ROUTING']['OUTPUT_INTERVAL'],
                     outlet_lon=outlet_lon,
                     outlet_lat=outlet_lat,
                     outlet_x_ind=outlet_x_ind,
                     outlet_y_ind=outlet_y_ind,
                     outlet_decomp_ind=outlet_decomp_ind,
                     outlet_number=outlet_number,
                     outlet_mask=outlet_mask,
                     source_lon=source_lon,
                     source_lat=source_lat,
                     source_x_ind=source_x_ind,
                     source_y_ind=source_y_ind,
                     source_decomp_ind=source_decomp_ind,
                     source_time_offset=source_time_offset,
                     source2outlet_ind=source2outlet_ind,
                     unit_hydrograph=unit_hydrograph)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write a summary of what was done to the log file.
    log.info('Parameter file includes %i outlets' % (len(outlets)))
    log.info('Parameter file includes %i Source Points' % (len(source_lon)))
    # ---------------------------------------------------------------- #

    return param_file, today
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def adjust_fractions(outlets, dom_fractions):
    """
    Constrain the fractions in the outles.
    The basic idea is that the sum of fraction from the outlets should not
    exceed the domain fractions.
    """

    log.debug('In adjust fractions now')

    # ---------------------------------------------------------------- #
    # Aggregate the fractions
    fractions = np.zeros(dom_fractions.shape)
    for cell_id, outlet in outlets.iteritems():
        y = outlet.y_source
        x = outlet.x_source

        fractions[y, x] += outlet.fractions
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Determin how to adjust the fractions
    diff_fractions = fractions - dom_fractions
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Only adjust fractions where the fractions are gt the domain fractions
    yi, xi = np.nonzero(fractions > dom_fractions)
    ratio_fraction = np.ones(diff_fractions.shape)
    ratio_fraction[yi, xi] = dom_fractions[yi, xi]/fractions[yi, xi]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust fracs based on ratio_fraction
    for cell_id, outlet in outlets.iteritems():
        y = outlet.y_source
        x = outlet.x_source
        outlets[cell_id].fractions *= ratio_fraction[y, x]
    # ---------------------------------------------------------------- #

    return outlets

# -------------------------------------------------------------------- #
# Shorten the unit hydrograph
def subset(outlets, subset_length=None, threshold=PRECISION):
    """ Shorten the Unit Hydrograph"""
    if not threshold:
        threshold = PRECISION

    for i, (cell_id, outlet) in enumerate(outlets.iteritems()):
        if i == 0:
            full_time_length = outlet.unit_hydrograph.shape[0]
            if not subset_length:
                subset_length = full_time_length

        outlets[cell_id].offset = np.empty(outlet.unit_hydrograph.shape[1], dtype=int)
        out_uh = np.zeros((subset_length, outlet.unit_hydrograph.shape[1]))
        for i in xrange(outlet.unit_hydrograph.shape[1]):
            # find position of first index > threshold
            inds = np.nonzero(outlet.unit_hydrograph[:, i] > threshold)
            if len(inds[0]) > 0:
                start = inds[0][0]
            else:
                start = 0
            # find end point
            end = np.minimum(full_time_length, start + subset_length)
            start = np.minimum(start, full_time_length - subset_length)
            outlets[cell_id].offset[i] = start
            # clip and normalize
            out_uh[:, i] = outlet.unit_hydrograph[start:end, i] / outlet.unit_hydrograph[start:end, i].sum()
        outlets[cell_id].unit_hydrograph = out_uh

    return outlets, full_time_length
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def group(outlets):
    """
    group the outlets into one set of arrays
    """

    for i, (cell_id, outlet) in enumerate(outlets.iteritems()):

        y = outlet.y_source
        x = outlet.x_source

        if i == 0:
            # -------------------------------------------------------- #
            # Source specific values
            unit_hydrograph = outlet.unit_hydrograph
            frac_sources = outlet.fractions
            source_lon = outlet.lon_source
            source_lat = outlet.lat_source
            source_x_ind = x
            source_y_ind = y
            source_decomp_ind = outlet.cell_id_source
            source_time_offset = outlet.offset
            source2outlet_ind = np.zeros(len(outlet.fractions))
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_lon = np.array(outlet.lon)
            outlet_lat = np.array(outlet.lat)
            outlet_x_ind = np.array(outlet.x)
            outlet_y_ind = np.array(outlet.y)
            outlet_decomp_ind = np.array(cell_id)
            outlet_number = np.array(i)
        else:
            # -------------------------------------------------------- #
            # Point specific values
            unit_hydrograph = np.append(unit_hydrograph, outlet.unit_hydrograph, axis=1)
            frac_sources = np.append(frac_sources, outlet.fractions)
            source_lon = np.append(source_lon, outlet.lon_source)
            source_lat = np.append(source_lat, outlet.lat_source)
            source_x_ind = np.append(source_x_ind, x)
            source_y_ind = np.append(source_y_ind, y)
            source_decomp_ind = np.append(source_decomp_ind, outlet.cell_id_source)
            source_time_offset = np.append(source_time_offset, outlet.offset)
            source2outlet_ind = np.append(source2outlet_ind, np.zeros_like(outlet.offset) + i)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_lon = np.append(outlet_lon, outlet.lon)
            outlet_lat = np.append(outlet_lat, outlet.lat)
            outlet_x_ind = np.append(outlet_x_ind, outlet.x)
            outlet_y_ind = np.append(outlet_y_ind, outlet.y)
            outlet_decomp_ind = np.append(outlet_decomp_ind, cell_id)
            outlet_number = np.append(outlet_number, i)
            # -------------------------------------------------------- #

    return unit_hydrograph, frac_sources, source_lon, source_lat, source_x_ind, \
    source_y_ind, source_decomp_ind, source_time_offset, source2outlet_ind, \
    outlet_lon, outlet_lat, outlet_x_ind, outlet_y_ind, outlet_decomp_ind, \
    outlet_number
# -------------------------------------------------------------------- #
