"""
param_file.py
"""
import numpy as np
import logging
from log import LOG_NAME
from write import write_param_file
from share import NcGlobals, SECSPERDAY
import os
from datetime import date
import plots


# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Wrap up functiions to finish the parameter file
def finish_params(outlets, dom_data, config_dict, directories):
    """
    Adjust the unit hydrographs and pack for parameter file
    """
    options = config_dict['OPTIONS']
    routing = config_dict['ROUTING']
    domain = config_dict['DOMAIN']
    dom_area = domain['AREA_VAR']
    dom_frac = domain['FRACTION_VAR']

    # ------------------------------------------------------------ #
    # netCDF variable options
    ncvaropts = {}
    if 'NETCDF_ZLIB' in options:
        ncvaropts['zlib'] = options['NETCDF_ZLIB']
    if 'NETCDF_COMPLEVEL' in options:
        ncvaropts['complevel'] = options['NETCDF_COMPLEVEL']
    if 'NETCDF_SIGFIGS' in options:
        ncvaropts['least_significant_digit'] = options['NETCDF_SIGFIGS']
    # ------------------------------------------------------------ #

    # ---------------------------------------------------------------- #
    # subset (shorten time base)
    if options['SUBSET_DAYS'] and \
            options['SUBSET_DAYS'] < routing['BASIN_FLOWDAYS']:
        subset_length = options['SUBSET_DAYS']*SECSPERDAY/routing['OUTPUT_INTERVAL']
        outlets, full_time_length, \
            before, after = subset(outlets, subset_length=subset_length)

        log.debug('plotting unit hydrograph timeseries now for before'
                  ' / after subseting')

        title = 'UHS before subset'
        pfname = plots.uhs(before, title, options['CASEID'],
                           directories['plots'])
        log.info('%s Plot:  %s', title, pfname)

        title = 'UHS after subset'
        pfname = plots.uhs(after, title, options['CASEID'],
                           directories['plots'])
        log.info('%s Plot:  %s', title, pfname)
    else:
        subset_length = routing['BASIN_FLOWDAYS']*SECSPERDAY/routing['OUTPUT_INTERVAL']
        log.info('Not subsetting because either SUBSET_DAYS is null or '
                 'SUBSET_DAYS<BASIN_FLOWDAYS')
        for i, (cell_id, outlet) in enumerate(outlets.iteritems()):
            outlets[cell_id].offset = np.zeros(outlet.unit_hydrograph.shape[1],
                                               dtype=np.int16)
        full_time_length = outlet.unit_hydrograph.shape[0]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # adjust fractions
    if options['CONSTRAIN_FRACTIONS']:
        adjust = True
        log.info('Adjusting Fractions to be less than or equal to '
                 'domain fractions')
    else:
        adjust = False
    outlets, plot_dict = adjust_fractions(outlets,
                                          dom_data[domain['FRACTION_VAR']],
                                          adjust=adjust)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Calculate the upstream area and upstream grid cells
    # The upstream_area must be calculated after adjust_fractions
    for i, outlet in outlets.iteritems():
        outlet.upstream_gridcells = len(outlet.y_source)
        outlet.upstream_area = np.sum(dom_data[dom_area][outlet.y_source, outlet.x_source] *
                                      dom_data[dom_frac][outlet.y_source, outlet.x_source])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Group
    grouped_data = group(outlets)

    # unpack grouped data
    unit_hydrograph = grouped_data['unit_hydrograph']
    frac_sources = grouped_data['frac_sources']
    source_lon = grouped_data['source_lon']
    source_lat = grouped_data['source_lat']
    source_x_ind = grouped_data['source_x_ind']
    source_y_ind = grouped_data['source_y_ind']
    source_decomp_ind = grouped_data['source_decomp_ind']
    source_time_offset = grouped_data['source_time_offset']
    source2outlet_ind = grouped_data['source2outlet_ind']
    outlet_lon = grouped_data['outlet_lon']
    outlet_lat = grouped_data['outlet_lat']
    outlet_x_ind = grouped_data['outlet_x_ind']
    outlet_y_ind = grouped_data['outlet_y_ind']
    outlet_decomp_ind = grouped_data['outlet_decomp_ind']
    outlet_number = grouped_data['outlet_number']
    outlet_name = grouped_data['outlet_name']
    outlet_upstream_area = grouped_data['outlet_upstream_area']
    outlet_upstream_gridcells = grouped_data['outlet_upstream_gridcells']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust Unit Hydrographs for differences in source/outlet areas and
    # fractions
    area = dom_data[domain['AREA_VAR']]

    for source, outlet in enumerate(source2outlet_ind):
        if outlet_y_ind.ndim == 0 or outlet_x_ind.ndim == 0:
            unit_hydrograph[:, source] *= area[source_y_ind[source],
                                               source_x_ind[source]]
            unit_hydrograph[:, source] /= area[outlet_y_ind[()],
                                               outlet_x_ind[()]]
            unit_hydrograph[:, source] *= frac_sources[source]
        else:
            unit_hydrograph[:, source] *= area[source_y_ind[source],
                                               source_x_ind[source]]
            unit_hydrograph[:, source] /= area[outlet_y_ind[outlet],
                                               outlet_x_ind[outlet]]
            unit_hydrograph[:, source] *= frac_sources[source]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make diagnostic plots
    sum_after = np.zeros(dom_data[domain['FRACTION_VAR']].shape)
    for i, (y, x) in enumerate(zip(source_y_ind, source_x_ind)):
        sum_after[y, x] += unit_hydrograph[:, i].sum()

    plot_dict['Sum UH Final'] = sum_after

    dom_y = dom_data[domain['LATITUDE_VAR']]
    dom_x = dom_data[domain['LONGITUDE_VAR']]

    for title, data in plot_dict.iteritems():
        pfname = plots.fractions(data, dom_x, dom_y, title, options['CASEID'],
                                 directories['plots'])
        log.info('%s Plot: %s', title, pfname)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # fill in some misc arrays
    if outlet_y_ind.ndim == 0:
        numoutlets = 1
    else:
        numoutlets = len(outlet_lon)
    outlet_mask = np.zeros(numoutlets)
    newshape = unit_hydrograph.shape + (1, )
    unit_hydrograph = unit_hydrograph.reshape(newshape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Write parameter file
    today = date.today().strftime('%Y%m%d')
    param_file = os.path.join(directories['params'],
                              '{0}.rvic.prm.{1}.{2}.'
                              'nc'.format(options['CASEID'],
                                          options['GRIDID'],
                                          today))

    if 'NEW_DOMAIN' in config_dict.keys():
        dom_file_name = config_dict['NEW_DOMAIN']['FILE_NAME']
    else:
        dom_file_name = config_dict['DOMAIN']['FILE_NAME']

    glob_atts = NcGlobals(title='RVIC parameter file',
                          RvicPourPointsFile=os.path.split(config_dict['POUR_POINTS']['FILE_NAME'])[1],
                          RvicUHFile=os.path.split(config_dict['UH_BOX']['FILE_NAME'])[1],
                          RvicFdrFile=os.path.split(routing['FILE_NAME'])[1],
                          RvicDomainFile=os.path.split(dom_file_name)[1])

    write_param_file(param_file,
                     nc_format=options['NETCDF_FORMAT'],
                     glob_atts=glob_atts,
                     full_time_length=full_time_length,
                     subset_length=subset_length,
                     unit_hydrograph_dt=routing['OUTPUT_INTERVAL'],
                     outlet_lon=outlet_lon,
                     outlet_lat=outlet_lat,
                     outlet_x_ind=outlet_x_ind,
                     outlet_y_ind=outlet_y_ind,
                     outlet_decomp_ind=outlet_decomp_ind,
                     outlet_number=outlet_number,
                     outlet_mask=outlet_mask,
                     outlet_name=outlet_name,
                     outlet_upstream_gridcells=outlet_upstream_gridcells,
                     outlet_upstream_area=outlet_upstream_area,
                     source_lon=source_lon,
                     source_lat=source_lat,
                     source_x_ind=source_x_ind,
                     source_y_ind=source_y_ind,
                     source_decomp_ind=source_decomp_ind,
                     source_time_offset=source_time_offset,
                     source2outlet_ind=source2outlet_ind,
                     unit_hydrograph=unit_hydrograph,
                     **ncvaropts)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write a summary of what was done to the log file.
    log.info('Parameter file includes %i outlets', len(outlets))
    log.info('Parameter file includes %i Source Points', len(source_lon))
    # ---------------------------------------------------------------- #

    return param_file, today
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def adjust_fractions(outlets, dom_fractions, adjust=True):
    """
    Constrain the fractions in the outles.
    The basic idea is that the sum of fraction from the outlets should not
    exceed the domain fractions.
    """

    log.info('Adjusting fractions now')

    fractions = np.zeros(dom_fractions.shape, dtype=np.float64)
    ratio_fraction = np.ones(fractions.shape, dtype=np.float64)
    adjusted_fractions = np.zeros(dom_fractions.shape, dtype=np.float64)
    sum_uh_fractions = np.zeros(dom_fractions.shape, dtype=np.float64)

    # ---------------------------------------------------------------- #
    # Aggregate the fractions
    for cell_id, outlet in outlets.iteritems():
        y = outlet.y_source
        x = outlet.x_source

        fractions[y, x] += outlet.fractions
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # First set fractions to zero where there is no land in the domain
    yd, xd = np.nonzero(dom_fractions == 0.0)
    fractions[yd, xd] = 0.0
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Only adjust fractions where the aggregated fractions are gt the domain
    # fractions
    yi, xi = np.nonzero(fractions > dom_fractions)
    log.info('Adjust fractions for %s grid cells', len(yi))
    ratio_fraction[yi, xi] = dom_fractions[yi, xi]/fractions[yi, xi]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust fracs based on ratio_fraction
    for cell_id, outlet in outlets.iteritems():
        y = outlet.y_source
        x = outlet.x_source
        if adjust:
            outlets[cell_id].fractions *= ratio_fraction[y, x]
        # For Diagnostics only
        adjusted_fractions[y, x] += outlets[cell_id].fractions
        sum_uh_fractions[y, x] += outlets[cell_id].unit_hydrograph.sum(axis=0)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make Fractions Dict for plotting
    plot_dict = {'Domain Fractions': dom_fractions,
                 'Aggregated Fractions': fractions,
                 'Ratio Fractions': ratio_fraction,
                 'Adjusted Fractions': adjusted_fractions,
                 'Sum UH Before': sum_uh_fractions}
    # ---------------------------------------------------------------- #
    return outlets, plot_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Shorten the unit hydrograph
def subset(outlets, subset_length=None):
    """ Shorten the Unit Hydrograph"""

    log.info('subsetting unit-hydrographs now...')
    log.debug('Subset Length:  %s' % subset_length)

    for i, (cell_id, outlet) in enumerate(outlets.iteritems()):
        if i == 0:
            full_time_length = outlet.unit_hydrograph.shape[0]
            log.debug('Subset Length:  %s' % subset_length)
            if not subset_length:
                subset_length = full_time_length
                log.debug('No subset_length provided, using full_time_length')
            before = outlets[cell_id].unit_hydrograph
        else:
            before = np.append(before, outlets[cell_id].unit_hydrograph,
                               axis=1)

        outlets[cell_id].offset = np.empty(outlet.unit_hydrograph.shape[1],
                                           dtype=np.int16)
        out_uh = np.zeros((subset_length, outlet.unit_hydrograph.shape[1]),
                          dtype=np.float64)

        d_left = -1*subset_length/2
        d_right = subset_length/2

        for j in xrange(outlet.unit_hydrograph.shape[1]):
            # find index position of maximum
            maxind = np.argmax(outlet.unit_hydrograph[:, j])

            # find bounds
            left = maxind + d_left
            right = maxind + d_right

            # make sure left and right fit in unit hydrograph array,
            # if not adjust
            if left < 0:
                left = 0
                right = subset_length
            if right > full_time_length:
                right = full_time_length
                left = full_time_length - subset_length

                log.warning('Subset centered on UH max extends beyond length '
                            'of unit hydrograph.')
                log.warning('--> Outlet %s', outlet)
                log.warning('----> Max Index is %s', maxind)
                log.warning('----> Last value in subset '
                            'is %s', outlet.unit_hydrograph[-1, j])
                if maxind == full_time_length:
                    log.warning('maxind == full_time_length, not able to '
                                'resolve unithydrograph')
            if left < 0 or right > full_time_length:
                raise ValueError('Subsetting failed left:{0} or right {1} does'
                                 ' not fit inside bounds'.format(left, right))

            outlets[cell_id].offset[j] = left

            # clip and normalize
            tot = outlet.unit_hydrograph[left:right, j].sum()
            out_uh[:, j] = outlet.unit_hydrograph[left:right, j] / tot

        outlets[cell_id].unit_hydrograph = out_uh

        if i == 0:
            after = outlets[cell_id].unit_hydrograph
        else:
            after = np.append(after, outlets[cell_id].unit_hydrograph, axis=1)

    log.info('Done subsetting')

    return outlets, full_time_length, before, after
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
            outlet_lon = np.array(outlet.lon, dtype=np.float64)
            outlet_lat = np.array(outlet.lat, dtype=np.float64)
            outlet_x_ind = np.array(outlet.domx, dtype=np.int16)
            outlet_y_ind = np.array(outlet.domy, dtype=np.int16)
            outlet_decomp_ind = np.array(cell_id, dtype=np.int16)
            outlet_number = np.array(i, dtype=np.int16)
            outlet_name = np.array(outlet.name)
            outlet_upstream_gridcells = np.array(outlet.upstream_gridcells,
                                                 dtype=np.int16)
            outlet_upstream_area = np.array(outlet.upstream_area,
                                            dtype=np.float64)
        else:
            # -------------------------------------------------------- #
            # Point specific values
            unit_hydrograph = np.append(unit_hydrograph,
                                        outlet.unit_hydrograph, axis=1)
            frac_sources = np.append(frac_sources, outlet.fractions)
            source_lon = np.append(source_lon, outlet.lon_source)
            source_lat = np.append(source_lat, outlet.lat_source)
            source_x_ind = np.append(source_x_ind, x)
            source_y_ind = np.append(source_y_ind, y)
            source_decomp_ind = np.append(source_decomp_ind,
                                          outlet.cell_id_source)
            source_time_offset = np.append(source_time_offset, outlet.offset)
            source2outlet_ind = np.append(source2outlet_ind,
                                          np.zeros_like(outlet.offset) + i)
            # -------------------------------------------------------- #

            # -------------------------------------------------------- #
            # outlet specific inputs
            outlet_lon = np.append(outlet_lon, outlet.lon)
            outlet_lat = np.append(outlet_lat, outlet.lat)
            outlet_x_ind = np.append(outlet_x_ind, outlet.domx)
            outlet_y_ind = np.append(outlet_y_ind, outlet.domy)
            outlet_decomp_ind = np.append(outlet_decomp_ind, cell_id)
            outlet_number = np.append(outlet_number, i)
            outlet_name = np.append(outlet_name, outlet.name)
            outlet_upstream_gridcells = np.append(outlet_upstream_gridcells,
                                                  outlet.upstream_gridcells)
            outlet_upstream_area = np.append(outlet_upstream_area,
                                             outlet.upstream_area)
            # -------------------------------------------------------- #

    grouped_data = {'unit_hydrograph': unit_hydrograph,
                    'frac_sources': frac_sources,
                    'source_lon': source_lon,
                    'source_lat': source_lat,
                    'source_x_ind': source_x_ind,
                    'source_y_ind': source_y_ind,
                    'source_decomp_ind': source_decomp_ind,
                    'source_time_offset': source_time_offset,
                    'source2outlet_ind': source2outlet_ind,
                    'outlet_lon': outlet_lon,
                    'outlet_lat': outlet_lat,
                    'outlet_x_ind': outlet_x_ind,
                    'outlet_y_ind': outlet_y_ind,
                    'outlet_decomp_ind': outlet_decomp_ind,
                    'outlet_number': outlet_number,
                    'outlet_name': outlet_name,
                    'outlet_upstream_gridcells': outlet_upstream_gridcells,
                    'outlet_upstream_area': outlet_upstream_area}

    return grouped_data
# -------------------------------------------------------------------- #
