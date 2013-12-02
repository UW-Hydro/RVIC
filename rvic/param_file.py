"""
param_file.py
"""
import numpy as np
import logging
from log import LOG_NAME
from write import write_param_file
from share import NcGlobals, PRECISION, SECSPERDAY
import os
from datetime import date
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

    # ---------------------------------------------------------------- #
    # subset (shorten time base)
    if options['SUBSET_DAYS'] and options['SUBSET_DAYS']<config_dict['ROUTING']['BASIN_FLOWDAYS']:
        subset_length = options['SUBSET_DAYS']*SECSPERDAY/config_dict['ROUTING']['OUTPUT_INTERVAL']
        outlets, full_time_length, before, after = subset(outlets,
                                           subset_length=subset_length)

        log.debug('plotting unit hydrograph timeseries now for before / after subseting')

        title = 'UHS before subset'
        pfname = plot_uhs(before, title, options['CASEID'], directories['plots'])
        log.info('%s Plot:  %s', title, pfname)

        title = 'UHS after subset'
        pfname = plot_uhs(after, title, options['CASEID'], directories['plots'])
        log.info('%s Plot:  %s', title, pfname)
    else:
        subset_length = config_dict['ROUTING']['BASIN_FLOWDAYS']*SECSPERDAY/config_dict['ROUTING']['OUTPUT_INTERVAL']
        log.info('Not subsetting because either SUBSET_DAYS is null or SUBSET_DAYS<BASIN_FLOWDAYS')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # adjust fractions
    if options['CONSTRAIN_FRACTIONS']:
        adjust=True
        log.info('Adjusting Fractions to be less than or equal to domain fractions')
    else:
        adjust=False
    outlets, plot_dict = adjust_fractions(outlets,
                                          dom_data[config_dict['DOMAIN']['FRACTION_VAR']],
                                          adjust=adjust)
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
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust Unit Hydrographs for differences in source/outlet areas and fractions
    area = dom_data[config_dict['DOMAIN']['AREA_VAR']]

    for source, outlet in enumerate(source2outlet_ind):
        unit_hydrograph[:, source] *= area[source_y_ind[source], source_x_ind[source]]
        unit_hydrograph[:, source] /= area[outlet_y_ind[outlet], outlet_x_ind[outlet]]
        unit_hydrograph[:, source] *= frac_sources[source]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make diagnostic plots
    sum_after = np.zeros(dom_data[config_dict['DOMAIN']['FRACTION_VAR']].shape)
    for i, (y, x) in enumerate(zip(source_y_ind, source_x_ind)):
        sum_after[y, x] += unit_hydrograph[:, i].sum()

    plot_dict['Sum UH Final'] = sum_after

    for title, data in plot_dict.iteritems():
        pfname = plot_fractions(data, title, options['CASEID'], directories['plots'])
        log.info('%s Plot:  %s', title, pfname)
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

    if 'NEW_DOMAIN' in config_dict.keys():
        dom_file_name = config_dict['NEW_DOMAIN']['FILE_NAME']
    else:
        dom_file_name = config_dict['DOMAIN']['FILE_NAME']

    write_param_file(param_file,
                     nc_format=options['NETCDF_FORMAT'],
                     glob_atts=NcGlobals(title='RVIC parameter file',
                                         RvicPourPointsFile=os.path.split(config_dict['POUR_POINTS']['FILE_NAME'])[1],
                                         RvicUHFile=os.path.split(config_dict['UH_BOX']['FILE_NAME'])[1],
                                         RvicFdrFile=os.path.split(config_dict['ROUTING']['FILE_NAME'])[1],
                                         RvicDomainFile=os.path.split(dom_file_name)[1]),
                     full_time_length=full_time_length,
                     subset_length=subset_length,
                     unit_hydrograph_dt=config_dict['ROUTING']['OUTPUT_INTERVAL'],
                     outlet_lon=outlet_lon,
                     outlet_lat=outlet_lat,
                     outlet_x_ind=outlet_x_ind,
                     outlet_y_ind=outlet_y_ind,
                     outlet_decomp_ind=outlet_decomp_ind,
                     outlet_number=outlet_number,
                     outlet_mask=outlet_mask,
                     outlet_name=outlet_name,
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
def adjust_fractions(outlets, dom_fractions, adjust=True):
    """
    Constrain the fractions in the outles.
    The basic idea is that the sum of fraction from the outlets should not
    exceed the domain fractions.
    """

    log.info('Adjusting fractions now')

    fractions = np.zeros(dom_fractions.shape)
    ratio_fraction = np.ones(fractions.shape)
    adjusted_fractions = np.zeros(dom_fractions.shape)
    sum_uh_fractions = np.zeros(dom_fractions.shape)

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
    # Only adjust fractions where the aggregated fractions are gt the domain fractions
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
    plot_dict= {'Domain Fractions':dom_fractions,
                'Aggregated Fractions':fractions,
                'Ratio Fractions':ratio_fraction,
                'Adjusted Fractions':adjusted_fractions,
                'Sum UH Before':sum_uh_fractions}
    # ---------------------------------------------------------------- #

    return outlets, plot_dict

def plot_uhs(data, title, case_id, plot_dir):
    """
    Plot diagnostic plot showing all unit hydrographs
    """
    today = date.today().strftime('%Y%m%d')
    file_name = "{}_{}_{}.png".format(title.lower().replace(" ", "_"),
                                      case_id.lower().replace(" ", "_"),
                                      today)
    pfname = os.path.join(plot_dir, file_name)

    fig = plt.figure()
    plt.plot(data)
    plt.title(title)
    plt.xlabel('timesteps')
    plt.ylabel('unit-hydrograph')
    fig.savefig(pfname)

    return pfname

def plot_fractions(data, title, case_id, plot_dir):
    """
    Plot diagnostic plots of fraction variables
    """
    # ---------------------------------------------------------------- #
    # Plot Fractions
    today = date.today().strftime('%Y%m%d')
    file_name = "{}_{}_{}.png".format(title.lower().replace(" ", "_"),
                                      case_id.lower().replace(" ", "_"),
                                      today)
    pfname = os.path.join(plot_dir, file_name)

    fig = plt.figure()
    plt.pcolormesh(data)
    plt.autoscale(tight=True)
    plt.axis('tight')
    plt.colorbar()
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    fig.savefig(pfname)

    return pfname
    # ---------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Shorten the unit hydrograph
def subset(outlets, subset_length=None):
    """ Shorten the Unit Hydrograph"""

    log.info('subsetting unit-hydrographs now...')
    log.debug('Subset Length:  %s' % subset_length)

    for i, (cell_id, outlet) in enumerate(outlets.iteritems()):
        if i == 0:
            full_time_length = outlet.unit_hydrograph.shape[0]
            if not subset_length:
                subset_length = full_time_length
                log.debug('No subset_length provided, using full_time_length')
            before = outlets[cell_id].unit_hydrograph
        else:
            before = np.append(before, outlets[cell_id].unit_hydrograph, axis=1)

        outlets[cell_id].offset = np.empty(outlet.unit_hydrograph.shape[1], dtype=int)
        out_uh = np.zeros((subset_length, outlet.unit_hydrograph.shape[1]))

        d_left = -1*subset_length/2
        d_right = subset_length/2


        for j in xrange(outlet.unit_hydrograph.shape[1]):
            # find index position of maximum
            maxind = np.argmax(outlet.unit_hydrograph[:, j])

            # find bounds
            left = maxind + d_left
            right = maxind + d_right

            # make sure left and right fit in unit hydrograph array, if not adjust
            if left < 0:
                left = 0
                right = subset_length
            if right > full_time_length:
                right = full_time_length
                left = full_time_length - subset_length

                log.warning('Subset centered on UH max extends beyond length of unit hydrograph.')
                log.warning('--> Outlet %s' % outlet)
                log.warning('----> Max Index is %s' % maxind)
                log.warning('----> Last value in subset is %s' % outlet.unit_hydrograph[-1, j])
                if maxind == full_time_length:
                    log.warning('maxind == full_time_length, not able to resolve unithydrograph')

            outlets[cell_id].offset[j] = left

            # clip and normalize
            out_uh[:, j] = outlet.unit_hydrograph[left:right, j] / outlet.unit_hydrograph[left:right, j].sum()

        outlets[cell_id].unit_hydrograph = out_uh

        if i == 0:
            after = outlets[cell_id].unit_hydrograph
        else:
            after = np.append(after, outlets[cell_id].unit_hydrograph, axis=1)


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
            outlet_lon = np.array(outlet.lon)
            outlet_lat = np.array(outlet.lat)
            outlet_x_ind = np.array(outlet.gridx)
            outlet_y_ind = np.array(outlet.gridy)
            outlet_decomp_ind = np.array(cell_id)
            outlet_number = np.array(i)
            outlet_name = np.array(outlet.name)
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
            outlet_x_ind = np.append(outlet_x_ind, outlet.gridx)
            outlet_y_ind = np.append(outlet_y_ind, outlet.gridy)
            outlet_decomp_ind = np.append(outlet_decomp_ind, cell_id)
            outlet_number = np.append(outlet_number, i)
            outlet_name = np.append(outlet_name, outlet.name)
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
                    'outlet_name': outlet_name}

    return grouped_data
# -------------------------------------------------------------------- #
