"""
convert.py
"""
import os
import re
import numpy as np
import logging
from log import LOG_NAME
from variables import Point
from share import FILLVALUE_I

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(LOG_NAME)
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# read station file
def read_station_file(file_name, dom_data, config_dict):
    """
    Read a standard routing station file
    http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Documentation/Routing/StationLocation.shtml
    """
    outlets = {}

    f = open(file_name, 'r')
    i = -1
    while True:
        i += 1
        line1 = re.sub(' +',' ',f.readline())
        if not line1:
            break
        active, name, x, y, area = line1.split(' ')
        uhs_file = f.readline().strip()

        # make sure files exist
        log.info('On station: %s, active: %s' %(name, active))
        if active == '1':
            if os.path.isfile(uhs_file):
                active = True
            elif os.path.isfile(os.path.join(config_dict['UHS_FILES']['ROUT_DIR'], name+'.uh_s2')):
                uhs_file = os.path.join(config_dict['UHS_FILES']['ROUT_DIR'], name+'.uh_s2')
            else:
                raise ValueError('missing uhs_file: (%s or %s)' %(uhs_file, os.path.join(config_dict['UHS_FILES']['ROUT_DIR'], name+'.uh_s2')))
        else:
            active = False

        if active:
            outlets[i] = Point(x=int(x), y=int(y))
            outlets[i].name = name
            outlets[i].area = float(area)
            outlets[i].uhs_file = uhs_file
            outlets[i].cell_id = i
            outlets[i].outlet_decomp_ind = dom_data['cell_ids'][y, x]
            outlets[i].lon = dom_data[config_dict['DOMAIN']['LONGITUDE_VAR']][y, x]
            outlets[i].lat = dom_data[config_dict['DOMAIN']['LATITUDE_VAR']][y, x]
        else:
            log.info('%s not active... skipping' %name)
    f.close()
    return outlets
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Read uhs files
def read_uhs_files(outlets, dom_data, config_dict):
    """
    Read a standard routing uh_s file
    Format:
    line0: num_sources
    line1: active lon lat fraction x y
    line2: unit_hydrograph_time_series
    line1: ...
    line2: ...
    """
    if config_dict['UHS_FILES']['ROUT_PROGRAM'] == 'C':
        for cell_id, outlet in outlets.iteritems():
            log.info('Reading outlet %i: %s' %(cell_id, outlet.name))
            log.debug(outlet.uhs_file)
            f = open(outlet.uhs_file, 'r')
            num_sources = int(f.readline())
            log.debug('Number of sources in file: %i' %num_sources)
            # setup some empty arrays
            outlets[cell_id].lon_source = np.empty(num_sources)
            outlets[cell_id].lat_source = np.empty(num_sources)
            outlets[cell_id].fractions = np.empty(num_sources)
            outlets[cell_id].x_source = np.empty(num_sources, dtype=int)
            outlets[cell_id].y_source = np.empty(num_sources, dtype=int)

            uh = []

            # loop over the source points
            for j in xrange(num_sources):
                line = re.sub(' +',' ',f.readline())
                lon, lat, fracs, x, y = line.split()
                outlets[cell_id].lon_source[j] = float(lon)
                outlets[cell_id].lat_source[j] = float(lat)
                outlets[cell_id].fractions[j] = float(fracs)
                outlets[cell_id].x_source[j] = int(x)
                outlets[cell_id].y_source[j] = int(y)
                line = re.sub(' +',' ',f.readline())
                uh.append(map(float, line.split()))

            outlets[cell_id].unit_hydrograph = np.rot90(np.array(uh), k=-1)
            outlets[cell_id].source_decomp_ind = dom_data['cell_ids'][outlets[cell_id].y_source, outlets[cell_id].x_source]
            f.close()

            outlets[cell_id].cell_id_source = dom_data['cell_ids'][outlets[cell_id].y_source, outlets[cell_id].x_source]

    elif config_dict['UHS_FILES']['ROUT_PROGRAM'] == 'Fortran':
        raise ValueError('Fortran conversion not working...')
        # log.info('parsing fortran uhs files')
        # # setup for finding x, y inds
        # dom_lon = dom_data[config_dict['DOMAIN']['LONGITUDE_VAR']]
        # dom_lat = dom_data[config_dict['DOMAIN']['LATITUDE_VAR']]
        # combined = np.dstack(([dom_lat.ravel(), dom_lon.ravel()]))[0]
        # mytree = cKDTree(combined)

        # for cell_id, outlet in outlets.iteritems():
        #     # read lons, lats, and unit hydrographs
        #     log.debug('reading %s' %outlet.ll_file)
        #     lons, lats = np.genfromtxt(outlet.ll_file, dtype="f8", unpack=True)

        #     log.debug('reading %s' %outlet.uhs_file)
        #     uh = np.genfromtxt(outlet.uhs_file, dtype="f8")
        #     outlets[cell_id].unit_hydrograph = np.rot90(uh, k=-1)
        #     if len(lons) != uh.shape[0]:
        #         raise ValueError('length mismatch in ll file and uhs file %s %s' %(lons.shape, uh.shape))

        #     # now find the y_source and x_sources
        #     points = list(np.vstack((np.array(lats), np.array(lons))).transpose())
        #     dist, indexes = mytree.query(points, k=1)
        #     yinds, xinds = np.unravel_index(np.array(indexes), dom_lat.shape)

        #     # finally, get the fractions and the source_decomp_ind from the domain data
        #     outlets[cell_id].cell_id_source = dom_data['cell_ids'][yinds, xinds]
        #     outlets[cell_id].fractions = dom_data[config_dict['DOMAIN']['FRACTION_VAR']][yinds, xinds]
        #     print 'fractions', outlets[cell_id].fractions

        #     # now store all the data
        #     outlets[cell_id].lat_source = lats
        #     outlets[cell_id].lon_source = lons
        #     outlets[cell_id].y_source = yinds
        #     outlets[cell_id].x_source = xinds
    else:
        raise ValueError('UHS_FILES[ROUT_PROGRAM] must be either C or Fortran')

    return outlets
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Adjust Domain Bounds
def move_domain(dom_data, new_dom_data, outlets):

    # ---------------------------------------------------------------- #
    # Create lookup arrays (size of dom_data but contains mapping to new_dom)
    if dom_data['cord_lons'].ndim == 1:
        new_y = np.zeros(len(dom_data['cord_lats']), dtype=np.int) - FILLVALUE_I
        new_x = np.zeros(len(dom_data['cord_lons']), dtype=np.int) - FILLVALUE_I
    else:
        raise ValueError('Grids must be regular and coordinate variables must have only 1 dimension to move domain to smaller size')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Check that lons/lats in new_domain match those in original_domain
    xi = np.searchsorted(dom_data['cord_lons'], new_dom_data['cord_lons'])
    yi = np.searchsorted(dom_data['cord_lats'], new_dom_data['cord_lats'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Fill in the lookup arrays with the correct value
    new_y[yi] = np.arange(len(yi), dtype=np.int)
    new_x[xi] = np.arange(len(xi), dtype=np.int)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust locations
    for cell_id, outlet in outlets.iteritems():

        outlets[cell_id].y = new_y[outlet.y]
        outlets[cell_id].x = new_x[outlet.x]

        outlets[cell_id].y_source = new_y[outlet.y_source]
        outlets[cell_id].x_source = new_x[outlet.x_source]
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Return outles
    log.info('Finished moving indicies to new domain')
    return outlets
    # ---------------------------------------------------------------- #
# -------------------------------------------------------------------- #
