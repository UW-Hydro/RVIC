# -*- coding: utf-8 -*-
'''
utilities.py
'''
import os
import tarfile
from scipy.spatial import cKDTree
import numpy as np
from shutil import rmtree, copyfile
from .pycompat import iteritems, SafeConfigParser
from netCDF4 import Dataset
from logging import getLogger
from .log import LOG_NAME
from .share import TIMESTAMPFORM, RPOINTER, EARTHRADIUS, METERSPERMILE
from .share import METERS2PERACRE, METERSPERKM, VALID_CHARS
from .config import read_config
from .pycompat import pyzip

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# find x y coordinates
def latlon2yx(plats, plons, glats, glons):
    '''find y x coordinates '''

    # use astronomical conventions for longitude
    # (i.e. negative longitudes to the east of 0)
    if (glons.max() > 180):
        posinds = np.nonzero(glons > 180)
        glons[posinds] -= 360
        log.info('adjusted grid lon to astronomical conventions')
    if (plons.max() > 180):
        posinds = np.nonzero(plons > 180)
        plons[posinds] -= 360
        log.info('adjusted point lon to astronomical conventions')

    if glons.ndim == 1 or glats.ndim == 1:
        glons, glats = np.meshgrid(glons, glats)

    combined = np.dstack(([glats.ravel(), glons.ravel()]))[0]
    points = list(np.vstack((np.array(plats), np.array(plons))).transpose())

    mytree = cKDTree(combined)
    indexes = mytree.query(points, k=1)[1]
    y, x = np.unravel_index(indexes, glons.shape)
    return y, x

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Search neighboring grid cells for channel
def search_for_channel(source_area, routys, routxs, search=1, tol=10):
    '''Search neighboring grid cells for channel'''

    log.debug('serching for channel, tol: %f, search: %i', tol, search)

    new_ys = np.copy(routys)
    new_xs = np.copy(routxs)

    ysize, xsize = source_area.shape

    for i, (y, x) in enumerate(pyzip(routys, routxs)):
        area0 = source_area[y, x]

        for j in range(search + 1):
            ymin = np.clip(y - j, 0, ysize)
            ymax = np.clip(y + j + 1, 0, ysize)
            xmin = np.clip(x - j, 0, xsize)
            xmax = np.clip(x + j + 1, 0, xsize)

            search_area = source_area[ymin:ymax, xmin:xmax]

            if np.any(search_area / area0 > tol):
                sy, sx = np.unravel_index(search_area.argmax(),
                                          search_area.shape)

                new_ys[i] = np.clip(y + sy - j, 0, ysize)
                new_xs[i] = np.clip(x + sx - j, 0, xsize)

                log.debug('Moving pour point to channel y: %s->%s, x: %s->%s',
                          y, new_ys[i], x, new_xs[i])
                log.debug('Source Area has increased from %s to %s',
                          area0, source_area[new_ys[i], new_xs[i]])

                break

    return new_ys, new_xs

# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write rpointer file
def write_rpointer(restart_dir, restart_file, timestamp):
    ''' Write a configuration file with restart file and time '''
    rpointer_file = os.path.join(restart_dir, RPOINTER)

    config = SafeConfigParser()
    config.optionxform = str

    time_str = timestamp.strftime(TIMESTAMPFORM)

    config.add_section('RESTART')
    config.set('RESTART', 'FILE_NAME', os.path.join(restart_dir, restart_file))
    config.set('RESTART', 'TIMESTAMP', time_str)

    with open(rpointer_file, 'w') as configfile:
        config.write(configfile)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# A helper function to read a netcdf file
def read_netcdf(nc_file, variables=None, coords=None):
    '''
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named
    by variable
    '''

    f = Dataset(nc_file, 'r')

    if not variables:
        variables = list(f.variables.keys())
    if not coords:
        coords = slice(None)

    log.debug('Reading input data variables: %s, from file: %s', variables,
              nc_file)

    d = {}
    a = {}
    g = {}

    for var in variables:
        d[var] = f.variables[var][coords]
        a[var] = f.variables[var].__dict__

    for attr in f.ncattrs():
        g[attr] = getattr(f, attr)

    f.close()

    return d, a, g
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Check to make sure all the expected variables are present in the dictionary
def check_ncvars(config_section, nckeys):
    '''
    Make sure the variables listed in the config file are present in the netcdf
    '''
    for key, value in iteritems(config_section):
        if key.endswith('var'):
            if value not in nckeys:
                log.error('%s (%s) not in %s', value, key,
                          config_section['FILE_NAME'])
                raise NameError('Check netcdf that netcdf variable names match'
                                ' those in the configuration file')
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the index of the the nearest value
def find_nearest(array, value):
    ''' Find the index location in (array) with value nearest to (value)'''
    return np.abs(array - value).argmin()
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Delete all the files in a directory
def clean_dir(directory):
    ''' Clean all files in a directory'''
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception:
            log.exception('Error cleaning file: %s', file_path)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Delete a particular file
def clean_file(file_name):
    ''' Delete the file'''
    try:
        if os.path.isfile(file_name):
            os.unlink(file_name)
    except Exception:
        log.exception('Error cleaning file: %s', file_name)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make a set of directories
def make_directories(rundir, subdir_names):
    '''Make rvic directory structure'''
    if not os.path.exists(rundir):
        os.makedirs(rundir)

    paths = {}
    for s in subdir_names:
        paths[s] = os.path.join(rundir, s)
        if not os.path.exists(paths[s]):
            os.makedirs(paths[s])
    return paths
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Move all the input files to a central location
def copy_inputs(config_file, inputs_dir):

    config_dict = read_config(config_file)

    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)

    new_config = os.path.join(inputs_dir, os.path.split(config_file)[1])

    # ---------------------------------------------------------------- #
    # copy the inputs
    for key, section in iteritems(config_dict):
        if 'FILE_NAME' in list(section.keys()):
            new_file_name = os.path.join(
                inputs_dir, os.path.split(section['FILE_NAME'])[1])

            copyfile(section['FILE_NAME'], new_file_name)

            # update the config file for an easy restart
            config.set(key, 'FILE_NAME',
                       os.path.join(inputs_dir,
                                    os.path.split(section['FILE_NAME'])[1]))

            # update the config_dict with the new value
            config_dict[key]['FILE_NAME'] = new_file_name
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write the new configuration file
    with open(new_config, 'w') as configfile:
        config.write(configfile)
    # ---------------------------------------------------------------- #

    return config_dict
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def tar_inputs(inputs, suffix='', tar_type='tar'):
    ''' Tar the inputss directory or file at the end of a run'''
    # ---------------------------------------------------------------- #
    # Make the TarFile
    if tar_type == 'tar':
        end = '.tar'
        mode = 'w:'
    elif tar_type in ['tgz', 'tar.gz', 'gunzip']:
        end = '.tgz'
        mode = 'w:'
    else:
        log.warning('Unknown tar_type: %s, proceeding with gunzipped mode',
                    tar_type)
        end = '.tgz'
        mode = 'w:'

    tar_file = inputs + suffix + end
    log.info('tarfile: %s', tar_file)

    if os.path.isdir(inputs):
        arcname = os.path.basename(os.path.normpath(inputs))
    else:
        arcname = os.path.split(inputs)[1]

    with tarfile.open(tar_file, mode) as tar:
        tar.add(inputs, arcname=arcname)

    # ---------------------------------------------------------------- #
    # Check to make sure the TarFile exists before deleting the sources
    if os.path.exists(tar_file):
        # ------------------------------------------------------------ #
        # Remove the inputs
        if os.path.isdir(inputs):
            rmtree(inputs)
        elif os.path.isfile(inputs):
            os.unlink(inputs)
        # ------------------------------------------------------------ #
    else:
        log.error('Problem removing inputs: %s', inputs)
    # ---------------------------------------------------------------- #
    return tar_file
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the domain
def read_domain(domain_dict, lat0_is_min=False):
    '''
    Read the domain file and return all the variables and attributes.
    Area is returned in m2
    '''
    dom_data, dom_vatts, dom_gatts = read_netcdf(domain_dict['FILE_NAME'])

    check_ncvars(domain_dict, list(dom_data.keys()))

    # ---------------------------------------------------------------- #
    # Create the cell_ids variable
    dom_mask = domain_dict['LAND_MASK_VAR']
    temp = np.arange(dom_data[dom_mask].size)
    dom_data['cell_ids'] = temp.reshape(dom_data[dom_mask].shape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make sure the longitude / latitude vars are 2d
    dom_lat = domain_dict['LATITUDE_VAR']
    dom_lon = domain_dict['LONGITUDE_VAR']

    dom_data['cord_lons'] = dom_data[dom_lon][:]
    dom_data['cord_lats'] = dom_data[dom_lat][:]

    if dom_data[dom_lon].ndim == 1:
        # ------------------------------------------------------------- #
        # Check latitude order, flip if necessary.
        if (dom_data[dom_lat][-1] > dom_data[dom_lat][0]) != lat0_is_min:
            log.debug('Domain Inputs came in upside down, flipping everything '
                      'now.')
            var_list = list(dom_data.keys())
            var_list.remove(dom_lon)
            for var in var_list:
                dom_data[var] = np.flipud(dom_data[var])
        # ------------------------------------------------------------ #

        # ------------------------------------------------------------- #
        # Make 2d coordinate vars
        dom_data[dom_lon], dom_data[dom_lat] = np.meshgrid(dom_data[dom_lon],
                                                           dom_data[dom_lat])
        # ------------------------------------------------------------- #
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make sure the area is in m2
    dom_area = domain_dict['AREA_VAR']
    area_units = dom_vatts[dom_area]['units']

    if area_units in ['rad2', 'radians2', 'radian2', 'radian^2', 'rad^2',
                      'radians^2', 'rads^2', 'radians squared',
                      'square-radians']:
        dom_data[dom_area] = dom_data[dom_area] * EARTHRADIUS * EARTHRADIUS
    elif area_units in ['m2', 'm^2', 'meters^2', 'meters2', 'square-meters',
                        'meters squared']:
        dom_data[dom_area] = dom_data[dom_area]
    elif area_units in ['km2', 'km^2', 'kilometers^2', 'kilometers2',
                        'square-kilometers', 'kilometers squared']:
        dom_data[dom_area] = dom_data[dom_area] * METERSPERKM * METERSPERKM
    elif area_units in ['mi2', 'mi^2', 'miles^2', 'miles', 'square-miles',
                        'miles squared']:
        dom_data[dom_area] = dom_data[dom_area] * METERSPERMILE * METERSPERMILE
    elif area_units in ['acres', 'ac', 'ac.']:
        dom_data[dom_area] = dom_data[dom_area] * METERS2PERACRE
    else:
        log.warning('WARNING: UNKNOWN AREA units (%s), ASSUMING THEY ARE IN '
                    'SQUARE METERS',
                    dom_data[domain_dict['AREA_VAR']]['units'])
    # ---------------------------------------------------------------- #
    return dom_data, dom_vatts, dom_gatts
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def strip_non_ascii(in_string):
    ''' Returns the string without non ASCII characters'''
    stripped = (c for c in in_string if 0 < ord(c) < 127)
    return ''.join(stripped)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def strip_invalid_char(in_string):
    ''' Returns the string without invalid characters for filenames'''
    return ''.join(c for c in in_string if c in VALID_CHARS)
# -------------------------------------------------------------------- #
