"""
utilities.py
"""
import os
import tarfile
from scipy.spatial import cKDTree
import numpy as np
from shutil import rmtree, copyfile
from ConfigParser import SafeConfigParser
from netCDF4 import Dataset
from cdo import Cdo
cdo = Cdo()
from logging import getLogger
from log import LOG_NAME
from share import TIMESTAMPFORM, RPOINTER, EARTHRADIUS, METERSPERMILE, METERS2PERACRE, METERSPERKM

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# find x y coordinates
def latlon2yx(plats, plons, glats, glons):
    """find y x coordinates """

    if glons.ndim == 1 or glats.ndim == 1:
        glons, glats = np.meshgrid(glons, glats)

    combined = np.dstack(([glats.ravel(), glons.ravel()]))[0]
    points = list(np.vstack((np.array(plats), np.array(plons))).transpose())

    mytree = cKDTree(combined)
    dist, indexes = mytree.query(points, k=1)
    y, x = np.unravel_index(indexes, glons.shape)
    return y, x

# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Write rpointer file
def write_rpointer(restart_dir, restart_file, timestamp):
    """ Write a configuration file with restart file and time """
    rpointer_file = os.path.join(restart_dir, RPOINTER)

    config = SafeConfigParser()
    config.optionxform=str

    time_str = timestamp.strftime(TIMESTAMPFORM)

    config.add_section('RESTART')
    config.set('RESTART', 'FILE_NAME', os.path.join(restart_dir, restart_file))
    config.set('RESTART', 'TIMESTAMP', time_str)

    with open(rpointer_file, 'w') as configfile:
        config.write(configfile)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# A helper function to read a netcdf file
def read_netcdf(nc_file, variables=None, coords=None):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """

    f = Dataset(nc_file, 'r')

    if not variables:
        variables = f.variables.keys()
    if not coords:
        coords = slice(None)

    log.debug('Reading input data variables: %s, from file: %s' % (variables, nc_file))

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
    """Make sure the variables listed in the config file are present in the netcdf"""
    for key, value in config_section.iteritems():
        if key.endswith('var'):
            if value not in nckeys:
                log.error('%s (%s) not in %s' % (value, key, config_section['FILE_NAME']))
                raise NameError('Check netcdf that netcdf variable names match those in the configuration file')
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the Configuration File
def read_config(config_file):
    """
    Return a dictionary with subdictionaries of all configFile options/values
    """
    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    sections = config.sections()
    dict1 = {}
    for section in sections:
        options = config.options(section)
        dict2 = {}
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        dict1[section] = dict2
    return dict1
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the type of the config options
def config_type(value):
    """
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    """
    val_list = [x.strip() for x in value.split(',')]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE', 'T']:
            return True
        elif value in ['false', 'False', 'FALSE', 'F']:
            return False
        elif value in ['none', 'None', 'NONE', '']:
            return None
        else:
            try:
                return float(value)
            except:
                return value
    else:
        try:
            return map(float, val_list)
        except:
            return val_list
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the index of the the nearest value
def find_nearest(array, value):
    """ Find the index location in (array) with value nearest to (value)"""
    return np.abs(array-value).argmin()
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Delete all the files in a directory
def clean_dir(directory):
    """ Clean all files in a directory"""
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except:
            log.exception('Error cleaning file: %s' % file_path)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Delete a particular file
def clean_file(file_name):
    """ Delete the file"""
    try:
        if os.path.isfile(file_name):
            os.unlink(file_name)
    except:
        log.exception('Error cleaning file: %s' % file_name)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Remap a file using CDO
def remap(grid_file, in_file, out_file, operator='remapcon',
          remap_options=None):
    """Remap infile using cdo"""

    remap_method = getattr(cdo, operator)

    if remap_options:
        remap_method(grid_file, input=in_file, output=out_file, options=remap_options)
    else:
        remap_method(grid_file, input=in_file, output=out_file)

    log.debug('remapped to out file: %s' % os.path.split(out_file)[1])

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make a set of directories
def make_directories(rundir, subdir_names):
    """Make rvic directory structure"""
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
def copy_inputs(config_file, InputsDir):

    config_dict = read_config(config_file)

    config = SafeConfigParser()
    config.read(config_file)

    new_config = os.path.join(InputsDir, os.path.split(config_file)[1])

    # ---------------------------------------------------------------- #
    # copy the inputs
    for key, section in config_dict.iteritems():
        if 'FILE_NAME' in section.keys():
            new_file_name = os.path.join(InputsDir,
                                         os.path.split(section['FILE_NAME'])[1])

            copyfile(section['FILE_NAME'], new_file_name)

            # update the config file for an easy restart
            config.set(key, 'FILE_NAME',
                       os.path.join(InputsDir,
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


def tar_inputs(inputs, suffix=''):
    """ Tar the inputss directory or file at the end of a run"""
    # ---------------------------------------------------------------- #
    # Make the TarFile
    tar_file = inputs + suffix + '.tar.gz'
    print 'tarfile: %s' %tar_file

    if os.path.isdir(inputs):
        arcname = os.path.basename(os.path.normpath(inputs))
    else:
        arcname = os.path.split(inputs)[1]

    with tarfile.open(tar_file, "w:gz") as tar:
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
        log.error('Problem removing inputs: %s' % inputs)
    # ---------------------------------------------------------------- #
    return tar_file
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the domain
def read_domain(domain_dict):
    """
    Read the domain file and return all the variables and attributes.
    Area is returned in m2
    """
    dom_data, dom_vatts, dom_gatts = read_netcdf(domain_dict['FILE_NAME'])

    check_ncvars(domain_dict, dom_data.keys())

    # ---------------------------------------------------------------- #
    # Create the cell_ids variable
    dom_data['cell_ids'] = np.arange(dom_data[domain_dict['LAND_MASK_VAR']].size).reshape(dom_data[domain_dict['LAND_MASK_VAR']].shape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make sure the longitude / latitude vars are 2d
    dom_data['cord_lons'] = dom_data[domain_dict['LONGITUDE_VAR']]
    dom_data['cord_lats'] = dom_data[domain_dict['LATITUDE_VAR']]
    if dom_data[domain_dict['LONGITUDE_VAR']].ndim == 1:
        dom_data[domain_dict['LONGITUDE_VAR']], \
        dom_data[domain_dict['LATITUDE_VAR']] = np.meshgrid(dom_data[domain_dict['LONGITUDE_VAR']],
                                                            dom_data[domain_dict['LATITUDE_VAR']])

    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make sure the area is in m2
    area_units = dom_vatts[domain_dict['AREA_VAR']]['units']

    if area_units in ["rad2", "radians2", "radian2", "radian^2", "rad^2",
                      "radians^2", "rads^2", "radians squared", "square-radians"]:
        dom_data[domain_dict['AREA_VAR']] = dom_data[domain_dict['AREA_VAR']]*EARTHRADIUS*EARTHRADIUS
    elif area_units in ["m2", "m^2", "meters^2", "meters2", "square-meters",
                        "meters squared"]:
        dom_data[domain_dict['AREA_VAR']] = dom_data[domain_dict['AREA_VAR']]
    elif area_units in ["km2", "km^2", "kilometers^2", "kilometers2",
                        "square-kilometers", "kilometers squared"]:
        dom_data[domain_dict['AREA_VAR']] = dom_data[domain_dict['AREA_VAR']]*METERSPERKM*METERSPERKM
    elif area_units in ["mi2", "mi^2", "miles^2", "miles", "square-miles",
                        "miles squared"]:
        dom_data[domain_dict['AREA_VAR']] = dom_data[domain_dict['AREA_VAR']]*METERSPERMILE*METERSPERMILE
    elif area_units in ["acres", "ac", "ac."]:
        dom_data[domain_dict['AREA_VAR']] = dom_data[domain_dict['AREA_VAR']]*METERS2PERACRE
    else:
        log.warning("WARNING: UNKNOWN AREA units (%s), ASSUMING THEY ARE IN "
                    "SQUARE METERS" % dom_data[domain_dict['AREA_VAR']]['units'])
    # ---------------------------------------------------------------- #
    return dom_data, dom_vatts, dom_gatts
# -------------------------------------------------------------------- #
