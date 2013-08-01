"""
utilities.py
"""
import os
import shutil
from ConfigParser import SafeConfigParser
import numpy as np
from netCDF4 import Dataset
from cdo import *
cdo = Cdo()
from rvic.log import log_name
import logging
from share import timeStampForm, rpointer, earthRadius
import tarfile

# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(log_name)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Write rpointer file
def write_rpointer(restart_dir, restart_file, timestamp):
    """ Write a configuration file with restart file and time """
    rpointer_file = os.path.join(restart_dir, rpointer)

    config = SafeConfigParser()

    time_str = timestamp.strftime(timeStampForm)

    config.add_section('restart')
    config.set('restart', 'file_name', os.path.join(restart_dir, restart_file))
    config.set('restart', 'timestamp', time_str)

    with open(rpointer_file, 'a') as configfile:
        config.write(configfile)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# A helper function to read a netcdf file
def ReadNetcdf(ncFile, variables=None, coords=None):
    """
    Read data from input netCDF. Will read all variables if none provided.
    Will also return all variable attributes.
    Both variables (data and attributes) are returned as dictionaries named by variable
    """

    f = Dataset(ncFile, 'r+')

    if not variables:
        variables = f.variables.keys()
    if not coords:
        coords = slice(None)

    log.info('Reading input data variables: %s, from file: %s' % (variables, ncFile))

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
def CheckNcvars(configSection, nckeys):
    """Make sure the variables listed in the config file are present in the netcdf"""
    for key, value in configSection.iteritems():
        if key.endswith('var'):
            if value not in nckeys:
                log.error('%s (%s) not in %s' % (value, key, configSection['file_name']))
                raise NameError('Check netcdf that netcdf variable names match those in the configuration file')
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the Configuration File
def ReadConfig(ConfigFile):
    """
    Return a dictionary with subdictionaries of all configFile options/values
    """
    Config = SafeConfigParser()
    Config.optionxform = str
    Config.read(ConfigFile)
    sections = Config.sections()
    dict1 = {}
    for section in sections:
        options = Config.options(section)
        dict2 = {}
        for option in options:
            dict2[option] = configType(Config.get(section, option))
        dict1[section] = dict2
    return dict1
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the type of the config options
def configType(value):
    """
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    """
    val_list = [x.strip() for x in value.split(',')]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE']:
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
    idx = (np.abs(array-value)).argmin()
    return idx
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Delete all the files in a directory
def clean_dir(dir):
    """ Clean all files in a directory"""
    for file in os.listdir(dir):
        file_path = os.path.join(dir, file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except:
            log.exception('Error cleaning file: %s' % file_path)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Delete a particular file
def clean_file(file_path):
    """ Delete the file"""
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
    except:
        log.exception('Error cleaning file: %s' % file_path)
    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Remap a file using CDO
def remap(gridFile, inFile, remapFile, operator='remapcon',
          remap_options=None):
    """Remap infile using cdo"""

    remap_method = getattr(cdo, operator)

    if remap_options:
        remap_method(gridFile, input=inFile, output=remapFile, options=remap_options)
    else:
        remap_method(gridFile, input=inFile, output=remapFile)

    log.info('remapped to out file: %s' % os.path.split(remapFile)[1])

    return
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Make a set of directories
def MakeDirs(rundir, subdir_names):
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
def CopyInputs(ConfigFile, InputsDir):

    ConfigDict = ReadConfig(ConfigFile)

    Config = SafeConfigParser()
    Config.read(ConfigFile)

    new_config = os.path.join(InputsDir, os.path.split(ConfigFile)[1])

    # ---------------------------------------------------------------- #
    # copy the inputs
    for key, section in ConfigDict.iteritems():
        if 'file_name' in section.keys():
            new_file_name = os.path.join(InputsDir,
                                         os.path.split(section['file_name'])[1])

            shutil.copyfile(section['file_name'], new_file_name)

            # update the config file for an easy restart
            Config.set(key, 'file_name',
                       os.path.join(InputsDir,
                                    os.path.split(section['file_name'])[1]))

            # update the ConfigDict with the new value
            ConfigDict[key]['file_name'] = new_file_name
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # write the new configuration file
    with open(new_config, 'w') as configfile:
        Config.write(configfile)
    # ---------------------------------------------------------------- #

    return ConfigDict
# -------------------------------------------------------------------- #


def TarInputs(Input, suffix=''):
    """ Tar the inputs directory or file at the end of a run"""

    # ---------------------------------------------------------------- #
    # Make the TarFile
    tar_file = Input + suffix + '.tar.gz'
    with tarfile.open(tar_file, "w:gz") as tar:
        tar.add(Input)

    # ---------------------------------------------------------------- #
    # Check to make sure the TarFile exists before deleting the sources
    if os.path.exists(tar_file):
        # ------------------------------------------------------------ #
        # Remove the Input
        if os.path.isdir(Input):
            shutil.rmtree(Input)
        elif os.path.isfile(Input):
            os.unlink(Input)
        # ------------------------------------------------------------ #
    else:
        log.error('Problem removing Input: %s' % Input)
    # ---------------------------------------------------------------- #
    return tar_file
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Shorten the unit hydrograph
def subset(uh, subset, threshold):
    """ Shorten the Unit Hydrograph"""
    full_length = uh.shape[0]
    offset = np.empty(uh.shape[1], dtype=int)
    out_uh = np.zeros((subset, uh.shape[1]))

    for i in xrange(uh.shape[1]):
        offset[i] = np.nonzero(uh[:, i] > threshold)[0][0]       # find position of first index > threshold
        end = np.minimum(uh.shape[0], offset[i]+subset)          # find end point
        out_uh[:end-offset[i], i] = uh[offset[i]:end, i] / uh[offset[i]:end, i].sum()  # clip and normalize

    return offset, out_uh, full_length
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Read the domain
def ReadDomain(DomainDict):
    """
    Read the domain file and return all the variables and attributes.
    Area is returned in m2
    """
    DomData, DomVats, DomGats = ReadNetcdf(DomainDict['file_name'])

    CheckNcvars(DomainDict, DomData.keys())

    # ---------------------------------------------------------------- #
    # Create the cell_ids variable
    DomData['cell_ids'] = np.arange(DomData[DomainDict['land_mask_var']].size).reshape(DomData[DomainDict['land_mask_var']].shape)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Make sure the area is in m2
    if DomVats[DomainDict['area_var']]['units'] in ["rad2", "radians2", "radian2", "rad^2",
                                                    "radians^2", "rads^2", "radians squared",
                                                    "square-radians"]:
        DomData[DomainDict['area_var']] = DomData[DomainDict['area_var']] * earthRadius * earthRadius
    elif DomVats[DomainDict['area_var']]['units'] in ["m2", "m^2", "meters^2", "meters2",
                                                      "square-meters", "meters squared"]:
        DomData[DomainDict['area_var']] = DomData[DomainDict['area_var']]
    elif DomVats[DomainDict['area_var']]['units'] in ["km2", "km^2", "kilometers^2",
                                                      "kilometers2", "square-kilometers",
                                                      "kilometers squared"]:
        DomData[DomainDict['area_var']] = DomData[DomainDict['area_var']] * metersPerKm * metersPerKm
    elif DomVats[DomainDict['area_var']]['units'] in ["mi2", "mi^2", "miles^2", "miles",
                                                      "square-miles", "miles squared"]:
        DomData[DomainDict['area_var']] = DomData[DomainDict['area_var']] * metersPerMile * metersPerMile
    elif DomVats[DomainDict['area_var']]['units'] in ["acres", "ac", "ac."]:
        DomData[DomainDict['area_var']] = DomData[DomainDict['area_var']] * meters2PerAcre
    else:
        log.warning("WARNING: UNKNOWN AREA UNITS (%s), ASSUMING THEY ARE IN "
                    "SQUARE METERS" % DomData[DomainDict['area_var']]['units'])
    # ---------------------------------------------------------------- #
    return DomData, DomVats, DomGats
# -------------------------------------------------------------------- #
