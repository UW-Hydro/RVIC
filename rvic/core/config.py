# -*- coding: utf-8 -*-
'''
config.py

'''

import os
from .pycompat import OrderedDict, SafeConfigParser


class Config(object):
    def __init__(self, **kwargs):
        pass


class ConvertConfig(Config):
    pass


class ConvolutionConfig(Config):
    pass


class ParametersConfig(Config):
    pass


# -------------------------------------------------------------------- #
# Read the Configuration File
def read_config(config_file):
    '''
    Return a dictionary with subdictionaries of all configFile options/values
    '''
    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    sections = config.sections()
    dict1 = OrderedDict()
    for section in sections:
        options = config.options(section)
        dict2 = OrderedDict()
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        dict1[section] = dict2
    return dict1
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Find the type of the config options
def config_type(value):
    '''
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    '''
    val_list = [x.strip() for x in value.split(',')]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE', 'T']:
            return True
        elif value in ['false', 'False', 'FALSE', 'F']:
            return False
        elif value in ['none', 'None', 'NONE', '']:
            return None
        elif isint(value):
            return int(value)
        elif isfloat(value):
            return float(value)
        else:
            return os.path.expandvars(value)
    else:
        try:
            return list(map(float, val_list))
        except (TypeError, ValueError):
            pass
        try:
            return list(map(int, val_list))
        except (TypeError, ValueError):
            return val_list
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def isfloat(x):
    '''Test if value is a float'''
    try:
        float(x)
    except ValueError:
        return False
    else:
        return True
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def isint(x):
    '''Test if value is an integer'''
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b
# -------------------------------------------------------------------- #
