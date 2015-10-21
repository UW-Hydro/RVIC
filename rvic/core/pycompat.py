# -*- coding: utf-8 -*-
'''module to support compatability between Python 2 and Python 3 and to handle
default and optional packages'''

import sys

py3 = sys.version_info[0] >= 3

if py3:
    basestring = str
    unicode_type = str
    bytes_type = bytes

    def iteritems(d):
        return iter(d.items())

    def itervalues(d):
        return iter(d.values())
    pyrange = range
    pyzip = zip
    from configparser import SafeConfigParser
else:
    # Python 2
    basestring = basestring
    unicode_type = unicode
    bytes_type = str

    def iteritems(d):
        return d.iteritems()

    def itervalues(d):
        return d.itervalues()
    pyrange = xrange
    from itertools import izip as pyzip
    from ConfigParser import SafeConfigParser
try:
    from cyordereddict import OrderedDict
except ImportError:
    try:
        from collections import OrderedDict
    except ImportError:
        from ordereddict import OrderedDict
