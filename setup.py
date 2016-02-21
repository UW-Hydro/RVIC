#!/usr/bin/env python
import os
import re
import sys
import warnings
try:
    from setuptools import setup
    from setuptools.extension import Extension
except:
    from distutils.core import setup
    from distutils.extension import Extension

MAJOR = 1
MINOR = 1
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
QUALIFIER = 'dev'

FULLVERSION = VERSION
write_version = True


# -------------------------------------------------------------------- #
def write_version_py(filename=None):
    cnt = '''\
version = '%s'
short_version = '%s'
'''
    if not filename:
        filename = os.path.join(
            os.path.dirname(__file__), 'rvic', 'version.py')

    a = open(filename, 'w')
    try:
        a.write(cnt % (FULLVERSION, VERSION))
    finally:
        a.close()

    return
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Get version string
if ISRELEASED:
    FULLVERSION += QUALIFIER
else:
    import subprocess
    FULLVERSION += '.dev'

    pipe = None
    for cmd in ['git', 'git.cmd']:
        try:
            pipe = subprocess.Popen(
                [cmd, 'describe', '--always', '--match', 'v[0-9]*'],
                stdout=subprocess.PIPE)
            (so, serr) = pipe.communicate()
            if pipe.returncode == 0:
                break
        except:
            pass

    if pipe is None or pipe.returncode != 0:
        # no git, or not in git dir
        if os.path.exists('rvic/version.py'):
            warnings.warn("WARNING: Couldn't get git revision, using existing "
                          "rvic/version.py")
            write_version = False
        else:
            warnings.warn("WARNING: Couldn't get git revision, using generic "
                          "version string")
    else:
        # have git, in git dir, but may have used a shallow clone
        # (travis does this)
        rev = so.strip()
        # makes distutils blow up on Python 2.7
        if sys.version_info[0] >= 3:
            rev = rev.decode('ascii')

        if not rev.startswith('v') and re.match('[a-zA-Z0-9]{7,9}', rev):
            # partial clone, manually construct version string
            # this is the format before we started using git-describe
            # to get an ordering on dev version strings.
            rev = 'v%s.dev-%s' % (VERSION, rev)

        # Strip leading v from tags format "vx.y.z" to get th version string
        FULLVERSION = rev.lstrip('v')
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
if write_version:
    write_version_py()
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Run Setup
setup(name='rvic',
      version=FULLVERSION,
      description='The RVIC Streamflow Routing Model',
      long_description='The RVIC streamflow routing model is based on the '
                       'original model of Lohmann, et al., 1996, Tellus, '
                       '48(A), 708-721',
      author='Joe Hamman',
      author_email='jhamman1@uw.edu',
      classifiers=['Development Status :: 4 - Beta',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Operating System :: OS Independent',
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5',
                   'Topic :: Scientific/Engineering'],
      install_requires=['scipy >= 0.13', 'numpy >= 1.8',
                        'netCDF4 >= 1.0.6', 'matplotlib >= 1.3.1',
                        'pandas >= 0.15.1'],
      tests_require=['pytest >= 2.5.2'],
      url='https://github.com/UW-Hydro/RVIC',
      packages=['rvic', 'rvic.core'],
      py_modules=['rvic.parameters', 'rvic.convolution', 'rvic.convert'],
      scripts=['scripts/rvic', 'tools/find_pour_points.py',
               'tools/fraction2domain.bash'],
      ext_modules=[Extension('rvic_convolution',
                             sources=['rvic/clib/rvic_convolution.c'])])
# -------------------------------------------------------------------- #
