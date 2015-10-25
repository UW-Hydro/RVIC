# -*- coding: utf-8 -*-
'''
plots.py
'''
import os
import logging
from .log import LOG_NAME
import numpy as np
from datetime import date
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    matplotlib_available = True
    try:
        from mpl_toolkits.basemap import Basemap
        basemap_available = False
    except ImportError:
        basemap_available = False
except ImportError:
    matplotlib_available = False


# -------------------------------------------------------------------- #
# create logger
log = logging.getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def uhs(data, title, case_id, plot_dir):
    '''
    Plot diagnostic plot showing all unit hydrographs
    '''
    pfname = _make_filename(title, case_id, plot_dir)

    fig = plt.figure()
    plt.plot(data)
    plt.title(title)
    plt.xlabel('timesteps')
    plt.ylabel('unit-hydrograph')
    fig.savefig(pfname)

    return pfname
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def _fractions_grid(data, dom_x, dom_y, title, case_id, plot_dir):
    '''
    Plot diagnostic plots of fraction variables
    '''
    # ---------------------------------------------------------------- #
    # Plot Fractions
    pfname = _make_filename(title, case_id, plot_dir)

    mask = data <= 0.0
    data = np.ma.array(data, mask=mask)

    cmap = matplotlib.cm.cool
    cmap.set_bad(color='w')

    fig = plt.figure()
    plt.pcolormesh(data, cmap=cmap)
    plt.autoscale(tight=True)
    plt.axis('tight')
    plt.colorbar()
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.ylim([0, dom_y.shape[0]])
    plt.xlim([0, dom_x.shape[1]])
    fig.savefig(pfname)
    # ---------------------------------------------------------------- #
    return pfname
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def _fractions_map(data, dom_x, dom_y, title, case_id, plot_dir):
    '''
    Plot diagnostic plots of fraction variables using Basemap
    '''
    # ---------------------------------------------------------------- #
    # Plot Fractions
    pfname = _make_filename(title, case_id, plot_dir)

    fig = plt.figure(figsize=(8, 8))
    fig.add_axes([0.1, 0.1, 0.8, 0.8])

    dom_x[dom_x < 0] += 360.0

    mask = data <= 0.0
    data = np.ma.array(data, mask=mask)

    cmap = matplotlib.cm.cool
    cmap.set_bad(color='w')

    # define projection
    midx = int(dom_x.shape[1] / 2)
    midy = int(dom_x.shape[0] / 2)

    projection = {'projection': 'stere',
                  'lon_0': dom_x[-1, midx],
                  'lat_0': dom_y[midy, midx],
                  'llcrnrlat': dom_y[-1, 0],
                  'urcrnrlat': dom_y[0, -1],
                  'llcrnrlon': dom_x[-1, 0],
                  'urcrnrlon': dom_x[0, -1],
                  'resolution': 'l'}

    log.debug('Projection: %s', projection)

    m = Basemap(**projection)

    m.drawcoastlines()
    m.drawcountries()

    # draw parallels.
    parallels = np.arange(-90., 90, 10.)
    m.drawparallels(parallels, labels=[1, 0, 0, 0])

    # draw meridians
    meridians = np.arange(-180., 180., 10.)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1])

    x, y = m(dom_x, dom_y)  # compute map proj coordinates.
    cs = m.pcolormesh(x, y, data, cmap=cmap)
    m.colorbar(cs, location='right', pad='5%')
    plt.title(title)
    fig.savefig(pfname)
    # ---------------------------------------------------------------- #
    return pfname
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def _make_filename(title, case_id, plot_dir):
    today = date.today().strftime('%Y%m%d')
    file_name = '{0}_{1}_{2}.png'.format(title.lower().replace(' ', '_'),
                                         case_id.lower().replace(' ', '_'),
                                         today)
    pfname = os.path.join(plot_dir, file_name)
    return pfname
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
def _fractions_dummy(*args):
    '''
    Pass on plotting
    '''
    return 'None <-- could not import matplotlib'
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Set function aliases
if matplotlib_available and basemap_available:
    fractions = _fractions_map
elif matplotlib_available and not basemap_available:
    fractions = _fractions_grid
elif not matplotlib_available or not basemap_available:
    fractions = _fractions_dummy
# -------------------------------------------------------------------- #
