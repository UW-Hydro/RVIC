# -*- coding: utf-8 -*-
'''
remap.py
'''

import os
from cdo import Cdo
cdo = Cdo()
from logging import getLogger
from .log import LOG_NAME

# -------------------------------------------------------------------- #
# create logger
log = getLogger(LOG_NAME)
# -------------------------------------------------------------------- #


# -------------------------------------------------------------------- #
# Remap a file using CDO
def remap(grid_file, in_file, out_file, operator='remapcon',
          remap_options=None):
    '''Remap infile using cdo'''

    log.info('Remapping %s to %s', in_file, out_file)

    remap_method = getattr(cdo, operator)

    if remap_options:
        remap_method(grid_file, input=in_file, output=out_file,
                     options=remap_options)
    else:
        remap_method(grid_file, input=in_file, output=out_file)

    log.debug('remapped to out file: %s', os.path.split(out_file)[1])

    return
# -------------------------------------------------------------------- #
