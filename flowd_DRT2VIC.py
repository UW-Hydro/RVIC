#!/usr/bin/python

import argparse
from numpy import *
from numpy import copy

# This script opens a DRT flow direction file and converts it to VIC flow direction format
# The conversion is as follows
# 	Direction	DRT	VIC
#	 east		 1	 3
#	 s.east		 2	 4
#	 south	 	 4	 5
# 	 s.west		 8	 6
# 	 west		 16	 7
#	 n.west		 32 	 8
#	 north		 64	 1
#	 n.east		 128	 2

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("inpath", help="path where DRT grid is found")
parser.add_argument("outpath", help="path where VIC grid should go")
parser.add_argument("nodata", help="nodata value",type=int)
parser.add_argument("res", help="input grid resolution, i.e. 16th, this can take multiple entries",nargs='+' )
args = parser.parse_args()

#inputs = ['10th', '12th', '16th', '1', '2', '8th', 'half', 'qd']
input= args.res

for i in input:
  invar = 'DRT_' + i + '_FDR_globe.asc'
  outvar = 'DRT_' + i + '_FDR_globe_vic.asc'

  #inpath = '/raid/jhamman/DRT/upscaled_global_hydrography/by_HydroSHEDS_Hydro1k/flow_direction/'
  inpath = arg.inpath
  #outpath = '/raid/jhamman/DRT/upscaled_global_hydrography/by_HydroSHEDS_Hydro1k/flow_direction/VIC_flow_direction/'
  outpath =arg.outpath

  d = {1: 3, 2: 4, 4: 5, 8: 6, 16: 7, 32: 8, 64: 1, 128: 2, -9999: arg.nodata}

  # open the input file (DRT flow direction)
  f1 = inpath+invar
  drtdata = loadtxt(f1, skiprows = 6, dtype = int)

  # replace DRT values with VIC values
  vicdata = copy(drtdata)
  for k, v in d.iteritems(): vicdata[drtdata==k] = v

  # save converted data to tmp file
  with open(outpath+'temp/tmp', 'w') as ftemp:
    savetxt(outpath+'temp/tmp', vicdata, fmt='%1.1i', delimiter=' ', newline='\n')

  # get original header (minus no data line (6))
  with open(f1) as myfile:
      head=[myfile.next() for x in xrange(5)]

  # open output file
  f2 = open(outpath+outvar, 'a+')
  tmp = open(outpath+'temp/tmp', 'r')

  #write header data
  for hh in head:
    f2.write(hh)
  #write nodata line
  f2.write('NODATA_value  arg.nodata\n')
  #write coverted data
  for line in tmp:
    f2.write(line)

  f2.close()
  tmp.close()
