#!/usr/bin/python

import numpy as np
import time
from osgeo import gdal


in_path = "/raid2/jhamman/route_rasm/inputs/Arctic_basins/"
# Load global inputs
f1 = gdal.Open(in_path+"flowdir.asc")
f2 = gdal.Open(in_path+"flowdist.asc")
f3 = gdal.Open(in_path+"basins.asc")
FDR = np.array(f1.GetRasterBand(1).ReadAsArray())
DIST = np.array(f2.GetRasterBand(1).ReadAsArray())
BASINS = np.array(f3.GetRasterBand(1).ReadAsArray())
FDR = np.ma.masked_values(FDR,-9999)
DIST = np.ma.masked_values(DIST,-9999)
BASINS = np.ma.masked_values(BASINS,-9999)

basin_info = np.loadtxt(in_path+"arctic_pour_points.txt", dtype={'names': ('OID', 'longitude', 'latitude','basin_area','min_lon','min_lat','max_lon','max_lat'),
                                                                 ...'formats': ('i8', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4')})
#Extract raster information (dimensions, min,max, etc)
width = f1.RasterXSize
height = f1.RasterYSize
gt = f1.GetGeoTransform()
res = gt[1]
minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5]
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3]

#Make arrays of lats and lons
lons = np.linspace(minx,maxx-res,num=width)
lats = np.linspace(maxy-res,miny,num=height)
#Make arrays of same dimensions as input arrays of lat/lon values
x,y = np.meshgrid(lons,lats)


######################## Trim each basin ##############

basin_id = 57087
# Make unique 5 digit alpha-numeric id for each basin 
basin = "Basin_"+str(basin_id).zfill(6)







Station File
1       YENIG   4264    1981    -9999
NONE
1       OBSAL   3946    1962    -9999
NONE
1       PECUT   3717    1943    -9999
NONE
1       DVIUP   3550    1923    -9999
NONE


Input File

FLOW_DIREC_FILE  /raid2/jhamman/route_rasm/inputs/run5/flow_dir.asc
VELOCITY_FILE    1
DIFF_FILE        3000
XMASK_FILE       /raid2/jhamman/route_rasm/inputs/run5/flow_dist.asc
FRACTION_FILE    /raid2/jhamman/route_rasm/inputs/run5/land_mask.asc
STATION_FILE     /raid2/jhamman/route_rasm/inputs/run5/stnloc.arctic.a
INPUT_FILE_PREC  4
UNIT_HD_FILE     /raid2/jhamman/route_rasm/inputs/run2/UH1

#Write Routing Shell Script

#qsub routing shell script


