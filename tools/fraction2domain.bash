#!/bin/bash -x 

set -e

# This script uses gdal, nco and cdo to make a CESM type domain file from a ascii fraction file
# gdal --> http://www.gdal.org/
# nco ---> http://nco.sourceforge.net/
# cdo ---> https://code.zmaw.de/projects/cdo

# -------------------------------------------------------------------- #
# Get command line inputs
fraction_file=$1
grid_name=$2
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Local variables
temp_frac1=temp_frac1.nc
temp_frac=temp_frac.nc
temp_area=temp_area.nc
temp_mask1=temp_mask1.nc
temp_mask=temp_mask.nc

outfile="domain.rvic."$grid_name"."`eval date +%Y%m%d`".nc"
script=$0
PLANET_RADIUS=6371000
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Make fraction variable
echo 'making fraction variable'
gdal_translate -of netCDF $fraction_file $temp_frac1
ncap2 -s 'Band1=double(Band1)' $temp_frac1 $temp_frac
ncrename -v Band1,frac $temp_frac
ncatted -O -a long_name,frac,o,c,"fraction of grid cell that is active" $temp_frac
ncatted -O -a note,frac,a,c,"unitless" $temp_frac

ncks -v frac $temp_frac $outfile
echo 'done with fraction variable'
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Make area file (m2)
echo 'making area variable'
cdo gridarea $temp_frac $temp_area

ncrename -v cell_area,area $temp_area
ncks -A -v area $temp_area $outfile
echo 'done making area variable'
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# Make land mask file (1 = land, 0 = not land)
echo 'making mask variable'
ncap2 -O -v -s 'mask=frac' temp_frac.nc $temp_mask1
ncap2 -O -v -s 'where(mask>0.0000001) mask=1; elsewhere mask=0;' $temp_mask1 $temp_mask
ncatted -O -a long_name,mask,a,c,"domain mask" $temp_mask
ncatted -O -a note,mask,a,c,"unitless" $temp_mask
ncatted -O -a comment,mask,a,c,"0 value indicates cell is not active" $temp_mask

ncks -A -v mask $temp_mask $outfile
echo 'done making mask variable'
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# add global attributes
ncatted -O -a title,global,a,c,"RVIC domain data" $outfile
ncatted -O -a history,global,o,c,"created by $(whoami) "`eval date +%Y-%m-%d-%T`" from $fraction_file" $outfile
ncatted -O -a source,global,o,c,"$0" $outfile
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
# remove temporary files
rm $temp_frac1 $temp_frac $temp_mask1 $temp_mask $temp_area
# -------------------------------------------------------------------- #

