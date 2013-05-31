import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset

# open the UH netcdf file
basin_nc ='/raid/jhamman/temp_uh_files/run6/MACRR_UH.nc'
f1 = Dataset(basin_nc)
uh1 = f1.variables['unit_hydrograph'][:]
t_uh1 = f1.variables['time'][:]
lon1 = f1.variables['lon'][:]
lat1 = f1.variables['lat'][:]

x1,y1 = np.meshgrid(lon1,lat1)

# open the UH netcdf file
basin_nc ='/raid/jhamman/temp_uh_files/run6/MACRR_RASM_UH.nc'
f2 = Dataset(basin_nc)
uh2 = f2.variables['unit_hydrograph'][:]
t_uh2 = f2.variables['time'][:]
x2 = f2.variables['xc'][:]
y2 = f2.variables['yc'][:]

ii = [25]
for i in ii:
    fig=plt.figure()
    plt.subplots_adjust(bottom=0.15,right=0.8,top=0.9)

    # Define Basemap Properties
    m=Basemap(projection='stere', lat_0 = 70, lon_0=-120,lat_ts=80.5,\
              llcrnrlon = -135,llcrnrlat = 50,\
              urcrnrlon = -90,urcrnrlat= 70,\
              resolution='l')
    m.drawcoastlines()
    m.drawlsmask(land_color = 'white',ocean_color='0.8')
    # draw parallels
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels = [0,0,0,0],fontsize=8)
    # draw meridians
    meridians = np.arange(0,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    #Subplot 1 - DRT Grid
    plt.subplot(121)
    plt.suptitle('Impulse Response Functions, Mackenzie River at Arctic Red River',fontsize=24)
    xi1,yi1 = m(x1,y1)

    # Define Basemap Properties
    m=Basemap(projection='stere', lat_0 = 70, lon_0=-120,lat_ts=80.5,\
              llcrnrlon = -135,llcrnrlat = 50,\
              urcrnrlon = -90,urcrnrlat= 70,\
              resolution='l')
    m.drawcoastlines()
    m.drawlsmask(land_color = 'white',ocean_color='0.8')
    # draw parallels
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels = [0,0,0,0],fontsize=8)
    # draw meridians
    meridians = np.arange(0,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=8)

    #Plot
    m.pcolor(xi1,yi1,np.squeeze(uh1[i,:,:]),vmax=1,vmin=0)
    plt.xlabel("DRT GRID")

    #Subplot 2 - RASM Grid
    plt.subplot(122)
    #plt.title(('Travel Time = '+str(i)+' Days'),fontsize=12)

    xi2,yi2 = m(x2,y2)

    # Define Basemap Properties
    m=Basemap(projection='stere', lat_0 = 70, lon_0=-120,lat_ts=80.5,\
              llcrnrlon = -135,llcrnrlat = 50,\
              urcrnrlon = -90,urcrnrlat= 70,\
              resolution='l')
    m.drawcoastlines()
    m.drawlsmask(land_color = 'white',ocean_color='0.8')
    # draw parallels
    parallels = np.arange(0.,90,10.)
    m.drawparallels(parallels,labels = [0,0,0,0],fontsize=8)
    # draw meridians
    meridians = np.arange(0,360.,10.)
    m.drawmeridians(meridians,labels=[0,0,0,0],fontsize=8)

    #Plot
    m.pcolor(xi2,yi2,np.squeeze(uh2[i,:,:]),vmax=1,vmin=0)
    plt.xlabel("RASM GRID")
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    cbar = plt.colorbar(cax=cax,orientation='horizontal')
    plt.tight_layout()
    #plt.show()

    plt.savefig('/usr1/jhamman/RASM/rout_rasm/Figs/MACRR'+str(i).zfill(3)+'.png',dpi=
                300,format = 'png',transparent='True')
    print 'Done with step '+str(i)
    plt.close()

print 'Done with making figs'    

