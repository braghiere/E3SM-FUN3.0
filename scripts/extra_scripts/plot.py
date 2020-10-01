import matplotlib
from pylab import *
from netCDF4 import Dataset
import numpy as np


import cartopy.crs as ccrs
#squared
proj=ccrs.PlateCarree 
#proj=ccrs.Robinson
import xarray
from cartopy.util import add_cyclic_point
import os.path
import sys

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

datadir='.'

def letter_label(ax=None,xpos=0.03,ypos=0.92):
    if ax is None:
        ax=gca()
    from string import ascii_lowercase
    plotnum=ax.get_subplotspec().num1
    return text(xpos,ypos,'('+ascii_lowercase[plotnum]+')',transform=ax.transAxes)

### Stephanie Kivlin's figure of observation-based and model-based species distributions
#ECM_obs=xarray.open_rasterio(os.path.join(datadir,'ECMAll.grd'))
#ECM_mod=xarray.open_rasterio(os.path.join(datadir,'EMFBENModel.grd'))
#AM_obs=xarray.open_rasterio(os.path.join(datadir,'AMFAll.grd'))
#AM_mod=xarray.open_rasterio(os.path.join(datadir,'AMFBENModel.grd'))

f = Dataset('time_average_h0.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
gpp = f.variables['GPP'][:]
pactive = f.variables['PACTIVE'][:]

print(np.max(pactive))

upper = cm.jet(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


figure(6);clf()
levs=arange(0,1.0570833e-07,1e-8)
levs=[-0.125e-8,0.,0.125e-8,0.25e-8,0.625e-8,1e-8,2e-8,3e-8,4e-8,5e-8,6e-8,7e-8,7.5e-8,8e-8,9e-8,1.0570833e-07]
subplot(111,projection=proj())
pactive,lon=add_cyclic_point(pactive,lon)
#contourf(lon,lat,pactive,cmap=get_cmap('jet'),levels=levs,extend='both')
contourf(lon,lat,pactive,cmap=cmap,levels=levs)
cmap.set_under('w')
gca().coastlines('50m')
gca().xlocator = mticker.FixedLocator([-180,-45,0,45,180])
xformatter = LONGITUDE_FORMATTER
title('Observation-based species dist. model')
letter_label(ypos=1.05)
cb=colorbar()
cb.set_label('More AM $\longleftarrow \longrightarrow$ More ECM')

tight_layout()
savefig('test_sulman.png',dpi=300)
sys.exit()

figure(4);clf()
subplot(223,projection=proj())
Eobs,lon=add_cyclic_point(ECM_obs.isel(band=0).values,ECM_obs.x)
contourf(lon,ECM_obs.y,Eobs,vmin=0.0,vmax=1.0,cmap=get_cmap('YlGn'))
gca().coastlines()
title('Observation-based ECM')
colorbar()



subplot(221,projection=proj())
Aobs=add_cyclic_point(AM_obs.isel(band=0).values)
contourf(lon,AM_obs.y,Aobs,vmin=0.0,vmax=1.0,cmap=get_cmap('YlGn'))
gca().coastlines()
title('Observation-based AM')
colorbar()

subplot(224,projection=proj())
Emod=add_cyclic_point(ECM_mod.isel(band=0).values)
contourf(lon,ECM_obs.y,Emod,vmin=0.0,vmax=1.0,cmap=get_cmap('YlGn'))
gca().coastlines()
title('Model-based ECM')
colorbar()

subplot(222,projection=proj())
Amod=add_cyclic_point(AM_mod.isel(band=0).values)
contourf(lon,AM_obs.y,Amod,vmin=0.0,vmax=1.0,cmap=get_cmap('YlGn'))
gca().coastlines()
title('Model-based AM')
colorbar()

tight_layout()

figure(6);clf()
levs=arange(-.6,.61,0.1)
subplot(211,projection=proj())
contourf(lon,ECM_obs.y,Eobs-Aobs,cmap=get_cmap('RdBu'),levels=levs,extend='both')
gca().coastlines()
title('Observation-based species dist. model')
letter_label(ypos=1.05)
cb=colorbar()
cb.set_label('More AM $\longleftarrow \longrightarrow$ More ECM')

subplot(212,projection=proj())
contourf(lon,AM_obs.y,Emod-Amod,cmap=get_cmap('RdBu'),levels=levs,extend='both')
gca().coastlines()
title('Model-based species dist. model')
letter_label(ypos=1.05)
cb=colorbar()
cb.set_label('More AM $\longleftarrow \longrightarrow$ More ECM')


tight_layout()


show()
sys.exit()

############################Make netcdf############################
file = 'ECM_obs_Sulman.nc'
ncfile = Dataset(file, 'w', format='NETCDF4_CLASSIC')

#lat = ECM_obs.y.values[::-1]
lat = ECM_obs.y.values
lon = lon

#create dimensions
latitude = ncfile.createDimension('latitude',np.shape(ECM_obs.y)[0])
longitude = ncfile.createDimension('longitude',np.shape(lon)[0])

#define variables
latitudes = ncfile.createVariable('latitude','d',('latitude',))
longitudes = ncfile.createVariable('longitude','d',('longitude',))
layer = ncfile.createVariable('layer','d',('latitude','longitude'))

latitudes[:] = lat
longitudes[:] = lon
#layer[:,:] = Eobs[::-1,:]
layer[:,:] = Eobs[:,:]

#the coordinates
longitudes.units = 'degrees_east'
longitudes.long_name = 'longitude'
latitudes.units = 'degrees_north'
latitudes.long_name = 'latitude'
layer.fill_value = -3.4e+38
layer.missing_value = -3.4e+38
layer.long_name = 'layer'
layer.min = 7.2365492087556e-06
layer.max = 0.973557472229004


#Global Attributes 
ncfile.proj4 = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
ncfile.Conventions = 'CF-1.4'

#close ncfile
ncfile.close()

sys.exit()

##########################################################
