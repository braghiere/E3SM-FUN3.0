from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from netCDF4 import Dataset, date2index
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
import numpy.ma as ma
import numpy as np
from mpl_toolkits.basemap import maskoceans
from mpl_toolkits.basemap import interp
from matplotlib.colors import ListedColormap
import seaborn as sns
# High resolution cost lines
#http://introtopython.org/visualization_earthquakes.html
# High resolution cost lines
#http://basemaptutorial.readthedocs.io/en/latest/utilities.html
plt.style.use('ggplot')
SIZE = 48
plt.rc('font', size=SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=SIZE)  # fontsize of the x any y labels
plt.rc('xtick', labelsize=SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SIZE)  # legend fontsize
plt.rc('figure', titlesize=SIZE)  # # size of the figure title

#f = Dataset('clumping005.nc','r')
f = Dataset('/home/renato/Desktop/clm5.0_NOFUN.clm2.h0.time.0102-1012.nc','r')
f1 = Dataset('/home/renato/Desktop/clm5.0_FUN.clm2.h0.time.0102-1012.nc','r')
f2 = Dataset('/home/renato/Desktop/clm5.0_FUN_NEWPERECM.clm2.h0.time.0102-1012.nc','r')

#surf = Dataset('/home/renato/Steindinger/surfdata_1.9x2.5_hist_78pfts_CMIP6_simyr1850_c190304.nc','r')
surf = Dataset('/home/renato/datasets/surfacedata_clm/surfdata_1.9x2.5_hist_78pfts_CMIP6_simyr1850_c191015_steidinger_unmasked.nc','r')

#surf_new = Dataset('/home/renato/datasets/surfacedata_clm/surfdata_1.9x2.5_hist_78pfts_CMIP6_simyr1850_c191015_steidinger_unmasked.nc','r')
surf_new = Dataset('/home/renato/datasets/surfacedata_clm/surfdata_1.9x2.5_hist_78pfts_CMIP6_simyr1850_c190304_steidinger_rcp85.nc','r')


params = Dataset('/home/renato/Steindinger/clm5_params.c171117.nc','r')
params_new = Dataset('/home/renato/datasets/sulman_2019/sulman_modified_clm5_params.c171117.nc','r')


#f1 = Dataset('/home/renato/src/shi_pft_ecms.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
clump = f.variables['GPP'][:]
clump1 = f1.variables['GPP'][:]
clump2 = f2.variables['GPP'][:]


#perecm = params.variables['perecm'][:]


perecm = surf.variables['PERECM'][:]
perecm_new = surf_new.variables['PERECM'][:]




pft_map_nat = surf.variables['PCT_NAT_PFT'][:]
pft_map_crop = surf.variables['PCT_CFT'][:]

pct_natveg = surf.variables['PCT_NATVEG'][:]
pct_crop = surf.variables['PCT_CROP'][:]

#Declaring empty matrix
perecm_map_nat = np.zeros((96,144))
perecm_map_nat_new = np.zeros((96,144))

perecm_map_crop = np.zeros((96,144))
perecm_map_crop_new = np.zeros((96,144))

perecm_map_tot = np.zeros((96,144))
perecm_map_tot_new = np.zeros((96,144))

for i in xrange(14):
    print(i,perecm[i],perecm_new[i])
    perecm_map_nat = perecm_map_nat + (perecm[i]*pft_map_nat[i,:,:])*pct_natveg[:,:]/100.

    #perecm_map_nat_new = perecm_map_nat_new + (perecm_new[i]*pft_map_nat[i,:,:])*pct_natveg[:,:]/100.
    perecm_map_nat_new = perecm_map_nat_new + (perecm_new[i,:,:]*pft_map_nat[i,:,:])*pct_natveg[:,:]/100.

perecm_map_nat = ma.array(perecm_map_nat,mask=[pft_map_nat[0,:,:]>90.])

perecm_map_nat_new = ma.array(perecm_map_nat_new,mask=[pft_map_nat[0,:,:]>90.])



for i in xrange(64):
    print(i+14,perecm[i+14],perecm_new[i+14])
    perecm_map_crop = perecm_map_crop + (perecm[i+14]*pft_map_crop[i,:,:])*pct_crop[:,:]/100.
    
    #perecm_map_crop_new = perecm_map_crop_new + (perecm_new[i+14]*pft_map_crop[i,:,:])*pct_crop[:,:]/100.
    perecm_map_crop_new = perecm_map_crop_new + (perecm_new[i+14,:,:]*pft_map_crop[i,:,:])*pct_crop[:,:]/100.

perecm_map_crop = ma.array(perecm_map_crop,mask=[pft_map_nat[0,:,:]>90.])

perecm_map_crop_new = ma.array(perecm_map_crop_new,mask=[pft_map_nat[0,:,:]>90.])

#perecm_map_tot = perecm_map_nat + perecm_map_crop
#perecm_map_tot_new = perecm_map_nat_new + perecm_map_crop_new

#ONLY NATURAL VEGETATION
perecm_map_tot = perecm_map_nat 
perecm_map_tot_new = perecm_map_nat_new 

print lat,lon,clump
#lon = np.clip(lon, -180., 180.)
#lon, lat = np.meshgrid(np.linspace(-180, 180, 144), np.linspace(-90, 90, 96))

# shifting grid to run from -180 to 180 rather than 0-360
#clump,lon = shiftgrid(180., clump, lon, start=False)
perecm_map_nat,lon = shiftgrid(180., perecm_map_nat, lon, start=False)

lon = f.variables['lon'][:]
perecm_map_nat_new,lon = shiftgrid(180., perecm_map_nat_new, lon, start=False)

lon = f.variables['lon'][:]
perecm_map_tot,lon = shiftgrid(180., perecm_map_tot, lon, start=False)

lon = f.variables['lon'][:]
perecm_map_crop,lon = shiftgrid(180., perecm_map_crop, lon, start=False)

lon = f.variables['lon'][:]
perecm_map_crop_new,lon = shiftgrid(180., perecm_map_crop_new, lon, start=False)

lon = f.variables['lon'][:]
perecm_map_tot_new,lon = shiftgrid(180., perecm_map_tot_new, lon, start=False)


m = Basemap(projection='robin', lon_0=0.,resolution='l')


#x, y = m(lon, lat)
x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)

print(lon,lat)
print(len(lon),len(lat))



v = perecm_map_tot_new[:,:] - perecm_map_tot[:,:]
#v = perecm_map_tot_new[:,:]
#v = perecm_map_tot[:,:]

x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*20)
y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*30)

#x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*2)
#y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*3)

x2,y2 = np.meshgrid(x2,y2)
X2,Y2 = m(x2, y2)

data2 = interp(v,x[0],y[:,0],x2,y2,order=1)



mdata = maskoceans(x2, y2, data2,resolution='h',grid=1.25,inlands=True)
#mdata = maskoceans(x, y, v)


#lats = lat
#lons = lon
#lon_len = len(perecm_map_tot[0,:])
#lat_len = len(perecm_map_tot[:,0])
#root = Dataset('shi_ecm_pft_1p9x2p5.nc','w',format='NETCDF4')
#root.description = 'ECM percentage by PFT based on Shi et al. (2016)'
#root.createDimension('lon',lon_len)	
#root.createDimension('lat',lat_len)
#root.createDimension('z',18)
#latitudes = root.createVariable('lat', 'f4', ('lat',))	
#longitudes = root.createVariable('lon', 'f4', ('lon',))	
#var_shi = root.createVariable('perecm (fraction)', 'f4', ('z','lat','lon',),fill_value=-1e+20)
#var_shi = root.createVariable('perecm (fraction)', 'f4', ('lat','lon',),fill_value=-1e+20)
#latitudes[:] = lats[0:lat_len]	
#longitudes[:] = lons[0:lon_len]


#var_shi[:,:] = mdata[:,:]
#root.close()
#sys.exit()



for i in xrange(1):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.6', linewidth=0.5)
   m.drawcountries(color='0.6', linewidth=0.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=0.25, color='0.5',fontsize='x-large')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5',fontsize='x-large')

 
   #PLOT ABSOLUTE
   levels = np.linspace(0.,100.,21)
   vmin,vmax = (0,100)
   cs = m.contourf(X2,Y2,mdata,levels,vmin=vmin,vmax=vmax,cmap=plt.cm.jet,extend='both')
   
   #PLOT DIFFERENCE
   vmin,vmax = (-25,25)
   levels = np.linspace(vmin,vmax,16)   
   cs = m.contourf(X2,Y2,mdata,levels,vmin=vmin,vmax=vmax,cmap=plt.cm.bwr,extend='both')
   

   plt.tight_layout()
   cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=np.linspace(vmin,vmax,7),format='%.1f')
   #cbar.ax.get_yaxis().labelpad = 60
   cbar.ax.get_xaxis().labelpad = 60
   #cbar.ax.set_ylabel('EM (%)', rotation=270)
   cbar.ax.set_xlabel('Change in fraction of ECM (%)', rotation=0,color='black', size=78)
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   cbar.solids.set_edgecolor("face")
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'Sulman et al. (2019) - Shi et al. (2016)')
   cbar.ax.tick_params(labelsize='xx-large')
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('diff_em_steidinger_future_present_v2.png',bbox_inches="tight",dpi=300)
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   plt.show()


sys.exit()

air = f.variables['t1p5m_gb']
wlite = f.variables['wlitev_diag']
wcarb = f.variables['wcarbv_diag']
wexp = f.variables['wexpv_diag']
fapar = f.variables['faparv_diag']
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]


W = np.zeros((2,10,360,720))


a = np.minimum(wlite, wcarb)
W = np.minimum(a, wexp)

print W

a = W[:,:,:,:]/wlite[:,:,:,:] == 1 
b = W[:,:,:,:]/wcarb[:,:,:,:] == 1 
c = W[:,:,:,:]/wexp[:,:,:,:] == 1 

W = a + 2*b + 3*c
print W

cMap = ListedColormap([sns.xkcd_rgb["medium green"], sns.xkcd_rgb["denim blue"],sns.xkcd_rgb["pale red"]])
#cMap = ListedColormap(['faded green', 'denim blue','pale red'])



#for i in range(10):
m = Basemap(projection='robin', lon_0=0,resolution='l')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)

fig = plt.figure(figsize=(48, 48)) 
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

mdata = maskoceans(lon, lat, W[0,9,:,:])


cs = m.contourf(x,y,mdata,cmap=cMap)
#cs = m.contourf(x,y,air[0,:,:],20)

heatmap = plt.pcolor(mdata, cmap=cMap)
cbar = m.colorbar(heatmap)
#cbar.solids.set_edgecolor("face")

cbar.ax.get_yaxis().set_ticks([])
for j, lab in enumerate(['Light','Carbon','Export']):
    #cbar.ax.text(.5, (2 * j + 1) / 8.0, lab, ha='center', va='center')
     cbar.ax.text(.5, (3* j+ 1) / 8.0, lab, ha='center', va='center', rotation=270)     #0.5, 1/8, 0
                                                                          #0.5, 3/8, 1
                                                                          #0.5, 5/8, 2
cbar.ax.get_yaxis().labelpad = 60
cbar.ax.set_ylabel('Limiting regime', rotation=270)


#cbar.set_ylabel('Limiting regimes', rotation=270)
#cbar.set_ticks([0,100])
plt.title(r'01/01/2000 ')
#plt.savefig('Export_limited.png')
#plt.savefig('Limiting_regimes_bottom.png')
plt.show()

sys.exit()

m = Basemap(projection='robin', lon_0=0,resolution='c')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

cs = m.contourf(x,y,wlite[0,0,:,:]*1e6,20)
cbar = m.colorbar()
#cbar.solids.set_edgecolor("face")
#cbar.set_ticks([0,100])

plt.title(r'01/01/2000 - W$_{lite}$')
plt.savefig('wlite.png')
plt.show()

m = Basemap(projection='robin', lon_0=0,resolution='c')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

cs = m.contourf(x,y,wcarb[0,0,:,:]*1e6,20)
cbar = m.colorbar()
cbar.solids.set_edgecolor("face")
#cbar.set_ticks([0,100])
plt.title(r'01/01/2000 - W$_{carb}$')
plt.savefig('wcarb.png')
plt.show()

m = Basemap(projection='robin', lon_0=0,resolution='c')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

cs = m.contourf(x,y,wexp[0,0,:,:]*1e6,20)
cbar = m.colorbar()
cbar.solids.set_edgecolor("face")
#cbar.set_ticks([0,100])
plt.title(r'01/01/2000 - W$_{exp}$')
plt.savefig('wexp.png')
plt.show()

m = Basemap(projection='robin', lon_0=0,resolution='c')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

cs = m.contourf(x,y,fapar[0,0,:,:],20)
cbar = m.colorbar()
cbar.solids.set_edgecolor("face")
#cbar.set_ticks([0,100])
plt.title(r'01/01/2000 - fAPAR')
plt.savefig('fapar.png')
plt.show()

m = Basemap(projection='robin', lon_0=0,resolution='c')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

cs = m.contourf(x,y,W[0,9,:,:]/wlite[0,9,:,:]==1,20)
cbar = m.colorbar()
cbar.solids.set_edgecolor("face")
#cbar.set_ticks([0,100])
plt.title(r'01/01/2000 - Light limited')
#plt.savefig('Light_limited.png')
plt.show()

m = Basemap(projection='robin', lon_0=0,resolution='c')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

cs = m.contourf(x,y,W[0,9,:,:]/wcarb[0,9,:,:] == 1,20)
cbar = m.colorbar()
cbar.solids.set_edgecolor("face")
#cbar.set_ticks([0,100])
plt.title(r'01/01/2000 - Carbon limited')
#plt.savefig('Carbon_limited.png')
plt.show()

m = Basemap(projection='robin', lon_0=0,resolution='c')
  # resolution c, l, i, h, f in that order
x, y = m(lon, lat)
m.drawmapboundary(fill_color='white', zorder=-1)
m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
m.drawcoastlines(color='0.6', linewidth=0.5)
m.drawcountries(color='0.6', linewidth=0.5)
m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')
m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

cs = m.contourf(x,y,W[0,9,:,:]/wexp[0,9,:,:] == 1,20)
cbar = m.colorbar()
cbar.solids.set_edgecolor("face")
#cbar.set_ticks([0,100])
plt.title(r'01/01/2000 - Export limited')
plt.savefig('Export_limited.png')
#plt.savefig('Limiting_regimes.png')
plt.show()
