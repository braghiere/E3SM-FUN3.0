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

mpl.font_manager.findfont('Times')
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
plt.rcParams["font.family"] = "sans-serif"  # # font

pftname_fig = [ "not_vegetated",
  "needleleaf_evergreen_temperate",
  "needleleaf_evergreen_boreal",
  "needleleaf_deciduous_boreal",
  "broadleaf_evergreen_tropical",
  "broadleaf_evergreen_temperate",
  "broadleaf_deciduous_tropical",
  "broadleaf_deciduous_temperate",
  "broadleaf_deciduous_boreal",
  "evergreen_shrub",
  "deciduous_temperate_shrub",
  "deciduous_boreal_shrub",
  "c3_arctic_grass",
  "c3_non-arctic_grass",
  "c4_grass",
  "c3_crop",
  "c3_irrigated"]

pftname_fig_v1 = [ "broadleaf_deciduous",
  "grassland",
  "needleleaf_evergreen",
  "needleleaf_deciduous",
  "broadleaf_evergreen",
  "shrubland"]





#f = Dataset('files/time_average_h0.nc','r')
f1 = Dataset('files/fix_global_v6_funp_time_avg.nc','r')
lat = f1.variables['lat'][:]
lon = f1.variables['lon'][:]

f = Dataset('files/surfdata.nc','r')


pft = f.variables['natpft'][:]
#gpp = f.variables['GPP'][:]
pactive = f.variables['PCT_NAT_PFT'][:]




# shifting grid to run from -180 to 180 rather than 0-360
#clump,lon = shiftgrid(180., clump, lon, start=False)


#lon = f.variables['lon'][:]
pactive,lon = shiftgrid(180., pactive, lon, start=False)


m = Basemap(projection='robin', lon_0=0.,resolution='l')



x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)


pactive = [(pactive[6,:,:]+pactive[7,:,:]+pactive[8,:,:]), \
           (pactive[12,:,:]+pactive[13,:,:]+pactive[14,:,:]+ \
                                          pactive[15,:,:]), \
           (pactive[1,:,:]+pactive[2,:,:]), \
           (pactive[3,:,:]), \
           (pactive[4,:,:]+pactive[5,:,:]), \
           (pactive[10,:,:]+pactive[11,:,:])]

pactive=np.array(pactive)

print(np.shape(pactive))

v = pactive[:,:,:]


#x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*20)
#y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*30)

x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*40)
y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*60)

x2,y2 = np.meshgrid(x2,y2)
X2,Y2 = m(x2, y2)

data2 = []
for i in range(6):
  print("Interpolating %s ..." % pftname_fig_v1[i])
  #order=0 for nearest-neighbor, order=1 for bilinear, order=3 cubic
  data2.append(interp(v[i,:,:],x[0],y[:,0],x2,y2,order=1))

data2 = np.array(data2)


mdata = []

for i in range(6):
   print("Masking oceans of %s ..." % pftname_fig_v1[i])
   mdata.append(maskoceans(x2, y2, data2[i,:,:],resolution='l',grid=1.25,inlands=True))

mdata = np.array(mdata)
#My colorbar

#upper = plt.cm.jet(np.arange(256))
upper = plt.cm.viridis(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


#max_val = np.max(v)
max_val = 100.

levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]


for i in range(6):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1, linewidth=4.5)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.0', linewidth=4.5)
   #m.drawcountries(color='0.', linewidth=4.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')

 
   #PLOT ABSOLUTE
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap,extend='both')
   cs = m.contourf(X2,Y2,mdata[i,:,:],levels,cmap=cmap)
   vmin,vmax = (-0.01251*max_val,max_val)


   plt.tight_layout()
   #cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=np.linspace(vmin,vmax,7),format='%.1f')
   ticks = [0.,0.0251*max_val,.1*max_val,.3*max_val,.5*max_val,.7*max_val,.8*max_val,1.*max_val]
   cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=ticks,format='%.2f')
   #cbar.ax.get_yaxis().labelpad = 60
   cbar.ax.get_xaxis().labelpad = 45
   #cbar.ax.set_ylabel('EM (%)', rotation=270)
   cbar.ax.set_xlabel('Plant Functional Type (%)', rotation=0,color='black', size=78, fontname='Times')
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   #no coloredge
   #cbar.solids.set_edgecolor("face")
   cbar.solids.set_edgecolor("black")
   cbar.solids.set_linewidth(6)
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   plt.title(r'%s' % pftname_fig_v1[i], fontname='Times', fontsize=92,pad=26)
   cbar.ax.tick_params(labelsize='xx-large')
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('figures/PFT_sum_%s.png' % pftname_fig_v1[i],bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()


sys.exit()

