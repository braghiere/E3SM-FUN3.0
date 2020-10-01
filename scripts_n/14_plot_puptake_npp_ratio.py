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


#f = Dataset('files/time_average_h0.nc','r')
f = Dataset('files/fix_global_v7_funp_time_avg.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
gpp = f.variables['GPP'][:]
pactive = f.variables['NUPTAKE_NPP_FRACTION'][:]
leafc_to_litter = f.variables['LEAFC_TO_LITTER'][:]
leafp = f.variables['LEAFP'][:]
leafc = f.variables['LEAFC'][:]

retransp_to_ppool = f.variables['RETRANSP_TO_PPOOL'][:]
free_retransp_to_ppool = f.variables['FREE_RETRANSP_TO_PPOOL'][:]
# shifting grid to run from -180 to 180 rather than 0-360
#clump,lon = shiftgrid(180., clump, lon, start=False)

paid_retransp_to_ppool = retransp_to_ppool - free_retransp_to_ppool

#lon = f.variables['lon'][:]
pactive,lon = shiftgrid(180., pactive, lon, start=False)
lon = f.variables['lon'][:]
leafp,lon = shiftgrid(180., leafp, lon, start=False)
lon = f.variables['lon'][:]
leafc,lon = shiftgrid(180., leafc, lon, start=False)
lon = f.variables['lon'][:]
leafc_to_litter,lon = shiftgrid(180., leafc_to_litter, lon, start=False)
lon = f.variables['lon'][:]
retransp_to_ppool,lon = shiftgrid(180., retransp_to_ppool, lon, start=False)

m = Basemap(projection='robin', lon_0=0.,resolution='l')



x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)





leafp_to_litter = leafc_to_litter[:,:]*(leafp[:,:]/leafc[:,:])


#leafp_to_litter = np.ma.masked_where(leafp_to_litter < 10e-10, leafp_to_litter)



leafp_to_litter[leafp_to_litter<np.max(leafp_to_litter)/100.]=np.max(leafp_to_litter)/1000.

#leafp_to_litter = leafp_to_litter.filled(fill_value=0.0)





pactive = np.ma.masked_where(pactive > 1., pactive)
pactive = np.ma.masked_where(pactive < 0., pactive)
#pactive[pactive<np.max(pactive)/100.]=0.0

#plt.imshow(pactive[:,:])
#plt.show()
#sys.exit()

v = np.zeros((96,144))

for i in xrange(144):
  for j in xrange(96):
    if (leafp_to_litter[j,i] != 0.0):
      v[j,i] = 100.*(pactive[j,i]/leafp_to_litter[j,i])
    else:
      v[j,i] = 0.0

#plt.imshow(v[:,:])
#plt.show()
#sys.exit()


v = 100.*pactive
 

#x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*20)
#y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*30)

x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*40)
y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*60)

x2,y2 = np.meshgrid(x2,y2)
X2,Y2 = m(x2, y2)

#order=0 for nearest-neighbor, order=1 for bilinear, order=3 cubic
data2 = interp(v,x[0],y[:,0],x2,y2,order=1)



mdata = maskoceans(x2, y2, data2,resolution='l',grid=1.25,inlands=True)
#mdata = maskoceans(x, y, v)

#My colorbar

upper = plt.cm.jet(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


max_val = np.max(v)
print(max_val)
#sys.exit()

#levels=[-0.001251*max_val,0.,0.001251*max_val,0.00625*max_val,.01*max_val,.05*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.8*max_val,.9*max_val,1.*max_val]
levels=[-0.0051251*max_val,0.,0.0051251*max_val,0.01251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]


for i in xrange(1):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1, linewidth=4.5)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.0', linewidth=4.5)
   #m.drawcountries(color='0.', linewidth=4.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large')


   #mdata[mdata==np.nan]=10e-20
   #mdata = mdata.filled(fill_value=10e-20)
 
   #PLOT ABSOLUTE
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap,extend='both')
   cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap)
   vmin,vmax = (-0.01251*max_val,max_val)


   plt.tight_layout()

  #levels=[-0.0051251*max_val,0.,0.0051251*max_val,0.01251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

   ticks = [0.,0.01251*max_val,.1*max_val,.3*max_val,.5*max_val,.7*max_val,.8*max_val,1.*max_val]
   #ticks = levels[1:]
   cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=ticks,format='%.1f')
   #cbar.ax.text(0.98,1,r'$\times$10$^{-3}$',va='bottom',ha='right',size=78)
   #cbar.ax.get_yaxis().labelpad = 60
   cbar.ax.get_xaxis().labelpad = 45
   #cbar.ax.set_ylabel('EM (%)', rotation=270)
   #cbar.ax.set_xlabel('g P m$^{-2}$ yr$^{-1}$', rotation=0,color='black', size=78, fontname='Times')
   cbar.ax.set_xlabel('%', rotation=0,color='black', size=78)
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   #no coloredge
   #cbar.solids.set_edgecolor("face")
   cbar.solids.set_edgecolor("black")
   cbar.solids.set_linewidth(6)
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'Retranslocation', fontname='Times', fontsize=92,pad=26)
   cbar.ax.tick_params(labelsize='xx-large')
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('figures/Nitrogen/NUPTAKE_NPP_ratio.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()


sys.exit()

