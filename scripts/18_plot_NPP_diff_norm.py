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
f = Dataset('files/fix_global_v6_funp_time_avg.nc','r')
f1 = Dataset('../fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
gpp = f.variables['GPP'][:]
npp_elm = f1.variables['NPP'][:]
npp_fun = f.variables['NPP'][:]

solutionp = f.variables['SOLUTIONP'][:]
leafn = f.variables['LEAFN'][:]
leafp = f.variables['LEAFP'][:]

# shifting grid to run from -180 to 180 rather than 0-360
#clump,lon = shiftgrid(180., clump, lon, start=False)


#lon = f.variables['lon'][:]
npp_elm,lon = shiftgrid(180., npp_elm, lon, start=False)
lon = f.variables['lon'][:]
npp_fun,lon = shiftgrid(180., npp_fun, lon, start=False)
lon = f.variables['lon'][:]
solutionp,lon = shiftgrid(180., solutionp, lon, start=False)
lon = f.variables['lon'][:]
leafn,lon = shiftgrid(180., leafn, lon, start=False)
lon = f.variables['lon'][:]
leafp,lon = shiftgrid(180., leafp, lon, start=False)


m = Basemap(projection='robin', lon_0=0.,resolution='l')



x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)





npp_elm = npp_elm[:,:]*365*24*60*60
npp_fun = npp_fun[:,:]*365*24*60*60

npp_diff = (npp_fun - npp_elm)

npp_diff = npp_diff/(solutionp[:,:])

npp_diff = npp_diff

npp_diff.filled()

print(np.max(abs(npp_diff)))


npp_diff = npp_diff/np.max(abs(npp_diff))

npp_diff = np.ma.filled(npp_diff,0.0)

print(npp_diff)


#x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*20)
#y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*30)

x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*40)
y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*60)

x2,y2 = np.meshgrid(x2,y2)
X2,Y2 = m(x2, y2)

#order=0 for nearest-neighbor, order=1 for bilinear, order=3 cubic
data2 = interp(npp_diff,x[0],y[:,0],x2,y2,order=1)



mdata = maskoceans(x2, y2, data2,resolution='l',grid=1.25,inlands=True)
#mdata = maskoceans(x, y, v)

#My colorbar

#upper = plt.cm.jet(np.arange(256))
upper = plt.cm.bwr(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


max_val = np.max(abs(npp_diff)/2.)

#levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

levels=[-1.*max_val,-0.8*max_val,-0.6*max_val,-0.4*max_val,-.2*max_val,-.1*max_val,-.05*max_val,0.,.05*max_val,.1*max_val,.2*max_val,.4*max_val,.6*max_val,.8*max_val,1.*max_val]

print(levels)

for i in range(1):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1, linewidth=4.5)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.0', linewidth=4.5)
   #m.drawcountries(color='0.', linewidth=4.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')

 
   #PLOT ABSOLUTE
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap,extend='both')
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap)
   cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.bwr,extend='both')

   #vmin,vmax = (-0.01251*max_val,max_val)
   vmin,vmax = (-1.*max_val,1.*max_val)

   plt.tight_layout()
   #cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=np.linspace(vmin,vmax,7),format='%.1f')

   #levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

   #levels=[-1.*max_val,0.8*max_val,0.6*max_val,0.4*max_val,.2*max_val,0.,.2*max_val,.4*max_val,.6*max_val,.8*max_val,1.*max_val]

   #ticks = [0.,0.0251*max_val,.1*max_val,.3*max_val,.5*max_val,.7*max_val,.8*max_val,1.*max_val]

   ticks = [-1.*max_val,-.6*max_val,-.2*max_val,-.05*max_val,0.,.05*max_val,.2*max_val,.6*max_val,1.*max_val]


   cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=ticks,format='%.2f')
   #cbar.ax.get_yaxis().labelpad = 60
   cbar.ax.get_xaxis().labelpad = 45
   #cbar.ax.set_ylabel('EM (%)', rotation=270)
   cbar.ax.set_xlabel(r'$\Delta$ Net Primary Productivity Normalized by Soil P', rotation=0,color='black', size=78, fontname='Times')
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   #no coloredge
   #cbar.solids.set_edgecolor("face")
   cbar.solids.set_edgecolor("black")
   cbar.solids.set_linewidth(6)
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'$\Delta$ Net Primary Productivity Normalized by Soil P', fontname='Times', fontsize=92,pad=26)
   cbar.ax.tick_params(labelsize='xx-large')
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('figures/delta_NPP_norm.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()


   plt.close()


########################################Nitrogen############################

f = Dataset('files/fix_global_v6_funp_time_avg.nc','r')
f1 = Dataset('../fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
gpp = f.variables['GPP'][:]
npp_elm = f1.variables['NPP'][:]
npp_fun = f.variables['NPP'][:]

sminn = f.variables['SMINN'][:]

# shifting grid to run from -180 to 180 rather than 0-360
#clump,lon = shiftgrid(180., clump, lon, start=False)


#lon = f.variables['lon'][:]
npp_elm,lon = shiftgrid(180., npp_elm, lon, start=False)
lon = f.variables['lon'][:]
npp_fun,lon = shiftgrid(180., npp_fun, lon, start=False)
lon = f.variables['lon'][:]
sminn,lon = shiftgrid(180., sminn, lon, start=False)

m = Basemap(projection='robin', lon_0=0.,resolution='l')



x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)


npp_elm = npp_elm[:,:]*365*24*60*60
npp_fun = npp_fun[:,:]*365*24*60*60

npp_diff = (npp_fun - npp_elm)

npp_diff = npp_diff/(sminn[:,:])

npp_diff = npp_diff

npp_diff.filled()

print(np.max(abs(npp_diff)))


npp_diff = npp_diff/np.max(abs(npp_diff))

npp_diff = np.ma.filled(npp_diff,0.0)

print(npp_diff)


#x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*20)
#y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*30)

x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*40)
y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*60)

x2,y2 = np.meshgrid(x2,y2)
X2,Y2 = m(x2, y2)

#order=0 for nearest-neighbor, order=1 for bilinear, order=3 cubic
data2 = interp(npp_diff,x[0],y[:,0],x2,y2,order=1)



mdata = maskoceans(x2, y2, data2,resolution='l',grid=1.25,inlands=True)
#mdata = maskoceans(x, y, v)

#My colorbar

#upper = plt.cm.jet(np.arange(256))
upper = plt.cm.bwr(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


max_val = np.max(abs(npp_diff)/2.)

#levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

levels=[-1.*max_val,-0.8*max_val,-0.6*max_val,-0.4*max_val,-.2*max_val,-.1*max_val,-.05*max_val,0.,.05*max_val,.1*max_val,.2*max_val,.4*max_val,.6*max_val,.8*max_val,1.*max_val]

print(levels)

for i in range(1):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1, linewidth=4.5)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.0', linewidth=4.5)
   #m.drawcountries(color='0.', linewidth=4.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')

 
   #PLOT ABSOLUTE
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap,extend='both')
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap)
   cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.bwr_r,extend='both')

   #vmin,vmax = (-0.01251*max_val,max_val)
   vmin,vmax = (-1.*max_val,1.*max_val)

   plt.tight_layout()
   #cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=np.linspace(vmin,vmax,7),format='%.1f')

   #levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

   #levels=[-1.*max_val,0.8*max_val,0.6*max_val,0.4*max_val,.2*max_val,0.,.2*max_val,.4*max_val,.6*max_val,.8*max_val,1.*max_val]

   #ticks = [0.,0.0251*max_val,.1*max_val,.3*max_val,.5*max_val,.7*max_val,.8*max_val,1.*max_val]

   ticks = [-1.*max_val,-.6*max_val,-.2*max_val,-.05*max_val,0.,.05*max_val,.2*max_val,.6*max_val,1.*max_val]


   cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=ticks,format='%.2f')
   #cbar.ax.get_yaxis().labelpad = 60
   cbar.ax.get_xaxis().labelpad = 45
   #cbar.ax.set_ylabel('EM (%)', rotation=270)
   cbar.ax.set_xlabel(r'$\Delta$ Net Primary Productivity Normalized by Soil N', rotation=0,color='black', size=78, fontname='Times')
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   #no coloredge
   #cbar.solids.set_edgecolor("face")
   cbar.solids.set_edgecolor("black")
   cbar.solids.set_linewidth(6)
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'$\Delta$ Net Primary Productivity Normalized by Soil P', fontname='Times', fontsize=92,pad=26)
   cbar.ax.tick_params(labelsize='xx-large')
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('figures/delta_NPP_norm_Nitrogen.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()


   plt.close()



####################### N&P limitation #####################

npp_diff = (npp_fun - npp_elm)

npp_diff_n = npp_diff/(sminn[:,:])
npp_diff_p = npp_diff/(solutionp[:,:]*leafn[:,:]/leafp[:,:])


npp_diff_n = npp_diff_n/np.max(abs(npp_diff_n))
npp_diff_p = npp_diff_p/np.max(abs(npp_diff_p))


multiplier = npp_elm/np.max(abs(npp_elm))

count_n = 0
count_p = 0
count_np = 0

for i in range(len(npp_diff[:,0])):
  for j in range(len(npp_diff[0,:])):
     if( abs(npp_diff_n[i,j]) >= abs(npp_diff_p[i,j]) + 0.01*abs(npp_diff_p[i,j]) ):
        npp_diff[i,j] = -1.*npp_diff_n[i,j]*multiplier[i,j]/np.cos(np.deg2rad(lat[i]))
        count_n = count_n + 1
     elif( abs(npp_diff_n[i,j]) > abs(npp_diff_p[i,j]) ):
         npp_diff[i,j] = -1.*npp_diff_n[i,j]*multiplier[i,j]/np.cos(np.deg2rad(lat[i]))
         count_np = count_np + 1
     elif( abs(npp_diff_n[i,j]) < abs(npp_diff_p[i,j]) - 0.01*abs(npp_diff_p[i,j]) ):
         npp_diff[i,j] = 1.*npp_diff_p[i,j]*multiplier[i,j]/np.cos(np.deg2rad(lat[i]))
         count_np = count_np + 1
     else:
        npp_diff[i,j] = npp_diff_p[i,j]*multiplier[i,j]/np.cos(np.deg2rad(lat[i]))
        count_p = count_p + 1 




print('N lim =',100.*count_n/(count_n+count_p+count_np),'%')
print('P lim =',100.*count_p/(count_n+count_p+count_np),'%')
print('NP lim =',100.*count_np/(count_n+count_p+count_np),'%')

npp_diff.filled()

print(np.max(abs(npp_diff)))


npp_diff = npp_diff/np.max(abs(npp_diff))

npp_diff = np.ma.filled(npp_diff,0.0)

print(npp_diff)


#x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*20)
#y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*30)

x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*40)
y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*60)

x2,y2 = np.meshgrid(x2,y2)
X2,Y2 = m(x2, y2)

#order=0 for nearest-neighbor, order=1 for bilinear, order=3 cubic
data2 = interp(npp_diff,x[0],y[:,0],x2,y2,order=1)



mdata = maskoceans(x2, y2, data2,resolution='l',grid=1.25,inlands=True)
#mdata = maskoceans(x, y, v)

#My colorbar

#upper = plt.cm.jet(np.arange(256))
upper = plt.cm.bwr(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


max_val = np.max(abs(npp_diff)/2.)

#levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

levels=[-1.*max_val,-0.8*max_val,-0.6*max_val,-0.4*max_val,-.2*max_val,-.1*max_val,-.05*max_val,0.,.05*max_val,.1*max_val,.2*max_val,.4*max_val,.6*max_val,.8*max_val,1.*max_val]

print(levels)

for i in range(1):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1, linewidth=4.5)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.0', linewidth=4.5)
   #m.drawcountries(color='0.', linewidth=4.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=1.0, color='0.5',fontsize='x-large', fontname='Times')

 
   #PLOT ABSOLUTE
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap,extend='both')
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap)
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.bwr,extend='both')
   cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.coolwarm,extend='both')

   #vmin,vmax = (-0.01251*max_val,max_val)
   vmin,vmax = (-1.*max_val,1.*max_val)

   plt.tight_layout()
   #cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=np.linspace(vmin,vmax,7),format='%.1f')

   #levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

   #levels=[-1.*max_val,0.8*max_val,0.6*max_val,0.4*max_val,.2*max_val,0.,.2*max_val,.4*max_val,.6*max_val,.8*max_val,1.*max_val]

   #ticks = [0.,0.0251*max_val,.1*max_val,.3*max_val,.5*max_val,.7*max_val,.8*max_val,1.*max_val]

   ticks = [-1.*max_val,-.6*max_val,-.2*max_val,-.05*max_val,0.,.05*max_val,.2*max_val,.6*max_val,1.*max_val]


   cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=ticks,format='%.2f')
   #cbar.ax.get_yaxis().labelpad = 60
   cbar.ax.get_xaxis().labelpad = 45
   #cbar.ax.set_ylabel('EM (%)', rotation=270)
   cbar.ax.set_xlabel(r'More P limited $\longleftarrow$                    NP co-limited                    $\longrightarrow$ More N limited', rotation=0,color='black', size=78, fontname='Times')
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   #no coloredge
   #cbar.solids.set_edgecolor("face")
   cbar.solids.set_edgecolor("black")
   cbar.solids.set_linewidth(6)
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'$\Delta$ Net Primary Productivity Normalized by Soil P', fontname='Times', fontsize=92,pad=26)
   cbar.ax.tick_params(labelsize='xx-large')
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('figures/delta_NPP_norm_NP.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()


   plt.close()


