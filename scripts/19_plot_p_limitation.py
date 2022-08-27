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
f = Dataset('/home/renato/ELM_FUNP/updated_color_scheme/files/fix_global_v6_funp_time_avg.nc','r')
f1 = Dataset('/home/renato/ELM_FUNP/fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
gpp = f.variables['GPP'][:]
plant_pdemand = f.variables['PLANT_PDEMAND'][:]
puptake = f.variables['PUPTAKE'][:]
nuptake = f.variables['NUPTAKE'][:]

pretrans = f.variables['PRETRANS'][:]
nretrans = f.variables['NRETRANS'][:]

leafp = f.variables['LEAFP'][:]
leafn = f.variables['LEAFN'][:]

solutionp = f.variables['SOLUTIONP'][:]



# shifting grid to run from -180 to 180 rather than 0-360
#clump,lon = shiftgrid(180., clump, lon, start=False)


#lon = f.variables['lon'][:]
plant_pdemand,lon = shiftgrid(180., plant_pdemand, lon, start=False)
lon = f.variables['lon'][:]
puptake,lon = shiftgrid(180., puptake, lon, start=False)
lon = f.variables['lon'][:]
solutionp,lon = shiftgrid(180., solutionp, lon, start=False)
lon = f.variables['lon'][:]
leafn,lon = shiftgrid(180., leafn, lon, start=False)
lon = f.variables['lon'][:]
leafp,lon = shiftgrid(180., leafp, lon, start=False)
lon = f.variables['lon'][:]
nretrans,lon = shiftgrid(180., nretrans, lon, start=False)
lon = f.variables['lon'][:]
pretrans,lon = shiftgrid(180., pretrans, lon, start=False)

m = Basemap(projection='robin', lon_0=0.,resolution='l')



x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)





puptake = puptake[:,:]*365*24*60*60
plant_pdemand = plant_pdemand[:,:]*365*24*60*60

npp_diff = 100.*((puptake/plant_pdemand))

npp_diff[npp_diff>100.0] = 100.0
npp_diff[npp_diff<0] = 0.0

#npp_diff = npp_diff/(solutionp[:,:])

#npp_diff = npp_diff

#npp_diff.filled()

#print(np.max(abs(npp_diff)))

#npp_diff = npp_diff/np.max(abs(npp_diff))

npp_diff = np.ma.filled(npp_diff,0.0)

#print(npp_diff)


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
upper = plt.cm.YlGnBu(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


max_val = np.max(abs(npp_diff)/1.)

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
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.jet,extend='both')

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
   plt.savefig('Plimitation.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()


   plt.close()

#sys.exit()
########################################Nitrogen############################

f = Dataset('/home/renato/ELM_FUNP/updated_color_scheme/files/fix_global_v6_funp_time_avg.nc','r')
f1 = Dataset('/home/renato/ELM_FUNP/fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
gpp = f.variables['GPP'][:]
plant_ndemand = f.variables['PLANT_NDEMAND'][:]
nuptake = f.variables['NUPTAKE'][:]

sminn = f.variables['SMINN'][:]

# shifting grid to run from -180 to 180 rather than 0-360
#clump,lon = shiftgrid(180., clump, lon, start=False)


#lon = f.variables['lon'][:]
gpp,lon = shiftgrid(180., gpp, lon, start=False)
lon = f.variables['lon'][:]
plant_ndemand,lon = shiftgrid(180., plant_ndemand, lon, start=False)
lon = f.variables['lon'][:]
nuptake,lon = shiftgrid(180., nuptake, lon, start=False)
lon = f.variables['lon'][:]
sminn,lon = shiftgrid(180., sminn, lon, start=False)

m = Basemap(projection='robin', lon_0=0.,resolution='l')



x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)


plant_ndemand = plant_ndemand[:,:]*365*24*60*60
nuptake = nuptake[:,:]*365*24*60*60

npp_diff = 100.*(1.-(nuptake/plant_ndemand))

npp_diff[npp_diff>100.0] = 100.0
npp_diff[npp_diff<0] = 0.0

#npp_diff = npp_diff/(sminn[:,:])

#npp_diff = npp_diff

#npp_diff.filled()

#print(np.max(abs(npp_diff)))

npp_diff = npp_diff/np.max(abs(npp_diff))

npp_diff = np.ma.filled(npp_diff,0.0)

#print(npp_diff)


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


max_val = np.max(abs(npp_diff)/1.)

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
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.jet,extend='both')

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
   plt.savefig('Nlimitation.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()


   plt.close()


#######################N&PLIMITATIONN#####################################




nlim = (nuptake/plant_ndemand)
nlim = nlim/np.max(abs(nlim))

plim = (puptake/plant_pdemand)
plim = plim/np.max(abs(plim))

npp_diff = 100.*(nlim - plim) 

npp_diff[npp_diff>100.0] = 100.0
npp_diff[npp_diff<0] = 0.0

#################################

multiplier = gpp/np.max(abs(gpp))


leafnp = leafn/leafp
#leafnp = leafnp/np.max(abs(leafnp))
retransnp = nretrans/pretrans
#retransnp = retransnp/np.max(abs(retransnp))

limnp = leafnp/retransnp
limnp = limnp/np.max(abs(limnp))

limpn = -1.*retransnp/leafnp
limpn = limpn/np.max(abs(limpn))


count_n = 0
count_p = 0
count_np = 0


for i in range(len(limnp[:,0])):
  for j in range(len(limnp[0,:])):
     if(retransnp[i,j]>=(leafnp[i,j] + leafnp[i,j]*0.01)):
        limnp[i,j] = limnp[i,j]*np.cos(np.deg2rad(lat[i]))*multiplier[i,j]
        count_n = count_n + 1
     elif(retransnp[i,j]>leafnp[i,j]):
         limnp[i,j] = (limnp[i,j]/100.)*np.cos(np.deg2rad(lat[i]))*multiplier[i,j]
         count_np = count_np + 1
     elif(retransnp[i,j]<=(leafnp[i,j] - leafnp[i,j]*0.01)):
         limnp[i,j] = limpn[i,j]*np.cos(np.deg2rad(lat[i]))*multiplier[i,j]
         count_p = count_p + 1
     else:
         limnp[i,j] = (limpn[i,j]/100.)*np.cos(np.deg2rad(lat[i]))*multiplier[i,j]
         count_np = count_np + 1

print('N lim =',100.*count_n/(count_n+count_p+count_np),'%')
print('P lim =',100.*count_p/(count_n+count_p+count_np),'%')
print('NP lim =',100.*count_np/(count_n+count_p+count_np),'%')

#sys.exit()
npp_diff = limnp

#npp_diff = npp_diff/(sminn[:,:])

#npp_diff = npp_diff

#npp_diff.filled()

#print(np.max(abs(npp_diff)))


f = Dataset('/mnt/n_p_lim_0p5.nc','r')


lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
npp_diff = f.variables['Band1'][:]


lon = f.variables['lon'][:]
npp_diff,lon = shiftgrid(180., npp_diff, lon, start=False)

x1,y1 = np.meshgrid(lon, lat) 


#npp_diff = npp_diff/np.max(abs(npp_diff))

npp_diff = np.ma.filled(npp_diff,0.0)




#print(npp_diff)


x4 = np.linspace(x1[0][0],x1[0][-1],x1.shape[1]*8)
y4 = np.linspace(y1[0][0],y1[-1][0],y1.shape[0]*16)

#x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*40)
#y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*60)

x3,y3 = np.meshgrid(x4,y4)
X3,Y3 = m(x3, y3)

#order=0 for nearest-neighbor, order=1 for bilinear, order=3 cubic
data2 = interp(npp_diff,x1[0],y1[:,0],x3,y3,order=1)

data3 = interp(limnp,x[0],y[:,0],x2,y2,order=1)

print(np.shape(data2), np.shape(data3))

#plt.imshow(data3-data2)
#plt.show()

#sys.exit()




mdata = maskoceans(x3, y3, data2,resolution='l',grid=1.25,inlands=True)

conc_matrix = np.ones(np.shape(data3))

conc_matrix[(data3/data2) < 0.] = -1.

agree = np.zeros(np.shape(data3))

#print(x4,y4)
#print(np.shape(x4),np.shape(y4))


for i in range(len(data3[:,0])):
  for j in range(len(data3[0,:])):
     if(abs(data3[i,j])>=abs(data2[i,j]) and conc_matrix[i,j] == 1.):
       agree[i,j] = abs(data2[i,j])/abs(data3[i,j])
     if(abs(data3[i,j])<abs(data2[i,j]) and conc_matrix[i,j] == 1.):       
       agree[i,j] = abs(data3[i,j])/abs(data2[i,j])
     if(abs(data3[i,j])>=abs(data2[i,j]) and conc_matrix[i,j] == -1.):       
       agree[i,j] = -1.*abs(data2[i,j])/abs(data3[i,j])
     if(abs(data3[i,j])<abs(data2[i,j]) and conc_matrix[i,j] == -1.):       
       agree[i,j] = -1*abs(data3[i,j])/abs(data2[i,j])
     print('%2.4f done!, agree = %1.2f, data2 = %1.2f, data3 = %1.2f' % (float(100.*i/len(data3[:,0])),agree[i,j], data2[i,j], data3[i,j]))

#plt.imshow(agree)
#plt.show()
#sys.exit()
     
mdata = maskoceans(x3, y3, agree,resolution='l',grid=1.25,inlands=True)
#mdata = maskoceans(x, y, v)

###################wwriting to netcdf file################

import netCDF4 as nc

fn = 'agreementmap_elm_du2020.nc'
ds = nc.Dataset(fn, 'w', format='NETCDF4')

lat = ds.createDimension('lat', len(y4))
lon = ds.createDimension('lon', len(x4))

lats = ds.createVariable('lat', 'f4', ('lat',))
lons = ds.createVariable('lon', 'f4', ('lon',))
value = ds.createVariable('agree_index', 'f4', ('lat', 'lon',))
value.units = '-'

lats[:] = y4
lons[:] = x4

value[:, :] = mdata

ds.close()

#My colorbar

#upper = plt.cm.jet(np.arange(256))
upper = plt.cm.bwr(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


max_val = np.max(abs(npp_diff)/1.)
max_val = 1.0

#levels=[-0.01251*max_val,0.,0.01251*max_val,0.0251*max_val,0.0625*max_val,.1*max_val,.2*max_val,.3*max_val,.4*max_val,.5*max_val,.6*max_val,.7*max_val,.75*max_val,.8*max_val,.9*max_val,1.*max_val]

levels=[-1.*max_val,-0.8*max_val,-0.6*max_val,-0.4*max_val,-.2*max_val,-.1*max_val,-.05*max_val,0.,.05*max_val,.1*max_val,.2*max_val,.4*max_val,.6*max_val,.8*max_val,1.*max_val]

print(levels)

for i in range(1):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1, linewidth=4.5)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.0', linewidth=4.5)
   #m.drawcountries(color='0.', linewidth=4.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=1.0, color='0.5',fontsize='xx-large', fontname='Times')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=1.0, color='0.5',fontsize='xx-large', fontname='Times')

 
   #PLOT ABSOLUTE
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap,extend='both')
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cmap)
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.bwr,extend='both')


   import matplotlib.colors as mcol
   cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["royalblue","gainsboro","mediumpurple","gainsboro","crimson"])
   norm=plt.Normalize(-0.5,0.5)

   #cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.gist_stern,extend='both')
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=cm1,norm=norm,extend='both')
   cs = m.contourf(X3,Y3,mdata,levels,cmap=plt.cm.BrBG,norm=norm,extend='both')
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.jet,extend='both')
   #cs = m.contourf(X2,Y2,mdata,levels,cmap=plt.cm.coolwarm,extend='both')

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
   #cbar.ax.set_xlabel(r'More P limited $\longleftarrow \longrightarrow$ More N limited', rotation=0,color='black', size=78, fontname='Times')
   #cbar.ax.set_xlabel(r'More P limited $\longleftarrow$            NP co-limited            $\longrightarrow$ More N limited', rotation=0,color='black', size=78*1.5, fontname='Times')
   cbar.ax.set_xlabel(r'Low agreement $\longleftarrow$            Medium agreement            $\longrightarrow$ High agreement', rotation=0,color='black', size=78*1.5, fontname='Times')
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   #no coloredge
   #cbar.solids.set_edgecolor("face")
   cbar.solids.set_edgecolor("black")
   cbar.solids.set_linewidth(6)
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'$\Delta$ Net Primary Productivity Normalized by Soil P', fontname='Times', fontsize=92,pad=26)
   cbar.ax.tick_params(labelsize=92)
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('diff_NPlimitation_du2020_v4.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()

   plt.close()

