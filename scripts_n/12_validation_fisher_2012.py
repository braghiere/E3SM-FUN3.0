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

from scipy import stats

mpl.font_manager.findfont('Times')
# High resolution cost lines
#http://introtopython.org/visualization_earthquakes.html
# High resolution cost lines
#http://basemaptutorial.readthedocs.io/en/latest/utilities.html
#plt.style.use('ggplot')
plt.style.use('seaborn-white')
SIZE = 48
plt.rc('font', size=SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=SIZE)  # fontsize of the x any y labels
plt.rc('xtick', labelsize=SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SIZE)  # legend fontsize
plt.rc('figure', titlesize=SIZE)  # # size of the figure title

#f = Dataset('files/time_average_h0.nc','r')
f = Dataset('files/fix_global_v6_funp_time_avg.nc','r')
f11 = Dataset('files/fix_global_v6_fun_time_avg.nc','r')
f0 = Dataset('../fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')
f1 = Dataset('/home/renato/ELM_FUNP/gimms_nlim_1p875_2p5.nc','r')


lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
gpp = f.variables['GPP'][:]*365*24*60*60
gpp0 = f0.variables['GPP'][:]*365*24*60*60
pactive = f.variables['NUPTAKE_NPP_FRACTION'][:]
pactive_n = f11.variables['NUPTAKE_NPP_FRACTION'][:]
pactive_p = f.variables['PUPTAKE_NPP_FRACTION'][:]
leafc_to_litter = f.variables['LEAFC_TO_LITTER'][:]
leafp = f.variables['LEAFP'][:]
leafc = f.variables['LEAFC'][:]

nlim = f1.variables['Band1'][:]
lon1 = f1.variables['lon'][:]




from sklearn.linear_model import LinearRegression
model = LinearRegression(fit_intercept=True)

#nutri_lim = 1000.*pactive
nutri_lim_np = 100.*(pactive+pactive_p)/np.max((pactive+pactive_p))
nutri_lim_n = 100.*(pactive_n)/np.max((pactive_n))

nutri_lim_np,lon = shiftgrid(180., nutri_lim_np, lon, start=False)
lon = f.variables['lon'][:]
nutri_lim_n,lon = shiftgrid(180., nutri_lim_n, lon, start=False)
lon = f.variables['lon'][:]
gpp,lon = shiftgrid(180., gpp, lon, start=False)
lon = f.variables['lon'][:]
gpp0,lon = shiftgrid(180., gpp0, lon, start=False)
lon = f.variables['lon'][:]

rel_gpp=100.*gpp/gpp0
rel_gpp[rel_gpp<0.]=np.nan
rel_gpp[rel_gpp>100.]=np.nan

#plt.imshow(rel_gpp)
#plt.show()
#sys.exit()


nlim[nlim==0.]=np.nan
#nlim[nlim<1.]=np.nan
nlim[gpp<100.]=np.nan
nlim=np.ma.masked_array(nlim,np.isnan(nlim))

plt.imshow(nlim)
plt.show()
plt.close()

nlim = np.average(nlim,axis=1)

#nutri_lim = rel_gpp
#nutri_lim[nutri_lim<1.]=np.nan
nutri_lim_np[gpp<50.]=np.nan
nutri_lim_n[gpp<50.]=np.nan

nutri_lim_np=np.ma.masked_array(nutri_lim_np,np.isnan(nutri_lim_np))
nutri_lim_n=np.ma.masked_array(nutri_lim_n,np.isnan(nutri_lim_n))


nutri_lim_np = np.average(nutri_lim_np,axis=1)
nutri_lim_n = np.average(nutri_lim_n,axis=1)


plt.plot(nlim,label="Nlim")
plt.plot(nutri_lim_np,label="E3SM-FUN3")
plt.plot(nutri_lim_n,label="E3SM-FUN2")
plt.show()
plt.close()

#for i in xrange(len(nutri_lim[:,0])):
#   for j in xrange(len(nutri_lim[0,:])):
#       print nlim[i,j],nutri_lim[i,j] 
#sys.exit()

#plt.imshow(nutri_lim - nlim)
#plt.show()


#nlim[nlim<5.]=np.nan
#nutri_lim[nutri_lim<5.]=np.nan

#nlim=np.ma.masked_array(nlim,np.isnan(nlim))
#nutri_lim=np.ma.masked_array(nutri_lim,np.isnan(nutri_lim))

mask = ~np.isnan(nlim) & ~np.isnan(nutri_lim_np)
mask1 = ~np.isnan(nlim) & ~np.isnan(nutri_lim_n)

####### Scatter plot #####
xfit = np.arange(0,100.,1)

#model.fit(nlim[mask].flatten()[:, np.newaxis],nutri_lim[mask].flatten())
#print("Model slope: ", model.coef_[0])
#print("Model intercept: ", model.intercept_)

#r1 = model.score(nlim.flatten()[:, np.newaxis],nutri_lim.flatten())

#yfit1 = model.predict(xfit[:, np.newaxis])



slope1, intercept1, r1, p1, std1 = stats.linregress(nlim[mask].flatten(),nutri_lim_np[mask].flatten())

print(nlim[mask].flatten(),nutri_lim_np[mask].flatten())

yfit1= slope1*xfit + intercept1

print("intercept: ", intercept1)
print("slope: ", slope1)
print("R2: ", r1)

slope2, intercept2, r2, p1, std2 = stats.linregress(nlim[mask1].flatten(),nutri_lim_n[mask1].flatten())

print(nlim[mask1].flatten(),nutri_lim_n[mask1].flatten())

yfit2= slope2*xfit + intercept2

print("intercept: ", intercept2)
print("slope: ", slope2)
print("R2: ", r2)


nlim_nan = np.ma.filled(nlim,0)
nutri_lim_np_nan = np.ma.filled(nutri_lim_np,0)
nutri_lim_n_nan = np.ma.filled(nutri_lim_n,0)

yfit11= slope1*nlim_nan + intercept1

yfit22= slope2*nlim_nan + intercept2

#print(yfit)

rmse1 = np.sqrt(sum((yfit11 - nutri_lim_np_nan)**2)/len(nutri_lim_np_nan))
print("RMSE1: ", rmse1)

rmse2 = np.sqrt(sum((yfit22 - nutri_lim_n_nan)**2)/len(nutri_lim_n_nan))
print("RMSE2: ", rmse2)

fig = plt.figure(figsize=(36, 36)) 
        
#From gCm-2s-1 to TgCyr-1
const = (365*24*60*60)*1000
     
plt.scatter(nlim,nutri_lim_n,label=r'E3SM-FUN2 (R$^{2}$ = %.2f; RMSE = %.2f)' % (r2, rmse2),linewidth=9.3,s=600,color='darkorange')
plt.scatter(nlim,nutri_lim_np,label=r'E3SM-FUN3 (R$^{2}$ = %.2f; RMSE = %.2f)' % (r1, rmse1),linewidth=9.3,s=600,color='green')

plt.plot(np.arange(0,100,1),np.arange(0,100,1),label='1:1',color='k',linestyle='--',linewidth=9.3)
              
plt.plot(xfit,yfit1,linewidth=9.3,color='green')
plt.plot(xfit,yfit2,linewidth=9.3,color='darkorange')


 
max_val = np.max(30.)

plt.legend(ncol=1, prop={'size': 60})
axes = plt.gca()
axes.set_xlim([0.0,max_val + 0.1*max_val])
axes.set_ylim([0.0,max_val + 0.1*max_val])



plt.grid(True)

plt.tick_params(axis='both', which='major', pad=20)
plt.xlabel(r'Nutrient limitation from Fisher et al. (2012) (%)')
plt.ylabel(r'Normalized C use ratio for nutrient uptake (%)' )

#plt.show()
#plt.savefig('figures/Nitrogen/scatter_elm_fisher2012_noarid.png',bbox_inches="tight")

plt.close()



#plt.scatter(nlim,1000.*pactive)
#plt.plot(np.arange(0,100,1),np.arange(0,100,1),'k--')
#plt.show()
#sys.exit()




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


nlim = f1.variables['Band1'][:]
pactive = nlim



#pactive = np.ma.masked_where(pactive > 1., pactive)
pactive = np.ma.masked_where(pactive > 100., pactive)
pactive = np.ma.masked_where(pactive < 0., pactive)
#pactive[pactive<np.max(pactive)/100.]=0.0

#plt.imshow(pactive[:,:])
#plt.show()
#sys.exit()

#v = np.zeros((96,144))






#v = 100.*pactive
v = 1.*pactive
#v = nutri_lim-nlim
 


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

upper = plt.cm.viridis(np.arange(256))

lower = np.ones((int(256/4),4))

for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

cmap = np.vstack(( lower, upper))

cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


max_val = np.max(v)
max_val = 100.
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
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=1.0, color='0.5',fontsize='xx-large')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=1.0, color='0.5',fontsize='xx-large')


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
   #cbar.ax.set_xlabel('C use ratio (%)', rotation=0,color='black', size=78)
   cbar.ax.set_xlabel('Total Nutrient limitation (%)', rotation=0,color='black', size=78*1.5)
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   #no coloredge
   #cbar.solids.set_edgecolor("face")
   cbar.solids.set_edgecolor("black")
   cbar.solids.set_linewidth(6)
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'Retranslocation', fontname='Times', fontsize=92,pad=26)
   #cbar.ax.tick_params(labelsize='xx-large')
   cbar.ax.tick_params(labelsize=92)
   #plt.savefig('em_steindinger_tot.pdf',bbox_inches="tight",dpi=300)
   plt.savefig('figures/Nitrogen/Nlim_Fisher2012.png',bbox_inches="tight")
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   #plt.show()
   plt.close()



sys.exit()

