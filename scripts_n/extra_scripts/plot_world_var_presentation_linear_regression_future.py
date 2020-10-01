from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from netCDF4 import Dataset, date2index
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
import numpy.ma as ma
import numpy as np
from mpl_toolkits.basemap import maskoceans, interp
from matplotlib.colors import ListedColormap
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn import metrics 
import scipy.stats as stats

def movingaverage(interval, window_size):
    window = np.ones(int(window_size)/(window_size))
    return np.convolve(interval,window,'same')

def smooth1(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_padded = np.pad(y, (box_pts//2,box_pts-1-box_pts//2), mode='edge')
    y_smooth = np.convolve(y_padded, box, mode='valid')
    return y_smooth

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

#future map
f = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_all/ISSP585Clm50BgcCrop_all.NPP_NUPTAKE.nc','r')
fi = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_all/ISSP585Clm50BgcCrop_all.NUPTAKE.nc','r')
fj = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_all/ISSP585Clm50BgcCrop_all.PLANT_NDEMAND.nc','r')
fk = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_all/ISSP585Clm50BgcCrop_all.NPP.nc','r')
#present map
f1 = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_myc/ISSP585Clm50BgcCrop_myc.NPP_NUPTAKE.nc','r')
f1i = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_myc/ISSP585Clm50BgcCrop_myc.NUPTAKE.nc','r')
f1j = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_myc/ISSP585Clm50BgcCrop_myc.PLANT_NDEMAND.nc','r')
f1k = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/ISSP585Clm50BgcCrop_myc/ISSP585Clm50BgcCrop_myc.NPP.nc','r')
f2 = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSteidinger/IHistClm50BgcSteidinger.NPP_NUPTAKE.nc','r')
f2i = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSteidinger/IHistClm50BgcSteidinger.NUPTAKE.nc','r')
f2j = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSteidinger/IHistClm50BgcSteidinger.PLANT_NDEMAND.nc','r')
f2k = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSteidinger/IHistClm50BgcSteidinger.NPP.nc','r')
f3 = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSoudzi/IHistClm50BgcSoudzi.NPP_NUPTAKE.nc','r')
f3i = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSoudzi/IHistClm50BgcSoudzi.NUPTAKE.nc','r')
f3j = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSoudzi/IHistClm50BgcSoudzi.PLANT_NDEMAND.nc','r')
f3k = Dataset('/home/renato/datasets/IHistClm50Bgc_outputs/IHistClm50BgcSoudzi/IHistClm50BgcSoudzi.NPP.nc','r')

surf = Dataset('/home/renato/Steindinger/surfdata_1.9x2.5_hist_78pfts_CMIP6_simyr1850_c190304.nc','r')
params = Dataset('/home/renato/Steindinger/clm5_params.c171117.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
time = f.variables['time'][:]
var = 'NPP_NUPTAKE'
nee = f1.variables[var]
unit = nee.units
long_name = nee.long_name
nee_stein1 = f1.variables[var][:]
nee_stein2 = f2.variables[var][:]
nee_stein3 = f3.variables[var][:]
var = 'NPP_NUPTAKE'
nee = f.variables[var][:]
nee_stein1 = f1.variables[var][:]
var = 'NUPTAKE'
neei = fi.variables[var][:]
neei_stein1 = f1i.variables[var][:]
neei_stein2 = f2i.variables[var][:]
neei_stein3 = f3i.variables[var][:]
var = 'PLANT_NDEMAND'
neej = fj.variables[var][:]
neej_stein1 = f1j.variables[var][:]
neej_stein2 = f2j.variables[var][:]
neej_stein3 = f3j.variables[var][:]
var = 'NPP'
neek = fk.variables[var][:]
neek_stein1 = f1k.variables[var][:]
neek_stein2 = f2k.variables[var][:]
neek_stein3 = f3k.variables[var][:]

# Function to calculate the weights of latitude
#1.875 longitude X 1.25 latitude)
radius = 6367449
res_lon = 1.875
#res_lon = 1.90
res_lat = 1.25
m_pi = 3.14159265358979323846
dtor = 360./(2*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat/dtor
lat_u = (lat + res_lat)/dtor

weights = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))

average = np.ma.average(nee,axis=1,weights=weights)
average = np.ma.average(average,axis=1)
average_stein1 = np.ma.average(nee_stein1,axis=1,weights=weights)
average_stein1 = np.ma.average(average_stein1,axis=1)
average_stein2 = np.ma.average(nee_stein2,axis=1,weights=weights)
average_stein2 = np.ma.average(average_stein2,axis=1)
average_stein3 = np.ma.average(nee_stein3,axis=1,weights=weights)
average_stein3 = np.ma.average(average_stein3,axis=1)

averagei = np.ma.average(neei,axis=1,weights=weights)
averagei = np.ma.average(averagei,axis=1)
averagei_stein1 = np.ma.average(neei_stein1,axis=1,weights=weights)
averagei_stein1 = np.ma.average(averagei_stein1,axis=1)
averagei_stein2 = np.ma.average(neei_stein2,axis=1,weights=weights)
averagei_stein2 = np.ma.average(averagei_stein2,axis=1)
averagei_stein3 = np.ma.average(neei_stein3,axis=1,weights=weights)
averagei_stein3 = np.ma.average(averagei_stein3,axis=1)

averagej = np.ma.average(neej,axis=1,weights=weights)
averagej = np.ma.average(averagej,axis=1)
averagej_stein1 = np.ma.average(neej_stein1,axis=1,weights=weights)
averagej_stein1 = np.ma.average(averagej_stein1,axis=1)
averagej_stein2 = np.ma.average(neej_stein2,axis=1,weights=weights)
averagej_stein2 = np.ma.average(averagej_stein2,axis=1)
averagej_stein3 = np.ma.average(neej_stein3,axis=1,weights=weights)
averagej_stein3 = np.ma.average(averagej_stein3,axis=1)

averagek = np.ma.average(neek,axis=1,weights=weights)
averagek = np.ma.average(averagek,axis=1)
averagek_stein1 = np.ma.average(neek_stein1,axis=1,weights=weights)
averagek_stein1 = np.ma.average(averagek_stein1,axis=1)
averagek_stein2 = np.ma.average(neek_stein2,axis=1,weights=weights)
averagek_stein2 = np.ma.average(averagek_stein2,axis=1)
averagek_stein3 = np.ma.average(neek_stein3,axis=1,weights=weights)
averagek_stein3 = np.ma.average(averagek_stein3,axis=1)

#plt.plot(averagej,'k',label='PLANT_NDEMAND')
#plt.plot(averagei,'b',label='NUPTAKE')
#plt.legend(loc=4)
#plt.savefig('noffer_vs_demand.png')
#plt.show()
#sys.exit()


slope_map = np.zeros((96,144))
slopei_map = np.zeros((96,144))
slopej_map = np.zeros((96,144))
slopek_map = np.zeros((96,144))

nee_av = nee
nee_av_stein1 = nee_stein1
nee_av_stein2 = nee_stein2
nee_av_stein3 = nee_stein3

neei_av = neei
neei_av_stein1 = neei_stein1
neei_av_stein2 = neei_stein2
neei_av_stein3 = neei_stein3

neej_av = neej
neej_av_stein1 = neej_stein1
neej_av_stein2 = neej_stein2
neej_av_stein3 = neej_stein3

neek_av = neek
neek_av_stein1 = neek_stein1
neek_av_stein2 = neek_stein2
neek_av_stein3 = neek_stein3


const_TgN= (1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)
const_PgC= (365*24*60*60)*148847000*(1000*1000)/(10**15)

#smoothing the whole matrix 
for i in xrange(96):
  for j in xrange(144):
     nee_av[:,i,j] = smooth1(nee[:,i,j],12)*const_PgC
     nee_av_stein1[:,i,j] = smooth1(nee_stein1[:,i,j],12)*const_PgC
     nee_av_stein2[:,i,j] = smooth1(nee_stein2[:,i,j],12)*const_PgC
     nee_av_stein3[:,i,j] = smooth1(nee_stein3[:,i,j],12)*const_PgC

     neei_av[:,i,j] = smooth1(neei[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)
     neei_av_stein1[:,i,j] = smooth1(neei_stein1[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)
     neei_av_stein2[:,i,j] = smooth1(neei_stein2[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)
     neei_av_stein3[:,i,j] = smooth1(neei_stein3[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)

     neej_av[:,i,j] = smooth1(neej[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)
     neej_av_stein1[:,i,j] = smooth1(neej_stein1[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)
     neej_av_stein2[:,i,j] = smooth1(neej_stein2[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)
     neej_av_stein3[:,i,j] = smooth1(neej_stein3[:,i,j],12)*(1000)*(365*24*60*60)*148847000*(1000*1000)/(10**15)

     neek_av[:,i,j] = smooth1(neek[:,i,j],12)*const_PgC
     neek_av_stein1[:,i,j] = smooth1(neek_stein1[:,i,j],12)*const_PgC
     neek_av_stein2[:,i,j] = smooth1(neek_stein2[:,i,j],12)*const_PgC
     neek_av_stein3[:,i,j] = smooth1(neek_stein3[:,i,j],12)*const_PgC




for i in xrange(96):
  for j in xrange(144):
     #NPP_NUPTAKE
     slope, intercept, r_value, p_value, std_err = stats.linregress(time[:],nee_av[:,i,j])
     slope_stein1, intercept, r_value, p_value, std_err = stats.linregress(time[:],nee_av_stein1[:,i,j])
     #slope_stein2, intercept, r_value, p_value, std_err = stats.linregress(time[:],nee_av_stein2[:,i,j])
     #slope_stein3, intercept, r_value, p_value, std_err = stats.linregress(time[:],nee_av_stein3[:,i,j])
     #NUPTAKE
     slopei, intercept, r_value, p_value, std_err = stats.linregress(time[:],neei_av[:,i,j])
     slopei_stein1, intercept, r_value, p_value, std_err = stats.linregress(time[:],neei_av_stein1[:,i,j])
     #slopei_stein2, intercept, r_value, p_value, std_err = stats.linregress(time[:],neei_av_stein2[:,i,j])
     #slopei_stein3, intercept, r_value, p_value, std_err = stats.linregress(time[:],neei_av_stein3[:,i,j])
     #PLANT_NDEMAND
     slopej, intercept, r_value, p_value, std_err = stats.linregress(time[:],neej_av[:,i,j])
     slopej_stein1, intercept, r_value, p_value, std_err = stats.linregress(time[:],neej_av_stein1[:,i,j])
     #slopej_stein2, intercept, r_value, p_value, std_err = stats.linregress(time[:],neej_av_stein2[:,i,j])
     #slopej_stein3, intercept, r_value, p_value, std_err = stats.linregress(time[:],neej_av_stein3[:,i,j])
     #NPP
     slopek, intercept, r_value, p_value, std_err = stats.linregress(time[:],neek_av[:,i,j])
     slopek_stein1, intercept, r_value, p_value, std_err = stats.linregress(time[:],neek_av_stein1[:,i,j])
     #slopek_stein2, intercept, r_value, p_value, std_err = stats.linregress(time[:],neek_av_stein2[:,i,j])
     #slopek_stein3, intercept, r_value, p_value, std_err = stats.linregress(time[:],neek_av_stein3[:,i,j])



     if((slopej > 0.0 and slopei > 0.0 and abs(slopej) > abs(slopei))): 
        #slope_map[i,j] =  (slopei/slopej) - 1.
        slope_map[i,j] =  (slopei/slopej)
     #if((slopej < 0.0 and slopei < 0.0 and abs(slopej) > abs(slopei))): 
      #  slope_map[i,j] =   (slopei/slopej)
        #slope_map[i,j] =  1. - (slopei/slopej)
     #if(slope > 0.0):
     #   slope_map[i,j] = slope
     #if(slope < 0.0):
     #   slope_map[i,j] = slope
     #if(slope == 0.0):
     #   slope_map[i,j] = slope


slope_map = slope_map/abs(np.max(slope_map))
slope_map = 1. - slope_map

#NPP_NUPTAKE
#slope_map = nee_av[:,:,:] - nee_av_stein1[:,:,:]
#NPP
slope_map = 100.*(neek_av - neek_av_stein1)/neek_av_stein1
#slope_map = 1000.*(neek_av - neek_av_stein1)

slope_map = (neek_av - neek_av_stein1)*(365*24*60*60)/const_PgC

#LAST 10 YEARS AVERAGE
#slope_map = np.ma.average(slope_map[((12*1)-1)*51:((12*1)-1)*61,:,:],axis=0)
#ALL
slope_map = np.ma.average(slope_map[:,:,:],axis=0)

print(np.min(slope_map), np.max(slope_map))
print(np.mean(slope_map), np.std(slope_map))

lon = f.variables['lon'][:]
slope_map,lon = shiftgrid(180., slope_map, lon, start=False)


pft_map_nat = surf.variables['PCT_NAT_PFT'][:]
lon = f.variables['lon'][:]
pft_map_nat,lon = shiftgrid(180., pft_map_nat, lon, start=False)

slope_map = ma.array(slope_map,mask=[pft_map_nat[0,:,:]>90.])



m = Basemap(projection='robin', lon_0=0.,resolution='l')
x,y = np.meshgrid(lon, lat) 
X,Y = m(x, y)

x2 = np.linspace(x[0][0],x[0][-1],x.shape[1]*20)
y2 = np.linspace(y[0][0],y[-1][0],y.shape[0]*30)

x2,y2 = np.meshgrid(x2,y2)
X2,Y2 = m(x2, y2)

data2 = interp(slope_map,x[0],y[:,0],x2,y2,order=1)

mdata = maskoceans(x2, y2, data2,resolution='h',grid=1.25,inlands=True)



for i in xrange(1):
   fig = plt.figure(figsize=(48, 48)) 
   m.drawmapboundary(fill_color='white', zorder=-1)
   m.fillcontinents(color='0.8', lake_color='white', zorder=0)
 
   m.drawcoastlines(color='0.6', linewidth=0.5)
   m.drawcountries(color='0.6', linewidth=0.5)
   m.drawparallels(np.arange(-90.,91.,30.), labels=[1,0,0,1],    dashes=[1,1], linewidth=0.25, color='0.5',fontsize='x-large')
   m.drawmeridians(np.arange(0., 360., 60.), labels=[1,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5',fontsize='x-large')

 
   #PLOT ABSOLUTE
   levels = np.linspace(-8.,8.,21)
   vmin,vmax = (-8,8)
   cs = m.contourf(X2,Y2,mdata,levels,vmin=vmin,vmax=vmax,cmap=plt.cm.jet,extend='both')
   
   #PLOT DIFFERENCE
   #vmin,vmax = (0.25,1.)
   vmean = np.mean(slope_map)
   vstd = np.std(slope_map)
   #vmin,vmax = (-1000.,1000.)
   #if(abs(np.min(slope_map)) > abs(np.max(slope_map))):
   #   vmin,vmax = (np.mean(slope_map), -1.*np.mean(slope_map))
   #else:
   #   vmin,vmax = (-1.*np.mean(slope_map), np.mean(slope_map))
   levels = np.linspace(vmin,vmax,16,endpoint=True)   
   cs = m.contourf(X2,Y2,mdata,levels,vmin=vmin,vmax=vmax,cmap=plt.cm.Spectral,extend='both')
   

   plt.tight_layout()
   cbar = m.colorbar(cs,location='bottom',pad='10%',ticks=np.linspace(vmin,vmax,7),format='%.1f')
   #cbar.ax.set_xticklabels(['Low risk','Medium risk','High risk'])
   #cbar.ax.get_yaxis().labelpad = 60
   cbar.ax.get_xaxis().labelpad = 60
   #cbar.ax.set_ylabel('EM (%)', rotation=270)
   cbar.ax.set_xlabel(r'$\Delta$ NPP (gC.m$^{-2}$.yr$^{-1}$)', rotation=0,color='black', size=78)
   #cbar.ax.set_xlabel(r'$\Delta$ NPP (%)', rotation=0,color='black', size=78)
   #cbar.ax.set_xlabel('ECM tree basal area (%)', rotation=0,color='black', size=78)
   cbar.solids.set_edgecolor("face")
   #cbar.set_clim(0.0,100)
   cbar.set_clim(vmin,vmax)
   #plt.title(r'Sulman et al. (2019) - Shi et al. (2016)')
   cbar.ax.tick_params(labelsize='xx-large')
   #plt.annotate(r'$\Delta$ NPP = %.2f +- %.2f PgC.yr$^{-1}$' % ((np.mean(slope_map)) , (np.std(slope_map))), xy=(0.35, 0), xycoords='axes fraction')
   plt.savefig('delt_npp_future_present.png',bbox_inches="tight")
   #plt.savefig('diff_em_tot_3.pdf',bbox_inches="tight",dpi=300)
   #plt.savefig('ecm_orig_shi_1p9x2p5.png',bbox_inches="tight",dpi=300)
   plt.show()

sys.exit()

