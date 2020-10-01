from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
from netCDF4 import Dataset, date2index
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
import numpy.ma as ma
import numpy as np
from mpl_toolkits.basemap import maskoceans
from matplotlib.colors import ListedColormap
import seaborn as sns
from scipy import stats 

# High resolution cost lines
#http://introtopython.org/visualization_earthquakes.html
# High resolution cost lines
#http://basemaptutorial.readthedocs.io/en/latest/utilities.html
plt.style.use('seaborn-white')
SIZE = 78
plt.rc('font', size=SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=SIZE)  # fontsize of the x any y labels
plt.rc('xtick', labelsize=SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SIZE)  # legend fontsize
plt.rc('figure', titlesize=SIZE)  # # size of the figure title

import scipy.stats as st

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = a.count()
    m, se = np.mean(a,axis=1), scipy.stats.sem(a,axis=1)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h


# lots of tips in http://www.subcriticalflow.com/2015/11/computing-zonal-and-weighted-statistics.html



f = Dataset('files/elm_pft_monthly_mean.nc','r')


lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
pft = f.variables['pft'][:]
time = f.variables['time'][:]

#print(f)
#sys.exit()

puptake = f.variables['PNONMYC']
#unit = puptake.units
#long_name = puptake.long_name

puptake = f.variables['PNONMYC'][:]
pnonmyc = f.variables['PNONMYC'][:]
pactive = f.variables['PACTIVE'][:]
pretrans = f.variables['PRETRANS'][:]
#pam = f.variables['PAM'][:]
#pecm = f.variables['PECM'][:]


puptake,lon = shiftgrid(180., puptake, lon, start=False)
lon = f.variables['lon'][:]
pnonmyc,lon = shiftgrid(180., pnonmyc, lon, start=False)
lon = f.variables['lon'][:]
pactive,lon = shiftgrid(180., pactive, lon, start=False)
lon = f.variables['lon'][:]
pretrans,lon = shiftgrid(180., pretrans, lon, start=False)
#lon = f.variables['lon'][:]
#pam,lon = shiftgrid(180., pam, lon, start=False)
#lon = f.variables['lon'][:]
#pecm,lon = shiftgrid(180., pecm, lon, start=False)


#Function to calculate the weights of latitude
radius = 6367449
res_lat = 1.9
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat/dtor
lat_u = (lat + res_lat)/dtor
weights = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))




gpp_mte_ori_mean = []
gpp_mte_ori_mean_pnonmyc = []
gpp_mte_ori_mean_pactive = []
gpp_mte_ori_mean_pretrans = []
gpp_mte_ori_ci34_0 = []
gpp_mte_ori_ci34_1 = []

for i in xrange(1):

      #print np.shape(npp[i])
      #gpp_mte_ori_mean.append(np.mean(npp[i],axis = 1))
      gpp_mte_ori_ci34_0.append(st.t.interval(0.341, np.nanmean(puptake,axis = 1).count()-1, loc=np.nanmean(puptake,axis = 1), scale=np.std(puptake,axis = 1))[0])

      gpp_mte_ori_ci34_1.append(st.t.interval(0.341, np.nanmean(puptake,axis = 1).count()-1, loc=np.nanmean(puptake,axis = 1), scale=np.std(puptake,axis = 1))[1])
      gpp_mte_ori_mean.append(np.mean(puptake,axis = 1))
      gpp_mte_ori_mean_pnonmyc.append(np.mean(pnonmyc,axis = 1))
      gpp_mte_ori_mean_pactive.append(np.mean(pactive,axis = 1))
      gpp_mte_ori_mean_pretrans.append(np.mean(pretrans,axis = 1))


print(np.shape(gpp_mte_ori_mean_pnonmyc))

sys.exit()
my_pfts = ['All']
x = range(len(lat))

for i, pft in enumerate(my_pfts):

   	fig = plt.figure(figsize=(48, 24)) 
        
        #From gCm-2s-1 to TgCyr-1
        const = (365*24*60*60)*weights/(10**12)
        plt.stackplot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),np.array(gpp_mte_ori_mean_pnonmyc[0][:])*(const),np.array(gpp_mte_ori_mean_pam[0][:])*(const),np.array(gpp_mte_ori_mean_pecm[0][:])*(const),np.array(gpp_mte_ori_mean_pretrans[0][:])*(const),edgecolor='k',linewidth=3.3,colors=['chartreuse','b','r','yellow'],labels=['Direct root uptake','AM root uptake','ECM root uptake','Retranslocation']) 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),np.array(gpp_mte_ori_mean[0][:])*(const),'k',linewidth=15, label = 'Shi et al. (2016)') 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),( gpp_mte_ori_ci34_1[0][:])*(const),'k--',linewidth=7, label = r'$\pm 1. \sigma$ ') 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),(gpp_mte_ori_ci34_0[0][:])*(const),'k--',linewidth=7)
        #plt.fill_between(np.floor(((np.asarray(x) + 1)*res_lat)-90.), (gpp_mte_ori_ci34_0[0][:])*(const), (gpp_mte_ori_ci34_1[0][:])*(const) ,alpha=0.3, facecolor='k')



	plt.legend()
	axes = plt.gca()
	axes.set_xlim([-90.,90.])
	#axes.set_ylim([0.0,3500.0])

        plt.xticks(np.arange(min((np.asarray(x)*res_lat) -90.), max((np.asarray(x)*res_lat) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'Latitude ($\degree$)')
	plt.ylabel(r'Total annual P uptake (Tg P yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('figures/PUPTAKE_monthly_mean.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	plt.show()


        sys.exit()



