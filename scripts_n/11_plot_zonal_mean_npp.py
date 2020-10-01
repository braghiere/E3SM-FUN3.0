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



#f = Dataset('files/time_average_h0.nc','r')
#f1 = Dataset('files/E3SM_latest_time_average_1994_2005.nc','r')

f = Dataset('files/fix_global_v7_funp_time_avg.nc','r')
#f1 = Dataset('../fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')
#f2 = Dataset('../fix_global_v3_1994_2005/files/fix_global_v3_fun_time_avg.nc','r')

f1 = Dataset('../fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')
f2 = Dataset('files/fix_global_v7_fun_time_avg.nc','r')

lon1 = f1.variables['lon'][:]
npp_elm = f1.variables['NPP'][:]
npp_elm,lon1 = shiftgrid(180., npp_elm, lon1, start=False)

lon2 = f2.variables['lon'][:]
npp_elm_fun = f2.variables['NPP'][:]
npp_elm_fun,lon2 = shiftgrid(180., npp_elm_fun, lon2, start=False)

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]



puptake = f.variables['NUPTAKE']
unit = puptake.units
long_name = puptake.long_name

puptake = f.variables['NUPTAKE'][:]
pnonmyc = f.variables['NPP'][:]

pactive = f.variables['NACTIVE'][:]
pretrans = f.variables['NRETRANS'][:]
pam = f.variables['NAM'][:]
pecm = f.variables['NECM'][:]


lon = f.variables['lon'][:]
puptake,lon = shiftgrid(180., puptake, lon, start=False)
lon = f.variables['lon'][:]
pnonmyc,lon = shiftgrid(180., pnonmyc, lon, start=False)
lon = f.variables['lon'][:]
pactive,lon = shiftgrid(180., pactive, lon, start=False)
lon = f.variables['lon'][:]
pretrans,lon = shiftgrid(180., pretrans, lon, start=False)
lon = f.variables['lon'][:]
pam,lon = shiftgrid(180., pam, lon, start=False)
lon = f.variables['lon'][:]
pecm,lon = shiftgrid(180., pecm, lon, start=False)







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


#print(weights,np.shape(weights))
#sys.exit()

#PUPTAKE
puptake = np.array(puptake)
puptake[puptake==1e+36]=np.nan

puptake=np.ma.masked_array(puptake,np.isnan(puptake))
average=np.ma.average(puptake,axis=0,weights=weights)
variance=np.ma.average((puptake-average)**2,axis=0,weights=weights)

mean_corr = np.mean(np.ma.average(puptake,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(puptake))))

print 'Global puptake =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

#PNONMYC
pnonmyc = np.array(pnonmyc)
pnonmyc[pnonmyc==1e+36]=np.nan

pnonmyc=np.ma.masked_array(pnonmyc,np.isnan(pnonmyc))
average=np.ma.average(pnonmyc,axis=0,weights=weights)
variance=np.ma.average((pnonmyc-average)**2,axis=0,weights=weights)

mean_corr = np.mean(np.ma.average(pnonmyc,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pnonmyc))))

print 'Global pnonmyc =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 
#PECM
pecm = np.array(pecm)
pecm[pecm==1e+36]=np.nan

pecm=np.ma.masked_array(pecm,np.isnan(pecm))
average=np.ma.average(pecm,axis=0,weights=weights)
variance=np.ma.average((pecm-average)**2,axis=0,weights=weights)

mean_corr = np.mean(np.ma.average(pecm,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pecm))))

print 'Global pecm =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 
#PAM
pam = np.array(pam)
pam[pam==1e+36]=np.nan

pam=np.ma.masked_array(pam,np.isnan(pam))
average=np.ma.average(pam,axis=0,weights=weights)
variance=np.ma.average((pam-average)**2,axis=0,weights=weights)

mean_corr = np.mean(np.ma.average(pam,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pam))))

print 'Global pam =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

#PRETRANS
pretrans = np.array(pretrans)
pretrans[pretrans==1e+36]=np.nan

pretrans=np.ma.masked_array(pretrans,np.isnan(pretrans))
average=np.ma.average(pretrans,axis=0,weights=weights)
variance=np.ma.average((pretrans-average)**2,axis=0,weights=weights)

mean_corr = np.mean(np.ma.average(pretrans,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pretrans))))

print 'Global pretrans =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

#NPP
npp_elm = np.array(npp_elm)
npp_elm[npp_elm==1e+36]=np.nan

npp_elm=np.ma.masked_array(npp_elm,np.isnan(npp_elm))
average=np.ma.average(npp_elm,axis=0,weights=weights)
variance=np.ma.average((npp_elm-average)**2,axis=0,weights=weights)

mean_corr = np.mean(np.ma.average(npp_elm,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_elm))))

print 'Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

#NPP_FUN
npp_elm_fun = np.array(npp_elm_fun)
npp_elm_fun[npp_elm_fun==1e+36]=np.nan

npp_elm_fun=np.ma.masked_array(npp_elm_fun,np.isnan(npp_elm_fun))
average=np.ma.average(npp_elm_fun,axis=0,weights=weights)
variance=np.ma.average((npp_elm_fun-average)**2,axis=0,weights=weights)

mean_corr = np.mean(np.ma.average(npp_elm_fun,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_elm_fun))))

print 'Global NPP_FUN =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'



gpp_mte_ori_mean = []
gpp_mte_ori_mean_pnonmyc = []
gpp_mte_ori_mean_npp_elm = []
gpp_mte_ori_mean_pecm = []
gpp_mte_ori_mean_pam = []
gpp_mte_ori_mean_pretrans = []
gpp_mte_ori_ci34_0 = []
gpp_mte_ori_ci34_1 = []
gpp_mte_ori_mean_npp_elm_fun = []

for i in xrange(1):

      #print np.shape(npp[i])
      #gpp_mte_ori_mean.append(np.mean(npp[i],axis = 1))
      gpp_mte_ori_ci34_0.append(st.t.interval(0.341, np.nanmean(puptake,axis = 1).count()-1, loc=np.nanmean(puptake,axis = 1), scale=np.std(puptake,axis = 1))[0])

      gpp_mte_ori_ci34_1.append(st.t.interval(0.341, np.nanmean(puptake,axis = 1).count()-1, loc=np.nanmean(puptake,axis = 1), scale=np.std(puptake,axis = 1))[1])
      gpp_mte_ori_mean.append(np.mean(puptake,axis = 1))
      gpp_mte_ori_mean_pnonmyc.append(np.mean(pnonmyc,axis = 1))
      gpp_mte_ori_mean_npp_elm.append(np.mean(npp_elm,axis = 1))
      gpp_mte_ori_mean_pecm.append(np.mean(pecm,axis = 1))
      gpp_mte_ori_mean_pam.append(np.mean(pam,axis = 1))
      gpp_mte_ori_mean_pretrans.append(np.mean(pretrans,axis = 1))
      gpp_mte_ori_mean_npp_elm_fun.append(np.mean(npp_elm_fun,axis = 1))

#sys.exit()
my_pfts = ['All']
x = range(len(lat))

for i, pft in enumerate(my_pfts):

   	fig = plt.figure(figsize=(48, 24)) 
        
        #From gCm-2s-1 to TgCyr-1
        const = (365*24*60*60)
        plt.stackplot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),np.array(gpp_mte_ori_mean_pnonmyc[0][:])*(const),np.array(gpp_mte_ori_mean_npp_elm_fun[0][:])*(const)-np.array(gpp_mte_ori_mean_pnonmyc[0][:])*(const),np.array(gpp_mte_ori_mean_npp_elm[0][:])*(const)-np.array(gpp_mte_ori_mean_npp_elm_fun[0][:])*(const),edgecolor='k',linewidth=3.3,colors=['darkorange','b','r'],labels=['ELMv1 -FUN3.0','ELMv1-FUN2.0','ELMv1']) 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),np.array(gpp_mte_ori_mean[0][:])*(const),'k',linewidth=15, label = 'Shi et al. (2016)') 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),( gpp_mte_ori_ci34_1[0][:])*(const),'k--',linewidth=7, label = r'$\pm 1. \sigma$ ') 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),(gpp_mte_ori_ci34_0[0][:])*(const),'k--',linewidth=7)
        #plt.fill_between(np.floor(((np.asarray(x) + 1)*res_lat)-90.), (gpp_mte_ori_ci34_0[0][:])*(const), (gpp_mte_ori_ci34_1[0][:])*(const) ,alpha=0.3, facecolor='k')

        max_val = np.max(gpp_mte_ori_mean_npp_elm[0][:]*const)

	plt.legend()
	axes = plt.gca()
	axes.set_xlim([-90.,90.])
	axes.set_ylim([0.0,max_val + 0.1*max_val])

        plt.xticks(np.arange(min((np.asarray(x)*res_lat) -90.), max((np.asarray(x)*res_lat) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'Latitude ($\degree$)')
	plt.ylabel(r'NPP (g C m$^{-2}$ yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('figures/Nitrogen/NPP_zonal_mean.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()


        sys.exit()



