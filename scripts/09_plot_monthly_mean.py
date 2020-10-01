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

pftname = [ "not_vegetated                           ",
  "needleleaf_evergreen_temperate_tree     ",
  "needleleaf_evergreen_boreal_tree        ",
  "needleleaf_deciduous_boreal_tree        ",
  "broadleaf_evergreen_tropical_tree       ",
  "broadleaf_evergreen_temperate_tree      ",
  "broadleaf_deciduous_tropical_tree       ",
  "broadleaf_deciduous_temperate_tree      ",
  "broadleaf_deciduous_boreal_tree         ",
  "broadleaf_evergreen_shrub               ",
  "broadleaf_deciduous_temperate_shrub     ",
  "broadleaf_deciduous_boreal_shrub        ",
  "c3_arctic_grass                         ",
  "c3_non-arctic_grass                     ",
  "c4_grass                                ",
  "c3_crop                                 ",
  "c3_irrigated                            ",
  "corn                                    ",
  "irrigated_corn                          ",
  "spring_temperate_cereal                 ",
  "irrigated_spring_temperate_cereal       ",
  "winter_temperate_cereal                 ",
  "irrigated_winter_temperate_cereal       ",
  "soybean                                 ",
  "irrigated_soybean                       "]



#f = Dataset('time_average_elm_pft.nc','r')
#f = Dataset('files/elm_pft_monthly_mean.nc','r')
f = Dataset('files/fix_global_v7_funp_pft_mean.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
pft = f.variables['pft'][:]
time = f.variables['time'][:]

#unit = puptake.units
#long_name = puptake.long_name


pnonmyc_ori = f.variables['PNONMYC'][:]
pactive_ori = f.variables['PACTIVE'][:]
pretrans_ori = f.variables['PRETRANS'][:]


pnonmyc_ori,lon = shiftgrid(180., pnonmyc_ori, lon, start=False)
lon = f.variables['lon'][:]
pactive_ori,lon = shiftgrid(180., pactive_ori, lon, start=False)
lon = f.variables['lon'][:]
pretrans_ori,lon = shiftgrid(180., pretrans_ori, lon, start=False)



pactive_ori[pactive_ori > 1e30]=np.nan
pactive_ori[pnonmyc_ori > 1e30]=np.nan
pactive_ori[pretrans_ori > 1e30]=np.nan



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

time = ['JAN','FEB','MAR','APR','MAY','JUN', \
        'JUL','AUG','SEP','OCT','NOV','DEC']

for j in xrange(len(time)):
 for i in xrange(len(pft)):

   print('time=',time[j],'pft=',pft[i])
   #PACTIVE
   pactive = np.array(pactive_ori[j,i,:,:])
   pactive[pactive==1e+36]=np.nan

   pactive=np.ma.masked_array(pactive,np.isnan(pactive))
   average=np.ma.average(pactive,axis=0,weights=weights)
   variance=np.ma.average((pactive-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pactive,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pactive))))

   print 'Global pactive =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 
   #PNONMYC
   pnonmyc = np.array(pnonmyc_ori[j,i,:,:])
   pnonmyc[pnonmyc==1e+36]=np.nan

   pnonmyc=np.ma.masked_array(pnonmyc,np.isnan(pnonmyc))
   average=np.ma.average(pnonmyc,axis=0,weights=weights)
   variance=np.ma.average((pnonmyc-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pnonmyc,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pnonmyc))))

   print 'Global pnonmyc =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 

   #PRETRANS
   pretrans = np.array(pretrans_ori[j,i,:,:])
   pretrans[pretrans==1e+36]=np.nan

   pretrans=np.ma.masked_array(pretrans,np.isnan(pretrans))
   average=np.ma.average(pretrans,axis=0,weights=weights)
   variance=np.ma.average((pretrans-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pretrans,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pretrans))))

   print 'Global pretrans =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


N = len(pft)

#PNONMYC
pnonmyc_arr = []
pnonmyc_ave = []
pnonmyc_var = []
pnonmyc_std_corr = []
pnonmyc_h_dn_corr = []

#PACTIVE
pactive_arr = []
pactive_ave = []
pactive_var = []
pactive_std_corr = []
pactive_h_dn_corr = []

#PRETRANS
pretrans_arr = []
pretrans_ave = []
pretrans_var = []
pretrans_std_corr = []
pretrans_h_dn_corr = []


pactive = np.array(pactive_ori[:,:,:,:])
pactive[pactive==1e+36]=np.nan
pactive=np.ma.masked_array(pactive,np.isnan(pactive))

pnonmyc= np.array(pnonmyc_ori[:,:,:,:])
pnonmyc[pactive==1e+36]=np.nan
pnonmyc=np.ma.masked_array(pnonmyc,np.isnan(pnonmyc))

pretrans= np.array(pretrans_ori[:,:,:,:])
pretrans[pretrans==1e+36]=np.nan
pretrans=np.ma.masked_array(pretrans,np.isnan(pretrans))


pnonmyc_arr_time= np.zeros((len(time),len(pft)))
pactive_arr_time= np.zeros((len(time),len(pft)))
pretrans_arr_time= np.zeros((len(time),len(pft)))

for j in xrange(len(time)):
  for i in xrange(len(pft)):

   #PNONMYC
   pnonmyc_arr_time[j,i] = (np.mean(np.ma.average(pnonmyc[j,i,:,:], \
                    axis=0,weights=weights)))

   pnonmyc_arr.append(np.mean(np.ma.average(pnonmyc[j,i,:,:], \
                    axis=0,weights=weights)))

   pnonmyc_ave.append(np.ma.average(pnonmyc[j,i,:,:], \
                    axis=0,weights=weights))

   pnonmyc_var.append(np.ma.average((pnonmyc[j,i,:,:]-pnonmyc_ave[i])**2,axis=0,weights=weights))

   #pnonmyc_std_corr.append(np.mean(np.sqrt(pnonmyc_var[i])))

   #pnonmyc_h_dn_corr.append(stats.norm.interval(0.66, loc=pnonmyc_arr[i], scale= pnonmyc_std_corr[i]/np.sqrt(len(np.ma.array(pnonmyc[j,i,:,:]))))[0])
 
   #PACTIVE

   pactive_arr_time[j,i] = (np.mean(np.ma.average(pactive[j,i,:,:], \
                    axis=0,weights=weights)))

   pactive_arr.append(np.mean(np.ma.average(pactive[j,i,:,:], \
                    axis=0,weights=weights)))

   pactive_ave.append(np.ma.average(pactive[j,i,:,:], \
                    axis=0,weights=weights))

   pactive_var.append(np.ma.average((pactive[j,i,:,:]-pactive_ave[i])**2,axis=0,weights=weights))

   #pactive_std_corr.append(np.mean(np.sqrt(pactive_var[j,i])))

   #pactive_h_dn_corr.append(stats.norm.interval(0.66, loc=pactive_arr[j,i], scale= pactive_std_corr[j,i]/np.sqrt(len(np.ma.array(pactive[j,i,:,:]))))[0])

   #PRETRANS

   pretrans_arr_time[j,i] = (np.mean(np.ma.average(pretrans[j,i,:,:], \
                    axis=0,weights=weights)))

   pretrans_arr.append(np.mean(np.ma.average(pretrans[j,i,:,:], \
                     axis=0,weights=weights)))

   pretrans_ave.append(np.ma.average(pretrans[j,i,:,:], \
                    axis=0,weights=weights))

   pretrans_var.append(np.ma.average((pretrans[j,i,:,:]-pretrans_ave[i])**2,axis=0,weights=weights))

   #pretrans_std_corr.append(np.mean(np.sqrt(pretrans_var[j,i])))

   #pretrans_h_dn_corr.append(stats.norm.interval(0.66, loc=pretrans_arr[j,i], scale= pretrans_std_corr[j,i]/np.sqrt(len(np.ma.array(pretrans[j,i,:,:]))))[0])


#print(np.shape(pnonmyc_arr_time))

pnonmyc_arr = pnonmyc_arr_time
pactive_arr = pactive_arr_time
pretrans_arr = pretrans_arr_time

ind = np.arange(6)
width = 0.35

#pnonmyc_h_dn_corr = np.array(pnonmyc_arr) - np.array(pnonmyc_h_dn_corr)

pnonmyc_arr_glob = np.zeros((len(time),6))
pactive_arr_glob = np.zeros((len(time),6))
pretrans_arr_glob = np.zeros((len(time),6))

for j in xrange(len(time)): 
  pnonmyc_arr_glob[j] = [np.sum([pnonmyc_arr[j,6],pnonmyc_arr[j,7],pnonmyc_arr[j,8]]), \
                    np.sum([pnonmyc_arr[j,12],pnonmyc_arr[j,13],pnonmyc_arr[j,14], \
                                          pnonmyc_arr[j,15],pnonmyc_arr[j,16]]), \
                    np.sum([pnonmyc_arr[j,1],pnonmyc_arr[j,2]]), \
                    np.sum([pnonmyc_arr[j,3]]), \
                    np.sum([pnonmyc_arr[j,4],pnonmyc_arr[j,5]]), \
                    np.sum([pnonmyc_arr[j,10],pnonmyc_arr[j,11]])]

  pactive_arr_glob[j] = [np.sum([pactive_arr[j,6],pactive_arr[j,7],pactive_arr[j,8]]), \
                    np.sum([pactive_arr[j,12],pactive_arr[j,13],pactive_arr[j,14], \
                                          pactive_arr[j,15],pactive_arr[j,16]]), \
                    np.sum([pactive_arr[j,1],pnonmyc_arr[j,2]]), \
                    np.sum([pactive_arr[j,3]]), \
                    np.sum([pactive_arr[j,4],pactive_arr[j,5]]), \
                    np.sum([pactive_arr[j,10],pactive_arr[j,11]])]

  pretrans_arr_glob[j] = [np.sum([pretrans_arr[j,6],pretrans_arr[j,7],pretrans_arr[j,8]]), \
                    np.sum([pretrans_arr[j,12],pretrans_arr[j,13],pretrans_arr[j,14], \
                                          pretrans_arr[j,15],pretrans_arr[j,16]]), \
                    np.sum([pretrans_arr[j,1],pretrans_arr[j,2]]), \
                    np.sum([pretrans_arr[j,3]]), \
                    np.sum([pretrans_arr[j,4],pretrans_arr[j,5]]), \
                    np.sum([pretrans_arr[j,10],pretrans_arr[j,11]])]

#print(pnonmyc_arr_glob)

#reorganizing from most to least expensive:
pnonmyc_arr_glob = [pnonmyc_arr_glob[:,4],pnonmyc_arr_glob[:,2], \
                     pnonmyc_arr_glob[:,1],pnonmyc_arr_glob[:,0], \
                     pnonmyc_arr_glob[:,5],pnonmyc_arr_glob[:,3]]

pactive_arr_glob = [pactive_arr_glob[:,4],pactive_arr_glob[:,2], \
                     pactive_arr_glob[:,1],pactive_arr_glob[:,0], \
                     pactive_arr_glob[:,5],pactive_arr_glob[:,3]]

pretrans_arr_glob = [pretrans_arr_glob[:,4],pretrans_arr_glob[:,2], \
                     pretrans_arr_glob[:,1],pretrans_arr_glob[:,0], \
                     pretrans_arr_glob[:,5],pretrans_arr_glob[:,3]]

my_pfts = ['Evergreen broadleaf','Evergreen needleleaf','Grassland','Decidous broadleaf','Shrubland','Deciduos needleleaf']
x = range(len(time))



for i, pft in enumerate(my_pfts):

   	fig = plt.figure(figsize=(48, 24)) 
        
        #From gCm-2s-1 to TgCyr-1
        const = 30*24*60*60
        #const = (365*24*60*60)*weights/(10**12)
        plt.stackplot(time[:],np.array(pnonmyc_arr_glob[i][:])*(const),np.array(pactive_arr_glob[i][:])*(const),np.array(pretrans_arr_glob[i][:]/1.)*(const),edgecolor='k',linewidth=3.3,colors=['b','r','darkorange'],labels=['Direct root uptake','Mycorrhizal uptake','Retranslocation']) 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),np.array(gpp_mte_ori_mean[0][:])*(const),'k',linewidth=15, label = 'Shi et al. (2016)') 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),( gpp_mte_ori_ci34_1[0][:])*(const),'k--',linewidth=7, label = r'$\pm 1. \sigma$ ') 
        #plt.plot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),(gpp_mte_ori_ci34_0[0][:])*(const),'k--',linewidth=7)
        #plt.fill_between(np.floor(((np.asarray(x) + 1)*res_lat)-90.), (gpp_mte_ori_ci34_0[0][:])*(const), (gpp_mte_ori_ci34_1[0][:])*(const) ,alpha=0.3, facecolor='k')


        if(i==2):
	 plt.legend()
	axes = plt.gca()
	#axes.set_xlim([-90.,90.])
	#axes.set_ylim([0.0,3500.0])

        #plt.xticks(np.arange(min((np.asarray(x)*res_lat) -90.), max((np.asarray(x)*res_lat) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'Month')
	plt.ylabel(r'Plant P uptake (g P m$^{-2}$ month$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        ttl=plt.title(r' %s' % (pft) ,fontsize = 84)
        ttl.set_position([.5, 1.05])

        plt.savefig('figures/PUPTAKE_monthly_mean_%s.png' % (pft),bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()


        #sys.exit()







sys.exit()
pnonmyc_std_glob = [np.sqrt(np.sum([pnonmyc_h_dn_corr[6]**2,pnonmyc_h_dn_corr[7]**2,pnonmyc_h_dn_corr[8]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[12]**2,pnonmyc_h_dn_corr[13]**2,pnonmyc_h_dn_corr[14], \
                                          pnonmyc_h_dn_corr[15]**2,pnonmyc_h_dn_corr[16]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[1]**2,pnonmyc_h_dn_corr[2]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[3]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[4]**2,pnonmyc_h_dn_corr[5]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[10]**2,pnonmyc_h_dn_corr[11]**2]))]



##############################################################

pactive_h_dn_corr = np.array(pactive_arr) - np.array(pactive_h_dn_corr)


pactive_arr_glob = [np.sum([pactive_arr[6],pactive_arr[7],pactive_arr[8]]), \
                    np.sum([pactive_arr[12],pactive_arr[13],pactive_arr[14], \
                                          pactive_arr[15],pactive_arr[16]]), \
                    np.sum([pactive_arr[1],pactive_arr[2]]), \
                    np.sum([pactive_arr[3]]), \
                    np.sum([pactive_arr[4],pactive_arr[5]]), \
                    np.sum([pactive_arr[10],pactive_arr[11]])]

pactive_std_glob = [np.sqrt(np.sum([pactive_h_dn_corr[6]**2,pactive_h_dn_corr[7]**2,pactive_h_dn_corr[8]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[12]**2,pactive_h_dn_corr[13]**2,pactive_h_dn_corr[14], \
                                          pactive_h_dn_corr[15]**2,pactive_h_dn_corr[16]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[1]**2,pactive_h_dn_corr[2]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[3]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[4]**2,pactive_h_dn_corr[5]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[10]**2,pactive_h_dn_corr[11]**2]))]

##############################################################

pretrans_h_dn_corr = np.array(pretrans_arr) - np.array(pretrans_h_dn_corr)

pretrans_arr_glob = [np.sum([pretrans_arr[6],pretrans_arr[7],pretrans_arr[8]]), \
                    np.sum([pretrans_arr[12],pretrans_arr[13],pretrans_arr[14], \
                                          pretrans_arr[15],pretrans_arr[16]]), \
                    np.sum([pretrans_arr[1],pretrans_arr[2]]), \
                    np.sum([pretrans_arr[3]]), \
                    np.sum([pretrans_arr[4],pretrans_arr[5]]), \
                    np.sum([pretrans_arr[10],pretrans_arr[11]])]

pretrans_std_glob = [np.sqrt(np.sum([pretrans_h_dn_corr[6]**2,pretrans_h_dn_corr[7]**2,pretrans_h_dn_corr[8]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[12]**2,pretrans_h_dn_corr[13]**2,pretrans_h_dn_corr[14], \
                                          pretrans_h_dn_corr[15]**2,pretrans_h_dn_corr[16]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[1]**2,pretrans_h_dn_corr[2]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[3]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[4]**2,pretrans_h_dn_corr[5]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[10]**2,pretrans_h_dn_corr[11]**2]))]




#reorganizing from most to least expensive:
pnonmyc_arr_glob = [pnonmyc_arr_glob[4],pnonmyc_arr_glob[2], \
                     pnonmyc_arr_glob[1],pnonmyc_arr_glob[0], \
                     pnonmyc_arr_glob[5],pnonmyc_arr_glob[3]]

pnonmyc_std_glob = [pnonmyc_std_glob[4],pnonmyc_std_glob[2], \
                     pnonmyc_std_glob[1],pnonmyc_std_glob[0], \
                     pnonmyc_std_glob[5],pnonmyc_std_glob[3]]

pactive_arr_glob = [pactive_arr_glob[4],pactive_arr_glob[2], \
                     pactive_arr_glob[1],pactive_arr_glob[0], \
                     pactive_arr_glob[5],pactive_arr_glob[3]]

pactive_std_glob = [pactive_std_glob[4],pactive_std_glob[2], \
                     pactive_std_glob[1],pactive_std_glob[0], \
                     pactive_std_glob[5],pactive_std_glob[3]]

pretrans_arr_glob = [pretrans_arr_glob[4],pretrans_arr_glob[2], \
                     pretrans_arr_glob[1],pretrans_arr_glob[0], \
                     pretrans_arr_glob[5],pretrans_arr_glob[3]]

pretrans_std_glob = [pretrans_std_glob[4],pretrans_std_glob[2], \
                     pretrans_std_glob[1],pretrans_std_glob[0], \
                     pretrans_std_glob[5],pretrans_std_glob[3]]


ptotal_std_glob = np.sqrt(np.array(pnonmyc_std_glob)**2+np.array(pactive_std_glob)**2+(np.array(pretrans_std_glob)/1.)**2)



fig = plt.figure(figsize=(48, 24))
p1 = plt.bar(ind, pnonmyc_arr_glob, width, color = 'b',edgecolor='k',linewidth=3.3)
p2 = plt.bar(ind, pactive_arr_glob, width, color = 'r',edgecolor='k',linewidth=3.3,bottom=pnonmyc_arr_glob)
p3 = plt.bar(ind, np.array(pretrans_arr_glob)/1., width,color = 'darkorange',edgecolor='k',linewidth=3.3, bottom=pactive_arr_glob,yerr=ptotal_std_glob, error_kw=dict(lw=5, capsize =5, capthick=3))
plt.ylabel(r'Carbon use (g C m$^{-2}$ yr$^{-1}$)')
plt.xticks(ind,('Evergreen\n broadleaf','Evergreen\n needleleaf','Grassland','Decidous\n broadleaf','Shrubland','Deciduos\n needleleaf'))
plt.legend((p1[0],p2[0],p3[0]),(r'Direct root uptake','Mycorrhizal uptake','Retranslocation '))
#($\div$ 100)
plt.ticklabel_format(axis="y",style="sci",scilimits=(0,0))
plt.savefig('figures/COST_PUPTAKE_zonal_pft.png',bbox_inches="tight")
#plt.show()


sys.exit()

