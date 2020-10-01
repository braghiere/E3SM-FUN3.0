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

pftname = [ "not_vegetated",
  "needleleaf_evergreen_temperate_tree",
  "needleleaf_evergreen_boreal_tree",
  "needleleaf_deciduous_boreal_tree",
  "broadleaf_evergreen_tropical_tree",
  "broadleaf_evergreen_temperate_tree",
  "broadleaf_deciduous_tropical_tree",
  "broadleaf_deciduous_temperate_tree",
  "broadleaf_deciduous_boreal_tree",
  "broadleaf_evergreen_shrub",
  "broadleaf_deciduous_temperate_shrub",
  "broadleaf_deciduous_boreal_shrub",
  "c3_arctic_grass",
  "c3_non-arctic_grass",
  "c4_grass",
  "c3_crop",
  "c3_irrigated                            ",
  "corn                                    ",
  "irrigated_corn                          ",
  "spring_temperate_cereal                 ",
  "irrigated_spring_temperate_cereal       ",
  "winter_temperate_cereal                 ",
  "irrigated_winter_temperate_cereal       ",
  "soybean                                 ",
  "irrigated_soybean                       "]

pftname_fig = [ "not \n vegetated",
  "needleleaf_evergreen \n temperate",
  "needleleaf_evergreen \n boreal",
  "needleleaf_deciduous \n boreal",
  "broadleaf_evergreen \n tropical",
  "broadleaf_evergreen \n temperate",
  "broadleaf_deciduous \n tropical",
  "broadleaf_deciduous \n temperate",
  "broadleaf_deciduous \n boreal",
  "evergreen_shrub",
  "deciduous_temperate \n shrub",
  "deciduous_boreal \n shrub",
  "c3_arctic_grass",
  "c3_non-arctic \n grass",
  "c4_grass",
  "c3_crop",
  "c3_irrigated"]



#f = Dataset('files/time_average_elm_pft.nc','r')
f = Dataset('files/fix_global_v7_funp_pft_year.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
pft = f.variables['pft'][:]


#unit = puptake.units
#long_name = puptake.long_name


#pnonmyc_ori = f.variables['COST_PNONMYC'][:]
#pactive_ori = f.variables['COST_PACTIVE'][:]
#pretrans_ori = f.variables['COST_PRETRANS'][:]

pnonmyc_ori = f.variables['COST_PNONMYC'][:]
pactive_ori = f.variables['COST_PACTIVE'][:]
pretrans_ori = f.variables['COST_PRETRANS'][:]

#pnonmyc_ori = f.variables['NNONMYC'][:]/f.variables['COST_NNONMYC'][:]
#pactive_ori = f.variables['NACTIVE'][:]/f.variables['COST_NACTIVE'][:]
#pretrans_ori = f.variables['NRETRANS'][:]/f.variables['COST_NRETRANS'][:]


pnonmyc_ori,lon = shiftgrid(180., pnonmyc_ori, lon, start=False)
lon = f.variables['lon'][:]
pactive_ori,lon = shiftgrid(180., pactive_ori, lon, start=False)
lon = f.variables['lon'][:]
pretrans_ori,lon = shiftgrid(180., pretrans_ori, lon, start=False)



pactive_ori[pactive_ori >= 1.e+2]=np.nan
#pactive_ori[pactive_ori == 1]=np.nan
pactive_ori[pactive_ori == 0]=np.nan

pnonmyc_ori[pnonmyc_ori >= 1.e+2]=np.nan
#pnonmyc_ori[pnonmyc_ori == 1]=np.nan
pnonmyc_ori[pnonmyc_ori == 0]=np.nan

pretrans_ori[pretrans_ori >= 1.e+2]=np.nan
#pretrans_ori[pretrans_ori == 1]=np.nan
pretrans_ori[pretrans_ori == 0]=np.nan


#for i in xrange(16):
#   plt.imshow(pnonmyc_ori[i,:,:])
#   plt.show()
#sys.exit()


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

npp_pactive_tot = [] 
npp_pnonmyc_tot = [] 
npp_pretrans_tot = [] 

for i in xrange(len(pft)):

   print('pft=',pft[i])
   #PACTIVE
   pactive = np.array(pactive_ori[i,:,:])
   #pactive[pactive==1e+36]=np.nan

   pactive=np.ma.masked_array(pactive,np.isnan(pactive))
   average=np.ma.average(pactive,axis=0,weights=weights)
   variance=np.ma.average((pactive-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pactive,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pactive))))

   print 'Global pactive =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 
   npp_pactive_tot.append(mean_corr)
   #PNONMYC
   pnonmyc = np.array(pnonmyc_ori[i,:,:])
   #pnonmyc[pnonmyc==1e+36]=np.nan

   pnonmyc=np.ma.masked_array(pnonmyc,np.isnan(pnonmyc))
   average=np.ma.average(pnonmyc,axis=0,weights=weights)
   variance=np.ma.average((pnonmyc-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pnonmyc,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pnonmyc))))

   print 'Global pnonmyc =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 
   npp_pnonmyc_tot.append(mean_corr)

   #PRETRANS
   pretrans = np.array(pretrans_ori[i,:,:])
   #pretrans[pretrans==1e+36]=np.nan

   pretrans=np.ma.masked_array(pretrans,np.isnan(pretrans))
   average=np.ma.average(pretrans,axis=0,weights=weights)
   variance=np.ma.average((pretrans-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pretrans,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pretrans))))

   print 'Global pretrans =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

   npp_pretrans_tot.append(mean_corr)


print('Myc',np.sum(npp_pactive_tot))
print('Non-Myc',np.sum(npp_pnonmyc_tot))
print('Retrans',np.sum(npp_pretrans_tot))
print('Total',np.sum(npp_pretrans_tot+npp_pnonmyc_tot+npp_pactive_tot))
#sys.exit()

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


#const = 365*24*60*60
const = 1.

pactive = np.array(pactive_ori[:,:,:]*const)
pactive[pactive==0]=np.nan
pactive[pactive==1]=np.nan
pactive=np.ma.masked_array(pactive,np.isnan(pactive))

pnonmyc= np.array(pnonmyc_ori[:,:,:]*const)
pnonmyc[pnonmyc==0]=np.nan
pnonmyc[pnonmyc==1]=np.nan
pnonmyc=np.ma.masked_array(pnonmyc,np.isnan(pnonmyc))

pretrans= np.array(pretrans_ori[:,:,:]*const)
pretrans[pretrans==0]=np.nan
pretrans[pretrans==1]=np.nan
pretrans=np.ma.masked_array(pretrans,np.isnan(pretrans))

#for i in xrange(16):
#   plt.imshow(pactive[i,:,:])
#   plt.show()
#sys.exit()

pretrans_mean = np.nanmean(np.nanmean(pretrans,axis=1),axis=1)
pnonmyc_mean = np.nanmean(np.nanmean(pnonmyc,axis=1),axis=1)
pactive_mean = np.nanmean(np.nanmean(pactive,axis=1),axis=1)

pretrans_std = np.nanstd(np.nanstd(pretrans,axis=1),axis=1)
pnonmyc_std = np.nanstd(np.nanstd(pnonmyc,axis=1),axis=1)
pactive_std = np.nanstd(np.nanstd(pactive,axis=1),axis=1)

total_std = np.sqrt(pretrans_std**2 + pnonmyc_std**2 + pactive_std**2)

#plt.plot(pretrans_mean,label='pretrans')
#plt.plot(pactive_mean,label='pactive')
#plt.plot(pnonmyc_mean,label='pnonmyc')
#plt.legend()



pftname_fig_v1 = [ "broadleaf \n deciduous",
  "grassland",
  "needleleaf \n evergreen",
  "needleleaf \n deciduous",
  "broadleaf \n evergreen",
  "shrubland"]


#pretrans_mean = [np.sum([pretrans_mean[6],pretrans_mean[7],pretrans_mean[8]]), \
                    #np.sum([pretrans_mean[12],pretrans_mean[13],pretrans_mean[14], \
#                                          pretrans_mean[15]]), \
##                    np.sum([pretrans_mean[1],pretrans_mean[2]]), \
#                   np.sum([pretrans_mean[3]]), \
#                    np.sum([pretrans_mean[4],pretrans_mean[5]]), \
#                    np.sum([pretrans_mean[10],pretrans_mean[11]])]

#pnonmyc_mean = [np.sum([pnonmyc_mean[6],pnonmyc_mean[7],pnonmyc_mean[8]]), \
                    #np.sum([pnonmyc_mean[12],pnonmyc_mean[13],pnonmyc_mean[14], \
#                                          pnonmyc_mean[15]]), \
#                    np.sum([pnonmyc_mean[1],pnonmyc_mean[2]]), \
#                    np.sum([pnonmyc_mean[3]]), \
#                    np.sum([pnonmyc_mean[4],pnonmyc_mean[5]]), \
#                    np.sum([pnonmyc_mean[10],pnonmyc_mean[11]])]

#pactive_mean = [np.sum([pactive_mean[6],pactive_mean[7],pactive_mean[8]]), \
                    #np.sum([pactive_mean[12],pactive_mean[13],pactive_mean[14], \
#                                          pactive_mean[15]]), \
#                    np.sum([pactive_mean[1],pactive_mean[2]]), \
#                    np.sum([pactive_mean[3]]), \
#                    np.sum([pactive_mean[4],pactive_mean[5]]), \
#                    np.sum([pactive_mean[10],pactive_mean[11]])]


ind = np.arange(17)
width = 0.35

print(np.shape(pactive_mean))

fig = plt.figure(figsize=(48, 24))

p1 = plt.bar(ind, np.array(pactive_mean), width, color = 'b',edgecolor='k',linewidth=3.3)

p2 = plt.bar(ind, np.array(pnonmyc_mean), width, color = 'r',edgecolor='k',linewidth=3.3, bottom=np.array(pactive_mean))


p3 = plt.bar(ind, np.array(pretrans_mean), width,color = 'darkorange',edgecolor='k',linewidth=3.3,bottom=np.array(pactive_mean)+np.array(pnonmyc_mean))



plt.ylabel(r'Carbon cost (g P g C $^{-1}$)')
#plt.xticks(ind,('Evergreen\n broadleaf','Evergreen\n needleleaf','Decidous\n broadleaf','Grassland','Shrubland','Deciduos\n needleleaf'))
plt.xticks(ind,pftname_fig,rotation=45)
plt.legend((p1[0],p2[0],p3[0]),('Direct root uptake','Mycorrhizal uptake','Retranslocation'))
#($\div$ 100)
plt.ticklabel_format(axis="y",style="sci",scilimits=(0,0))
plt.savefig('figures/10_test_b.png')
plt.close()
#sys.exit()

for i in xrange(len(pft)):

 #PNONMYC
 pnonmyc_arr.append(np.mean(np.ma.average(pnonmyc[i,:,:], \
                    axis=0,weights=weights)))

 pnonmyc_ave.append(np.ma.average(pnonmyc[i,:,:], \
                    axis=0,weights=weights))

 pnonmyc_var.append(np.ma.average((pnonmyc[i,:,:]-pnonmyc_ave[i])**2,axis=0,weights=weights))

 pnonmyc_std_corr.append(np.mean(np.sqrt(pnonmyc_var[i])))

 pnonmyc_h_dn_corr.append(stats.norm.interval(0.66, loc=pnonmyc_arr[i], scale= pnonmyc_std_corr[i]/np.sqrt(len(np.ma.array(pnonmyc[i,:,:]))))[0])
 
 #PACTIVE

 pactive_arr.append(np.mean(np.ma.average(pactive[i,:,:], \
                    axis=0,weights=weights)))

 pactive_ave.append(np.ma.average(pactive[i,:,:], \
                    axis=0,weights=weights))

 pactive_var.append(np.ma.average((pactive[i,:,:]-pactive_ave[i])**2,axis=0,weights=weights))

 pactive_std_corr.append(np.mean(np.sqrt(pactive_var[i])))

 pactive_h_dn_corr.append(stats.norm.interval(0.66, loc=pactive_arr[i], scale= pactive_std_corr[i]/np.sqrt(len(np.ma.array(pactive[i,:,:]))))[0])

 #PRETRANS
 pretrans_arr.append(np.mean(np.ma.average(pretrans[i,:,:], \
                     axis=0,weights=weights)))

 pretrans_ave.append(np.ma.average(pretrans[i,:,:], \
                    axis=0,weights=weights))

 pretrans_var.append(np.ma.average((pretrans[i,:,:]-pretrans_ave[i])**2,axis=0,weights=weights))

 pretrans_std_corr.append(np.mean(np.sqrt(pretrans_var[i])))

 pretrans_h_dn_corr.append(stats.norm.interval(0.66, loc=pretrans_arr[i], scale= pretrans_std_corr[i]/np.sqrt(len(np.ma.array(pretrans[i,:,:]))))[0])

ind = np.arange(6)
width = 0.35


pnonmyc_h_dn_corr = np.array(pnonmyc_arr) - np.array(pnonmyc_h_dn_corr)

pnonmyc_arr_glob = [np.sum([pnonmyc_arr[6],pnonmyc_arr[7],pnonmyc_arr[8]]), \
                    np.sum([pnonmyc_arr[12],pnonmyc_arr[13], \
                                          pnonmyc_arr[15]]), \
                    np.sum([pnonmyc_arr[1],pnonmyc_arr[2]]), \
                    np.sum([pnonmyc_arr[3]]), \
                    np.sum([pnonmyc_arr[4],pnonmyc_arr[5]]), \
                    np.sum([pnonmyc_arr[10],pnonmyc_arr[11]])]

pnonmyc_std_glob = [np.sqrt(np.sum([pnonmyc_h_dn_corr[6]**2,pnonmyc_h_dn_corr[7]**2,pnonmyc_h_dn_corr[8]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[12]**2,pnonmyc_h_dn_corr[13]**2, \
                                          pnonmyc_h_dn_corr[15]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[1]**2,pnonmyc_h_dn_corr[2]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[3]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[4]**2,pnonmyc_h_dn_corr[5]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[10]**2,pnonmyc_h_dn_corr[11]**2]))]



##############################################################

pactive_h_dn_corr = np.array(pactive_arr) - np.array(pactive_h_dn_corr)


pactive_arr_glob = [np.sum([pactive_arr[6],pactive_arr[7],pactive_arr[8]]), \
                    np.sum([pactive_arr[12],pactive_arr[13], \
                                          pactive_arr[15]]), \
                    np.sum([pactive_arr[1],pactive_arr[2]]), \
                    np.sum([pactive_arr[3]]), \
                    np.sum([pactive_arr[4],pactive_arr[5]]), \
                    np.sum([pactive_arr[10],pactive_arr[11]])]

pactive_std_glob = [np.sqrt(np.sum([pactive_h_dn_corr[6]**2,pactive_h_dn_corr[7]**2,pactive_h_dn_corr[8]**2])), \
                    np.sqrt(np.sum([pnonmyc_h_dn_corr[12]**2,pactive_h_dn_corr[13]**2, \
                                          pactive_h_dn_corr[15]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[1]**2,pactive_h_dn_corr[2]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[3]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[4]**2,pactive_h_dn_corr[5]**2])), \
                    np.sqrt(np.sum([pactive_h_dn_corr[10]**2,pactive_h_dn_corr[11]**2]))]

##############################################################

pretrans_h_dn_corr = np.array(pretrans_arr) - np.array(pretrans_h_dn_corr)

pretrans_arr_glob = [np.sum([pretrans_arr[6],pretrans_arr[7],pretrans_arr[8]]), \
                    np.sum([pretrans_arr[12],pretrans_arr[13], \
                                          pretrans_arr[15]]), \
                    np.sum([pretrans_arr[1],pretrans_arr[2]]), \
                    np.sum([pretrans_arr[3]]), \
                    np.sum([pretrans_arr[4],pretrans_arr[5]]), \
                    np.sum([pretrans_arr[10],pretrans_arr[11]])]

pretrans_std_glob = [np.sqrt(np.sum([pretrans_h_dn_corr[6]**2,pretrans_h_dn_corr[7]**2,pretrans_h_dn_corr[8]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[12]**2,pretrans_h_dn_corr[13]**2, \
                                          pretrans_h_dn_corr[15]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[1]**2,pretrans_h_dn_corr[2]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[3]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[4]**2,pretrans_h_dn_corr[5]**2])), \
                    np.sqrt(np.sum([pretrans_h_dn_corr[10]**2,pretrans_h_dn_corr[11]**2]))]





#reorganizing from most to least expensive:
pnonmyc_arr_glob = [pnonmyc_arr_glob[1],pnonmyc_arr_glob[5], \
                     pnonmyc_arr_glob[0],pnonmyc_arr_glob[4], \
                     pnonmyc_arr_glob[2],pnonmyc_arr_glob[3]]

pnonmyc_std_glob = [pnonmyc_std_glob[1],pnonmyc_std_glob[5], \
                     pnonmyc_std_glob[0],pnonmyc_std_glob[4], \
                     pnonmyc_std_glob[2],pnonmyc_std_glob[3]]

pactive_arr_glob = [pactive_arr_glob[1],pactive_arr_glob[5], \
                     pactive_arr_glob[0],pactive_arr_glob[4], \
                     pactive_arr_glob[2],pactive_arr_glob[3]]

pactive_std_glob = [pactive_std_glob[1],pactive_std_glob[5], \
                     pactive_std_glob[0],pactive_std_glob[4], \
                     pactive_std_glob[2],pactive_std_glob[3]]

pretrans_arr_glob = [pretrans_arr_glob[1],pretrans_arr_glob[5], \
                     pretrans_arr_glob[0],pretrans_arr_glob[4], \
                     pretrans_arr_glob[2],pretrans_arr_glob[3]]

pretrans_std_glob = [pretrans_std_glob[1],pretrans_std_glob[5], \
                     pretrans_std_glob[0],pretrans_std_glob[4], \
                     pretrans_std_glob[2],pretrans_std_glob[3]]


ptotal_std_glob = np.sqrt(np.array(pnonmyc_std_glob)**2+np.array(pactive_std_glob)**2+(np.array(pretrans_std_glob))**2)

print('Grassland','Shrubland','Decidous\n broadleaf','Evergreen\n broadleaf','Evergreen\n needleleaf','Deciduos\n needleleaf')
tot = np.array(pretrans_arr_glob) + np.array(pnonmyc_arr_glob) + np.array(pactive_arr_glob)
print('ptotal',tot,ptotal_std_glob)
print('pretrans',pretrans_arr_glob,pretrans_arr_glob/tot)
print('pnonmyc',pnonmyc_arr_glob,pnonmyc_arr_glob/tot)
print('pactive',pactive_arr_glob,pactive_arr_glob/tot)


#sys.exit()

fig = plt.figure(figsize=(48, 24))

#const = 365*24*60*60
const = 1.

p1 = plt.bar(ind, np.array(pactive_arr_glob), width, color = 'b',edgecolor='k',linewidth=3.3)

p2 = plt.bar(ind, np.array(pnonmyc_arr_glob), width, color = 'r',edgecolor='k',linewidth=3.3, bottom=np.array(pactive_arr_glob))


p3 = plt.bar(ind, np.array(pretrans_arr_glob), width,color = 'darkorange',edgecolor='k',linewidth=3.3,yerr=ptotal_std_glob, error_kw=dict(lw=5, capsize =5, capthick=3), bottom=np.array(pactive_arr_glob)+np.array(pnonmyc_arr_glob))



plt.ylabel(r'Carbon cost (g P g C$^{-1}$)')
plt.xticks(ind,('Grassland','Shrubland','Decidous\n broadleaf','Evergreen\n broadleaf','Evergreen\n needleleaf','Deciduos\n needleleaf'), rotation=45)
plt.legend((p1[0],p2[0],p3[0]),('Mycorrhizal uptake','Direct root uptake','Retranslocation'))
#($\div$ 100)
plt.ticklabel_format(axis="y",style="sci",scilimits=(0,0))
plt.savefig('figures/COST_PUPTAKE_zonal_pft_b.png',bbox_inches="tight")
#plt.show()


sys.exit()

