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



#f = Dataset('files/time_average_elm_pft.nc','r')
f = Dataset('files/fix_global_v6_fun_pft_year.nc','r')

lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
pft = f.variables['pft'][:]


#unit = puptake.units
#long_name = puptake.long_name


pnonmyc_ori = f.variables['COST_NNONMYC'][:]
pactive_ori = f.variables['COST_NACTIVE'][:]
pretrans_ori = f.variables['COST_NRETRANS'][:]
pfix_ori = f.variables['COST_NFIX'][:]

pnonmyc_ori,lon = shiftgrid(180., pnonmyc_ori, lon, start=False)
lon = f.variables['lon'][:]
pactive_ori,lon = shiftgrid(180., pactive_ori, lon, start=False)
lon = f.variables['lon'][:]
pretrans_ori,lon = shiftgrid(180., pretrans_ori, lon, start=False)
lon = f.variables['lon'][:]
pfix_ori,lon = shiftgrid(180., pfix_ori, lon, start=False)

pactive_ori[pactive_ori > 1e30]=np.nan
pnonmyc_ori[pnonmyc_ori > 1e30]=np.nan
pretrans_ori[pretrans_ori > 1e30]=np.nan
pfix_ori[pfix_ori > 1e30]=np.nan


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



for i in xrange(len(pft)):

   print('pft=',pft[i])
   #PACTIVE
   pactive = np.array(pactive_ori[i,:,:])
   pactive[pactive==1e+36]=np.nan

   pactive=np.ma.masked_array(pactive,np.isnan(pactive))
   average=np.ma.average(pactive,axis=0,weights=weights)
   variance=np.ma.average((pactive-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pactive,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pactive))))

   print 'Global pactive =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 
   #PNONMYC
   pnonmyc = np.array(pnonmyc_ori[i,:,:])
   pnonmyc[pnonmyc==1e+36]=np.nan

   pnonmyc=np.ma.masked_array(pnonmyc,np.isnan(pnonmyc))
   average=np.ma.average(pnonmyc,axis=0,weights=weights)
   variance=np.ma.average((pnonmyc-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pnonmyc,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pnonmyc))))

   print 'Global pnonmyc =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'
 

   #PRETRANS
   pretrans = np.array(pretrans_ori[i,:,:])
   pretrans[pretrans==1e+36]=np.nan

   pretrans=np.ma.masked_array(pretrans,np.isnan(pretrans))
   average=np.ma.average(pretrans,axis=0,weights=weights)
   variance=np.ma.average((pretrans-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pretrans,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pretrans))))

   print 'Global pretrans =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

   #PFIX
   pfix = np.array(pfix_ori[i,:,:])
   pfix[pretrans==1e+36]=np.nan

   pfix=np.ma.masked_array(pfix,np.isnan(pfix))
   average=np.ma.average(pfix,axis=0,weights=weights)
   variance=np.ma.average((pfix-average)**2,axis=0,weights=weights)

   mean_corr = np.mean(np.ma.average(pfix,axis=0,weights=weights))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

   h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(pfix))))

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

#PFIX
pfix_arr = []
pfix_ave = []
pfix_var = []
pfix_std_corr = []
pfix_h_dn_corr = []


pactive = np.array(pactive_ori[:,:,:])
pactive[pactive==1e+36]=np.nan
pactive=np.ma.masked_array(pactive,np.isnan(pactive))

pnonmyc= np.array(pnonmyc_ori[:,:,:])
pnonmyc[pactive==1e+36]=np.nan
pnonmyc=np.ma.masked_array(pnonmyc,np.isnan(pnonmyc))

pretrans= np.array(pretrans_ori[:,:,:])
pretrans[pretrans==1e+36]=np.nan
pretrans=np.ma.masked_array(pretrans,np.isnan(pretrans))

pfix= np.array(pfix_ori[:,:,:])
pfix[pretrans==1e+36]=np.nan
pfix=np.ma.masked_array(pfix,np.isnan(pfix))

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

 #PFIX
 pfix_arr.append(np.mean(np.ma.average(pfix[i,:,:], \
                     axis=0,weights=weights)))

 pfix_ave.append(np.ma.average(pfix[i,:,:], \
                    axis=0,weights=weights))

 pfix_var.append(np.ma.average((pfix[i,:,:]-pfix_ave[i])**2,axis=0,weights=weights))

 pfix_std_corr.append(np.mean(np.sqrt(pfix_var[i])))

 pfix_h_dn_corr.append(stats.norm.interval(0.66, loc=pfix_arr[i], scale= pfix_std_corr[i]/np.sqrt(len(np.ma.array(pfix[i,:,:]))))[0])




ind = np.arange(6)
width = 0.35

pnonmyc_h_dn_corr = np.array(pnonmyc_arr) - np.array(pnonmyc_h_dn_corr)

pnonmyc_arr_glob = [np.sum([pnonmyc_arr[6],pnonmyc_arr[7],pnonmyc_arr[8]]), \
                    np.sum([pnonmyc_arr[12],pnonmyc_arr[13],pnonmyc_arr[14], \
                                          pnonmyc_arr[15],pnonmyc_arr[16]]), \
                    np.sum([pnonmyc_arr[1],pnonmyc_arr[2]]), \
                    np.sum([pnonmyc_arr[3]]), \
                    np.sum([pnonmyc_arr[4],pnonmyc_arr[5]]), \
                    np.sum([pnonmyc_arr[10],pnonmyc_arr[11]])]

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

##############################################################

pfix_h_dn_corr = np.array(pfix_arr) - np.array(pfix_h_dn_corr)

pfix_arr_glob = [np.sum([pfix_arr[6],pfix_arr[7],pfix_arr[8]]), \
                    np.sum([pfix_arr[12],pfix_arr[13],pfix_arr[14], \
                                          pfix_arr[15],pfix_arr[16]]), \
                    np.sum([pfix_arr[1],pfix_arr[2]]), \
                    np.sum([pfix_arr[3]]), \
                    np.sum([pfix_arr[4],pfix_arr[5]]), \
                    np.sum([pfix_arr[10],pfix_arr[11]])]

pfix_std_glob = [np.sqrt(np.sum([pfix_h_dn_corr[6]**2,pfix_h_dn_corr[7]**2,pfix_h_dn_corr[8]**2])), \
                    np.sqrt(np.sum([pfix_h_dn_corr[12]**2,pfix_h_dn_corr[13]**2,pfix_h_dn_corr[14], \
                                          pfix_h_dn_corr[15]**2,pfix_h_dn_corr[16]**2])), \
                    np.sqrt(np.sum([pfix_h_dn_corr[1]**2,pfix_h_dn_corr[2]**2])), \
                    np.sqrt(np.sum([pfix_h_dn_corr[3]**2])), \
                    np.sqrt(np.sum([pfix_h_dn_corr[4]**2,pfix_h_dn_corr[5]**2])), \
                    np.sqrt(np.sum([pfix_h_dn_corr[10]**2,pfix_h_dn_corr[11]**2]))]



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

pfix_arr_glob = [pfix_arr_glob[4],pfix_arr_glob[2], \
                     pfix_arr_glob[1],pfix_arr_glob[0], \
                     pfix_arr_glob[5],pfix_arr_glob[3]]

pfix_std_glob = [pfix_std_glob[4],pfix_std_glob[2], \
                     pfix_std_glob[1],pfix_std_glob[0], \
                     pfix_std_glob[5],pfix_std_glob[3]]

ptotal_std_glob = np.sqrt(np.array(pnonmyc_std_glob)**2+np.array(pactive_std_glob)**2+(np.array(pretrans_std_glob)/1.)**2 +(np.array(pfix_std_glob)/1.)**2)



fig = plt.figure(figsize=(48, 24))
p1 = plt.bar(ind, pnonmyc_arr_glob, width, color = 'b',edgecolor='k',linewidth=3.3)
p2 = plt.bar(ind, pactive_arr_glob, width, color = 'r',edgecolor='k',linewidth=3.3,bottom=pnonmyc_arr_glob)
p3 = plt.bar(ind, np.array(pretrans_arr_glob)/1., width,color = 'darkorange',edgecolor='k',linewidth=3.3, bottom=pactive_arr_glob)
p4 = plt.bar(ind, np.array(pfix_arr_glob)/1., width,color = 'orchid',edgecolor='k',linewidth=3.3, bottom=pretrans_arr_glob,yerr=ptotal_std_glob, error_kw=dict(lw=5, capsize =5, capthick=3))
plt.ylabel(r'Carbon use (g C m$^{-2}$ yr$^{-1}$)')
plt.xticks(ind,('Evergreen\n broadleaf','Evergreen\n needleleaf','Grassland','Decidous\n broadleaf','Shrubland','Deciduos\n needleleaf'))
plt.legend((p1[0],p2[0],p3[0],p4[0]),(r'Direct root uptake','Mycorrhizal uptake','Retranslocation','Symbiotic fixation'))
plt.ticklabel_format(axis="y",style="sci",scilimits=(0,0))
plt.savefig('figures/Nitrogen_fun/COST_NUPTAKE_zonal_pft.png',bbox_inches="tight")
#plt.show()


sys.exit()

