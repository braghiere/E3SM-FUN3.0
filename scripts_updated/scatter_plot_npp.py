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

from sklearn.linear_model import LinearRegression

from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

model = LinearRegression(fit_intercept=True)


e3sm_fun = Dataset('/home/renato/ELM_FUNP/fix_global_v6_1994_2005/files/fix_global_v6_fun_time_avg.nc','r')
e3sm_funp = Dataset('/home/renato/ELM_FUNP/fix_global_v6_1994_2005/files/fix_global_v6_funp_time_avg.nc','r')
e3sm = Dataset('/home/renato/ELM_FUNP/fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')



modis = Dataset('/home/renato/ELM_FUNP/npp_modis_reg_v1.nc','r')
igbp = Dataset('/home/renato/ELM_FUNP/npp_igbp_1p9x2p5.nc','r')


######################################E3SM #########################
lon = e3sm.variables['lon'][:]
npp_e3sm = e3sm.variables['NPP'][:]

npp_e3sm,lon = shiftgrid(180., npp_e3sm, lon, start=False)

npp_e3sm[npp_e3sm==1e+36]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

#npp_mpi=np.mean(npp_igbp[(4*12)-1:(15*12)-1,:,:],axis=0)



lat_e3sm = e3sm.variables['lat'][:]
#lon = igbp.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_e3sm = 1.9
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_e3sm/dtor
lat_u = (lat_e3sm + res_lat_e3sm)/dtor
weights_e3sm = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
npp_e3sm = np.array(npp_e3sm)
npp_e3sm[npp_e3sm==1e+36]=np.nan

npp_e3sm=np.ma.masked_array(npp_e3sm,np.isnan(npp_e3sm))
average=np.ma.average(npp_e3sm,axis=0,weights=weights_e3sm)
variance=np.ma.average((npp_e3sm-average)**2,axis=0,weights=weights_e3sm)


mean_corr = np.mean(np.ma.average(npp_e3sm,axis=0,weights=weights_e3sm))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_e3sm))))

print(np.shape(npp_e3sm))
print 'E3SM - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

######################################E3SM-FUNP #########################
lon = e3sm_funp.variables['lon'][:]
npp_e3sm_funp = e3sm_funp.variables['NPP'][:]

npp_e3sm_funp,lon = shiftgrid(180., npp_e3sm_funp, lon, start=False)

npp_e3sm_funp[npp_e3sm_funp==1e+36]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

#npp_mpi=np.mean(npp_igbp[(4*12)-1:(15*12)-1,:,:],axis=0)



lat_e3sm_funp = e3sm_funp.variables['lat'][:]
#lon = igbp.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_e3sm_funp = 1.9
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_e3sm_funp/dtor
lat_u = (lat_e3sm_funp + res_lat_e3sm_funp)/dtor
weights_e3sm_funp = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
npp_e3sm_funp = np.array(npp_e3sm_funp)
npp_e3sm_funp[npp_e3sm_funp==1e+36]=np.nan

npp_e3sm_funp=np.ma.masked_array(npp_e3sm_funp,np.isnan(npp_e3sm_funp))
average=np.ma.average(npp_e3sm_funp,axis=0,weights=weights_e3sm_funp)
variance=np.ma.average((npp_e3sm_funp-average)**2,axis=0,weights=weights_e3sm_funp)


mean_corr = np.mean(np.ma.average(npp_e3sm_funp,axis=0,weights=weights_e3sm_funp))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_e3sm_funp))))

print(np.shape(npp_e3sm_funp))
print 'E3SM-FUNP - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

######################################E3SM-FUN #########################
lon = e3sm_fun.variables['lon'][:]
npp_e3sm_fun = e3sm_fun.variables['NPP'][:]

npp_e3sm_fun,lon = shiftgrid(180., npp_e3sm_fun, lon, start=False)

npp_e3sm_fun[npp_e3sm_fun==1e+36]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

#npp_mpi=np.mean(npp_igbp[(4*12)-1:(15*12)-1,:,:],axis=0)



lat_e3sm_fun = e3sm_fun.variables['lat'][:]
#lon = igbp.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_e3sm_fun = 1.9
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_e3sm_fun/dtor
lat_u = (lat_e3sm_fun + res_lat_e3sm_fun)/dtor
weights_e3sm_fun = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
npp_e3sm_fun = np.array(npp_e3sm_fun)
npp_e3sm_fun[npp_e3sm_fun==1e+36]=np.nan

npp_e3sm_fun=np.ma.masked_array(npp_e3sm_fun,np.isnan(npp_e3sm_fun))
average=np.ma.average(npp_e3sm_fun,axis=0,weights=weights_e3sm_fun)
variance=np.ma.average((npp_e3sm_fun-average)**2,axis=0,weights=weights_e3sm_fun)


mean_corr = np.mean(np.ma.average(npp_e3sm_fun,axis=0,weights=weights_e3sm_fun))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*(365*24*60*60)*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_e3sm_fun))))

print(np.shape(npp_e3sm_fun))
print 'E3SM-FUN - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


######################################MODIS #########################
lon = modis.variables['lon'][:]
npp_modis = modis.variables['npp'][:]

#npp_modis,lon = shiftgrid(180., npp_modis, lon, start=False)

#npp_modis[npp_modis==np.nan]=0.0
npp_modis=npp_modis.filled()
npp_modis[npp_modis==1e+20]=0.0

#npp_mpi=np.mean(npp_modis[(4*12)-1:(15*12)-1,:,:],axis=0)

#plt.imshow(npp_modis)
#plt.show()

lat_modis = modis.variables['lat'][:]

x,y = np.meshgrid(lon, lat_modis) 

npp_modis = maskoceans(x, y, npp_modis,resolution='l',grid=1.25,inlands=True)

#plt.imshow(npp_modis)
#plt.show()

#sys.exit()

lon = modis.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_modis = 1.875
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_modis/dtor
lat_u = (lat_modis + res_lat_modis)/dtor
weights_modis = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
#npp_modis = np.array(npp_modis)
#npp_modis[npp_modis==0.0]=np.nan

#npp_modis=np.ma.masked_array(npp_modis,np.isnan(npp_modis))
average=np.ma.average(npp_modis,axis=0,weights=weights_modis)
variance=np.ma.average((npp_modis-average)**2,axis=0,weights=weights_modis)


mean_corr = np.mean(np.ma.average(npp_modis,axis=0,weights=weights_modis))*365*24*60*60*1000*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*365*24*60*60*1000*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_modis))))

print(np.shape(npp_modis))
print 'MODIS - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


######################################IGBP #########################
lon = igbp.variables['lon'][:]
npp_igbp = igbp.variables['npp'][:]

#npp_modis,lon = shiftgrid(180., npp_modis, lon, start=False)
npp_igbp[npp_igbp==-99.99]=np.nan

#npp_igbp=npp_igbp.filled()
#npp_igbp[npp_igbp==1e+20]=0.0

#npp_mpi=np.mean(npp_modis[(4*12)-1:(15*12)-1,:,:],axis=0)



lat_igbp = igbp.variables['lat'][:]


plt.imshow(npp_modis)
plt.show()

#sys.exit()


#Function to calculate the weights of latitude
radius = 6367449
res_lat_igbp = 1.875
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_igbp/dtor
lat_u = (lat_igbp + res_lat_igbp)/dtor
weights_igbp = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
#npp_igbp = np.array(npp_igbp)
#npp_modis[npp_modis==-99.99]=np.nan

#npp_igbp=np.ma.masked_array(npp_igbp,np.isnan(npp_igbp))
average=np.ma.average(npp_igbp,axis=0,weights=weights_igbp)
variance=np.ma.average((npp_igbp-average)**2,axis=0,weights=weights_igbp)


mean_corr = np.mean(np.ma.average(npp_igbp,axis=0,weights=weights_igbp))*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))**148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_igbp))))

print(np.shape(npp_modis))
print 'IGBP - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'



gpp_mte_ori_e3sm = []
gpp_mte_ori_e3sm_fun = []
gpp_mte_ori_e3sm_funp = []


gpp_mte_ori_modis = []
gpp_mte_ori_igbp = []


gpp_mte_ori_e3sm.append(np.ma.average(npp_e3sm,axis = 1))
gpp_mte_ori_e3sm_fun.append(np.ma.average(npp_e3sm_fun,axis = 1))
gpp_mte_ori_e3sm_funp.append(np.ma.average(npp_e3sm_funp,axis = 1))


gpp_mte_ori_modis.append(np.ma.average(npp_modis,axis = 1))
gpp_mte_ori_igbp.append(np.ma.average(npp_igbp,axis = 1))


my_pfts = ['All']

const = (365*24*60*60)*1000


xfit = np.arange(0,2000,1)

print(np.size(npp_modis.flatten()))
print(npp_e3sm.flatten())

model.fit(const*npp_modis.flatten()[:, np.newaxis],(const/1000.)*npp_e3sm_fun.flatten())

print("Model slope: ", model.coef_[0])
print("Model intercept: ", model.intercept_)

r1 = model.score((const*npp_modis.flatten())[:, np.newaxis],(const/1000.)*npp_e3sm_fun.flatten())

print("R2: ", r1)

yfit1 = model.predict(xfit[:, np.newaxis])


plt.scatter(const*npp_modis.flatten(),(const/1000.)*npp_e3sm.flatten())
plt.scatter(const*npp_modis.flatten(),(const/1000.)*npp_e3sm_fun.flatten())
plt.scatter(const*npp_modis.flatten(),(const/1000.)*npp_e3sm_funp.flatten())
plt.plot(xfit,xfit,'k--')
plt.show()
#sys.exit()


for i, pft in enumerate(my_pfts):

   	fig = plt.figure(figsize=(48, 24)) 
        
        #From gCm-2s-1 to TgCyr-1

     
	plt.plot(lat_e3sm,np.array(gpp_mte_ori_e3sm[0][:]*const/1000.),label='E3SM',linewidth=9.3)
	plt.plot(lat_e3sm_fun,np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.),label='E3SM-FUN2',linewidth=9.3)
	plt.plot(lat_e3sm_funp,np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.),label='E3SM-FUN3',linewidth=9.3)

     
        plt.plot(lat_modis,np.array(gpp_mte_ori_modis[0][:]*const), label='MODIS',color='k',linestyle='--',linewidth=7.3)

        #max_val = np.max(gpp_mte_ori_mean_pnonmyc[0][:]*const)
        max_val = np.max(gpp_mte_ori_e3sm[0][:]*const/1000.)

	plt.legend(ncol=2, prop={'size': 60})
	axes = plt.gca()
	axes.set_xlim([-90.,90.])
	axes.set_ylim([0.0,max_val + 0.1*max_val])

        #plt.xticks(np.arange(min((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.), max((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'Latitude ($\degree$)')
	plt.ylabel(r'NPP (g C m$^{-2}$ yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('scatter_modis.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()
        plt.close()




        ####### Scatter plot #####
        xfit = np.arange(0,1400,1)

        model.fit(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope1 = model.coef_[0]

        r1 = model.score(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))

        r1 = r2_score(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))
        print("R2: ", r1)

        yfit1 = model.predict(xfit[:, np.newaxis])

        yfit11= model.coef_[0]*np.array(gpp_mte_ori_modis[0][:]*const) + model.intercept_

        rmse1 = np.sqrt(sum((yfit11 - np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))**2)/len(np.array(gpp_mte_ori_e3sm[0][:]*const/1000.)))
        
        rmse1 = np.sqrt(mean_squared_error(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.)))
        
        print("RMSE1: ", rmse1)


        model.fit(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)

        slope2 = model.coef_[0]

        r2 = model.score(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))
        r2 = r2_score(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))

        print("R2: ", r2)

        yfit2 = model.predict(xfit[:, np.newaxis])

        yfit22= model.coef_[0]*np.array(gpp_mte_ori_modis[0][:]*const) + model.intercept_

        rmse2 = np.sqrt(sum((yfit22 - np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))**2)/len(np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.)))
        rmse2 = np.sqrt(mean_squared_error(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.)))
        print("RMSE2: ", rmse2)

        model.fit(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope3 = model.coef_[0]

        r3 = model.score(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))
        r3 = r2_score(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))


        print("R2: ", r3)

        yfit3 = model.predict(xfit[:, np.newaxis])

        yfit33= model.coef_[0]*np.array(gpp_mte_ori_modis[0][:]*const) + model.intercept_

        rmse3 = np.sqrt(sum((yfit33 - np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))**2)/len(np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.)))
        rmse3 = np.sqrt(mean_squared_error(np.array(gpp_mte_ori_modis[0][:]*const)[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.)))
        print("RMSE3: ", rmse3)





   	fig = plt.figure(figsize=(36, 36)) 
        
        #From gCm-2s-1 to TgCyr-1
        const = (365*24*60*60)*1000
     
	plt.scatter(np.array(gpp_mte_ori_modis[0][:]*const),np.array(gpp_mte_ori_e3sm[0][:]*const/1000.),label='E3SM (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r1,rmse1,slope1),linewidth=9.3,s=600)
	plt.scatter(np.array(gpp_mte_ori_modis[0][:]*const),np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.),label='E3SM-FUN2 (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r2,rmse2,slope2),linewidth=9.3,s=600)
	plt.scatter(np.array(gpp_mte_ori_modis[0][:]*const),np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.),label='E3SM-FUN3 (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r3,rmse3,slope3),linewidth=9.3,s=600)
	plt.plot(np.arange(0,1600,1),np.arange(0,1600,1),label='1:1',color='k',linestyle='--',linewidth=9.3)
              
        plt.plot(xfit,yfit1,color='blue', \
        linewidth=9.3)
        plt.plot(xfit,yfit2,color='orange', \
        linewidth=9.3)
        plt.plot(xfit,yfit3,color='green', \
        linewidth=9.3)

     
        #plt.plot(lat_modis,np.array(gpp_mte_ori_modis[0][:]*const), label='MODIS',color='k',linestyle='--',linewidth=7.3)

        #max_val = np.max(gpp_mte_ori_mean_pnonmyc[0][:]*const)
        max_val = np.max(gpp_mte_ori_e3sm[0][:]*const/1000.)

	plt.legend(ncol=1, prop={'size': 54})
	axes = plt.gca()
	axes.set_xlim([0.0,max_val + 0.1*max_val])
	axes.set_ylim([0.0,max_val + 0.1*max_val])

        #plt.xticks(np.arange(min((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.), max((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'NPP MODIS (g C m$^{-2}$ yr$^{-1}$)')
	plt.ylabel(r'NPP model (g C m$^{-2}$ yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('scatter_modis_v4.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()
        plt.close()



xfit = np.arange(0,2000,1)

print(np.size(npp_igbp.flatten()))
print(npp_e3sm.flatten())

model.fit(npp_igbp.flatten()[:, np.newaxis],(const/1000.)*npp_e3sm_fun.flatten())

print("Model slope: ", model.coef_[0])
print("Model intercept: ", model.intercept_)

r1 = model.score((npp_igbp.flatten())[:, np.newaxis],(const/1000.)*npp_e3sm_fun.flatten())

print("R2: ", r1)

yfit1 = model.predict(xfit[:, np.newaxis])


plt.scatter(npp_igbp.flatten(),(const/1000.)*npp_e3sm.flatten())
plt.scatter(npp_igbp.flatten(),(const/1000.)*npp_e3sm_fun.flatten())
plt.scatter(npp_igbp.flatten(),(const/1000.)*npp_e3sm_funp.flatten())
plt.plot(xfit,xfit,'k--')
plt.show()
#sys.exit()


for i, pft in enumerate(my_pfts):

   	fig = plt.figure(figsize=(48, 24)) 
        
        #From gCm-2s-1 to TgCyr-1

     
	plt.plot(lat_e3sm,np.array(gpp_mte_ori_e3sm[0][:]*const/1000.),label='E3SM',linewidth=9.3)
	plt.plot(lat_e3sm_fun,np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.),label='E3SM-FUN',linewidth=9.3)
	plt.plot(lat_e3sm_funp,np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.),label='E3SM-FUNP',linewidth=9.3)

     
        plt.plot(lat_modis,np.array(gpp_mte_ori_igbp[0][:]), label='MODIS',color='k',linestyle='--',linewidth=7.3)

        #max_val = np.max(gpp_mte_ori_mean_pnonmyc[0][:]*const)
        max_val = np.max(gpp_mte_ori_e3sm[0][:]*const/1000.)

	plt.legend(ncol=2, prop={'size': 60})
	axes = plt.gca()
	axes.set_xlim([-90.,90.])
	axes.set_ylim([0.0,max_val + 0.1*max_val])

        #plt.xticks(np.arange(min((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.), max((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'Latitude ($\degree$)')
	plt.ylabel(r'NPP (g C m$^{-2}$ yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('scatter_igbp.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()
        plt.close()




        ####### Scatter plot #####
        xfit = np.arange(0,1400,1)

        model.fit(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope1 = model.coef_[0]

        r1 = model.score(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))
        r1 = r2_score(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))

        print("R2: ", r1)

        yfit1 = model.predict(xfit[:, np.newaxis])

        yfit11= model.coef_[0]*np.array(gpp_mte_ori_igbp[0][:]) + model.intercept_

        rmse1 = np.sqrt(sum((yfit11 - np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))**2)/len(np.array(gpp_mte_ori_e3sm[0][:]*const/1000.)))
        rmse1 = np.sqrt(mean_squared_error(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm[0][:]*const/1000.)))
        print("RMSE1: ", rmse1)
        

        model.fit(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope2 = model.coef_[0]

        r2 = model.score(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))
        r2 =r2_score(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))

        print("R2: ", r2)

        yfit2 = model.predict(xfit[:, np.newaxis])

        yfit22= model.coef_[0]*np.array(gpp_mte_ori_igbp[0][:]) + model.intercept_

        rmse2 = np.sqrt(sum((yfit22 - np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))**2)/len(np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.)))
        rmse2 = np.sqrt(mean_squared_error(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.)))
        print("RMSE2: ", rmse2)

        model.fit(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope3 = model.coef_[0]

        r3 = model.score(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))
        r3 = r2_score(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))

        print("R2: ", r3)

        yfit3 = model.predict(xfit[:, np.newaxis])

        yfit33= model.coef_[0]*np.array(gpp_mte_ori_igbp[0][:]) + model.intercept_

        rmse3 = np.sqrt(sum((yfit33 - np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))**2)/len(np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.)))
        rmse3 = np.sqrt(mean_squared_error(np.array(gpp_mte_ori_igbp[0][:])[:, np.newaxis],np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.)))
        print("RMSE3: ", rmse3)






   	fig = plt.figure(figsize=(36, 36)) 
        
        #From gCm-2s-1 to TgCyr-1
        const = (365*24*60*60)*1000
     
	plt.scatter(np.array(gpp_mte_ori_igbp[0][:]),np.array(gpp_mte_ori_e3sm[0][:]*const/1000.),label='E3SM (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r1,rmse1,slope1),linewidth=9.3,s=600)
	plt.scatter(np.array(gpp_mte_ori_igbp[0][:]),np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.),label='E3SM-FUN2 (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r2,rmse2,slope2),linewidth=9.3,s=600)
	plt.scatter(np.array(gpp_mte_ori_igbp[0][:]),np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.),label='E3SM-FUN3 (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r3,rmse3,slope3),linewidth=9.3,s=600)
	plt.plot(np.arange(0,1600,1),np.arange(0,1600,1),label='1:1',color='k',linestyle='--',linewidth=9.3)
              
        plt.plot(xfit,yfit1,color='blue', \
        linewidth=9.3)
        plt.plot(xfit,yfit2,color='orange', \
        linewidth=9.3)
        plt.plot(xfit,yfit3,color='green', \
        linewidth=9.3)

     
        #plt.plot(lat_modis,np.array(gpp_mte_ori_modis[0][:]*const), label='MODIS',color='k',linestyle='--',linewidth=7.3)

        #max_val = np.max(gpp_mte_ori_mean_pnonmyc[0][:]*const)
        max_val = np.max(gpp_mte_ori_e3sm[0][:]*const/1000.)

	plt.legend(ncol=1, prop={'size': 54})
	axes = plt.gca()
	axes.set_xlim([0.0,max_val + 0.1*max_val])
	axes.set_ylim([0.0,max_val + 0.1*max_val])

        #plt.xticks(np.arange(min((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.), max((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'NPP IGBP (g C m$^{-2}$ yr$^{-1}$)')
	plt.ylabel(r'NPP model (g C m$^{-2}$ yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('scatter_igbp_v4.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()
        plt.close()






sys.exit()




