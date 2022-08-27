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

from scipy.interpolate import interp1d

from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score

model = LinearRegression(fit_intercept=True)

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


#cesm2 = Dataset('/mnt/CMIP6/NPP/npp_Lmon_CESM2_historical_r1i1p1f1_gn_185001-201412.nc','r')
cesm2 = Dataset('/home/renato/ELM_FUNP/out.nc','r')
#area_cesm2 = Dataset('/mnt/CMIP6/area/areacella_fx_CESM2_historical_r1i1p1f1_gn.nc','r')
area_cesm2 = Dataset('/home/renato/ELM_FUNP/out_area.nc','r')
ukesm = Dataset('/mnt/CMIP6/NPP/UKESM1-0-LL/npp_Lmon_UKESM1-0-LL_historical_r1i1p1f2_gn_195001-201412.nc','r')
area_ukesm = Dataset('/mnt/CMIP6/area/areacella_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc','r')
#ukesm = Dataset('/mnt/CMIP6/AMIP/npp_Lmon_CESM2_amip_r1i1p1f1_gn_195001-201412.nc','r')
mpi = Dataset('/mnt/CMIP6/NPP/MPI-ESM1-2-LR/npp_Lmon_MPI-ESM1-2-LR_historical_r2i1p1f1_gn_199001-201412.nc','r')
#mpi = Dataset('/mnt/CMIP5/NPP/npp_Lmon_MPI-ESM-LR_historical_r1i1p1_185001-200512.nc','r')
area_mpi = Dataset('/mnt/CMIP6/area/areacella_fx_MPI-ESM1-2-LR_historical_r1i1p1f1_gn.nc','r')
#mpi = Dataset('/mnt/CMIP6/1pctco2/npp_Lmon_MPI-ESM1-2-LR_1pctCO2_r1i1p1f1_gn_199001-201412.nc','r')
gfdl = Dataset('/mnt/CMIP6/NPP/GFDL-ESM4/npp_Lmon_GFDL-ESM4_esm-hist_r1i1p1f1_gr1_195001-201412.nc','r')
area_gfdl = Dataset('/mnt/CMIP6/area/areacella_fx_GFDL-ESM4_historical_r1i1p1f1_gr1.nc','r')

ipsl = Dataset('/mnt/CMIP6/NPP/npp_Lmon_IPSL-CM6A-LR_historical_r1i1p1f1_gr_185001-201412.nc','r')
area_ipsl = Dataset('/mnt/CMIP6/area/areacella_fx_IPSL-CM6A-LR_historical_r1i1p1f1_gr.nc','r')

miroc = Dataset('/mnt/CMIP6/NPP/npp_Lmon_MIROC-ES2L_historical_r1i1p1f2_gn_185001-201412.nc','r')
area_miroc = Dataset('/mnt/CMIP6/area/areacella_fx_MIROC-ES2L_historical_r1i1p1f2_gn.nc','r')

access = Dataset('/mnt/CMIP6/NPP/npp_Lmon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc','r')
area_access = Dataset('/mnt/CMIP6/area/areacella_fx_ACCESS-ESM1-5_historical_r1i1p1f1_gn.nc','r')

cnrm = Dataset('/mnt/CMIP6/NPP/npp_Lmon_CNRM-ESM2-1_historical_r1i1p1f2_gr_185001-201412.nc','r')
area_cnrm = Dataset('/mnt/CMIP6/area/areacella_fx_CNRM-ESM2-1_esm-hist_r1i1p1f2_gr.nc','r')

canesm = Dataset('/mnt/CMIP6/NPP/npp_Lmon_CanESM5_esm-hist_r1i1p1f1_gn_185001-201412.nc','r')
area_canesm = Dataset('/mnt/CMIP6/area/areacella_fx_CanESM5_esm-hist_r1i1p1f1_gn.nc','r')

bcccsm = Dataset('/mnt/CMIP6/NPP/npp_Lmon_BCC-CSM2-MR_esm-hist_r1i1p1f1_gn_185001-201412.nc','r')
area_bcccsm = Dataset('/mnt/CMIP6/area/areacella_fx_BCC-CSM2-MR_hist-resIPO_r1i1p1f1_gn.nc','r')

noresm = Dataset('/mnt/CMIP6/NPP/NorESM2-LM/npp_Lmon_NorESM2-LM_historical_r1i1p1f1_gn_199001-200912.nc','r')
area_noresm = Dataset('/mnt/CMIP6/area/areacella_fx_NorESM2-LM_historical_r1i1p1f1_gn.nc','r')

e3sm_fun = Dataset('/home/renato/ELM_FUNP/fix_global_v6_1994_2005/files/fix_global_v6_fun_time_avg.nc','r')
e3sm_funp = Dataset('/home/renato/ELM_FUNP/fix_global_v6_1994_2005/files/fix_global_v6_funp_time_avg.nc','r')
e3sm = Dataset('/home/renato/ELM_FUNP/fix_global_v3_1994_2005/files/fix_global_v3_time_avg.nc','r')


igbp = Dataset('/home/renato/ISLSCP_MODEL_NPP_1027/data/model_npp_1deg/npp_igbp.nc','r')

modis = Dataset('/home/renato/ELM_FUNP/modis_npp/npp_modis.nc','r')

######################################CESM2 #########################
lon = cesm2.variables['lon'][:]
npp_cesm2 = cesm2.variables['npp'][:]


lon_area_cesm2 = area_cesm2.variables['lon'][:]
areacella_cesm2 = area_cesm2.variables['areacella'][:]


#npp_cesm2,lon = shiftgrid(180., npp_cesm2, lon, start=False)
#areacella_cesm2,lon_area_cesm2 = shiftgrid(180., areacella_cesm2, lon_area_cesm2, start=False)

npp_cesm2[npp_cesm2==1.e+20]=np.nan
#areacella_cesm2[areacella_cesm2==1.e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_cesm2=np.std(npp_cesm2[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_cesm2
npp_cesm2=np.mean(npp_cesm2[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_cesm2


lat_cesm2 = cesm2.variables['lat'][:]
lon = cesm2.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_cesm2 = 0.9375
res_lon = 1.25
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_cesm2/dtor
lat_u = (lat_cesm2 + res_lat_cesm2)/dtor
weights_cesm2 = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
#npp_cesm2 = np.array(npp_cesm2)
#npp_cesm2[npp_cesm2==1e+20]=np.nan
#npp[npp==0.0]=np.nan


average=np.ma.average(npp_cesm2*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_cesm2*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_cesm2*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_cesm2*365*24*60*60*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_cesm2))))


print(np.shape(npp_cesm2))
print 'CESM2 - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

#sys.exit()
######################################UKESM #########################
lon = ukesm.variables['lon'][:]
npp_ukesm = ukesm.variables['npp'][:]

lon_area_ukesm = area_ukesm.variables['lon'][:]
areacella_ukesm = area_ukesm.variables['areacella'][:]


npp_ukesm,lon = shiftgrid(180., npp_ukesm, lon, start=False)
areacella_ukesm,lon_area_ukesm = shiftgrid(180., areacella_ukesm, lon_area_ukesm, start=False)

npp_ukesm[npp_ukesm==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan



std_npp_ukesm=np.std(npp_ukesm[(44*12)-1:(55*12)-1,:,:],axis=0)*areacella_ukesm
npp_ukesm=np.mean(npp_ukesm[(44*12)-1:(55*12)-1,:,:],axis=0)*areacella_ukesm



lat_ukesm = ukesm.variables['lat'][:]
lon = ukesm.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_ukesm = 1.25
res_lon = 1.875
#res_lat_ukesm = 0.9
#res_lon = 1.25
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_ukesm/dtor
lat_u = (lat_ukesm + res_lat_ukesm)/dtor
weights_ukesm = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
#npp_ukesm = np.array(npp_ukesm)
#npp_ukesm[npp_ukesm==1e+20]=np.nan
#npp[npp==0.0]=np.nan


average=np.ma.average(npp_ukesm*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_ukesm*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_ukesm*365*24*60*60*1000)/(10**15)

std_corr = np.sum(std_npp_ukesm*365*24*60*60*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_ukesm))))



print(np.shape(npp_ukesm))
print 'UKESM - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

#sys.exit()

######################################MPI #########################
lon = mpi.variables['lon'][:]
npp_mpi = mpi.variables['npp'][:]

lon_area_mpi = area_mpi.variables['lon'][:]
areacella_mpi = area_mpi.variables['areacella'][:]

npp_mpi,lon = shiftgrid(180., npp_mpi, lon, start=False)
areacella_mpi,lon_area_mpi = shiftgrid(180., areacella_mpi, lon_area_mpi, start=False)




npp_mpi[npp_mpi==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_mpi = np.std(npp_mpi[(4*12)-1:(15*12)-1,:,:],axis=0)*areacella_mpi
npp_mpi=np.mean(npp_mpi[(4*12)-1:(15*12)-1,:,:],axis=0)*areacella_mpi

#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_mpi = mpi.variables['lat'][:]
lon = mpi.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_mpi = 1.875
res_lon = 1.875
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_mpi/dtor
lat_u = (lat_mpi + res_lat_mpi)/dtor
weights_mpi = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



#npp_mpi=np.ma.masked_array(npp_mpi,np.isnan(npp_mpi))
average=np.ma.average(npp_mpi*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_mpi*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_mpi*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_mpi*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_mpi))))

print(np.shape(npp_mpi))
print 'MPI - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

######################################GFDL #########################
lon = gfdl.variables['lon'][:]
npp_gfdl = gfdl.variables['npp'][:]

lon_area_gfdl = area_gfdl.variables['lon'][:]
areacella_gfdl = area_gfdl.variables['areacella'][:]

npp_gfdl,lon = shiftgrid(180., npp_gfdl, lon, start=False)
areacella_gfdl,lon_area_gfdl = shiftgrid(180., areacella_gfdl, lon_area_gfdl, start=False)




npp_gfdl[npp_gfdl==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_gfdl = np.std(npp_gfdl[(44*12)-1:(55*12)-1,:,:],axis=0)*areacella_gfdl
npp_gfdl=np.mean(npp_gfdl[(44*12)-1:(55*12)-1,:,:],axis=0)*areacella_gfdl

npp_gfdl[npp_gfdl==0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_gfdl = gfdl.variables['lat'][:]
lon = gfdl.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_gfdl = 1.
res_lon = 1.25
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_gfdl/dtor
lat_u = (lat_gfdl + res_lat_gfdl)/dtor
weights_gfdl = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



npp_gfdl=np.ma.masked_array(npp_gfdl,np.isnan(npp_gfdl))
average=np.ma.average(npp_gfdl*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_gfdl*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_gfdl*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_gfdl*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_gfdl))))

print(np.shape(npp_gfdl))
print 'GFDL - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

######################################IPSL #########################
lon = ipsl.variables['lon'][:]
npp_ipsl = ipsl.variables['npp'][:]

lon_area_ipsl = area_ipsl.variables['lon'][:]
areacella_ipsl = area_ipsl.variables['areacella'][:]

npp_ipsl,lon = shiftgrid(180., npp_ipsl, lon, start=False)
areacella_ipsl,lon_area_ipsl = shiftgrid(180., areacella_ipsl, lon_area_ipsl, start=False)




npp_ipsl[npp_ipsl==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_ipsl = np.std(npp_ipsl[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_ipsl
npp_ipsl=np.mean(npp_ipsl[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_ipsl

#npp_ipsl[npp_ipsl==0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_ipsl = ipsl.variables['lat'][:]
lon = ipsl.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_ipsl = 1.259
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_ipsl/dtor
lat_u = (lat_ipsl + res_lat_ipsl)/dtor
weights_ipsl = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



npp_ipsl=np.ma.masked_array(npp_ipsl,np.isnan(npp_ipsl))
average=np.ma.average(npp_ipsl*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_ipsl*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_ipsl*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_ipsl*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_ipsl))))

print(np.shape(npp_ipsl))
print 'IPSL - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

######################################MIROC #########################
lon = miroc.variables['lon'][:]
npp_miroc = miroc.variables['npp'][:]

lon_area_miroc = area_miroc.variables['lon'][:]
areacella_miroc = area_miroc.variables['areacella'][:]

npp_miroc,lon = shiftgrid(180., npp_miroc, lon, start=False)
areacella_miroc,lon_area_ipsl = shiftgrid(180., areacella_miroc, lon_area_miroc, start=False)




npp_miroc[npp_miroc==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_miroc = np.std(npp_miroc[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_miroc
npp_miroc=np.mean(npp_miroc[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_miroc

#npp_ipsl[npp_ipsl==0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_miroc = miroc.variables['lat'][:]
lon = miroc.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_miroc = 1.259
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_miroc/dtor
lat_u = (lat_miroc + res_lat_miroc)/dtor
weights_ipsl = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



#npp_miroc=np.ma.masked_array(npp_miroc,np.isnan(npp_miroc))
average=np.ma.average(npp_miroc*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_miroc*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_miroc*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_miroc*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_miroc))))

print(np.shape(npp_miroc))
print 'MIROC - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()


######################################ACCESS #########################
lon = access.variables['lon'][:]
npp_access = access.variables['npp'][:]

lon_area_access = area_access.variables['lon'][:]
areacella_access = area_access.variables['areacella'][:]

npp_access,lon = shiftgrid(180., npp_access, lon, start=False)
areacella_access,lon_area_access = shiftgrid(180., areacella_access, lon_area_access, start=False)




npp_access[npp_access==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_access = np.std(npp_access[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_access
npp_access=np.mean(npp_access[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_access

#npp_ipsl[npp_ipsl==0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_access = access.variables['lat'][:]
lon = access.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_access = 1.259
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_access/dtor
lat_u = (lat_access + res_lat_access)/dtor
weights_ipsl = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



#npp_access=np.ma.masked_array(npp_access,np.isnan(npp_access))
average=np.ma.average(npp_access*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_access*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_access*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_access*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_access))))

print(np.shape(npp_access))
print 'ACCESS - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

######################################CNRM #########################
lon = cnrm.variables['lon'][:]
npp_cnrm = cnrm.variables['npp'][:]

lon_area_cnrm = area_cnrm.variables['lon'][:]
areacella_cnrm = area_cnrm.variables['areacella'][:]

npp_cnrm,lon = shiftgrid(180., npp_cnrm, lon, start=False)
areacella_cnrm,lon_area_cnrm = shiftgrid(180., areacella_cnrm, lon_area_cnrm, start=False)




npp_cnrm[npp_cnrm==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_cnrm = np.std(npp_cnrm[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_cnrm
npp_cnrm=np.mean(npp_cnrm[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_cnrm

#npp_ipsl[npp_ipsl==0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_cnrm = cnrm.variables['lat'][:]
lon = cnrm.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_cnrm = 1.259
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_cnrm/dtor
lat_u = (lat_cnrm + res_lat_cnrm)/dtor
weights_cnrm = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



#npp_access=np.ma.masked_array(npp_access,np.isnan(npp_access))
average=np.ma.average(npp_cnrm*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_cnrm*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_cnrm*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_cnrm*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_cnrm))))

print(np.shape(npp_cnrm))
print 'CNRM - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

######################################canESM #########################
lon = canesm.variables['lon'][:]
npp_canesm = canesm.variables['npp'][:]

lon_area_canesm = area_canesm.variables['lon'][:]
areacella_canesm = area_canesm.variables['areacella'][:]

npp_canesm,lon = shiftgrid(180., npp_canesm, lon, start=False)
areacella_canesm,lon_area_canesm = shiftgrid(180., areacella_canesm, lon_area_canesm, start=False)





npp_canesm[npp_canesm==1e+20]=np.nan

#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_canesm = np.std(npp_canesm[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_canesm
npp_canesm=np.mean(npp_canesm[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_canesm




#npp_canesm[npp_canesm<0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_canesm = canesm.variables['lat'][:]



x,y = np.meshgrid(lon, lat_canesm) 



npp_canesm = maskoceans(x, y, npp_canesm,resolution='l',grid=1.25,inlands=True)

#plt.imshow(npp_canesm)
#plt.show()
#sys.exit()


lon = canesm.variables['lon'][:]


#Function to calculate the weights of latitude
radius = 6367449
res_lat_canesm = 1.259
res_lon = 2.5
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_canesm/dtor
lat_u = (lat_canesm + res_lat_canesm)/dtor
weights_canesm = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



#npp_canesm=np.ma.masked_array(npp_canesm,np.isnan(npp_canesm))
average=np.ma.average(npp_canesm*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_canesm*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_canesm*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_canesm*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_canesm))))

print(np.shape(npp_ipsl))
print 'CanESM - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

######################################BCC-CSM2 #########################
lon = bcccsm.variables['lon'][:]
npp_bcccsm = bcccsm.variables['npp'][:]

lon_area_bcccsm = area_bcccsm.variables['lon'][:]
areacella_bcccsm = area_bcccsm.variables['areacella'][:]

npp_bcccsm,lon = shiftgrid(180., npp_bcccsm, lon, start=False)
areacella_bcccsm,lon_area_bcccsm = shiftgrid(180., areacella_bcccsm, lon_area_bcccsm, start=False)




npp_bcccsm[npp_bcccsm==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_bcccsm = np.std(npp_bcccsm[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_bcccsm
npp_bcccsm=np.mean(npp_bcccsm[(144*12)-1:(155*12)-1,:,:],axis=0)*areacella_bcccsm

#npp_ipsl[npp_ipsl==0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_bcccsm = bcccsm.variables['lat'][:]
lon = bcccsm.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_bcccsm = 1.
res_lon = 1.
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_bcccsm/dtor
lat_u = (lat_bcccsm + res_lat_ipsl)/dtor
weights_ipsl = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



#npp_bcccsm=np.ma.masked_array(npp_bcccsm,np.isnan(npp_bcccsm))
average=np.ma.average(npp_bcccsm*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_bcccsm*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_bcccsm*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_bcccsm*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_bcccsm))))

print(np.shape(npp_ipsl))
print 'BCC-CSM - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

######################################Nor-ESM2 #########################
lon = noresm.variables['lon'][:]
npp_noresm = noresm.variables['npp'][:]

#getting record 1
npp_noresm = npp_noresm[0,:,:,:]

lon_area_noresm = area_noresm.variables['lon'][:]
areacella_noresm = area_noresm.variables['areacella'][:]

npp_noresm,lon = shiftgrid(180., npp_noresm, lon, start=False)
areacella_noresm,lon_area_noresm = shiftgrid(180., areacella_noresm, lon_area_noresm, start=False)



npp_noresm[npp_noresm==1e+20]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

std_npp_noresm = np.std(npp_noresm[(4*12)-1:(15*12)-1,:,:],axis=0)*areacella_noresm
npp_noresm=np.mean(npp_noresm[(4*12)-1:(15*12)-1,:,:],axis=0)*areacella_noresm



npp_noresm[npp_noresm==0.0]=np.nan
#npp_mpi=np.mean(npp_mpi[(144*12)-1:(155*12)-1,:,:],axis=0)



lat_noresm = noresm.variables['lat'][:]
lon = noresm.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_noresm = 1.
res_lon = 1.
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_noresm/dtor
lat_u = (lat_noresm + res_lat_noresm)/dtor
weights_ipsl = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))


#PUPTAKE
#npp_mpi = np.array(npp_mpi)
#npp_mpi[npp_mpi==1e+20]=np.nan



npp_noresm=np.ma.masked_array(npp_noresm,np.isnan(npp_noresm))
average=np.ma.average(npp_noresm*(365*24*60*60)*1000,axis=0)
variance=np.ma.average((npp_noresm*(365*24*60*60)*1000-average)**2,axis=0)


mean_corr = np.sum(npp_noresm*(365*24*60*60*1000))/(10**15)

std_corr = np.sum(std_npp_noresm*(365*24*60*60*1000)/(10**15))

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_noresm))))

print(np.shape(npp_ipsl))
print 'Nor-ESM2 - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'


#sys.exit()

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


######################################IGBP #########################
lon = igbp.variables['lon'][:]
npp_igbp = igbp.variables['npp'][:]

#npp_igbp,lon = shiftgrid(180., npp_igbp, lon, start=False)

npp_igbp[npp_igbp==-99.99]=np.nan
#npp[npp<(np.max(abs(npp))*0.01)]=np.nan

#npp_mpi=np.mean(npp_igbp[(4*12)-1:(15*12)-1,:,:],axis=0)



lat_igbp = igbp.variables['lat'][:]
#lon = igbp.variables['lon'][:]

#Function to calculate the weights of latitude
radius = 6367449
res_lat_igbp = 1.
res_lon = 1.
m_pi = 3.14159265358979323846
dtor = 360./(2.*m_pi)

lon_u = res_lon/dtor
lon_l = 0.0/dtor

lat_l = lat_igbp/dtor
lat_u = (lat_igbp + res_lat_igbp)/dtor
weights_igbp = (radius*radius*(lon_u - lon_l)*(np.sin(lat_u)-np.sin(lat_l)))



#PUPTAKE
npp_igbp = np.array(npp_igbp)
npp_igbp[npp_igbp==-99.99]=np.nan

npp_igbp=np.ma.masked_array(npp_igbp,np.isnan(npp_igbp))
average=np.ma.average(npp_igbp,axis=0,weights=weights_igbp)
variance=np.ma.average((npp_igbp-average)**2,axis=0,weights=weights_igbp)


mean_corr = np.mean(np.ma.average(npp_igbp,axis=0,weights=weights_igbp))*148847000*(1000*1000)/(10**15)

std_corr = np.mean(np.sqrt(variance))*148847000*(1000*1000)/(10**15)

h_dn_corr, h_up_corr = stats.norm.interval(0.66, loc=mean_corr, scale= std_corr/np.sqrt(len(np.ma.array(npp_igbp))))

print(np.shape(npp_igbp))
print 'IGBP - Global NPP =', mean_corr, '+-', mean_corr-h_dn_corr, 'PgCyr-1'

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
res_lat_modis = 1.
res_lon = 1.
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



###################################### #########################

#plt.plot(lat_cesm2,weights_cesm2,label='cesm2')
#plt.plot(lat_ukesm,weights_ukesm,label='ukesm')
#plt.plot(lat_mpi,weights_mpi,label='mpi')
#plt.plot(lat_igbp,weights_igbp,'k--',label='mpi')
#plt.legend()
#plt.show()


#plt.plot(lat_cesm2,areacella_cesm2[:,0],label='cesm2')
#plt.plot(lat_ukesm,areacella_ukesm[:,0],label='ukesm')
#plt.plot(lat_mpi,areacella_mpi[:,0],label='mpi')
#plt.legend()
#plt.show()

#plt.plot(lat_cesm2,areacella_cesm2[:,0]/weights_cesm2,label='cesm2')
#plt.plot(lat_ukesm,areacella_ukesm[:,0]/weights_ukesm,label='ukesm')
#plt.plot(lat_mpi,areacella_mpi[:,0]/weights_mpi,label='mpi')
#plt.legend()
#plt.show()
#sys.exit()




gpp_mte_ori_cesm2 = []
gpp_mte_ori_ukesm = []
gpp_mte_ori_mpi = []
gpp_mte_ori_gfdl = []
gpp_mte_ori_ipsl = []
gpp_mte_ori_miroc = []
gpp_mte_ori_access = []
gpp_mte_ori_cnrm = []
gpp_mte_ori_canesm = []
gpp_mte_ori_bcccsm = []
gpp_mte_ori_noresm = []
gpp_mte_ori_e3sm = []
gpp_mte_ori_e3sm_fun = []
gpp_mte_ori_e3sm_funp = []

gpp_mte_ori_igbp = []
gpp_mte_ori_modis = []




#gpp_mte_ori_cesm2.append(np.nanmean(npp_cesm2,axis = 1))
#gpp_mte_ori_ukesm.append(np.mean(npp_ukesm,axis = 1))
#gpp_mte_ori_mpi.append(np.mean(npp_mpi,axis = 1))
#gpp_mte_ori_igbp.append(np.mean(npp_igbp,axis = 1))

gpp_mte_ori_cesm2.append(np.ma.average(npp_cesm2/areacella_cesm2,axis = 1))
gpp_mte_ori_ukesm.append(np.ma.average(npp_ukesm/areacella_ukesm,axis = 1))
gpp_mte_ori_mpi.append(np.ma.average(npp_mpi/areacella_mpi,axis = 1))
gpp_mte_ori_gfdl.append(np.ma.average(npp_gfdl/areacella_gfdl,axis = 1))
gpp_mte_ori_ipsl.append(np.ma.average(npp_ipsl/areacella_ipsl,axis = 1))
gpp_mte_ori_miroc.append(np.ma.average(npp_miroc/areacella_miroc,axis = 1))
gpp_mte_ori_access.append(np.ma.average(npp_access/areacella_access,axis = 1))
gpp_mte_ori_cnrm.append(np.ma.average(npp_cnrm/areacella_cnrm,axis = 1))
gpp_mte_ori_canesm.append(np.ma.average(npp_canesm/areacella_canesm,axis = 1))
gpp_mte_ori_bcccsm.append(np.ma.average(npp_bcccsm/areacella_bcccsm,axis = 1))
gpp_mte_ori_noresm.append(np.ma.average(npp_noresm/areacella_noresm,axis = 1))

gpp_mte_ori_e3sm.append(np.ma.average(npp_e3sm,axis = 1))
gpp_mte_ori_e3sm_fun.append(np.ma.average(npp_e3sm_fun,axis = 1))
gpp_mte_ori_e3sm_funp.append(np.ma.average(npp_e3sm_funp,axis = 1))

gpp_mte_ori_igbp.append(np.ma.average(npp_igbp,axis = 1))
gpp_mte_ori_modis.append(np.ma.average(npp_modis,axis = 1))

#sys.exit()
my_pfts = ['All']
x_cesm2_1 = range(len(lat_cesm2))
x_ukesm_1 = range(len(lat_ukesm))
x_mpi_1 = range(len(lat_mpi))
x_igbp_1 = range(len(lat_igbp))

# Interpolating data for ensemble member
const = (365*24*60*60)*1000
f1 = interp1d(lat_cesm2,np.array(gpp_mte_ori_cesm2[0][:]*const))    
f2 = interp1d(lat_ukesm,np.array(gpp_mte_ori_ukesm[0][:]*const))
f3 = interp1d(lat_mpi,np.array(gpp_mte_ori_mpi[0][:]*const))
f4 = interp1d(lat_gfdl,np.array(gpp_mte_ori_gfdl[0][:]*const))
f5 = interp1d(lat_ipsl,np.array(gpp_mte_ori_ipsl[0][:]*const))
f6 = interp1d(lat_miroc,np.array(gpp_mte_ori_miroc[0][:]*const))
f7 = interp1d(lat_access,np.array(gpp_mte_ori_access[0][:]*const))
f8 = interp1d(lat_cnrm,np.array(gpp_mte_ori_cnrm[0][:]*const))        
f9 = interp1d(lat_canesm,np.array(gpp_mte_ori_canesm[0][:]*const))        
f10 = interp1d(lat_bcccsm,np.array(gpp_mte_ori_bcccsm[0][:]*const))
f11 = interp1d(lat_noresm,np.array(gpp_mte_ori_noresm[0][:]*const))

f12 = interp1d(lat_e3sm,np.array(gpp_mte_ori_e3sm[0][:]*const/1000.))
f13 = interp1d(lat_e3sm_fun,np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.))
f14 = interp1d(lat_e3sm_funp,np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.))

x = np.linspace(-87.5,87.5,num=172,endpoint=True)

ensemble = np.zeros((11,172))

ensemble[0] = f1(x)
ensemble[1] = f2(x)
ensemble[2] = f3(x)
ensemble[3] = f4(x)
ensemble[4] = f5(x)
ensemble[5] = f6(x)
ensemble[6] = f7(x)
ensemble[7] = f8(x)
ensemble[8] = f9(x)
ensemble[9] = f10(x)
ensemble[10] = f11(x)

        
ensemble_member = np.mean(ensemble,axis=0)
ensemble_member_std = np.std(ensemble,axis=0)
ensemble_member_min = np.min(ensemble,axis=0)
ensemble_member_max = np.max(ensemble,axis=0)


for i, pft in enumerate(my_pfts):

   	fig = plt.figure(figsize=(48, 24)) 
        
        #From gCm-2s-1 to TgCyr-1
        const = (365*24*60*60)*1000
        x_cesm2 = np.floor(((np.asarray(x_cesm2_1) + 1)*res_lat_cesm2) -90.) 
        x_ukesm = np.floor(((np.asarray(x_ukesm_1) + 1)*res_lat_ukesm) -89.375) 
        #plt.stackplot(np.floor(((np.asarray(x) + 1)*res_lat) -90.),np.array(gpp_mte_ori_mean[0][:])*(const),edgecolor='k',linewidth=3.3,colors=['darkorange'],labels=['CESM2']) 
        plt.plot(lat_cesm2,np.array(gpp_mte_ori_cesm2[0][:]*const),label='CESM2',linewidth=3.3, alpha =0.5)    
        plt.plot(lat_ukesm,np.array(gpp_mte_ori_ukesm[0][:]*const),label='UKESM',linewidth=3.3, alpha =0.5)
        plt.plot(lat_mpi,np.array(gpp_mte_ori_mpi[0][:]*const),label='MPI',linewidth=3.3, alpha =0.5)
        plt.plot(lat_gfdl,np.array(gpp_mte_ori_gfdl[0][:]*const),label='GFDL',linewidth=3.3, alpha =0.5)
        plt.plot(lat_ipsl,np.array(gpp_mte_ori_ipsl[0][:]*const),label='IPSL',linewidth=3.3, alpha =0.5)
        plt.plot(lat_miroc,np.array(gpp_mte_ori_miroc[0][:]*const),label='MIROC',linewidth=3.3, alpha =0.5)
        plt.plot(lat_access,np.array(gpp_mte_ori_access[0][:]*const),label='ACCESS',linewidth=3.3, alpha =0.5)
        plt.plot(lat_cnrm,np.array(gpp_mte_ori_cnrm[0][:]*const),label='CNRM',linewidth=3.3, alpha =0.5)        
	plt.plot(lat_canesm,np.array(gpp_mte_ori_canesm[0][:]*const),label='CanESM',linewidth=3.3, alpha =0.5)        
	plt.plot(lat_bcccsm,np.array(gpp_mte_ori_bcccsm[0][:]*const),label='BCC-CSM',linewidth=3.3, alpha =0.5)
	plt.plot(lat_noresm,np.array(gpp_mte_ori_noresm[0][:]*const),label='Nor-ESM2',linewidth=3.3, alpha =0.5)
	plt.plot(lat_e3sm,np.array(gpp_mte_ori_e3sm[0][:]*const/1000.),label='E3SM',linewidth=9.3, color='blue',alpha =1.0)
	plt.plot(lat_e3sm_fun,np.array(gpp_mte_ori_e3sm_fun[0][:]*const/1000.),label='E3SM-FUN2',linewidth=9.3,  color='darkorange',alpha =1.0)
	plt.plot(lat_e3sm_funp,np.array(gpp_mte_ori_e3sm_funp[0][:]*const/1000.),label='E3SM-FUN3',linewidth=9.3,  color='green',alpha =1.0)

        plt.plot(lat_igbp,np.array(gpp_mte_ori_igbp[0][:]), label='IGBP',color='k',linestyle='-.',linewidth=9.3)
        plt.plot(lat_modis,np.array(gpp_mte_ori_modis[0][:]*const), label='MODIS',color='k',linestyle='--',linewidth=9.3)
        #plt.plot(x,ensemble_member, label='Ensemble',color='k',linestyle='dotted',linewidth=9.3)
        plt.fill_between(x,ensemble_member_min,ensemble_member_max, label='CMIP6\nrange',facecolor='gray',alpha =0.15)

        #max_val = np.max(gpp_mte_ori_mean_pnonmyc[0][:]*const)
        max_val = np.max(gpp_mte_ori_mpi[0][:]*const)

	plt.legend(ncol=2, prop={'size': 60},loc=0)
	axes = plt.gca()
	axes.set_xlim([-60.,80.])
	axes.set_ylim([0.0,max_val + 0.1*max_val])

        #plt.xticks(np.arange(min((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.), max((np.asarray(x_ukesm_1)*res_lat_ukesm) -90.) + 1., 30.))

	plt.grid(True)

	plt.tick_params(axis='both', which='major', pad=20)
	plt.xlabel(r'Latitude ($\degree$)')
	plt.ylabel(r'NPP (g C m$^{-2}$ yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('CMIP6_zonal_mean_v2.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()
        plt.close()




        #sys.exit()




        ####### Scatter plot #####
        xfit = np.arange(0,1400,1)

        model.fit(ensemble_member[:, np.newaxis],f12(x))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope1 = model.coef_[0]

        r1 = model.score(ensemble_member[:, np.newaxis],f12(x))
        r1 = r2_score(ensemble_member[:, np.newaxis],f12(x))
        print("R2: ", r1)
 

        yfit1 = model.predict(xfit[:, np.newaxis])
        yfit11 = model.predict(ensemble_member[:, np.newaxis])

        rmse1 = np.sqrt(sum((yfit11 - f12(x))**2)/len(f12(x)))
        rmse1 = np.sqrt(mean_squared_error(ensemble_member[:, np.newaxis],f12(x)))
        print("RMSE1: ", rmse1)


        model.fit(ensemble_member[:, np.newaxis],f13(x))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope2 = model.coef_[0]

        r2 = model.score(ensemble_member[:, np.newaxis],f13(x))
        r2 = r2_score(ensemble_member[:, np.newaxis],f13(x))
        print("R2: ", r2)

        yfit2 = model.predict(xfit[:, np.newaxis])
        yfit22 = model.predict(ensemble_member[:, np.newaxis])

        rmse2 = np.sqrt(sum((yfit22 - f13(x))**2)/len(f13(x)))
        rmse2 = np.sqrt(mean_squared_error(ensemble_member[:, np.newaxis],f13(x)))
        print("RMSE2: ", rmse2)

        model.fit(ensemble_member[:, np.newaxis],f14(x))

        print("Model slope: ", model.coef_[0])
        print("Model intercept: ", model.intercept_)
        
        slope3 = model.coef_[0]

        r3 = model.score(ensemble_member[:, np.newaxis],f14(x))
        r3 = r2_score(ensemble_member[:, np.newaxis],f14(x))
        print("R2: ", r3)

        yfit3 = model.predict(xfit[:, np.newaxis])
        yfit33 = model.predict(ensemble_member[:, np.newaxis])

        rmse3 = np.sqrt(sum((yfit33 - f14(x))**2)/len(f14(x)))
        rmse3 = np.sqrt(mean_squared_error(ensemble_member[:, np.newaxis],f14(x)))
        print("RMSE3: ", rmse3)

   	fig = plt.figure(figsize=(36, 36)) 
        
        #From gCm-2s-1 to TgCyr-1

     
	plt.scatter(ensemble_member,f12(x),label='E3SM (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r1,rmse1,slope1),linewidth=9.3,s=600)
	plt.scatter(ensemble_member,f13(x),label='E3SM-FUN2.0 (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r2,rmse2,slope2),linewidth=9.3,s=600)
	plt.scatter(ensemble_member,f14(x),label='E3SM-FUN3.0 (R$^{2}$ = %.2f;\nRMSE = %.2f;\nslope = %.2f)' % (r3,rmse3,slope3),linewidth=9.3,s=600)
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
	plt.xlabel(r'NPP CMIP6 Ensemble (g C m$^{-2}$ yr$^{-1}$)')
	plt.ylabel(r'NPP model (g C m$^{-2}$ yr$^{-1}$)' )
        #plt.ylabel('' +long_name+' ('+unit[0:2]+' m$^2$ yr$^{-1}$)' )

        #ttl=plt.title(r'ZONAL MEAN - %s' % (pft) ,fontsize = 84)
        #ttl.set_position([.5, 1.05])

        plt.savefig('scatter_ensemble_v4.png',bbox_inches="tight")
        #print 'figures/opt4_%s.png saved!'% (pft)
	#plt.show()
        plt.close()



sys.exit()



