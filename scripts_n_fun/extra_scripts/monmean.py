import xarray as xr
import nctoolkit as nc

#ds = xr.open_dataset('elm_pft.nc')
#monthly_data = ds.resample(freq='m',dim='time',how='mean')
#monthly_data = ds.resample(time='1D').mean()
#monthly_data = ds.resample(time='m').mean()

data = nc.open_data('elm_pft.nc','r')
df = data.monthly_mean()

#df = data.to_dataframe
print(df)

#df.to_netcdf(path='test.nc',mode='w')


