import sys
import os
import datetime

begin = datetime.datetime.now()

print(datetime.datetime.now())

print('Ploting Mycorrhizal uptake...')
os.system('python scripts_n/01_plot_pactive.py')
print('Done! Look for PACTIVE.png')

print('Ploting Direct root uptake...')
os.system('python scripts_n/02_plot_pnonmyc.py')
print('Done! Look for PNONMYC.png')

print('Ploting Retranslocation...')
os.system('python scripts_n/03_plot_pretrans.py')
print('Done! Look for PRETRANS.png')

print('Ploting AM root uptake...')
os.system('python scripts_n/04_plot_pam.py')
print('Done! Look for PAM.png')

print('Ploting ECM root uptake...')
os.system('python scripts_n/05_plot_pecm.py')
print('Done! Look for PECM.png')

print('Ploting PRetrans ratio...')
os.system('python scripts_n/06_plot_pretrans_ratio.py')
print('Done! Look for PRETRANS_ratio.png')

print('Ploting Zonal mean total...')
os.system('python scripts_n/07_plot_zonal_mean_total.py')
print('Done! Look for PUPTAKE_zonal_mean.png')

print('Ploting Zonal mean total...')
os.system('python scripts_n/08_plot_zonal_mean_ratio.py')
print('Done! Look for PUPTAKE_zonal_mean_ratio.png')

print('Ploting Monthly means per PFT...')
os.system('python scripts_n/09_plot_monthly_mean.py')
print('Done! Look for PUPTAKE_monthly_mean_PFT_names.png ')

print('Ploting Carbon use total...')
os.system('python scripts_n/10_plot_carbon_use_pft.py')
print('Done! Look for COST_PUPTAKE_zonal_pft.png')

print('Ploting zonal NPP...')
os.system('python scripts_n/11_plot_zonal_mean_npp.py')
print('Done! Look for NPP_zonal_mean.png')

print('Ploting NFIX...')
os.system('python scripts_n/12_plot_nfix.py')
print('Done! Look for NFIX.png')

print('Ploting Zonal mean GPP...')
os.system('python scripts_n/13_plot_zonal_mean_gpp.py')
print('Ploting Zonal mean GPP...')

print('Ploting PUPTAKE_NPP_FRACTION ratio...')
os.system('python scripts_n/14_plot_puptake_npp_ratio.py')
print('Done! Look for PUPTAKE_NPP_ratio.png')

print('This took this time to finish:')
print(datetime.datetime.now() - begin)
print('Done!')
