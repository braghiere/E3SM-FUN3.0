import sys
import os
import datetime
import numpy as np

begin = datetime.datetime.now()

print(datetime.datetime.now())

print('Ploting Mycorrhizal uptake...')
os.system('python scripts/01_plot_pactive.py')
print('Done! Look for PACTIVE.png')

print('Ploting Direct root uptake...')
os.system('python scripts/02_plot_pnonmyc.py')
print('Done! Look for PNONMYC.png')

print('Ploting Retranslocation...')
os.system('python scripts/03_plot_pretrans.py')
print('Done! Look for PRETRANS.png')

print('Ploting AM root uptake...')
os.system('python scripts/04_plot_pam.py')
print('Done! Look for PAM.png')

print('Ploting ECM root uptake...')
os.system('python scripts/05_plot_pecm.py')
print('Done! Look for PECM.png')

print('Ploting PRetrans ratio...')
os.system('python scripts/06_plot_pretrans_ratio.py')
print('Done! Look for PRETRANS_ratio.png')

print('Ploting Zonal mean total...')
os.system('python scripts/07_plot_zonal_mean_total.py')
print('Done! Look for PUPTAKE_zonal_mean.png')

print('Ploting Zonal mean total...')
os.system('python scripts/08_plot_zonal_mean_ratio.py')
print('Done! Look for PUPTAKE_zonal_mean_ratio.png')

print('Ploting Monthly means per PFT...')
os.system('python scripts/09_plot_monthly_mean.py')
print('Done! Look for PUPTAKE_monthly_mean_PFT_names.png ')

print('Ploting Carbon use total...')
os.system('python scripts/10_plot_carbon_use_pft.py')
print('Done! Look for COST_PUPTAKE_zonal_pft.png')

print('Ploting zonal NPP...')
os.system('python scripts/11_plot_zonal_mean_npp.py')
print('Done! Look for NPP_zonal_mean.png')

print('Ploting PUPTAKE_NPP_FRACTION ratio...')
os.system('python scripts/12_plot_puptake_npp_ratio.py')
print('Done! Look for PUPTAKE_NPP_ratio.png')

print('Ploting PUPTAKE_NPP_FRACTION ratio...')
os.system('python scripts/13_plot_puptake_npp_zones.py')
print('Done! Look for PUPTAKE_NPP_ratio.png')

print('Ploting Mycorrhizal uptake...')
os.system('python scripts/14_plot_cost_pactive.py')
print('Done! Look for PACTIVE.png')

print('Ploting Mycorrhizal uptake...')
os.system('python scripts/15_plot_cost_pnonmyc.py')
print('Done! Look for PACTIVE.png')

print('Ploting Mycorrhizal uptake...')
os.system('python scripts/16_plot_cost_pretrans.py')
print('Done! Look for PACTIVE.png')

print('Ploting Solution P...')
os.system('python scripts/17_plot_*.py')
print('Done! Look for SolutionP.png')

print('Ploting Mycorrhizal uptake...')
os.system('python scripts/18_plot_NPP_diff_norm.py')
print('Done! Look for PACTIVE.png')

print('Ploting Mycorrhizal uptake...')
os.system('python scripts/19_plot_p_limitation.py')
print('Done! Look for PACTIVE.png')

for i in [20,21,22,23,24,25]:
   print('Ploting Solution P...')
   os.system('python scripts/%s_plot_*.py' % (int(i)))
   print('Done! Plotted figure %s' % (int(i)))

print('This took this time to finish:')
print(datetime.datetime.now() - begin)
print('Done!')
