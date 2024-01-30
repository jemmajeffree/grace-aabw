import numpy as np
import xarray as xr
import sys
sys.path.append('/g/data/x77/jj8842/')
import cosimagrace as cg

titlefontsize=8

ACCESS_OM2_01 = '01deg_jra55v140_iaf_cycle3'
ACCESS_OM2_01_ryf = '01deg_jra55v13_ryf9091'
ACCESS_OM2_025 = '025deg_jra55_iaf_omip2_cycle5'
ACCESS_OM2 = '1deg_jra55_iaf_omip2_cycle2'
MOM6_025 = 'OM4_025.JRA_RYF'

ocean_cutoff=0.99
NS_extent = (-10,10)
remove_climatology=True
fi = 3

add_SLR = (lambda x: x, lambda x : x+1)
do_nothing2 = (lambda x: x, lambda x : x)
do_nothing = (lambda x: x,)

def trim_add_SLR(t):
    return (lambda x: x[t:], lambda x : x[t:]+1)
def trim(t,n):
    return [lambda x: x[t:]]*n

JPL_sigma = 0.015
ANU_sigma = 0.021
GSFC_trend_sigma = 0.0025
add_JPL_noise = (lambda x : x+np.random.normal(0,JPL_sigma,x.shape),)
add_ANU_noise = (lambda x : x+np.random.normal(0,0.022,x.shape),)

saveloc = '/g/data/x77/jj8842/'
transect = xr.load_dataset('/g/data/x77/jj8842/data_v4/transect_definitions.nc')

t_names = ('West Pacific','East Pacific', 'West Atlantic', 'East Atlantic', 'West Indian', 'East Indian')
t_ids = (14,15,11,12,2,3)

wb_names = ('Pacific','Atlantic','Indian')
wb_ids = (39,40,41)
markers = ('o','s','^')

filters = (lambda x: cg.rollave(x,12),
           lambda x: cg.low_pass(x,6),
           lambda x: cg.low_pass(x,12),
           lambda x: cg.low_pass(x,24),
           lambda x: cg.low_pass(x,36),
           lambda x: cg.low_pass(x,48),
          lambda x: cg.low_pass(x,72),
          lambda x: cg.low_pass(x,96),)

window_func = [np.ones(12),
                        np.exp(-np.linspace(-3,3,6)**2/2)/(np.sqrt(2)*np.pi),
                        np.exp(-np.linspace(-3,3,12)**2/2)/(np.sqrt(2)*np.pi),
                        np.exp(-np.linspace(-3,3,24)**2/2)/(np.sqrt(2)*np.pi),
                        np.exp(-np.linspace(-3,3,36)**2/2)/(np.sqrt(2)*np.pi),
                        np.exp(-np.linspace(-3,3,48)**2/2)/(np.sqrt(2)*np.pi),
               np.exp(-np.linspace(-3,3,72)**2/2)/(np.sqrt(2)*np.pi),
               np.exp(-np.linspace(-3,3,96)**2/2)/(np.sqrt(2)*np.pi),
              ]
for i,w in enumerate(window_func):
    window_func[i] = w/np.sum(w)
                        

filter_names = ('12m roll_ave',
                '6m LP',
                '12m LP',
                '24m LP',
                '36m LP',
                '48m LP',
                '72m LP',
                '96m LP',
                'raw',
               )