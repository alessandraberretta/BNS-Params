import pandas as pd
import numpy as np
from scipy.stats import cauchy
import matplotlib.pyplot as plt
import cosmcalc_ale

rng = np.random.default_rng()

ns = 10000
one = np.ones((ns,))

thetaobs = np.linspace(0, 0.4, ns)

mu_thetaobs = 0.05
gamma_thetaobs = 0.15
fs_thetaobs = 3
dist_thetaobs = cauchy(mu_thetaobs, gamma_thetaobs)

z = np.linspace(0, 4, ns)
mu_z = 0.3
gamma_z = 0.8
fs_z = 34
dist_z = cauchy(mu_z, gamma_z)

scale_array = one*1/fs_z
loc_array = one*0

dl = np.linspace(200, 40000, ns)

mu_dl = 1000
gamma_dl = 6500
fs_dl = 370000
dist_dl = cauchy(mu_dl, gamma_dl)

def compute_bns_merger_angle():
    lor_rad = values_thetaobs
    lor_rad = np.array(lor_rad)
    lor_rad_deg = np.rad2deg(lor_rad)
    lor_rad_deg_complementary = 180 - lor_rad_deg
    lor_rad_final = np.deg2rad(lor_rad_deg_complementary)
    return lor_rad_final


values_z = []
values_dl = []
min_range = 0
max_range = 4
while len(values_z)<ns:
  tmp = dist_z.rvs(1)[0]
  if min_range<tmp<max_range:
    values_z.append(tmp)
    values_dl.append(cosmcalc_ale.compute_dl_Mpc(tmp))

values_thetaobs = []

min_range = 0
max_range = 0.4
while len(values_thetaobs)<ns:
  tmp = dist_thetaobs.rvs(1)[0]
  if min_range<tmp<max_range:
    values_thetaobs.append(tmp)

'''
min_range = 200
max_range = 40000
while len(values_dl)<ns:
  tmp = dist_dl.rvs(1)[0]
  if min_range<tmp<max_range:
    values_dl.append(tmp)
'''
distlum = []
zdist = []

for elm in rng.uniform(0, 4, size=(ns,)):
  zdist.append(elm)
  distlum.append(cosmcalc_ale.compute_dl_Mpc(elm))


parameters = pd.DataFrame.from_dict({
    'redshift': values_z,
    'luminosity_distance': values_dl, 
    'mass_1': rng.uniform(1, 2., size=(ns,)),
    'mass_2': rng.uniform(1, 2., size=(ns,)),
    # 'theta_jn': 0.2, 
    'theta_jn': compute_bns_merger_angle(),
    # 'theta_jn': np.arccos(rng.uniform(-1., 1., size=(ns,))), 
    'dec': np.arccos(rng.uniform(-1., 1., size=(ns,))) - np.pi / 2.,
    'ra': rng.uniform(0, 2. * np.pi, size=(ns,)),
    # 'dec': -10.8,
    # ra': 66.6,
    'psi': rng.uniform(0, 2. * np.pi, size=(ns,)),
    'phase': rng.uniform(0, 2. * np.pi, size=(ns,)),
    'lambda_1': rng.uniform(0, 1000, size=(ns,)),
    'lambda_2': rng.uniform(0, 1000, size=(ns,)),
    'geocent_time': rng.uniform(1735257618, 1766793618, size=(ns,)) # full year 2035
})

print(parameters['redshift'])
print(parameters['luminosity_distance'])
print(min(parameters['theta_jn']))
parameters.to_hdf('population_z_dl_thetaobs_10000.hdf5', mode='w', key='root')