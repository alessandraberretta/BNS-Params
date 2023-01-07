import pandas as pd
import numpy as np

rng = np.random.default_rng()

ns = 100000
one = np.ones((ns,))
parameters = pd.DataFrame.from_dict({
    'redshift': 0.1*one,
    'luminosity_distance': 383*one, 
    'mass_1': rng.uniform(1, 2., size=(ns,)),
    'mass_2': rng.uniform(1, 2., size=(ns,)),
    'theta_jn': np.arccos(rng.uniform(-1., 1., size=(ns,))),
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
parameters.to_hdf('100628A_population_100000.hdf5', mode='w', key='root')