from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from argparse import ArgumentParser
from scipy.stats import cauchy

import warnings
warnings.filterwarnings("ignore")

rng = np.random.default_rng()

L = rng.standard_cauchy(600)
L = L[(L>0) & (L<1)]

parser = ArgumentParser(description='Plot for BNS-GRB analysis with GWFish in the ET context')

parser.add_argument("-vartoplot", "--vartoplot", type=str,
                    dest="vartoplot", help="Make plot of redshift or observational angle")
args = parser.parse_args()

def read_data_file(file, sep='\t'):
    return pd.read_csv(file, sep)

path_dirs = '/Users/alessandraberretta/'
dirs_list = [path_dirs +
                 dir for dir in os.listdir(path_dirs) if dir.startswith('20')]
summary_list = []

GRB = ['111117A', '060801', '090510', '101219A', '060313', '090515', 
    '061201', '100625A', '061217', '100628A', '090417A', '070209', 
    '100117A', '100206A', '060502B', '140903A', '131105A', '081203A', 
    '091127', '120909A', '120804A', '140930B', '150821A', '061006', 
    '130603B', '050318']

for elm in dirs_list: 
    # if elm.endswith('results'):
    if 'results' in elm: 
        for file in [elm + '/' + file for file in os.listdir(elm) if file.startswith('new_sum_') and file.endswith('.csv')]:
            if file.split('_')[6][0:-4] in GRB: 
                print(file.split('_')[6][0:-4])
            # print(file.split('_')[6][0:-4])
                summary_list.append(file)
    if 'simbad_res' in elm: 
        for file in [elm + '/' + file for file in os.listdir(elm) if file.startswith('new_sum_') and file.endswith('.csv')]:
            if file.split('_')[7][0:-4] in GRB:
                print(file.split('_')[7][0:-4])
            # print(file.split('_')[6][0:-4])
                summary_list.append(file)

theta_obs = []
redshift = []
for idx, elm in enumerate(summary_list):
    DF = read_data_file(elm)
    theta_obs.append(float(DF['Values'][2]))
    redshift.append(float(DF['Values'][8]))

print(theta_obs)

z = [2.3, 1.1, 0.9, 0.7, 0.1, 0.4, 0.1, 0.4, 0.8, 0.1, 0.08, 1.1, 0.9, 0.4, 1.5, 0.3, 1.6, 2.1, 0.4, 3.9, 1.3, 4, 
    0.7, 0.4, 0.3, 1.4]

if args.vartoplot == 'z':
    plt.hist(redshift, bins=5, histtype='step', color='darkorange', lw=2)
    plt.xlabel('redshift')
    plt.ylabel('counts')
    plt.show()

if args.vartoplot == 'theta_obs':
    plt.hist(theta_obs, bins=3, histtype='step', color='darkred', lw=2)
    # plt.hist(L, bins=20, histtype='step', color='black', lw=2)
    plt.xlabel(r'$\theta_{obs}$ [rad]')
    plt.ylabel('counts')
    plt.show()

# df = pd.DataFrame(dict(thetaobs=theta_obs, redshift=redshift))
# df.to_csv('theta_obs_z.csv', sep='\t')
# Define the distribution parameters to be plotted
# gamma_values = [0.5, 1.0, 2.0]
'''
linestyles = ['-', '--', ':']
mu = 0.05
gamma = 0.15
x = np.linspace(0, 0.4, 1000)

#------------------------------------------------------------
# plot the distributions
# fig, ax = plt.subplots(figsize=(5, 3.75))

#for gamma, ls in zip(gamma_values, linestyles):
dist = cauchy(mu, gamma)

plt.plot(x, dist.pdf(x)*3, color='black',
             label=r'$\mu=%i,\ \gamma=%.1f$' % (mu, gamma))

plt.xlim(0, 0.6)
plt.ylim(0, 10)

plt.xlabel('$x$')
plt.ylabel(r'$p(x|\mu,\gamma)$')
plt.title('Cauchy Distribution')

plt.legend()
plt.show()
'''