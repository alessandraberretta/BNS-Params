import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np 
import sys

path_dirs_thetajn = '/Users/alessandraberretta/bns_analysis_res/thetajn/'
path_dirs_no_thetajn = '/Users/alessandraberretta/bns_analysis_res/no_thetajn/' 

thetajn_list = [path_dirs_thetajn + file for file in os.listdir(path_dirs_thetajn) if file.startswith('Errors_ET_') and file.endswith('1e5.txt')]
no_thetajn_list = [path_dirs_no_thetajn + file for file in os.listdir(path_dirs_no_thetajn) if file.startswith('Errors_ET_')]
print(no_thetajn_list)

rate_weight = []
rate = []

names = []
dl = []
z = []
rate_list = []
err_rate_list= []
thetajn = []

namesnotheta = []
dlnotheta = []
znotheta = []
ratenotheta = []
err_rate_listnotheta = []

for file_thetajn in thetajn_list:
    df = pd.read_csv(file_thetajn, sep=' ')
    names.append(file_thetajn.split('_')[5])
    dL = df['luminosity_distance']
    dl.append(dL[0])
    z.append(df['redshift'][0])
    thetajn.append(df['theta_jn'][0])
    rate_list.append(len(df))
    err_rate_list.append(np.sqrt(1e5*((len(df)/100000))*(1-(len(df)/100000))))
    rate_weight.append((len(df)/100000)*(365/dL[0]))
    rate.append((len(df)/100000))

for file_no_thetajn in no_thetajn_list:
    df = pd.read_csv(file_no_thetajn, sep=' ')
    dL = df['luminosity_distance']
    namesnotheta.append(file_no_thetajn.split('_')[6])
    dlnotheta.append(dL[0])
    znotheta.append(df['redshift'][0])
    ratenotheta.append((len(df)))
    err_rate_listnotheta.append(np.sqrt(1e5*((len(df)/100000))*(1-(len(df)/100000))))
    print(len(df))
    print(err_rate_listnotheta)
    print('no_thetajn')
    rate_weight.append((len(df)/100000)*(365/dL[0]))
    rate.append((len(df)/100000))

rate_list_err = []
for id, elm in enumerate(rate_list):
    rate_list_err.append(f'${elm} \pm {int(err_rate_list[id])}$')

rate_list_errnotheta = []
for id, elm in enumerate(ratenotheta):
    rate_list_errnotheta.append(f'${elm} \pm {int(err_rate_listnotheta[id])}$')

# dfthetajn =  pd.DataFrame(dict(GRB=names, d_L=dl, z=z, thetajn = thetajn, rate=rate_list_err))
# print(dfthetajn.to_latex(index=False, longtable=True, escape=False)) 
dfnothetajn =  pd.DataFrame(dict(GRB=namesnotheta, d_L=dlnotheta, z=znotheta, rate=rate_list_errnotheta))
print(dfnothetajn.to_latex(index=False, longtable=True, escape=False)) 

err_GRB_zdist = np.sqrt(26)/(14000*16)
err_GW_zdist = np.sqrt(1e4*(6142/10000)*(1-(6142/10000)))
# err_tot = 14000*np.sqrt(np.power(1e-4*err_GW/10000,2)+np.power(0.61*err_GRB,2))
err_tot = np.sqrt(np.power(err_GRB_zdist/1e-4,2)+np.power(err_GW_zdist/(10000*0.61),2))*0.86
print('total rate:', (6142/10000)*(1e-4)*14000)
print('err total rate:', err_tot)

# print('total rate (simple mean)', np.mean(rate)*1e-4*14000)

err_fraction = ((np.sqrt(1e5*(np.mean(rate))*(1-(np.mean(rate))))/100000)/np.sqrt(len(rate)))
print(err_fraction)
err_GRB = np.sqrt(26)/(14000*16)
tot_err_rate_mean = np.sqrt(np.power(err_GRB/1e-4,2)+np.power(err_fraction/np.mean(rate),2))*np.mean(rate)*1e-4*14000
print('total rate (simple mean)', np.mean(rate)*1e-4*14000)
print('error total rate (simple mean)', tot_err_rate_mean)

err_fraction_weight = ((np.sqrt(1e5*(np.mean(rate_weight))*(1-(np.mean(rate_weight))))/100000)/np.sqrt(len(rate_weight)))
err_GRB =  np.sqrt(26)/(14000*16)
tot_err_rate_weight = np.sqrt(np.power(err_GRB/1e-4,2)+np.power(err_fraction_weight/np.mean(rate_weight),2))*np.mean(rate_weight)*1e-4*14000
print('total rate (weighted mean)', np.mean(rate_weight)*1e-4*14000)
print('error total rate (weighted mean)', tot_err_rate_weight)
'''
redred=[]
for id, elm in enumerate(redshift):
    if skyloc_error[id] < 1e3:
        redred.append(elm)

redred2=[]
for id, elm in enumerate(redshift):
    if skyloc_error[id] < 1e2:
        redred2.append(elm)

redred3=[]
for id, elm in enumerate(redshift):
    if skyloc_error[id] < 1e1:
        redred3.append(elm)

snr_red = []
for id, elm in enumerate(snr):
    if skyloc_error[id] < 1e4:
        snr_red.append(elm)
         
skylocerr_red = []
for elm in skyloc_error:
    if elm < 1e4:
        skylocerr_red.append(elm)

red_thetajn = []
for elm in theta_jn:
    if elm < 3:
        red_thetajn.append(elm)

weigths_red = []
for elm in theta_jn:
    if elm < 3:
        red_thetajn.append(elm)

ones0 = np.ones(len(redshift))
weight0 = (1e-4*14000)/10000*ones0
ones1 = np.ones(len(redred))
weight1 = (1e-4*14000)/10000*ones1
ones2 = np.ones(len(redred2))
weight2 = (1e-4*14000)/10000*ones2
ones3 = np.ones(len(redred3))
weight3 = (1e-4*14000)/10000*ones3
n0, bins0, _ = plt.hist(redshift, bins = 15, weights=weight0, histtype='step', lw=2, label='all')
n1,bins1, _ = plt.hist(redred, bins = 15, weights=weight1, histtype='step', lw=2, color='red', label=r'$\Delta\Omega < 1e3 deg^{2}$')
n2, bins2, _ = plt.hist(redred2, bins = 15, weights=weight2, histtype='step', lw=2, color='green', label=r'$\Delta\Omega < 1e2 deg^{2}$')
n3, bins3, _ = plt.hist(redred3, bins = 15, weights=weight3, histtype='step', lw=2, color='orange', label=r'$\Delta\Omega < 10 deg^{2}$')
plt.xlabel('redshift')
plt.ylabel('Number of of joint detections/year')
plt.legend()
plt.show()
print('total rate:', (6142/10000)*(1e-4*14000))

err_GRB = np.sqrt(26)/14000
err_GW = np.sqrt(1e4*(6142/10000)*(1-(6142/10000)))
# err_tot = 14000*np.sqrt(np.power(1e-4*err_GW/10000,2)+np.power(0.61*err_GRB,2))
err_tot = np.sqrt(np.power(1e-4/err_GRB,2)+np.power(0.61*err_GW/10000,2))*0.86
print(err_tot)
'''





'''
ones = np.ones(307.1)
counts, _ = np.histogram(redshift, bins=20)
weight = counts*1e-4*14000/1e4 
'''
# plt.hist(redred, bins = 20, histtype='step', color='red', label=r'$\Delta\Omega < 1e3 deg^{2}$')
# plt.hist(redred2, bins = 20, histtype='step', color='green', label=r'$\Delta\Omega < 1e2 deg^{2}$')
# plt.hist(redred3, bins = 20, histtype='step', color='orange', label=r'$\Delta\Omega < 10 deg^{2}$')
# plt.legend()
'''
plt.scatter(snr_red, skylocerr_red)
plt.ylabel('sky localisation error [deg]^2')
plt.xlabel('network_SNR')
'''
'''
if args.plot == 'hist':
    if args.par == 'dec': 
        plt.hist(np.sin(df[args.par]), bins=100, color='blue', alpha= 0.9, label='detected')
        plt.hist(np.sin(df_sim[args.par]), bins=100, color='blue', alpha= 0.6, label='simulated')
        plt.ylabel('counts')
        plt.xlabel(r'$\sin(Dec)$')
        plt.legend()
        if args.save == True:
            plt.savefig(f'hist_{args.grb_name}_{df[args.par]}.png', dpi=300)
        plt.show()
    elif args.par == 'theta_jn': 
        plt.hist(np.cos(df[args.par]), bins=100, color='blue', alpha= 0.9, label='detected')
        plt.hist(np.cos(df_sim[args.par]), bins=100, color='blue', alpha = 0.6, label='simulated')
        plt.ylabel('counts')
        plt.xlabel(r'$\cos(\theta_{JN})$')
        plt.legend()
        if args.save == True:
            plt.savefig(f'hist_{args.grb_name}_{df[args.par]}.png', dpi=300)
        plt.show()
    elif args.par == 'ET_SNR':
        plt.hist(df_sim[args.par], bins=50, color='red', alpha=0.4)
        plt.axvline(x=12, color='darkred', linestyle='dashed', label='SNR threshold')
        plt.title(f'GRB{args.grb_name}')
        # plt.text(10, 110, 'SNR = 9', fontsize=12, verticalalignment='top', color='darkred', alpha=1)
        plt.ylabel('counts')
        plt.xlabel(args.par)
        plt.legend()
        if args.save == True:
            plt.savefig(f'hist_{args.grb_name}_{df[args.par]}.png', dpi=300)
        plt.show()
    elif args.par == 'err_sky_location':
        plt.hist(df[args.par]*(180/np.pi)**2, bins=50, color='red', alpha=0.5)
        # plt.axvline(x=12, color='darkred', linestyle='dashed')
        # plt.text(10, 110, 'SNR = 9', fontsize=12, verticalalignment='top', color='darkred', alpha=1)
        plt.ylabel('counts')
        plt.xlabel(args.par)
        if args.save == True:
            plt.savefig(f'hist_{args.grb_name}_{df[args.par]}.png', dpi=300)
        plt.show()
    else: 
        plt.hist(df[args.par], bins=100, color='blue', alpha= 0.9, label='detected')
        plt.hist(df_sim[args.par], bins=100, color='blue', alpha= 0.6, label='simulated')
        # plt.axvline(x=12, color='darkred', linestyle='dashed')
        # plt.text(10, 110, 'SNR = 9', fontsize=12, verticalalignment='top', color='darkred', alpha=1)
        plt.xlabel(r'$\Lambda_{2}$')
        plt.legend()
        plt.ylabel('counts')
        plt.show()

if args.plot == 'scatter':
    if args.scat1 == 'err_sky_location':
        df[args.scat1] = df[args.scat1]*(180/np.pi)**2
    if args.scat2 == 'err_sky_location':
        df[args.scat2] = df[args.scat2]*(180/np.pi)**2
    if args.colormap == True: 
        if args.parcm == 'mchirp':
            mchirp = conversions.mchirp_from_mass1_mass2(df['mass_1'], df['mass_2'])
            plt.scatter(df[args.scat1], df[args.scat2], c=(mchirp), cmap='viridis',  s=10)
            cb = plt.colorbar()
            cb.set_label(r'M_chirp (M$_{\odot}$)')
            # plt.title(f'GRB{args.grb_name} - z={redshift[0]} - {args.numsim} sim - {len(df)} > SNR=12')
            plt.title(f'GRB{args.grb_name}')
            plt.xlabel(args.scat1)
            plt.ylabel(f'{args.scat2} [deg]$^{2}$')
            if args.save == True:
                plt.savefig(f'scatmchirp_{args.grb_name}.png', dpi=300)
            plt.show()
        if args.parcm == 'lambdatilde':
            lambda_tilde = conversions.lambda_tilde(df['mass_1'], df['mass_2'], df['lambda_1'], df['lambda_2'])
            plt.scatter(df[args.scat1], df[args.scat2], c=(lambda_tilde), cmap='plasma',  s=10)
            plt.title(f'GRB{args.grb_name}')
            cb = plt.colorbar()
            cb.set_label(r'$\tilde{\Lambda}$')
            plt.xlabel(args.scat1)
            plt.ylabel(f'{args.scat2} [deg]$^{2}$')
            if args.save == True:
                plt.savefig(f'scatlambdatilde_{args.grb_name}.png', dpi=300)
            plt.show()
        if args.parcm == 'snr':
            plt.scatter(np.cos(df[args.scat1]), df[args.scat2], c=(df['network_SNR']), cmap='inferno',  s=10)
            plt.title(f'GRB{args.grb_name}')
            cb = plt.colorbar()
            cb.set_label('SNR')
            plt.xlabel(r'$\cos(\theta_{JN})$')
            plt.ylabel(f'{args.scat2} [deg]$^{2}$')
            if args.save == True:
                plt.savefig(f'scatsnr_{args.grb_name}.png', dpi=300)
            plt.show()
        if args.parcm == 'theta_jn':
            plt.scatter(df[args.scat1], df[args.scat2], c=(np.cos(df['theta_jn'])), cmap='inferno',  s=10)
            plt.title(f'GRB{args.grb_name}')
            cb = plt.colorbar()
            cb.set_label(r'$\cos(\theta_{JN})$')
            # plt.gca().xaxis.set_major_formatter(lambda x, pos: f'${np.exp(x):.2g}$')
            # plt.gca().yaxis.set_major_formatter(lambda x, pos: f'${np.exp(x):.2g}$')
            plt.xlabel(args.scat1)
            plt.ylabel(f'{args.scat2} [deg]$^{2}$')
            if args.save == True:
                plt.savefig(f'scatthetajn_{args.grb_name}.png', dpi=300)
            plt.show()
    else:
        if args.scat1 == 'theta_jn':
            plt.scatter(np.cos(df[args.scat1]), df[args.scat2], s=2, c='darkred')
            plt.xlabel(r'$\cos(\theta_{JN})$')
            plt.ylabel(args.scat2)
            if args.save == True:
                plt.savefig(f'scat_{args.grb_name}.png', dpi=300)
            plt.show()
        elif args.scat2 == 'theta_jn': 
            plt.scatter(df[args.scat1], np.cos(df[args.scat2]), s=2, c='darkred')
            plt.xlabel(args.scat1)
            plt.ylabel(r'$\cos(\theta_{JN})$')
            if args.save == True:
                plt.savefig(f'scat_{args.grb_name}.png', dpi=300)
            plt.show()
        else: 
            plt.scatter(df[args.scat1], df[args.scat2], s=2, c='indianred')
            plt.xlabel(args.scat1)
            plt.ylabel(args.scat2)
            if args.save == True:
                plt.savefig(f'scat_{args.grb_name}.png', dpi=300)
            plt.show()
'''