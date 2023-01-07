import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
from pycbc import conversions
import sys
from argparse import ArgumentParser

parser = ArgumentParser(description='Plot for BNS-GRB analysis with GWFish in the ET context')

parser.add_argument("-grb", "--grb", type=str,
                    dest="grb_name", help="GRB name")
parser.add_argument("-plot", "--plot", type=str,
                    dest="plot", help="Choose between scatter or hist")
parser.add_argument("-par", "--par", type=str,
                    dest="par", help="Choose a GW parameter") 
parser.add_argument("-parscat1", "--parscat1", type=str,
                    dest="scat1", help="First parameter for the scatter plot")  
parser.add_argument("-parscat2", "--parscat2", type=str,
                    dest="scat2", help="Second parameter for the scatter plot")  
parser.add_argument("-colormap", "--colormap", type=bool, default=False, 
                    dest="colormap", help="True if you want to visualise a color map based on the parcm")
parser.add_argument("-parcm", "--parcm", type=str,
                    dest="parcm", help="Parameter to use for the color map")
parser.add_argument("-save", "--save", type=bool,
                    dest="save", help="Save or not the figure") 
parser.add_argument("-sim", "--sim", type=bool,
                    dest="sim", help="If sim==True use the simulated files, otherwise the deteted")
parser.add_argument("-numsim", "--numsim", type=str,
                    dest="numsim", help="Number of simulations")

args = parser.parse_args()

df = pd.read_csv(f'Errors_ET_BNS_{args.grb_name}_SNR12.0.txt', sep=' ')
# _{args.numsim}
df_sim = pd.read_csv(f'Signals_BNS_{args.grb_name}.txt', sep=' ')

skyloc_error = df['err_sky_location']* (180/np.pi)**2
redshift = df['redshift']
ra = df['err_sky_location']
mass_1 = df['mass_1']
mass_2 = df['mass_2']
theta_jn = df['theta_jn']
ra = df['ra']
dec = df['dec']
'''
lambda1 = df['lambda_1']
lambda2 = df['lambda_2']
'''
snr = df['network_SNR']

print(len(df)/100000)

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
n1,bins1, _ = plt.hist(redred, bins = 15, weights=weight1, histtype='step', lw=2, color='red', label=r'$\Delta\Omega < 10^{3} deg^{2}$')
n2, bins2, _ = plt.hist(redred2, bins = 15, weights=weight2, histtype='step', lw=2, color='green', label=r'$\Delta\Omega < 10^{2} deg^{2}$')
n3, bins3, _ = plt.hist(redred3, bins = 15, weights=weight3, histtype='step', lw=2, color='orange', label=r'$\Delta\Omega < 10 deg^{2}$')
plt.xlabel('redshift', fontsize=20)
plt.ylabel('Number of of joint detections/year', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend()
plt.show()
print('total rate:', (6142/10000)*(1e-4)*14000)

# sbagliato
# err_GRB = np.sqrt(26)/14000
# err_GW = np.sqrt(1e4*(6142/10000)*(1-(6142/10000)))
# err_tot = 14000*np.sqrt(np.power(1e-4*err_GW/10000,2)+np.power(0.61*err_GRB,2))
# err_tot = np.sqrt(np.power(1e-4/err_GRB,2)+np.power(0.61*err_GW/10000,2))*0.86
# print(err_tot)

# giusto
err_GRB_zdist = np.sqrt(26)/(14000*16)
err_GW_zdist = np.sqrt(1e4*(6142/10000)*(1-(6142/10000)))
err_tot = np.sqrt(np.power(err_GRB_zdist/1e-4,2)+np.power(err_GW_zdist/(10000*0.61),2))*0.86
print('total rate:', (6142/10000)*(1e-4)*14000)
print('err total rate:', err_tot)





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
if args.plot == 'hist':
    if args.par == 'dec': 
        plt.hist(np.sin(df[args.par]), bins=100, color='blue', alpha= 0.9, label='detected')
        plt.hist(np.sin(df_sim[args.par]), bins=100, color='blue', alpha= 0.6, label='simulated')
        plt.ylabel('counts', fontsize=25)
        plt.xlabel(r'$\sin(Dec)$', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend()
        if args.save == True:
            plt.savefig(f'hist_{args.grb_name}_{df[args.par]}.png', dpi=300)
        plt.show()
    elif args.par == 'theta_jn': 
        plt.hist(np.cos(df[args.par]), bins=100, color='blue', alpha= 0.9, label='detected')
        plt.hist(np.cos(df_sim[args.par]), bins=100, color='blue', alpha = 0.6, label='simulated')
        plt.ylabel('counts', fontsize=25)
        plt.xlabel(r'$\cos(\theta_{JN})$', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
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
        plt.xlabel('RA [rad]', fontsize=20)
        plt.ylabel('counts', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        plt.legend()

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