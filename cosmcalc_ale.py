import sys
from math import *
import matplotlib.pyplot as plt

H0 =  69.6
WM = 0.286
WV = 0.714

#H0 = 75                         # Hubble constant
#WM = 0.3                        # Omega(matter)
#WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

WR = 0.        # Omega(radiation)
WK = 0.        # Omega curvaturve = 1-Omega(total)
c = 299792.458 # velocity of light in km/sec
Tyr = 977.8    # coefficent for converting 1/H into Gyr
DTT = 0.5      # time from z to now in units of 1/H0
DTT_Gyr = 0.0  # value of DTT in Gyr
age = 0.5      # age of Universe in units of 1/H0
age_Gyr = 0.0  # value of age in Gyr
zage = 0.1     # age of Universe at redshift z in units of 1/H0
zage_Gyr = 0.0 # value of zage in Gyr
DCMR = 0.0     # comoving radial distance in units of c/H0
DCMR_Mpc = 0.0 
DCMR_Gyr = 0.0
DA = 0.0       # angular size distance
DA_Mpc = 0.0
DA_Gyr = 0.0
kpc_DA = 0.0
DL = 0.0       # luminosity distance
DL_Mpc = 0.0
DL_Gyr = 0.0   # DL in units of billions of light years
V_Gpc = 0.0
a = 1.0        # 1/(1+z), the scale factor of the Universe
az = 0.5       # 1/(1+z(object))

h = H0/100.
WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
WK = 1-WM-WR-WV
age = 0.
n=1000         # number of points in integrals


# z = [2.3, 1.1, 0.9, 0.7, 0.1, 0.4, 0.1, 0.4, 0.8, 0.1, 0.08, 1.1, 0.9, 0.4, 1.5, 0.3, 1.6, 2.1, 0.4, 3.9, 1.3, 4, 
    # 0.7, 0.4, 0.3, 1.4]
# DL_Mpc_list = []

def compute_dl_Mpc(z):

    az = 1.0/(1+1.0*z)
    DTT = 0.5
    DCMR = 0.0     
    
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

    DTT = (1.-az)*DTT/n
    DCMR = (1.-az)*DCMR/n
    
    ratio = 1.00
    x = sqrt(abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio =  0.5*(exp(x)-exp(-x))/x 
        else:
            ratio = sin(x)/x
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
    
    DCMT = ratio*DCMR
    DA = az*DCMT
    DL = DA/(az*az)
    DL_Mpc = (c/H0)*DL
    # DL_Mpc_list.append(DL_Mpc)
    return DL_Mpc
