import sys
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
import numpy as np
from pylab import *
from astropy.io import fits
from scipy.optimize import curve_fit

folder_RM = '/home/eli/Documentos/Astronomia/proyectos/DES-Y1/catalogs/'
folder_Mice = '/home/eli/Documentos/Astronomia/proyectos/MICEv1.0/catalogs/'

clusters = fits.open(folder_RM+'redmapper_y1a1_public_v6.4_catalog.fits')[1].data
members  = fits.open(folder_RM+'redmapper_y1a1_public_v6.4_members.fits')[1].data
halos    = fits.open(folder_Mice+'MICEv1.0_halo_cat.fits')[1].data

# MICE PARAMETERS

mg_h = halos.g_des_true
mr_h = halos.r_des_true
mi_h = halos.i_des_true
mz_h = halos.z_des_true
  
z_h = halos.z_v
lMH  = halos.log_m

# redMaPPer PARAMETERS

IDc = clusters.ID
z_c = clusters.Z_LAMBDA
mg_c = clusters.MAG_AUTO_G
mr_c = clusters.MAG_AUTO_R
mi_c = clusters.MAG_AUTO_I
mz_c = clusters.MAG_AUTO_Z
Lambda = clusters.LAMBDA

## mass-richness relation from McClintock et al. 2018 (arXiv:1805.00039v2)
# Table 4
M0  = 10**(14.489)
Fl  = 1.356
Gz  = -0.30
# Eq. 64
lMC = np.log10(M0 * (Lambda/40.)**Fl * ((1. + z_c)/(1. + 0.35))**Gz)

# DISTRIBUTIONS

mc = (z_c > 0.2)*(z_c < 0.65)
mh = (z_h > 0.2)*(z_h < 0.65)*(lMH > 14.05)

mz35_c = (z_c >= 0.35)
mz35_h = (z_h >= 0.35)

color_c         = mg_c - mr_c
color_c[mz35_c] = (mr_c - mi_c)[mz35_c]

color_h         = mg_h - mr_h
color_h[mz35_h] = (mr_h - mi_h)[mz35_h]


plt.hist(lMC[mc],np.linspace(14,15.5,20),histtype='step', label = 'redMaPPer')  
plt.hist(lMH[mh],np.linspace(14,15.5,20),histtype='step',label = 'MICE')
plt.legend()
plt.xlabel('log(M)')  

plt.figure()
plt.hist(z_c[mc],np.linspace(0.2,0.65,20),histtype='step', label = 'redMaPPer')  
plt.hist(z_h[mh],np.linspace(0.2,0.65,20),histtype='step',label = 'MICE')
plt.legend()
plt.xlabel('z')  


for j in range(45):

    zmin = 0.2+0.01*j
    print(str('%.2f' % zmin))
    mz_c = (z_c >= zmin)*(z_c < (zmin + 0.01))
    mz_h = (z_h >= zmin)*(z_h < (zmin + 0.01))
   
    
    plt.figure()
    plt.plot(mr_h[mz_h],color_h[mz_h],'C7.',label = 'MICE')
    plt.plot(mr_c[mz_c],color_c[mz_c],'.', label = 'redMaPPer')  
    plt.plot(mr_h[mz_h*mh],color_h[mz_h*mh],'.',label = 'MICE - log(MH) > 14')
    plt.legend()
    plt.xlabel('mag_r')
    if zmin < 0.35:
        plt.ylabel('g-r')
    else:
        plt.ylabel('r-i')
    plt.savefig('/home/eli/Documentos/Astronomia/proyectos/MICEv1.0/RS_test_plots/rs_z'+str('%.2f' % zmin)+'_true.png')


def gauss(x, a,mu,sigma):
    return a*np.exp(-(x-mu)**2/(2.*sigma**2))

mz_h = (z_h > 0.2)*(z_h < 0.21)
mz_c = (z_c > 0.2)*(z_c < 0.21)
    
color_hf = color_h[mz_h]
n,c      = np.histogram(color_hf,15)      
c        = (c+(c[1]-c[0])*0.5)[:-1]

err   = np.ones(len(c))

fit_gauss = curve_fit(gauss,c,n,sigma=err,absolute_sigma=True)
a         = fit_gauss[0][0]
mu        = fit_gauss[0][1]
sigma     = abs(fit_gauss[0][2])

color_cf = color_c[mz_c]
n_c,c_c    = np.histogram(color_cf,15)      
c_c      = (c_c+(c_c[1]-c_c[0])*0.5)[:-1]

plt.figure()
plt.plot(mr_h[mz_h],color_hf,'C7.')
plt.plot(mr_c[mz_c*mc],color_cf,'.')



plt.figure()
plt.plot(c,n)
plt.plot(c_c,n_c)
plt.plot(c,gauss(c,a,mu,sigma))

