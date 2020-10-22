import sys
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
import numpy as np
from pylab import *
from astropy.io import fits
from scipy.optimize import curve_fit
from astropy.cosmology import LambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u 
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

"""
MICE QUERY

SELECT `id`, `ra`, `dec`, `z_v`, `g_des_true`, `r_des_true`, 
`i_des_true`, `z_des_true`, `g_des_realization`, 
`r_des_realization`, `i_des_realization`, `z_des_realization`, 
`log_m`, `nsat` FROM micecatv1_0_hpix 
WHERE `dec` < 56.5 AND `log_m` > 13 AND `flag` == 0 

"""
# ----------------------
# READ CATALOGS
# ----------------------

folder_RM = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/DES-Y1/catalogs/'
folder_Mice = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICEv1.0/catalogs/'

clusters = fits.open(folder_RM+'redmapper_y1a1_public_v6.4_catalog.fits')[1].data
halos    = fits.open(folder_Mice+'MICEv1.0_halo_cat.fits')[1].data

# ----------------------
# SELECT HALO MICE IN SPECIFIC REGION
# ----------------------

mregion = ((halos.dec < 1.5) + (halos.dec > 40.))*(halos.ra < 80.)

plt.figure()
plt.plot(halos.ra,halos.dec,'C7,',alpha=0.1)
plt.plot(halos.ra[mregion],halos.dec[mregion],'C3,',alpha=0.1)
plt.plot(clusters.RA,clusters.DEC,'C0,',alpha=1.)

halos = halos[mregion]

# ----------------------
# MICE PARAMETERS
# ----------------------

z_h = halos.z_v

mg_h = halos.g_des_true - 0.8 * (np.arctan(1.5 * z_h) - 0.1489)
mr_h = halos.r_des_true - 0.8 * (np.arctan(1.5 * z_h) - 0.1489)
mi_h = halos.i_des_true - 0.8 * (np.arctan(1.5 * z_h) - 0.1489)
mz_h = halos.z_des_true - 0.8 * (np.arctan(1.5 * z_h) - 0.1489)
  
lMH  = np.log10(10**(halos.log_m) * 0.7)

Dl_h = np.array(cosmo.luminosity_distance(z_h).value)*1.e6 

Mg_h=mg_h+5.0-5.0*np.log10(Dl_h)
Mr_h=mr_h+5.0-5.0*np.log10(Dl_h)
Mi_h=mi_h+5.0-5.0*np.log10(Dl_h)
Mz_h=mz_h+5.0-5.0*np.log10(Dl_h)

mh = (z_h > 0.2)*(z_h < 0.65)*(lMH > 13.77) # SELECT HALOS TO MIMIC redMaPPer

Dcosmo = cosmo.comoving_distance(z_h[mh])
catalog = SkyCoord(ra=halos.ra[mh]*u.degree, dec=halos.dec[mh]*u.degree, distance=Dcosmo)
idx, Rprox_h2, Rprox_h = catalog.match_to_catalog_3d(catalog, nthneighbor=2)
Rprox_h = np.array(Rprox_h.value)
Rprox_h2 = np.array(Rprox_h2.value)

# ----------------------
# redMaPPer PARAMETERS
# ----------------------

IDc = clusters.ID
z_c = clusters.Z_LAMBDA
mg_c = clusters.MAG_AUTO_G
mr_c = clusters.MAG_AUTO_R
mi_c = clusters.MAG_AUTO_I
mz_c = clusters.MAG_AUTO_Z
Lambda = clusters.LAMBDA

Dl_c = np.array(cosmo.luminosity_distance(z_c).value)*1.e6 

Mg_c = mg_c+5.0-5.0*np.log10(Dl_c)
Mr_c = mr_c+5.0-5.0*np.log10(Dl_c)
Mi_c = mi_c+5.0-5.0*np.log10(Dl_c)
Mz_c = mz_c+5.0-5.0*np.log10(Dl_c)

mc = (z_c > 0.2)*(z_c < 0.65)

Dcosmo = cosmo.comoving_distance(z_c[mc])
catalog = SkyCoord(ra=clusters.RA[mc]*u.degree, dec=clusters.DEC[mc]*u.degree, distance=Dcosmo)
idx, Rprox_c2, Rprox_c = catalog.match_to_catalog_3d(catalog, nthneighbor=2)
Rprox_c = np.array(Rprox_c.value)
Rprox_c2 = np.array(Rprox_c2.value)

## mass-richness relation from McClintock et al. 2018 (arXiv:1805.00039v2)
# Table 4
M0  = 10**(14.489)
Fl  = 1.356
Gz  = -0.30
# Eq. 64
lMC = np.log10(0.76 * 0.7 * M0 * (Lambda/40.)**Fl * ((1. + z_c)/(1. + 0.35))**Gz)

# ----------------------
# DISTRIBUTIONS
# ----------------------

mz35_c = (z_c >= 0.35)
mz35_h = (z_h >= 0.35)

color_c         = mg_c - mr_c
color_c[mz35_c] = (mr_c - mi_c)[mz35_c]

color_h         = mg_h - mr_h
color_h[mz35_h] = (mr_h - mi_h)[mz35_h]

plt.figure()
plt.hist(Rprox_c,np.linspace(0,100,50),histtype='step', label = 'redMaPPer')  
plt.hist(Rprox_h,np.linspace(0,100,50),histtype='step',label = 'MICE')
plt.legend(loc=1)
plt.xlabel('Rprox')  

plt.figure()
plt.hist(lMC[mc],np.linspace(13.8,15.5,20),histtype='step', label = 'redMaPPer')  
plt.hist(lMH[mh],np.linspace(13.8,15.5,20),histtype='step',label = 'MICE')
plt.legend(loc=1)
plt.xlabel('log(M)')  

plt.figure()
plt.hist(z_c[mc],np.linspace(0.2,0.65,20),histtype='step', label = 'redMaPPer')  
plt.hist(z_h[mh],np.linspace(0.2,0.65,20),histtype='step',label = 'MICE')
plt.legend(loc=2)
plt.xlabel('z')  

# ----------------------
# PLOT RED-SEQUENCE IN REDSHIFT BINS
# ----------------------

for j in range(45):

    zmin = 0.2+0.01*j
    print(str('%.2f' % zmin))
    mz_c = (z_c >= zmin)*(z_c < (zmin + 0.01))
    mz_h = (z_h >= zmin)*(z_h < (zmin + 0.01))
   
    
    plt.figure()
    plt.plot(Mr_h[mz_h],color_h[mz_h],'C7.',label = 'MICE')
    plt.plot(Mr_c[mz_c],color_c[mz_c],'.', label = 'redMaPPer')  
    plt.plot(Mr_h[mz_h*mh],color_h[mz_h*mh],'.',label = 'MICE - log(MH) > 14')
    plt.legend()
    plt.xlabel('mag_r')
    if zmin < 0.35:
        plt.ylabel('g-r')
    else:
        plt.ylabel('r-i')
    plt.savefig('/home/eli/Documentos/Astronomia/proyectos/MICEv1.0/RS_test_plots/rs_z'+str('%.2f' % zmin)+'_true.png')

# ----------------------
# FIT RED-SEQUENCE
# ----------------------

def gauss(x, a,mu,sigma):
    return a*np.exp(-(x-mu)**2/(2.*sigma**2))

mz_h = (z_h > 0.25)*(z_h < 0.26)
mz_c = (z_c > 0.25)*(z_c < 0.26)
    
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

mask   = (color_hf < (mu+sigma))*(color_hf > (mu-sigma))

plt.figure()
plt.plot(mr_h[mz_h],color_hf,'C7.')
plt.plot(mr_c[mz_c*mc],color_cf,'.')
plt.plot(mr_h[mz_h][mask],color_hf[mask],'.')

plt.figure()
plt.plot(c,n)
plt.plot(c_c,n_c)
plt.plot(c,gauss(c,a,mu,sigma))

