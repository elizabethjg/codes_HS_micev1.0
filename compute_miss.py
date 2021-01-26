import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from profiles_fit import *
from fit_profiles_curvefit import *
from astropy.constants import G,c,M_sun, pc

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)


    
folder = '/home/elizabeth/MICEv2.0/'


p_name = 'profiles/profile_ebin_142.fits'
m_name = 'mapas/mapa_bin_144.fits'

profile = fits.open(folder+p_name)
mapa = fits.open(folder+m_name)[1].data
fitmiss = fits.open(folder+'profiles/fitresults_mono_Rayleigh_0_2500_profile_ebin_142.fits')[0].header

print(p_name)

# '''
h   = profile[0].header
p   = profile[1].data
cov = profile[2].data

lM200_miss = fitmiss['lm200']
c200_miss = fitmiss['c200']
soff = fitmiss['soff']

cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)


zmean = h['z_mean']
q  = h['q2d_mean']
qr = h['q2dr_mean']

e = (1-q)/(1+q)
er = (1-qr)/(1+qr)

y = mapa.ympc
x = mapa.xmpc

theta  = np.arctan2(y,x)
r = np.sqrt(x**2 + y**2)


DS = Delta_Sigma_NFW_miss_parallel(r,zmean,10**lM200_miss,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=40)
gt,gx = GAMMA_components_miss_parallel(r,zmean,10**lM200_miss,ellip=e,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=40)

table = [fits.Column(name='xmpc', format='E', array=x),
            fits.Column(name='ympc', format='E', array=y),
            fits.Column(name='DS', format='E', array=DS),
            fits.Column(name='Gt', format='E', array=gt),
            fits.Column(name='Gx', format='E', array=gx)]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))


h = fits.Header()

primary_hdu = fits.PrimaryHDU(header=h)

hdul = fits.HDUList([primary_hdu, tbhdu])

hdul.writeto(folder+'mapas/mapa_bin_142_miss.fits',overwrite=True)
