import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
sys.path.append('/mnt/projects/lensing/lens_codes_v3.7')
from models_profiles import *
from astropy.constants import G,c,M_sun, pc

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)


    
folder = '/mnt/projects/lensing/HALO_SHAPE/MICEv2.0/'


p_name = 'profiles/profile_ebin_134.fits'
m_name = 'mapas/mapa_bin_134.fits'

profile = fits.open(folder+p_name)
mapa = fits.open(folder+m_name)[1].data
# fitmiss = fits.open(folder+'profiles/fitresults_fullmodel_allmis_0_2000_profile_ebin_134.fits')[0].header
fitmiss = fits.open(folder+'profiles/fitresults_mono_Rayleigh_0_2500_profile_ebin_134.fits')[0].header

print(p_name)

# '''
h   = profile[0].header
p   = profile[1].data
cov = profile[2].data

lM200_miss = fitmiss['lm200']
c200_miss = fitmiss['c200']
soff = 0.08

cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)


zmean = h['z_mean']
q  = fitmiss['c200']
qr = h['q2dr_mean']

e = (1-q)/(1+q)
er = (1-qr)/(1+qr)

y = mapa.ympc
x = mapa.xmpc

theta  = np.arctan2(y,x)
r = np.sqrt(x**2 + y**2)

R = r*np.sqrt(q*(np.cos(theta))**2 + (np.sin(theta))**2 / q)

print('Computing S0...')
S0 = Sigma_NFW_miss_parallel(r,zmean,10**lM200_miss,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=56)
print('Computing DS0...')
DS0 = Delta_Sigma_NFW_miss_parallel(r,zmean,10**lM200_miss,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=56)
print('Computing S...')
S = Sigma_NFW_miss_parallel(R,zmean,10**lM200_miss,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=56)
print('Computing DS...')
DS = Delta_Sigma_NFW_miss_parallel(R,zmean,10**lM200_miss,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=56)
print('Computing gt,gx...')
gt,gx,s2 = GAMMA_components_miss_parallel(r,zmean,10**lM200_miss,ellip=e,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo, return_S2 = True, ncores=56)


table = [fits.Column(name='rmpc', format='E', array=r),
            fits.Column(name='xmpc', format='E', array=x),
            fits.Column(name='ympc', format='E', array=y),
            fits.Column(name='S0', format='E', array=S0),
            fits.Column(name='DS0', format='E', array=DS0),
            fits.Column(name='S', format='E', array=S),
            fits.Column(name='DS', format='E', array=DS),
            fits.Column(name='Gt', format='E', array=gt),
            fits.Column(name='S2', format='E', array=s2),
            fits.Column(name='Gx', format='E', array=gx)]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))


h = fits.Header()

primary_hdu = fits.PrimaryHDU(header=h)

hdul = fits.HDUList([primary_hdu, tbhdu])

hdul.writeto(folder+'mapas/mapa_bin_134_miss.fits',overwrite=True)
