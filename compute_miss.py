import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
from profiles_fit import *
from models_profiles import *
from astropy.constants import G,c,M_sun, pc

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)


    
folder = '/mnt/clemente/lensing/HALO_SHAPE/MICE_v2.0/catalogs/'


p_name = 'profile_ebin_142.fits'
m_name = 'mapa_bin_142.fits'

profile = fits.open(folder+p_name)
mapa = fits.open(folder+m_name)[1].data
fitmiss = fits.open(folder+'fitresults_mono_Rayleigh_0_2500_profile_ebin_142.fits')[0].header

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

print('Computing S...')
S = Sigma_NFW_miss_parallel(r,zmean,10**lM200_miss,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=56)
print('Computing DS...')
DS = Delta_Sigma_NFW_miss_parallel(r,zmean,10**lM200_miss,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=56)
print('Computing gt,gx...')
gt,gx = GAMMA_components_miss_parallel(r,zmean,10**lM200_miss,ellip=e,s_off = soff,c200 = c200_miss, P_Roff= Rayleigh, cosmo=cosmo,ncores=56)

j = np.argsort(r)

table = [fits.Column(name='rmpc', format='E', array=r[j]),
            fits.Column(name='xmpc', format='E', array=x[j]),
            fits.Column(name='ympc', format='E', array=y[j]),
            fits.Column(name='S', format='E', array=S[j]),
            fits.Column(name='DS', format='E', array=DS[j]),
            fits.Column(name='Gt', format='E', array=gt[j]),
            fits.Column(name='Gx', format='E', array=gx[j])]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))


h = fits.Header()

primary_hdu = fits.PrimaryHDU(header=h)

hdul = fits.HDUList([primary_hdu, tbhdu])

hdul.writeto(folder+'mapa_bin_142_miss.fits',overwrite=True)
