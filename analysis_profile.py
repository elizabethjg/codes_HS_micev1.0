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

h = 1.0
cosmo = LambdaCDM(H0=100*h, Om0=0.3, Ode0=0.7)


p_name = 'profile_pru.fits'
profile = fits.open('../profiles/'+p_name)

h = profile[1].header
p = profile[1].data

zmean = h['z_mean']

H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
roc_mpc  = roc*((pc*1.0e6)**3.0)

nfw        = NFW_stack_fit(p.Rp,p.DSigma_T,p.error_DSigma_T,zmean,roc,h=1)

M200_NFW   = (800.0*np.pi*roc_mpc*(nfw[0]**3))/(3.0*Msun)

nfw2 = Delta_Sigma_fit(p.Rp,p.DSigma_T,p.error_DSigma_T,zmean,cosmo)

ndots = p.shape[0]

C_matrix = np.zeros(ndots,ndots)

for i in range(ndots):
    for j in range(ndots):
