import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
from profiles_fit import *
from fit_profiles_curvefit import *
from multipoles_shear import *
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


print(p_name)

# '''
h   = profile[0].header
p   = profile[1].data
cov = profile[2].data

CovDS  = cov.COV_ST.reshape(len(p),len(p))

cosmo  = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)


zmean = h['z_mean']
q  = h['q2d_mean']
qr = h['q2dr_mean']

e = (1-q)/(1+q)
er = (1-qr)/(1+qr)

nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo,True)

y = mapa.ympc
x = mapa.xmpc

theta  = np.arctan2(y,x)
r = np.sqrt(x**2 + y**2)

out = multipole_shear_parallel(r,M200=nfw.M200,ellip=0.25,z=0.2,
							 h=h['hcosmo'],misscentred=True,
							 s_off=0.2,components = ['t0','t','tcos','xsin'],
							 verbose = True, Yanmiss = False, ncores=35)

table = [fits.Column(name='xmpc', format='E', array=x),
            fits.Column(name='ympc', format='E', array=y),
            fits.Column(name='Gt0', format='E', array=out['Gt0']),
            fits.Column(name='Gt2', format='E', array=out['Gt2']),
            fits.Column(name='Gx2', format='E', array=out['Gx2']),
            fits.Column(name='Gt0_off', format='E', array=out['Gt0_off']),
            fits.Column(name='Gt_off', format='E', array=out['Gt_off']),
            fits.Column(name='Gt_off_cos', format='E', array=out['Gt_off_cos']),
            fits.Column(name='Gx_off_sin', format='E', array=out['Gx_off_sin'])]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))


h = fits.Header()

primary_hdu = fits.PrimaryHDU(header=h)

hdul = fits.HDUList([primary_hdu, tbhdu])

hdul.writeto(folder+'mapa_bin_142_missx.fits',overwrite=True)
