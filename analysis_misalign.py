import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from astropy.constants import G,c,M_sun, pc
from fit_models_colossus import Delta_Sigma_fit

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)
import corner
folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'


def qratios(samp):

    RIN = 250
    ROUT = 5000
    
    p_name = 'profile_'+samp+'.fits'  
    
    q1h_2h   = fits.open(folder+'fitresults_2h_'+str(int(RIN))+'_'+str(int(ROUT))+'_profile_'+samp+'.fits')[0].header['q']
    q1h_2h2q = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_profile_'+samp+'.fits')[0].header['q']
    q1h_q    = fits.open(folder+'fitresults_onlyq_'+str(int(RIN))+'_2000_profile_'+samp+'.fits')[0].header['q']
    q1h_miss = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_profile_'+samp+'_misalign.fits')[0].header['q']

    q2h_2h2q = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_profile_'+samp+'.fits')[0].header['q2h']
    q2h_miss = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_profile_'+samp+'_misalign.fits')[0].header['q2h']

    q1hr_2h   = fits.open(folder+'fitresults_2h_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_profile_'+samp+'.fits')[0].header['q']
    q1hr_2h2q = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_profile_'+samp+'.fits')[0].header['q']
    q1hr_q    = fits.open(folder+'fitresults_onlyq_'+str(int(RIN))+'_2000_reduced_profile_'+samp+'.fits')[0].header['q']
    q1hr_miss = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_profile_'+samp+'_misalign.fits')[0].header['q']

    q2hr_2h2q = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_profile_'+samp+'.fits')[0].header['q2h']
    q2hr_miss = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_profile_'+samp+'_misalign.fits')[0].header['q2h']
    
    print('q1h ratio standard')
    print(q1h_q/q1h_miss,q1h_2h/q1h_miss,q1h_2h2q/q1h_miss)
    print('q2h ratio standard')
    print(q2h_2h2q/q2h_miss)
    print('q1h ratio reduced')
    print(q1hr_q/q1hr_miss,q1hr_2h/q1hr_miss,q1hr_2h2q/q1hr_miss)
    print('q2h ratio reduced')
    print(q2hr_2h2q/q2hr_miss)
