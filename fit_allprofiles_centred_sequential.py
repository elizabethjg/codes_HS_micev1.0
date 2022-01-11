import sys
sys.path.append('/mnt/projects/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from models_profiles import *
from multiprocessing import Pool
from multiprocessing import Process
import argparse
from astropy.constants import G,c,M_sun,pc
import emcee
from models_profiles import *
from astropy.cosmology import LambdaCDM
from fit_profiles_curvefit import *
# import corner
import os
from colossus.cosmology import cosmology  
params = {'flat': True, 'H0': 70.0, 'Om0': 0.25, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.95}
cosmology.addCosmology('MICE', params)
cosmo = cosmology.setCosmology('MICE')
from colossus.halo import concentration
cmodel = 'diemer19'


'''
folder = '/home/elizabeth/Documentos/proyectos/HALO-SHAPE/MICEv2.0/profiles/'
cont = False
file_name = 'profile_bin_140.fits'
angle = 'standard'
ncores = 3
nit = 250
RIN = 0.
ROUT =5000.
# '''

parser = argparse.ArgumentParser()
#parser.add_argument('-folder', action='store', dest='folder', default='/mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/catalogs/')
parser.add_argument('-folder', action='store', dest='folder', default='/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles/')
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ang', action='store', dest='angle', default='standard')
parser.add_argument('-ncores', action='store', dest='ncores', default=2)
parser.add_argument('-RIN', action='store', dest='RIN', default=0)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=2500)
parser.add_argument('-nit', action='store', dest='nit', default=250)
parser.add_argument('-continue', action='store', dest='cont', default='False')
parser.add_argument('-components', action='store', dest='comp', default='all')
args = parser.parse_args()


folder    = args.folder
file_name = args.file_name
angle     = args.angle

if 'True' in args.cont:
	cont      = True
elif 'False' in args.cont:
	cont      = False

	
nit       = int(args.nit)
ncores    = args.ncores
ncores    = int(ncores)
RIN       = float(args.RIN)
ROUT      = float(args.ROUT)
if angle == 'standard':
    ang = ''
elif angle == 'reduced':
    ang = '_reduced'

if args.comp == 'all':
    outfile     = 'fitresults_onlyq_'+str(int(RIN))+'_'+str(int(ROUT))+ang+'_'+file_name
else:
    outfile     = 'fitresults_onlyq_'+args.comp+'_'+str(int(RIN))+'_'+str(int(ROUT))+ang+'_'+file_name
backup      = folder+'backup_'+outfile
plot_folder = folder+'plots_mcmc/'


print('fitting profiles')
print(folder)
print(file_name)
print(angle)
print('ncores = ',ncores)
print('RIN ',RIN)
print('ROUT ',ROUT)
print('nit', nit)
# print('continue',cont)
print('outfile',outfile)
print('fitting components ',args.comp)


profile = fits.open(folder+file_name)
h       = profile[0].header
p       = profile[1].data
cov     = profile[2].data
zmean   = h['Z_MEAN'] 

CovDS  = cov.COV_ST.reshape(p.shape[0],p.shape[0])

cosmo_as = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
nfw     = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo_as,True)
   
lM200 = np.log10(nfw.M200)
c200  = nfw.c200



def log_likelihood(data_model, R, profiles, iCOV):
    
    q = data_model
    
    e = (1.-q)/(1.+q)

    gt, gx = profiles
    iCgt, iCgx = iCOV 

    GT,GX   = GAMMA_components(R,zmean,ellip=e,M200 = 10**lM200,c200=c200,cosmo=cosmo_as)

    L_GT = -np.dot((gt-GT),np.dot(iCgt,(gt-GT)))/2.0
    L_GX = -np.dot((gx-GX),np.dot(iCgx,(gx-GX)))/2.0
    
    if args.comp == 'all':
        L = L_GT + L_GX
    elif args.comp == 'tangential':
        L = L_GT
    elif args.comp == 'cross':
        L = L_GX
    
    return L
    

def log_probability(data_model, R, profiles, iCOV):
    
    q = data_model
    
    if 0.2 < q < 1.0:
        return log_likelihood(data_model, R, profiles, iCOV)
        
    return -np.inf

# initializing

pos = np.array([np.random.uniform(0.2,0.9,15)]).T

nwalkers, ndim = pos.shape

#-------------------
# running emcee

maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))

mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]

CovGT  = cov['COV_GT'+ang].reshape(len(p.Rp),len(p.Rp))[mr]
CovGX  = cov['COV_GX'+ang].reshape(len(p.Rp),len(p.Rp))[mr]

p  = p[maskr]

GT  = p['GAMMA_Tcos'+ang]
GX  = p['GAMMA_Xsin'+ang]

CovGT  = CovGT.reshape(maskr.sum(),maskr.sum())
CovGX  = CovGX.reshape(maskr.sum(),maskr.sum())

profiles = [GT,GX]
iCov     = [np.linalg.inv(CovGT),np.linalg.inv(CovGX)]

t1 = time.time()
# backend = emcee.backends.HDFBackend(backup)
# if not cont:
    # backend.reset(nwalkers, ndim)


pool = Pool(processes=(ncores))    
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(p.Rp,profiles,iCov),pool = pool)
				


if cont:                                
    sampler.run_mcmc(None, nit, progress=True)
else:
    sampler.run_mcmc(pos, nit, progress=True)
pool.terminate()
    
    
print('TOTAL TIME FIT')    
print((time.time()-t1)/60.)

#-------------------
# saving mcmc out

mcmc_out = sampler.get_chain(flat=True)

table = [fits.Column(name='q', format='E', array=mcmc_out)]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))

lM   = lM200
q    = np.percentile(mcmc_out[1500:], [16, 50, 84])
c200 = c200



h = fits.Header()
h.append(('lM200',np.round(lM200,4)))
h.append(('c200',np.round(c200,4)))

h.append(('q',np.round(q[1],4)))
h.append(('eqM',np.round(np.diff(q)[0],4)))
h.append(('eqm',np.round(np.diff(q)[1],4)))


primary_hdu = fits.PrimaryHDU(header=h)

hdul = fits.HDUList([primary_hdu, tbhdu])

hdul.writeto(folder+outfile,overwrite=True)

print('SAVED FILE')

