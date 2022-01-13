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
folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles/'
cont = False
file_name = 'profile_LM_Lz_relaxed.fits'
angle = 'standard'
ncores = 32
nit = 10
RIN = 250.
ROUT =2500.
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
    outfile     = 'fitresults_2h_'+str(int(RIN))+'_'+str(int(ROUT))+ang+'_'+file_name
else:
    outfile     = 'fitresults_2h_'+args.comp+'_'+str(int(RIN))+'_'+str(int(ROUT))+ang+'_'+file_name
backup      = folder+'backup_'+outfile



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

# extracting data from profile
profile = fits.open(folder+file_name)
h       = profile[0].header
p       = profile[1].data
cov     = profile[2].data
zmean   = h['Z_MEAN'] 

maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))

mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]

CovDS  = cov['COV_ST'].reshape(len(p.Rp),len(p.Rp))[mr]
CovGT  = cov['COV_GT'+ang].reshape(len(p.Rp),len(p.Rp))[mr]
CovGX  = cov['COV_GX'+ang].reshape(len(p.Rp),len(p.Rp))[mr]

p  = p[maskr]


DSt = p.DSigma_T
GT  = p['GAMMA_Tcos'+ang]
GX  = p['GAMMA_Xsin'+ang]


CovDS  = CovDS.reshape(maskr.sum(),maskr.sum())
CovGT  = CovGT.reshape(maskr.sum(),maskr.sum())
CovGX  = CovGX.reshape(maskr.sum(),maskr.sum())

profiles = [GT,GX]
iCov     = [np.linalg.inv(CovGT),np.linalg.inv(CovGX)]
iCds     =  np.linalg.inv(CovDS)
  
# First running for DS

def log_likelihood_DS(data_model, R, ds, iCds):
    
    lM200, c200 = data_model
    
    DS   = Delta_Sigma_NFW_2h(R,zmean,M200 = 10**lM200,c200=c200,cosmo_params=params)

    L_DS = -np.dot((ds-DS),np.dot(iCds,(ds-DS)))/2.0
        
    return L_DS
    

def log_probability_DS(data_model, R, profiles, iCOV):
    
    lM200,c200 = data_model
    
    if 12.5 < lM200 < 16.0 and 1 < c200 < 7:
        return log_likelihood_DS(data_model, R, profiles, iCOV)
        
    return -np.inf

# initializing

pos = np.array([np.random.uniform(12.5,15.5,15),
                np.random.uniform(1,5,15)]).T

nwalkers, ndim = pos.shape

t1 = time.time()

pool = Pool(processes=(ncores))    
sampler_DS = emcee.EnsembleSampler(nwalkers, ndim, log_probability_DS, 
                                args=(p.Rp,DSt,iCds),pool = pool)

sampler_DS.run_mcmc(pos, nit, progress=True)
pool.terminate()

t2 = time.time()

print('TIME DS')    
print((t2-t1)/60.)

mcmc_out_DS = sampler_DS.get_chain(flat=True).T
lM     = np.percentile(mcmc_out_DS[0][1500:], [16, 50, 84])
c200   = np.percentile(mcmc_out_DS[1][1500:], [16, 50, 84])

# NOW FIT q with Gamma components
# initializing
def log_likelihood(data_model, R, profiles, iCOV):
    
    q = data_model
    
    e = (1.-q)/(1.+q)

    gt, gx = profiles
    iCgt, iCgx = iCOV 

    GT,GX   = GAMMA_components_2h(R,zmean,ellip=e,M200 = 10**lM[1],c200=c200[1],cosmo_params=params)

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

pos = np.array([np.random.uniform(0.2,0.9,15)]).T

nwalkers, ndim = pos.shape

pool = Pool(processes=(ncores))    
sampler_GC = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(p.Rp,profiles,iCov),pool = pool)
				
sampler_GC.run_mcmc(pos, nit, progress=True)
pool.terminate()

t3 = time.time()

print('TIME G components')    
print((t3-t2)/60.)

mcmc_out_GC = sampler_GC.get_chain(flat=True).T
q      = np.percentile(mcmc_out_GC[0][1500:], [16, 50, 84])
    
print('TOTAL TIME FIT')    
print((time.time()-t1)/60.)

#-------------------
# saving mcmc out

table = [fits.Column(name='lM200', format='E', array=mcmc_out_DS[0]),
            fits.Column(name='q', format='E', array=mcmc_out_GC[0]),
            fits.Column(name='c200', format='E', array=mcmc_out_DS[1])]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))

h = fits.Header()
h.append(('lM200',np.round(lM[1],4)))
h.append(('elM200M',np.round(np.diff(lM)[0],4)))
h.append(('elM200m',np.round(np.diff(lM)[1],4)))

h.append(('c200',np.round(c200[1],4)))
h.append(('ec200M',np.round(np.diff(c200)[0],4)))
h.append(('ec200m',np.round(np.diff(c200)[1],4)))

h.append(('q',np.round(q[1],4)))
h.append(('eqM',np.round(np.diff(q)[0],4)))
h.append(('eqm',np.round(np.diff(q)[1],4)))


primary_hdu = fits.PrimaryHDU(header=h)

hdul = fits.HDUList([primary_hdu, tbhdu])

hdul.writeto(folder+outfile,overwrite=True)

print('SAVED FILE')
