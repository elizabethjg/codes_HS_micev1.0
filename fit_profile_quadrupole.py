import sys
sys.path.append('/mnt/clemente/lensing')
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
import time
import numpy as np
from astropy.io import fits
from astropy.cosmology import LambdaCDM
from models_profile import *
from multiprocessing import Pool
from multiprocessing import Process
import argparse
from astropy.constants import G,c,M_sun,pc
import emcee
from models_profiles import *
from astropy.cosmology import LambdaCDM

cosmo = LambdaCDM(H0=100, Om0=0.25, Ode0=0.75)


parser = argparse.ArgumentParser()
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ang', action='store', dest='angle', default='standard')
parser.add_argument('-ncores', action='store', dest='ncores', default=4)
parser.add_argument('-RIN', action='store', dest='RIN', default=0)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=5000)
parser.add_argument('-nit', action='store', dest='nit', default=250)
parser.add_argument('-continue', action='store', dest='cont', default='False')
args = parser.parse_args()

folder    = '../profiles/'
file_name = args.file_name
angle     = args.angle

if 'True' in args.cont:
	cont      = True
elif 'False' in args.cont:
	cont      = False

	
component = args.component
nit       = int(args.nit)
ncores    = args.ncores
ncores    = int(ncores)
rin       = float(args.RIN)
rout      = float(args.ROUT)

outfile = folder+'fitted_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.fits'
backup  = folder+'backup_'+file_name[:-4]+'_'+angle+'_'+str(int(rin))+'_'+str(int(rout))+'.fits'



print('fitting profiles')
print(folder)
print(file_name)
print(angle)
print(component)
print('ncores = ',ncores)
print('RIN ',rin)
print('ROUT ',rout)
print('nit', nit)
print('continue',cont)
print('outfile',outfile)


profile = fits.open(folder+file_name)
h       = profile[1].header
p       = profile[1].data
zmean   = h['Z_MEAN']    
lMguess = np.log10(h['M200'])
cguess  = h['c200']


def log_likelihood(data_model, r, profiles, iCOV):
    lM200,c200,q = data_model

    e = (1.-q)/(1.+q)

    ds, gt, gx = profiles
    iCds, iCgt, iCgx = iCOV 

    DS      = Delta_Sigma_NFW(R,z,10**lM200,c200 = c200,cosmo=cosmo)
    GT,GX   = GAMMA_components(r,zmean,ellip=e,M200 = 10**lM200,c200 = c200,cosmo=cosmo)

    L_DS = -np.dot((ds-DS),np.dot(iCds,(ds-DS)))/2.0
    L_GT = -np.dot((gt-GT),np.dot(iCgt,(gt-GT)))/2.0
    L_GX = -np.dot((gx-GX),np.dot(iCgx,(gx-GX)))/2.0

    return L_DS + L_GT + L_GX
    

def log_probability(data_model, r, profiles, iCOV):
    
    lM200,c200,q = data_model
    
    if 0.2 < q < 1.0 and 12.5 < lM200 < 16.0 and 1 < c200 10:
        return log_likelihood(data_model, r, profiles, iCOV)
    return -np.inf

# initializing

pos = np.array([np.random.uniform(12.5,15.5,15),
                np.random.normal(cguess,0.5,15),
                np.random.uniform(0.2,0.9,15)]).T

qdist = pos[:,2]                
pos[qdist > 1.,1] = 1.

nwalkers, ndim = pos.shape

#-------------------
# running emcee

profile = np.loadtxt(folder+file_name[:-4]+'_'+angle+'.cat').T
maskr   = (profile[0]>(rin/1000.))*(profile[0]<(rout/1000.))
profile = profile[:,maskr]

t1 = time.time()

backend = emcee.backends.HDFBackend(backup)
if not cont:
    backend.reset(nwalkers, ndim)
    


if component == 'tcos':
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(profile[0],profile[1],profile[2]),
				backend=backend)
				
elif component == 'xsin':                                
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(profile[0],profile[3],profile[4]),
				backend=backend)
elif component == 'both':                                
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(np.append(profile[0],profile[0]),
				np.append(profile[1],profile[3]),
				np.append(profile[2],profile[4])),
				backend=backend)


if cont:                                
    sampler.run_mcmc(None, nit, progress=True)
else:
    sampler.run_mcmc(pos, nit, progress=True)
    
print (time.time()-t1)/60.

#-------------------
# saving mcmc out

mcmc_out = sampler.get_chain(flat=True)


f1=open(outfile,'w')
f1.write('# ellip \n')
np.savetxt(f1,mcmc_out,fmt = ['%12.6f'])
f1.close()
