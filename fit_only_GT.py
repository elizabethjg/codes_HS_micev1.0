import sys
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
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
import corner
import os


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
parser.add_argument('-folder', action='store', dest='folder', default='/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICEv2.0/profiles/')
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ang', action='store', dest='angle', default='standard')
parser.add_argument('-ncores', action='store', dest='ncores', default=40)
parser.add_argument('-RIN', action='store', dest='RIN', default=0)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=2500)
parser.add_argument('-nit', action='store', dest='nit', default=250)
parser.add_argument('-continue', action='store', dest='cont', default='False')
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

outfile     = 'fitresults_onlyGT_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+file_name
backup      = folder+'backup_'+outfile
plot_folder = folder+'plots_mcmc/'

os.system('mkdir '+plot_folder)

print('fitting profiles')
print(folder)
print(file_name)
print(angle)
print('ncores = ',ncores)
print('RIN ',RIN)
print('ROUT ',ROUT)
print('nit', nit)
print('continue',cont)
print('outfile',outfile)


profile = fits.open(folder+file_name)
h       = profile[0].header
p       = profile[1].data
cov     = profile[2].data
zmean   = h['Z_MEAN']    

cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)

def log_likelihood(data_model, R, gt, iCgt):
    lM200,q,c200 = data_model

    e = (1.-q)/(1.+q)

    GT,GX   = GAMMA_components(R,zmean,ellip=e,M200 = 10**lM200,c200 = c200,cosmo=cosmo)

    L_GT = -np.dot((gt-GT),np.dot(iCgt,(gt-GT)))/2.0

    return L_GT 
    

def log_probability(data_model, R, profiles, iCOV):
    
    lM200,q,c200 = data_model
    
    if 0.2 < q < 1.0 and 12.5 < lM200 < 16.0 and 1 < c200 < 7:
        return log_likelihood(data_model, R, profiles, iCOV)
        
    return -np.inf

# initializing

pos = np.array([np.random.uniform(12.5,15.5,15),
                np.random.uniform(0.2,0.9,15),
                np.random.uniform(1,5,15)]).T

qdist = pos[:,1]                
pos[qdist > 1.,1] = 1.

nwalkers, ndim = pos.shape


#-------------------
# running emcee

maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))
mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]

CovGT  = cov.COV_GT.reshape(len(p.Rp),len(p.Rp))[mr]
CovGT  = CovGT.reshape(maskr.sum(),maskr.sum())

p  = p[maskr]

GT  = p.GAMMA_Tcos


iCgt   = np.linalg.inv(CovGT)

t1 = time.time()
backend = emcee.backends.HDFBackend(backup)
if not cont:
    backend.reset(nwalkers, ndim)


pool = Pool(processes=(ncores))    
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(p.Rp,GT,iCgt),
                                backend=backend,pool = pool)
				


if cont:                                
    sampler.run_mcmc(None, nit, progress=True)
else:
    sampler.run_mcmc(pos, nit, progress=True)
pool.terminate()
    
    
print('TOTAL TIME FIT')    
print((time.time()-t1)/60.)

#-------------------
# saving mcmc out

mcmc_out = sampler.get_chain(flat=True).T

table = [fits.Column(name='lM200', format='E', array=mcmc_out[0]),
            fits.Column(name='q', format='E', array=mcmc_out[1]),
            fits.Column(name='c200', format='E', array=mcmc_out[2])]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))

lM     = np.percentile(mcmc_out[0][1500:], [16, 50, 84])
q      = np.percentile(mcmc_out[1][1500:], [16, 50, 84])
c200   = np.percentile(mcmc_out[2][1500:], [16, 50, 84])


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
