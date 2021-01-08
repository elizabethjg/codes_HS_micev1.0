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
parser.add_argument('-folder', action='store', dest='folder', default='/mnt/clemente/lensing/HALO_SHAPE/MICE_v2.0/catalogs/')
parser.add_argument('-file', action='store', dest='file_name', default='profile.cat')
parser.add_argument('-ncores', action='store', dest='ncores', default=40)
parser.add_argument('-RIN', action='store', dest='RIN', default=0)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=5000)
parser.add_argument('-nit', action='store', dest='nit', default=250)
parser.add_argument('-continue', action='store', dest='cont', default='False')
args = parser.parse_args()


folder    = args.folder
file_name = args.file_name

if 'True' in args.cont:
	cont      = True
elif 'False' in args.cont:
	cont      = False

	
nit       = int(args.nit)
ncores    = args.ncores
ncores    = int(ncores)
RIN       = float(args.RIN)
ROUT      = float(args.ROUT)



outfile     = 'fitresults_all_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+file_name
backup      = folder+'backup_'+outfile
plot_folder = folder+'plots_mcmc/'

os.system('mkdir '+plot_folder)

print('fitting profiles')
print(folder)
print(file_name)
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
lMguess = np.log10(h['M200'])
cguess  = h['c200']

cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)

def log_likelihood(data_model, R, ds, iCds):
    lM200,s_off = data_model

    DS      = Delta_Sigma_NFW_miss(R,zmean,10**lM200,s_off=s_off,cosmo=cosmo)

    L_DS = -np.dot((ds-DS),np.dot(iCds,(ds-DS)))/2.0

    return L_DS 
    

def log_probability(data_model, r, profiles, iCOV):
    
    lM200,s_off = data_model
    
    if 0.0 < s_off < 0.5 and 12.5 < lM200 < 16.0:
        return log_likelihood(data_model, r, profiles, iCOV)
        
    return -np.inf

# initializing

pos = np.array([np.random.uniform(12.5,15.5,15),
                np.random.uniform(0.,0.2,15)]).T

qdist = pos[:,1]                
pos[qdist > 1.,1] = 1.

nwalkers, ndim = pos.shape

#-------------------
# running emcee

maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))

mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]

CovDS  = cov.COV_ST.reshape(len(p.Rp),len(p.Rp))[mr]

p  = p[maskr]

t1 = time.time()

DSt = p.DSigma_T

CovDS  = CovDS.reshape(maskr.sum(),maskr.sum())
iCov     = np.linalg.inv(CovDS)

backend = emcee.backends.HDFBackend(backup)
if not cont:
    backend.reset(nwalkers, ndim)
    
pool = Pool(processes=(ncores))
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, 
                                args=(p.Rp,DSt,iCov),backend=backend,pool = pool)


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

table = [fits.Column(name='logM200', format='E', array=mcmc_out[0]),
            fits.Column(name='s_off', format='E', array=mcmc_out[1])]

tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))

lMout  = np.percentile(mcmc_out[0][1500:], [16, 50, 84])
s_off  = np.percentile(mcmc_out[1][1500:], [16, 50, 84])

h = fits.Header()
h.append(('lM200',np.round(lMout[1],4)))
h.append(('elM200M',np.round(np.diff(lMout)[0],4)))
h.append(('elM200m',np.round(np.diff(lMout)[1],4)))

h.append(('soff',np.round(s_off[1],4)))
h.append(('esoffM',np.round(np.diff(s_off)[0],4)))
h.append(('esoffm',np.round(np.diff(s_off)[1],4)))

primary_hdu = fits.PrimaryHDU(header=h)

hdul = fits.HDUList([primary_hdu, tbhdu])

hdul.writeto(folder+outfile,overwrite=True)

fig = corner.corner(mcmc_out.T, labels=['lM200','s_off'])
plt.savefig(plot_folder+'corner_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+file_name[8:-4]+'png')

f, ax = plt.subplots(2, 1, figsize=(6,3))

ax[0].plot(mcmc_out[0],'k.',alpha=0.3)
ax[0].axvline(1500)
ax[0].axhline(lMout[1])
ax[0].axhline(lMout[1] - np.diff(lMout)[0],ls='--')
ax[0].axhline(lMout[1] + np.diff(lMout)[0],ls='--')

ax[1].plot(mcmc_out[1],'k.',alpha=0.3)
ax[1].axvline(1500)
ax[1].axhline(qout[1])
ax[1].axhline(qout[1] - np.diff(qout)[0],ls='--')
ax[1].axhline(qout[1] + np.diff(qout)[0],ls='--')


f.subplots_adjust(hspace=0,wspace=0)
plt.savefig(plot_folder+'walk_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+file_name[8:-4]+'png')
