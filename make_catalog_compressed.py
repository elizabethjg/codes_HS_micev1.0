import sys
sys.path.append('/mnt/clemente/lensing')
import time
import numpy as np
import pandas as pd
from lensing37 import LensCat, gentools
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM
from mice import MICE
from multiprocessing import Pool
from multiprocessing import Process
from lensing37.LensCat.main import CompressedCatalog
import itertools
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

R_Mpc = 10.
ncores = 56

folder = '/mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/catalogs/'

f = fits.open(folder+'MICEv1.0_halo_cat.fits')
ft = Table(f[1].data)
df = ft.to_pandas()
mregion = ((df.dec < 1.5) | (df.dec > 40.))&(df.ra < 80.)   
df = df[mregion]

z_min = df['z_v'].min()


R_deg = gentools.Mpc2deg(R_Mpc=R_Mpc, z=df['z_v'], cosmo=cosmo)
Nlenses = df.shape[0]

t0 = time.time()
print( '*** {} ***'.format('MICE'))
print( 'Number of Lenses to Search: {}'.format(Nlenses))
# Search neighbors around each lens
MICE.load()
print('masking...')
mask = (MICE.data['z_v']>z_min)&(MICE.data['z_v']<1.3)
mregion = ((MICE.data['dec'] < 1.5) | (MICE.data['dec'] > 40.))&(MICE.data['ra'] < 80.)   
MICE.data = MICE.data[mask & mregion]
print('selecting neighbors...')

ra    = MICE.data.ra
dec   = MICE.data.dec
catid = MICE.data.CATID

def find_patch(lpar):
        RA0,DEC0,delta = lpar
        mask = (ra < (RA0+delta))&(ra > (RA0-delta))&(dec > (DEC0-delta))&(dec < (DEC0+delta))
        return np.array(catid[mask])

# SPLIT LENSING CAT

lbins = int(round(Nlenses/float(ncores), 0))
slices = ((np.arange(lbins)+1)*ncores).astype(int)
slices = slices[(slices < Nlenses)]
RA    = np.split(df.ra[:],slices)
DEC   = np.split(df.dec[:],slices)
delta = np.split(R_deg,slices)

ii = np.array([])
tslice = np.array([])

for l in range(len(RA)):
        print('########################')
        print('RUN ',l+1,' OF ',len(RA))

        t1 = time.time()

        num = len(RA[l])
        print(num)
        if num == 1:
                entrada = [RA[0],DEC[0],delta[0]]
                
                salida = [find_patch(entrada)]
        else:          
                entrada = np.array([RA[l],DEC[l],delta[l]]).T
                
                pool = Pool(processes=(num))
                salida = pool.map(find_patch, entrada)
                pool.terminate()
        
        t10 = time.time()
        ii = np.append(ii,salida)
        t20 = time.time()
        ts = (t20-t10)/60.
        print('t per add',ts)        
        
        t2 = time.time()
        ts = (t2-t1)/60.
        tslice = np.append(tslice,ts)
        print('------------------------')
        print('TIME SLICE')
        print(ts)
        print('Estimated ramaining time - minutes')
        print(np.mean(tslice)*(len(RA)-(l+1)))
        print('Estimated ramaining time - hours')
        print(np.mean(tslice)*(1./60.)*(len(RA)-(l+1)))
        print('########################')

print('MAKING CAT...')

# One catalog, two data frames, one for galaxies and one for groups
cat = CompressedCatalog(name = MICE.name, LensID='id')

print('Adding lenses')

# Lenses data
src_per_lens = np.fromiter(map(len, ii), dtype=np.int32)
mask_nsrc = src_per_lens>0
cat_ids = [list(MICE.data['CATID'][_]) for _ in ii]

dic ={'CATNAME': np.tile([MICE.name], len(ii)),
	'N_SOURCES': src_per_lens, 
	'CATID': np.array(cat_ids, dtype=np.object) }
ii_data = pd.DataFrame(dic)

df.reset_index(drop=True, inplace=True)
cat.data_L = pd.concat([df, ii_data], axis=1)[mask_nsrc]

print('Adding sources')

# Sources data
ii = list(itertools.chain.from_iterable(ii))
iu, im = np.unique(ii, return_counts=True)
extra_data = pd.DataFrame({'CATNAME': np.tile([MICE.name], len(iu)),
	'MULTIPLICITY': im})
mid = np.in1d(MICE.data.CATID,iu)
iu_data = MICE.data[mid].reset_index(drop=True)
cat.data_S = pd.concat([iu_data, extra_data], axis=1).reset_index(drop=True)

print('saving cat...')
cat.write_to('{}/gx_{}_matched_sources.fits'.format(folder, 'MICE'), overwrite=True)

print( 'Tiempo', (time.time() - t0)/60.)
print('C\'est fini :)')

