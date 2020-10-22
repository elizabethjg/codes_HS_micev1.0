import sys
sys.path.append('/mnt/clemente/lensing')
import time
import numpy as np
import pandas as pd
from lensing37 import LensCat, gentools
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM
import mice
cosmo = LambdaCDM(H0=100, Om0=0.3, Ode0=0.7)

f = fits.open('/mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/catalogs/MICEv1.0_halo_cat.fits')
ft = Table(f[1].data)
df = ft.to_pandas()
      
z_min = df['z_v'].min()

R_Mpc = 10.
R_deg = gentools.Mpc2deg(R_Mpc=R_Mpc, z=df['z_v'], cosmo=cosmo)

t0 = time.time()
print( '*** {} ***'.format('MICE'))
print( 'Number of Lenses to Search: {}'.format(df.shape[0]))
# Search neighbors around each lens
MICE.load()
print('masking...')
mask = (MICE.data['z_v']>z_min)*(MICE.data['z_v']<1.3)
MICE.data = MICE.data[mask]
print('selecting neighbors...')
cat = MICE.find_neighbors(centre=df[['ra','dec']],upper_radii=R_deg, append_data=df, compressed=True, njobs=1)
print('saving cat...')
cat.write_to('{}/gx_{}_{}.fits'.format(folder, cat.name, lensname), overwrite=True)

print( 'Tiempo', (time.time() - t0)/60.)
print('C\'est fini :)')

