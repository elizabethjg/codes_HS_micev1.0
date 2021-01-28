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


    
folder = '../../MICEv2.0/'

for j in range(200):

    m_name = 'mapas/mapa_ind_'+str(j)+'.fits'

    mapa = fits.open(folder+m_name)[1].data

    x = mapa.xmpc
    y = mapa.ympc

    f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)
    
    ax[0].scatter(x,y,c=mapa.K,vmin=-0.02,vmax=0.02)
    ax[1].scatter(x,y,c=mapa.GT,vmin=-50,vmax=50.)
    ax[2].scatter(x,y,c=mapa.GX,vmin=-50,vmax=50.)
    
    ax[0].set_title('Kappa')
    ax[1].set_title('GT')
    ax[2].set_title('GX')
    
    f.savefig(folder+'mapas/mapas_ind/mapa'+str(j)+'.png')
