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


    
folder = '../../MICEv2.0/mapas/'

ides = np.loadtxt(folder+'list_mpru').astype(int)

for j in ides:

    m_name = 'mapa_ind_pru_'+str(j)+'.fits'

    mapa = fits.open(folder+m_name)[1].data

    x = mapa.xmpc
    y = mapa.ympc

    
    plt.scatter(x,y,c=mapa.K,vmin=-0.005,vmax=0.005)
    
    plt.savefig(folder+'pru/mapas_ind/mapa'+str(j)+'.png')
