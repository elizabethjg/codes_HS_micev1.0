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

bnum = '136_145_selected'
    
folder = '../../MICEv2.0/'

ides = np.loadtxt(folder+'profiles/136_145_ides.list').astype(int)

for j in ides:

    m_name = 'mapa_ind_139_145_'+str(j)+'.fits'

    try:
        mapa = fits.open(folder+m_name)[1].data
    except:
        continue

    x = mapa.xmpc
    y = mapa.ympc

    
    plt.scatter(x,y,c=mapa.K,vmin=-0.005,vmax=0.005)
    
    plt.savefig(folder+'mapas/'+bnum+'/ind_selected/mapa'+str(j)+'.png')
