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

bnum = 'filtered'
    
folder = '../../MICEv2.0/'

ides = np.loadtxt(folder+'profiles/sel_ides_discard.list').astype(int)
allids = np.loadtxt(folder+'mapas/'+bnum+'/list').astype(int)

for j in allids:

    m_name = 'mapas/'+bnum+'/mapa_ind_selec_'+str(j)+'.fits'

    try:
        mapa = fits.open(folder+m_name)[1].data
    except:
        continue

    x = mapa.xmpc
    y = mapa.ympc

    plt.figure()
    plt.scatter(x,y,c=mapa.K,vmin=-0.015,vmax=0.015)
    
    if np.in1d(j,ides)[0]:
        plt.savefig(folder+'mapas/'+bnum+'/ind_discarded/mapa'+str(j)+'.png')
    else:
        plt.savefig(folder+'mapas/'+bnum+'/ind_selected/mapa'+str(j)+'.png')
