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

h = 1.0
cosmo = LambdaCDM(H0=100*h, Om0=0.3, Ode0=0.7)


p_name = 'profile_bin14.fits'
profile = fits.open('../profiles/'+p_name)

h = profile[1].header
p = profile[1].data

zmean = h['z_mean']

H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
roc_mpc  = roc*((pc*1.0e6)**3.0)


ndots = p.shape[0]


GT  = p.GAMMA_Tcos
GTr = p.GAMMA_Tcos_reduced
GTc = p.GAMMA_Tcos_control

GX  = p.GAMMA_Xsin
GXr = p.GAMMA_Xsin_reduced
GXc = p.GAMMA_Xsin_control

DS_k     = np.zeros((100,ndots))

GT_k     = np.zeros((100,ndots))
GTr_k    = np.zeros((100,ndots))
GTc_k    = np.zeros((100,ndots))

GX_k     = np.zeros((100,ndots))
GXr_k    = np.zeros((100,ndots))
GXc_k    = np.zeros((100,ndots))


for k in range(100):
    
    DS_k[k,:]     = p['DSigma_T_K'+str(k+1)]
    
    GT_k[k,:]     = p['GAMMA_Tcos_K'+str(k+1)]
    GTr_k[k,:]    = p['GAMMA_Tcos_K_reduced'+str(k+1)]
    GTc_k[k,:]    = p['GAMMA_Tcos_control_K'+str(k+1)]

    GX_k[k,:]     = p['GAMMA_Tcos_K'+str(k+1)]
    GXr_k[k,:]    = p['GAMMA_Tcos_K_reduced'+str(k+1)]
    GXc_k[k,:]    = p['GAMMA_Tcos_control_K'+str(k+1)]



DS_kmean = np.mean(DS_k,axis=0)
DS_kstd  = np.std(DS_k,axis=0)

GT_kmean  = np.mean(GT_k,axis=0)
GTr_kmean = np.mean(GTr_k,axis=0)
GTc_kmean = np.mean(GTc_k,axis=0)

GX_kmean  = np.mean(GX_k,axis=0)
GXr_kmean = np.mean(GXr_k,axis=0)
GXc_kmean = np.mean(GXc_k,axis=0)

CovDS  = np.zeros((ndots,ndots))

CovGT   = np.zeros((ndots,ndots))
CovGTr  = np.zeros((ndots,ndots))
CovGTc  = np.zeros((ndots,ndots))

CovGX   = np.zeros((ndots,ndots))
CovGXr  = np.zeros((ndots,ndots))
CovGXc  = np.zeros((ndots,ndots))

Corr = np.zeros((ndots,ndots))

for k in range(100):
    
    dif = (p['DSigma_T_K'+str(k+1)]-DS_kmean)
    CovDS += np.outer(dif,dif)        

    dif = (p['GAMMA_Tcos_K'+str(k+1)]-GT_kmean)
    CovGT += np.outer(dif,dif)        
    dif = (p['GAMMA_Tcos_K_reduced'+str(k+1)]-GTr_kmean)
    CovGTr += np.outer(dif,dif)        
    dif = (p['GAMMA_Tcos_control_K'+str(k+1)]-GTc_kmean)
    CovGTc += np.outer(dif,dif)        
    
    dif = (p['GAMMA_Tcos_K'+str(k+1)]-GT_kmean)
    CovGX += np.outer(dif,dif)        
    dif = (p['GAMMA_Tcos_K_reduced'+str(k+1)]-GTr_kmean)
    CovGXr += np.outer(dif,dif)        
    dif = (p['GAMMA_Tcos_control_K'+str(k+1)]-GTc_kmean)
    CovGXc += np.outer(dif,dif)        
        
    difw = (p['DSigma_T_K'+str(k+1)]-DS_kmean)/DS_kstd
    Corr += np.outer(difw,difw)        

CovDS  *= 99/100.

CovGT  *= 99/100.
CovGTr *= 99/100.
CovGTc *= 99/100.

CovGX  *= 99/100.
CovGXr *= 99/100.
CovGXc *= 99/100.

Corr /= 100.

nfw = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo)

gt,gx = GAMMA_components(p.Rp,zmean,ellip=0.25,M200 =nfw.M200,c200 = None,cosmo=cosmo)
