import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from astropy.constants import G,c,M_sun, pc

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)

mv = 2
samp = '136_145'
    
folder = '../../MICEv'+str(mv)+'.0/profiles/'


p_name = 'profile_'+samp+'.fits'

p_name_ind1 = 'profile_'+samp+'_individual_part1.fits'
p_name_ind2 = 'profile_'+samp+'_individual_part2.fits'

profile = fits.open(folder+p_name)
profile_ind1 = fits.open(folder+p_name_ind1)
profile_ind2 = fits.open(folder+p_name_ind2)

print(p_name)
    
# '''
h   = profile[0].header
p   = profile[1].data
p_ind1 = profile_ind1[1].data
p_ind2 = profile_ind2[1].data
cov = profile[2].data

cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
'''

h = profile[1].header
p = profile[1].data
'''
micev = str(h['MICE version'])

nlens = h['N_LENSES']

zmean = h['z_mean']
q  = h['q2d_mean']
qr = h['q2dr_mean']

e = (1-q)/(1+q)
er = (1-qr)/(1+qr)

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

# '''
CovDS  = cov.COV_ST.reshape(len(GT),len(GT))

CovGT  = cov.COV_GT.reshape(len(GT),len(GT))
CovGTr = cov.COV_GT_reduced.reshape(len(GT),len(GT))
CovGTc = cov.COV_GT_control.reshape(len(GT),len(GT))

CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
CovGXr = cov.COV_GX_reduced.reshape(len(GT),len(GT))
CovGXc = cov.COV_GX_control.reshape(len(GT),len(GT))

# FIT MONOPOLE
rplot = np.arange(0.1,5,0.05)

nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo,True)
gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
gtr,gxr = GAMMA_components(rplot,zmean,ellip=er,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)

mass = str(np.round(np.log10(nfw.M200),1))



mask = p.Rp < 1.

pall = np.zeros((nlens,mask.sum()))

ides = np.zeros((nlens))

for j in range(nlens):
    
    try:
    
        plt.plot(p.Rp,p_ind1['DS'+str(j)][1:],'k')
        pall[j,:] = p_ind1['DS'+str(j)][1:][mask]
        ides[j] = p_ind1['DS'+str(j)][0]
    except:
        plt.plot(p.Rp,p_ind2['DS'+str(j)][1:],'k')
        pall[j,:] = p_ind2['DS'+str(j)][1:][mask]
        ides[j] = p_ind2['DS'+str(j)][0]

    


m   = np.tile(np.median(pall,axis=0),(nlens,1))
q10 = np.tile(np.quantile(pall,0.1,axis=0),(nlens,1))
q95 = np.tile(np.quantile(pall,0.97,axis=0),(nlens,1))

conserve = np.all((pall > q10)*(pall < q95),axis=1)

for j in range(nlens):
    
    if conserve[j]:
        try:
            plt.plot(p.Rp,p_ind1['DS'+str(j)][1:],'C0')
        except:
            plt.plot(p.Rp,p_ind2['DS'+str(j)][1:],'C0')
        




plt.plot(p.Rp,p.DSigma_T,'C3')

    
np.savetxt(folder+samp+'_ides.list', ides[conserve],fmt='%100i')
np.savetxt(folder+samp+'_ides_discard.list', ides[~conserve],fmt='%100i')
