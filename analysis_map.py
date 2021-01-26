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


p_name = 'profiles/profile_ebin_142.fits'
m_name = 'mapas/mapa_bin_142.fits'

profile = fits.open(folder+p_name)
mapa = fits.open(folder+m_name)[1].data
fitmiss = fits.open(folder+'profiles/fitresults_mono_Rayleigh_0_2500_profile_ebin_140.fits')[0].header

print(p_name)

# '''
h   = profile[0].header
p   = profile[1].data
cov = profile[2].data

lM200_miss = fitmiss['lm200']
c200_miss = fitmiss['c200']
soff = fitmiss['soff']

cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)


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

f, ax = plt.subplots() 
f1, ax1 = plt.subplots() 
f2, ax2 = plt.subplots() 
f3, ax3 = plt.subplots() 

ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)


ax1.legend()
ax2.legend()


ax.plot(p.Rp,p.DSigma_T,'C1')
ax.plot(nfw.xplot,nfw.yplot,'C3',label='fited nfw')
ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.2)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r'$\Delta\Sigma$')
ax.set_xlabel('r [$h^{-1}$ Mpc]')
ax.set_ylim(2,200)
ax.set_xlim(0.1,10)
ax.xaxis.set_ticks([0.1,1,5,7])
ax.set_xticklabels([0.1,1,5,7])
ax.yaxis.set_ticks([5,10,100])
ax.set_yticklabels([5,10,100])
ax.legend()

# ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
# ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')

ax1.plot(p.Rp,GT,'C4',label = 'standard')
ax1.plot(p.Rp,GTr,'C0--',label = 'reduced')
ax1.plot(rplot,gt,'C3')
ax1.plot(rplot,gtr,'C3--')

ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C4',alpha=0.2)
ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C0',alpha=0.2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('r [$h^{-1}$ Mpc]')
ax1.set_ylabel(r'$\Gamma_T$')
ax1.set_ylim(1,100)
ax1.set_xlim(0.1,10)
ax1.xaxis.set_ticks([0.1,1,5,7])
ax1.set_xticklabels([0.1,1,5,7])
ax1.yaxis.set_ticks([0.3,10,100])
ax1.set_yticklabels([0.3,10,100])

# ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
# ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')

ax2.plot([0,5],[0,0],'C7')
ax2.plot(p.Rp,GX,'C2')
ax2.plot(p.Rp,GXr,'C5--')
ax2.plot(rplot,gx,'C3')
ax2.plot(rplot,gxr,'C3--')

ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C2',alpha=0.2)
ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C5',alpha=0.2)
ax2.set_xlabel('r [$h^{-1}$ Mpc]')
ax2.set_ylabel(r'$\Gamma_\times$')
ax2.set_xscale('log')
ax2.set_xlim(0.1,10)
ax2.set_ylim(-20,20)
ax2.xaxis.set_ticks([0.1,1,5,7])
ax2.set_xticklabels([0.1,1,5,7])


ax3.plot([0,5],[0,0],'C7')
ax3.plot(p.Rp,GTc,'k', label = 'GT control')
ax3.plot(p.Rp,GXc,'C8--', label = 'GX control')
ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C8',alpha=0.2)
ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C7',alpha=0.2)
ax3.set_xlabel('r [$h^{-1}$ Mpc]')
ax3.set_xscale('log')
ax3.set_xlim(0.1,10)
ax3.set_ylim(-10,7)
ax3.xaxis.set_ticks([0.1,1,5,7])
ax3.set_xticklabels([0.1,1,5,7])
ax3.legend()

y = mapa.ympc
x = mapa.xmpc

theta  = np.arctan2(y,x)
r = np.sqrt(x**2 + y**2)

gt0 = Delta_Sigma_NFW(r,zmean,nfw.M200,c200 = nfw.c200,cosmo=cosmo)
gtc,gxs  = GAMMA_components(r,zmean,ellip=e,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)

GT = gt0 + gtc*np.cos(2*theta)
GX = gxs*np.sin(2*theta)

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(x,y,c=mapa.GT_control,vmin=0,vmax=40.)
ax[1].scatter(x,y,c=gt0,vmin=0,vmax=40.)
ax[2].scatter(x,y,c=mapa.GT_control-gt0,vmin=-10,vmax=10.)

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(x,y,c=mapa.GT-gt0,vmin=-20,vmax=20.)
ax[1].scatter(x,y,c=GT-gt0,vmin=-20,vmax=20.)
ax[2].scatter(x,y,c=mapa.GT-GT,vmin=-20,vmax=20.)

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].scatter(x,y,c=mapa.GX,vmin=-10,vmax=10.)
ax[1].scatter(x,y,c=GX,vmin=-10,vmax=10.)
ax[2].scatter(x,y,c=mapa.GX-GX,vmin=-10,vmax=10.)
