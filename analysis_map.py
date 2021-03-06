import sys
import os
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from astropy.constants import G,c,M_sun, pc
cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)

folder = '../../MICEv2.0/'


bnum = '142'

p_name = 'profiles/profile_ebin_'+bnum+'.fits'
m_name = 'mapas/mapa_bin_'+bnum+'.fits'
m_name_miss = 'mapas/mapa_bin_'+bnum+'_miss.fits'


mapmodel_folder = 'mapas/'+bnum+'/map_models/'
pfolder = 'mapas/'+bnum+'/profiles/'
map_folder = 'mapas/'+bnum+'/maps_compare/'

os.system('mkdir '+folder+'mapas/'+bnum)
os.system('mkdir '+folder+pfolder)
os.system('mkdir '+folder+map_folder)
os.system('mkdir '+folder+mapmodel_folder)

fitmiss = fits.open(folder+'profiles/fitresults_mono_Rayleigh_0_2500_profile_ebin_'+bnum+'.fits')[0].header

profile = fits.open(folder+p_name)
mapa = fits.open(folder+m_name)[1].data
miss = fits.open(folder+m_name_miss)[1].data

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

y = mapa.ympc
x = mapa.xmpc

theta  = np.arctan2(y,x)
r = np.sqrt(x**2 + y**2)


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

ax1.legend()
ax2.legend()

o = np.argsort(r)

ax.plot(p.Rp,p.DSigma_T,'C1')
ax.plot(nfw.xplot,nfw.yplot,'C3',label = 'lM200='+mass+',c200='+str(nfw.c200),alpha=0.5) 
# ax.plot(r[o],miss.DS0[o],'C3',label = 'lM200='+str(lM200_miss)+',c200='+str(c200_miss)+',soff='+str(soff)) 
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
f.savefig(folder+pfolder+'profile_DS.png')

# ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
# ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')

ax1.plot(p.Rp,GT,'C4',label = 'standard')
# ax1.plot(p.Rp,GTr,'C0--',label = 'reduced')
ax1.plot(rplot,gt,'C3',alpha= 0.5)
# ax1.plot(rplot,gtr,'C3--')
ax1.plot(r[o],miss.Gt[o],'C3') 
ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C4',alpha=0.2)
# ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C0',alpha=0.2)
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
f1.savefig(folder+pfolder+'profile_GT.png')

# ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
# ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')

ax2.plot([0,5],[0,0],'C7')
ax2.plot(p.Rp,GX,'C2')
# ax2.plot(p.Rp,GXr,'C5--')
ax2.plot(rplot,gx,'C3',alpha=0.5)
# ax2.plot(rplot,gxr,'C3--')
ax2.plot(r[o],miss.Gx[o],'C3') 

ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C2',alpha=0.2)
# ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C5',alpha=0.2)
ax2.set_xlabel('r [$h^{-1}$ Mpc]')
ax2.set_ylabel(r'$\Gamma_\times$')
ax2.set_xscale('log')
ax2.set_xlim(0.1,10)
ax2.set_ylim(-20,20)
ax2.xaxis.set_ticks([0.1,1,5,7])
ax2.set_xticklabels([0.1,1,5,7])
f2.savefig(folder+pfolder+'profile_GX.png')

R = r*np.sqrt(q*(np.cos(theta))**2 + (np.sin(theta))**2 / q)


S0 = Sigma_NFW(r,zmean,nfw.M200,c200 = nfw.c200,cosmo=cosmo)
S  = Sigma_NFW(R,zmean,nfw.M200,c200 = nfw.c200,cosmo=cosmo)
S2  = quadrupole(r,zmean,nfw.M200,c200 = nfw.c200,cosmo=cosmo)

gt = Delta_Sigma_NFW(R,zmean,nfw.M200,c200 = nfw.c200,cosmo=cosmo)
gt0 = Delta_Sigma_NFW(r,zmean,nfw.M200,c200 = nfw.c200,cosmo=cosmo)
gtc,gxs  = GAMMA_components(r,zmean,ellip=e,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)

GT = gt0 + gtc*np.cos(2*theta)
GX = gxs*np.sin(2*theta)

GTmiss = miss.DS0 + miss.Gt*np.cos(2*theta)
GXmiss = miss.Gx*np.sin(2*theta)

Sr = S0+e*S2*np.cos(2*theta)
Sr_miss = (miss.S0-e*miss.S2*np.cos(2*theta))

# '''
# SMODEL

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=S,vmin=-10,vmax=50.)
ax[1].scatter(x,y,c=S0,vmin=-10,vmax=50.)
ax[2].scatter(x,y,c=S - S0,vmin=-10,vmax=50.)

ax[0].set_title('S(R)_model')
ax[1].set_title('S(r)_model')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+mapmodel_folder+'S_model_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=miss.S,vmin=-10,vmax=50.)
ax[1].scatter(x,y,c=miss.S0,vmin=-10,vmax=50.)
ax[2].scatter(x,y,c=miss.S - miss.S0,vmin=-10,vmax=50.)

ax[0].set_title('S(R)_model')
ax[1].set_title('S(r)_model')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+mapmodel_folder+'S_model_miscentred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=S,vmin=-10,vmax=50.)
ax[1].scatter(x,y,c=Sr,vmin=-10,vmax=50.)
ax[2].scatter(x,y,c=S-Sr,vmin=-10,vmax=50.)

ax[0].set_title('S(R)')
ax[1].set_title('S(r) + e S2(r) cos(2t)')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+mapmodel_folder+'S_model_vs_approx_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=miss.S,vmin=-10,vmax=50.)
ax[1].scatter(x,y,c=Sr_miss,vmin=-10,vmax=50.)
ax[2].scatter(x,y,c=miss.S-Sr_miss,vmin=-10,vmax=50.)

ax[0].set_title('S(R)')
ax[1].set_title('S(r) + e S2(r) cos(2t)')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+mapmodel_folder+'S_model_vs_approx_miscentred.png')


#### COMPARISON S

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=S0,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.K_control,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=S0 - mapa.K_control,vmin=-50,vmax=50.)

ax[0].set_title('S(r)_model')
ax[1].set_title('S_control')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'S0_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=S,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.K,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=S - mapa.K,vmin=-50,vmax=50.)

ax[0].set_title('S(R)_model')
ax[1].set_title('S')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'S_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=Sr,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.K,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=Sr - mapa.K,vmin=-50,vmax=50.)

ax[0].set_title('S(r) + e S2(r) cos(2t)')
ax[1].set_title('S')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'Sr_centred.png')


f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=Sr_miss,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.K,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=Sr_miss - mapa.K,vmin=-50,vmax=50.)

ax[0].set_title('S(r) + e S2(r) cos(2t)')
ax[1].set_title('S')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'Sr_miscentred.png')


f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=miss.S0,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.K_control,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=miss.S0 - mapa.K_control,vmin=-50,vmax=50.)

ax[0].set_title('S(r)_model')
ax[1].set_title('S_control')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'S0_miscentred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=miss.S,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.K,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=miss.S - mapa.K,vmin=-50,vmax=50.)

ax[0].set_title('S(R)_model')
ax[1].set_title('S')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'S_miscentred.png')

#COMPARISON control

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=mapa.GT_control,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.GT_control - gt0,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=mapa.GX_control,vmin=-50,vmax=50.)

ax[0].set_title('GT_control')
ax[1].set_title('GT_control - GT0_model')
ax[2].set_title('GX_control')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'control_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=mapa.GT_control,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.GT_control - miss.DS0,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=mapa.GX_control,vmin=-50,vmax=50.)

ax[0].set_title('GT_control')
ax[1].set_title('GT_control - GT0_model')
ax[2].set_title('GX_control')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'control_miscentred.png')


#COMPARISON monopole

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=gt0,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.GT_control,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=gt0-mapa.GT_control,vmin=-50,vmax=50.)

ax[0].set_title('GT0_model')
ax[1].set_title('GT_control')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'monopole_centred.png')


f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=miss.DS0,vmin=-50,vmax=50.)
ax[1].scatter(x,y,c=mapa.GT_control,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=miss.DS0-mapa.GT_control,vmin=-50,vmax=50.)

ax[0].set_title('GT0_model')
ax[1].set_title('GT_control')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'monopole_miscentred.png')

#COMPARISON GT

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=GT,vmin=-50,vmax=50.)  
ax[1].scatter(x,y,c=mapa.GT,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=GT-mapa.GT,vmin=-50,vmax=50.)

ax[0].set_title('(GT0+GT2)_model')
ax[1].set_title('GT')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'GT_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=GTmiss,vmin=-50,vmax=50.)  
ax[1].scatter(x,y,c=mapa.GT,vmin=-50,vmax=50.)
ax[2].scatter(x,y,c=GTmiss-mapa.GT,vmin=-50,vmax=50.)

ax[0].set_title('(GT0+GT2)_model')
ax[1].set_title('GT')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'GT_miscentred.png')

#COMPARISON GT2

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=GT-gt0,vmin=-10,vmax=10.)
ax[1].scatter(x,y,c=mapa.GT-gt0,vmin=-10,vmax=10.)
ax[2].scatter(x,y,c=GT-mapa.GT,vmin=-10,vmax=10.)

ax[0].set_title('GT2_model')
ax[1].set_title('GT - GT0_model')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'GT2_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=GTmiss-miss.DS0,vmin=-10,vmax=10.)
ax[1].scatter(x,y,c=mapa.GT-miss.DS0,vmin=-10,vmax=10.)
ax[2].scatter(x,y,c=GTmiss-mapa.GT,vmin=-10,vmax=10.)

ax[0].set_title('GT2_model')
ax[1].set_title('GT - GT0_model')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'GT2_miscentred.png')

#COMPARISON GX2

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=GX,vmin=-10,vmax=10.)
ax[1].scatter(x,y,c=mapa.GX,vmin=-10,vmax=10.)
ax[2].scatter(x,y,c=GX-mapa.GX,vmin=-10,vmax=10.)

ax[0].set_title('GX2_model')
ax[1].set_title('GX')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'GX2_centred.png')

f, ax = plt.subplots(1,3, figsize=(14,5), sharex=True, sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

im0 = ax[0].scatter(x,y,c=GXmiss,vmin=-10,vmax=10.)
ax[1].scatter(x,y,c=mapa.GX,vmin=-10,vmax=10.)
ax[2].scatter(x,y,c=GXmiss-mapa.GX,vmin=-10,vmax=10.)

ax[0].set_title('GX2_model')
ax[1].set_title('GX')
ax[2].set_title('Difference')
f.colorbar(im0, ax=ax, orientation='horizontal', fraction=.05)
f.savefig(folder+map_folder+'GX2_miscentred.png')
# '''
