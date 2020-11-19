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
# cosmo = LambdaCDM(H0=100*h, Om0=0.3, Ode0=0.7)
cosmo = LambdaCDM(H0=70, Om0=0.25, Ode0=0.75)


RM = np.loadtxt('/home/elizabeth/Documentos/posdoc/halo-elongation/redMapper/member_distribution/profiles/profile_total.cat').T
RMt  = np.loadtxt('/home/elizabeth/Documentos/posdoc/halo-elongation/redMapper/member_distribution/profiles/profile_total_t.cat').T
RMtp = np.loadtxt('/home/elizabeth/Documentos/posdoc/halo-elongation/redMapper/member_distribution/profiles/profile_total_tp.cat').T

def plt_profile(mbin,ax,ax1,ax2,ax3):

    p_name = 'profile_bin_'+str(int(mbin))+'.fits'
    profile = fits.open('../profiles/'+p_name)

    print(p_name)
    
    h = profile[1].header
    p = profile[1].data
    
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
    
    
    CovDS  = p.COV_ST.reshape(len(GT),len(GT))
    
    CovGT  = p.COV_GT.reshape(len(GT),len(GT))
    CovGTr = p.COV_GT_reduced.reshape(len(GT),len(GT))
    CovGTc = p.COV_GT_control.reshape(len(GT),len(GT))
    
    CovGX  = p.COV_GX.reshape(len(GT),len(GT))
    CovGXr = p.COV_GX_reduced.reshape(len(GT),len(GT))
    CovGXc = p.COV_GX_control.reshape(len(GT),len(GT))
    
    
    nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo,True)
    
    mass = str(np.round(np.log10(nfw.M200),1))
    
    rplot = np.arange(0.1,5,0.1)
    
    gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    gtr,gxr = GAMMA_components(rplot,zmean,ellip=er,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    
    gt_RM,gx_RM   = GAMMA_components(RM[0],z=0.25,ellip=0.22,M200 =2.e14,c200 = None,cosmo=cosmo)
    
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

    ax.legend()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    
    ax.plot(1000,1000,'w.',label='$M_{200}$'+mass)
    ax.plot(p.Rp,p.DSigma_T,'C1',label = 'standard')
    ax.plot(nfw.xplot,nfw.yplot,'C3')
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(2,200)
    ax.set_xlim(0.4,5)
    ax.xaxis.set_ticks([0.4,1,3])
    ax.set_xticklabels([0.4,1,3])
    ax.yaxis.set_ticks([5,10,100])
    ax.set_yticklabels([5,10,100])

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
    ax1.set_ylim(1,100)
    ax1.set_xlim(0.4,5)
    ax1.xaxis.set_ticks([0.4,1,3])
    ax1.set_xticklabels([0.4,1,3])
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
    ax2.set_xscale('log')
    ax2.set_xlim(0.4,5)
    ax2.set_ylim(-10,10)
    ax2.xaxis.set_ticks([0.4,1,3])
    ax2.set_xticklabels([0.4,1,3])
    
    
    ax3.plot([0,5],[0,0],'C7')
    ax3.plot(p.Rp,GTc,'k', label = 'tangential')
    ax3.plot(p.Rp,GXc,'C8--', label = 'cross')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C8',alpha=0.2)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C7',alpha=0.2)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.4,5)
    ax3.set_ylim(-10,7)
    ax3.xaxis.set_ticks([0.4,1,3])
    ax3.set_xticklabels([0.4,1,3])



f, ax = plt.subplots()
f1, ax1 = plt.subplots()
f2, ax2 = plt.subplots()
f3, ax3 = plt.subplots()



ft, axt  = plt.subplots(2,4, figsize=(12,6), sharey=True)
ft1, axt1 = plt.subplots(2,4, figsize=(12,6), sharey=True)
ft2, axt2 = plt.subplots(2,4, figsize=(12,6), sharey=True)
ft3, axt3 = plt.subplots(2,4, figsize=(12,6), sharey=True)

ft.subplots_adjust(hspace=0,wspace=0)
ft1.subplots_adjust(hspace=0,wspace=0)
ft2.subplots_adjust(hspace=0,wspace=0)
ft3.subplots_adjust(hspace=0,wspace=0)

ax  = axt.flatten()
ax1 = axt1.flatten()
ax2 = axt2.flatten()
ax3 = axt3.flatten()

for j in range(8):
    
    plt_profile(130+2*j,ax[j],ax1[j],ax2[j],ax3[j])
