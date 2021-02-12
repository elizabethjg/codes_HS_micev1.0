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

def plt_profile_compare(samp1,samp2,ax,ax1,ax2,ax3,RIN,ROUT,mv1,mv2):
    
    folder = '../../MICEv'+str(mv1)+'.0/profiles/'


    p_name = 'profile_'+samp1+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''
    micev = str(h['MICE version'])
    
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
    
    
    
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

    
    ax1.legend()
    ax2.legend()
    
    
    ax.plot(p.Rp,p.DSigma_T,'C0')
    ax.plot(nfw.xplot,nfw.yplot,'C0',label=samp1)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C0',alpha=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(2,200)
    ax.set_xlim(0.5,10)
    ax.xaxis.set_ticks([0.1,1,5,7])
    ax.set_xticklabels([0.1,1,5,7])
    ax.yaxis.set_ticks([1,10,100])
    ax.set_yticklabels([1,10,100])
    ax.legend()
    
    # ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
    # ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')
    
    ax1.plot(p.Rp,GT,'C0',label = samp1)
    ax1.plot(rplot,gt,'C0')

    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C0',alpha=0.2)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [$h^{-1}$ Mpc]')
    ax1.set_ylim(1,100)
    ax1.set_xlim(0.1,10)
    ax1.xaxis.set_ticks([0.1,1,5,7])
    ax1.set_xticklabels([0.1,1,5,7])
    ax1.yaxis.set_ticks([0.1,10,100])
    ax1.set_yticklabels([0.1,10,100])
    
    # ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
    # ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')
    
    ax2.plot([0,5],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C0')
    ax2.plot(rplot,gx,'C0')
    
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C0',alpha=0.2)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-20,20)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])
    
    
    ax3.plot([0,5],[0,0],'C7')
    ax3.plot(p.Rp,GTc,'C0', label = 'GT control')
    ax3.plot(p.Rp,GXc,'C0--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C0',alpha=0.2)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C0',alpha=0.2)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,10)
    ax3.set_ylim(-10,7)
    ax3.xaxis.set_ticks([0.1,1,5,7])
    ax3.set_xticklabels([0.1,1,5,7])
    ax3.legend()


    folder = '../../MICEv'+str(mv2)+'.0/profiles/'


    p_name = 'profile_'+samp2+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''
    micev = str(h['MICE version'])
    
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
    
    
    
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

    
    
    
    ax.plot(p.Rp,p.DSigma_T,'C1')
    ax.plot(nfw.xplot,nfw.yplot,'C1',label=samp2)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylabel(r'$\Delta\Sigma$')
    ax.set_ylim(0.5,200)
    ax.set_xlim(0.1,10)
    ax.xaxis.set_ticks([0.1,1,5,7])
    ax.set_xticklabels([0.1,1,5,7])
    ax.yaxis.set_ticks([1,10,100])
    ax.set_yticklabels([1,10,100])
    ax.legend()
    
    # ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
    # ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')
    
    ax1.plot(p.Rp,GT,'C1',label = samp2)
    ax1.plot(rplot,gt,'C1')

    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C1',alpha=0.2)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [$h^{-1}$ Mpc]')
    ax1.set_ylabel(r'$\Gamma_T$')
    ax1.set_ylim(1,100)
    ax1.set_xlim(0.1,10)
    ax1.xaxis.set_ticks([0.1,1,5,7])
    ax1.set_xticklabels([0.1,1,5,7])
    ax1.yaxis.set_ticks([0.1,10,100])
    ax1.set_yticklabels([0.1,10,100])


    
    # ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
    # ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')
    
    ax2.plot([0,5],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C1')
    ax2.plot(rplot,gx,'C1')
    
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C1',alpha=0.2)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-20,20)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])

    ax1.legend()
    ax2.legend()
    
    
    ax3.plot([0,5],[0,0],'C7')
    ax3.plot(p.Rp,GTc,'C1', label = 'GT control')
    ax3.plot(p.Rp,GXc,'C1--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C1',alpha=0.2)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C1',alpha=0.2)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,10)
    ax3.set_ylim(-10,7)
    ax3.xaxis.set_ticks([0.1,1,5,7])
    ax3.set_xticklabels([0.1,1,5,7])
    ax3.legend()

def plt_profile_wofit(samp,ax,ax1,ax2,ax3,RIN,ROUT,mv):
    
    folder = '../../MICEv'+str(mv)+'.0/profiles/'


    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''
    micev = str(h['MICE version'])
    
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
    
    
    
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

    
    ax1.legend()
    ax2.legend()
    
    ax.set_title('MICE v'+micev+'.0')
    ax1.set_title('MICE v'+micev+'.0')
    ax2.set_title('MICE v'+micev+'.0')
    ax3.set_title('MICE v'+micev+'.0')
    
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


def plt_profile_wofit_onlystandard(samp,ax,ax1,ax2,ax3,RIN,ROUT,mv):
    
    folder = '../../MICEv'+str(mv)+'.0/profiles/'


    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''
    micev = str(h['MICE version'])
    
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
    rplot = np.arange(0.1,10,0.1)
    
    nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo,True)
    gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    gtr,gxr = GAMMA_components(rplot,zmean,ellip=er,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    
    mass = str(np.round(np.log10(nfw.M200),1))
    con = str(np.round(nfw.c200,1))
    
        
    ax.plot(p.Rp,p.DSigma_T,'C1')
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.2)
    ax.plot(nfw.xplot,nfw.yplot,'C3',label='lM200 = '+mass+',c200 = '+con,alpha=0.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Delta\Sigma$')
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(0.5,200)
    ax.set_xlim(0.1,10)
    ax.xaxis.set_ticks([0.1,1,5,7])
    ax.set_xticklabels([0.1,1,5,7])
    ax.yaxis.set_ticks([1,10,100])
    ax.set_yticklabels([1,10,100])

    
    # ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
    # ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')
    
    ax1.plot(p.Rp,GT,'C4')
    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C4',alpha=0.2)
    ax1.plot(rplot,gt,'C3',alpha=0.5)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [$h^{-1}$ Mpc]')
    ax1.set_ylabel(r'$\Gamma_T$')
    ax1.set_ylim(0.5,200)
    ax1.set_xlim(0.1,10)
    ax1.xaxis.set_ticks([0.1,1,5,7])
    ax1.set_xticklabels([0.1,1,5,7])
    ax1.yaxis.set_ticks([1,10,100])
    ax1.set_yticklabels([1,10,100])
    
    # ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
    # ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')
    
    ax2.plot([0,5],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C2')
    ax2.plot(rplot,gx,'C3',alpha=0.5)
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C2',alpha=0.2)
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




def plt_profile(samp,ax,ax1,ax2,ax3,RIN,ROUT,mv):
    
    folder = '../../MICEv'+str(mv)+'.0/profiles/'

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''
    micev = str(h['MICE version'])
    
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
    rplot = np.arange(0.1,5,0.1)
    
    nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo,True)
    gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    gtr,gxr = GAMMA_components(rplot,zmean,ellip=er,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    
    mass = str(np.round(np.log10(nfw.M200),1))
    
    
    # MCMC results

    fall = fits.open(folder+'fitresults_all_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fq   = fits.open(folder+'fitresults_onlyq_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fGX  = fits.open(folder+'fitresults_onlyGX_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fGT  = fits.open(folder+'fitresults_onlyGT_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
  
    eall = (1. - fall['q']) / (1. + fall['q'])
    eq   = (1. - fq['q']) / (1. + fq['q'])
    eGX  = (1. - fGX['q']) / (1. + fGX['q'])
    eGT  = (1. - fGT['q']) / (1. + fGT['q'])
  
    DS_all = Delta_Sigma_NFW(rplot,zmean,10**fall['lM200'],cosmo=cosmo)
    DS_q   = Delta_Sigma_NFW(rplot,zmean,10**fq['lM200'],cosmo=cosmo)
    DS_GX  = Delta_Sigma_NFW(rplot,zmean,10**fGX['lM200'],cosmo=cosmo)
    DS_GT  = Delta_Sigma_NFW(rplot,zmean,10**fGT['lM200'],cosmo=cosmo)
    
    gt_all,gx_all = GAMMA_components(rplot,zmean,ellip=eall,M200 =10**fall['lM200'],cosmo=cosmo)
    gt_q  ,gx_q   = GAMMA_components(rplot,zmean,ellip=eq,M200 =10**fq['lM200'],cosmo=cosmo)
    gt_GX ,gx_GX  = GAMMA_components(rplot,zmean,ellip=eGX,M200 =10**fGX['lM200'],cosmo=cosmo)
    gt_GT ,gx_GT  = GAMMA_components(rplot,zmean,ellip=eGT,M200 =10**fGT['lM200'],cosmo=cosmo)
    

    
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

    
    ax1.legend()
    ax2.legend()
    
    ax.set_title('MICE v'+micev+'.0')
    ax1.set_title('MICE v'+micev+'.0')
    ax2.set_title('MICE v'+micev+'.0')
    ax3.set_title('MICE v'+micev+'.0')
    
    ax.plot(p.Rp,p.DSigma_T,'C1')
    ax.plot(nfw.xplot,nfw.yplot,'C3',label='fited nfw')
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.2)
    ax.plot(rplot,DS_all,'C1--',label='mcmc all')
    ax.plot(rplot,DS_q,'C0--',label='mcmc q')
    ax.plot(rplot,DS_GT,'C4--',label='mcmc gt')
    ax.plot(rplot,DS_GX,'C2--',label='mcmc gx')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Delta\Sigma$')
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(2,200)
    ax.set_xlim(0.1,5)
    ax.xaxis.set_ticks([0.1,1,3])
    ax.set_xticklabels([0.1,1,3])
    ax.yaxis.set_ticks([5,10,100])
    ax.set_yticklabels([5,10,100])
    ax.legend()
    
    # ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
    # ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')
    
    ax1.plot(p.Rp,GT,'C4',label = 'standard')
    ax1.plot(p.Rp,GTr,'C0--',label = 'reduced')
    ax1.plot(rplot,gt,'C3')
    ax1.plot(rplot,gtr,'C3--')

    ax1.plot(rplot,gt_all,'C1--')
    ax1.plot(rplot,gt_q,'C0--')
    ax1.plot(rplot,gt_GT,'C4--')
    ax1.plot(rplot,gt_GX,'C2--')


    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C4',alpha=0.2)
    ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C0',alpha=0.2)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [$h^{-1}$ Mpc]')
    ax1.set_ylabel(r'$\Gamma_T$')
    ax1.set_ylim(1,100)
    ax1.set_xlim(0.1,5)
    ax1.xaxis.set_ticks([0.1,1,3])
    ax1.set_xticklabels([0.1,1,3])
    ax1.yaxis.set_ticks([0.3,10,100])
    ax1.set_yticklabels([0.3,10,100])
    
    # ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
    # ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')
    
    ax2.plot([0,5],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C2')
    ax2.plot(p.Rp,GXr,'C5--')
    ax2.plot(rplot,gx,'C3')
    ax2.plot(rplot,gxr,'C3--')
    
    ax2.plot(rplot,gx_all,'C1--')
    ax2.plot(rplot,gx_q,'C0--')
    ax2.plot(rplot,gx_GT,'C4--')
    ax2.plot(rplot,gx_GX,'C2--')
    
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C2',alpha=0.2)
    ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C5',alpha=0.2)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,5)
    ax2.set_ylim(-10,10)
    ax2.xaxis.set_ticks([0.1,1,3])
    ax2.set_xticklabels([0.1,1,3])
    
    
    ax3.plot([0,5],[0,0],'C7')
    ax3.plot(p.Rp,GTc,'k', label = 'GT control')
    ax3.plot(p.Rp,GXc,'C8--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C8',alpha=0.2)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C7',alpha=0.2)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,5)
    ax3.set_ylim(-10,7)
    ax3.xaxis.set_ticks([0.1,1,3])
    ax3.set_xticklabels([0.1,1,3])
    ax3.legend()
    
    
'''
folder = '../../MICEv2.0/profiles/'

f, ax = plt.subplots()
f1, ax1 = plt.subplots()
f2, ax2 = plt.subplots()
f3, ax3 = plt.subplots()

RIN = 0
ROUT = 5000

# plt_profile_micec(samp,ax,ax1,ax2,ax3,RIN,ROUT)
plt_profile_compare('140','140',ax,ax1,ax2,ax3,RIN,ROUT,1,2)

samp = '140'

f.savefig(folder+'comparison-MICE/'+'DSigma_'+samp+'_'+str(RIN)+'_'+str(ROUT)+'.png')
f1.savefig(folder+'comparison-MICE/'+'GT_'+samp+'_'+str(RIN)+'_'+str(ROUT)+'.png')
f2.savefig(folder+'comparison-MICE/'+'GX_'+samp+'_'+str(RIN)+'_'+str(ROUT)+'.png')
f3.savefig(folder+'comparison-MICE/'+'control'+samp+'.png')





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

for j in range(7):
    
    plt_profile(136+2*j,ax[j],ax1[j],ax2[j],ax3[j])
'''
