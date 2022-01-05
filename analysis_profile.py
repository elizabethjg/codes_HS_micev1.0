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

folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles/'

def plt_profile_compare(samp1,samp2):
    
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
        
    ax.plot(nfw.xplot,nfw.yplot,'C0',label=samp1)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C0',alpha=0.2)
    
    ax1.plot(rplot,gt,'C0')
    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C0',alpha=0.2)
    
    ax2.plot(rplot,gx,'C0')    
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C0',alpha=0.2)        

    ax3.plot(p.Rp,GTc,'C0', label = 'GT control')
    ax3.plot(p.Rp,GXc,'C0--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C0',alpha=0.2)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C0',alpha=0.2)
    
    ##### SECOND PROFILE
    
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
        
    ax.plot(nfw.xplot,nfw.yplot,'C1',label=samp2)
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
    ax1.yaxis.set_ticks([0.3,10,100])
    ax1.set_yticklabels([0.3,10,100])
        
    ax2.plot([0,5],[0,0],'C7')
    ax2.plot(rplot,gx,'C1')
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C1',alpha=0.2)
    
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-20,20)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])
    
    
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

def plt_profile_wofit(samp):
    

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
    CovS   = cov.COV_S.reshape(len(GT),len(GT))
    
    CovGT  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGTr = cov.COV_GT_reduced.reshape(len(GT),len(GT))
    CovGTc = cov.COV_GT_control.reshape(len(GT),len(GT))
    
    CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
    CovGXr = cov.COV_GX_reduced.reshape(len(GT),len(GT))
    CovGXc = cov.COV_GX_control.reshape(len(GT),len(GT))

    # FIT MONOPOLE
    rplot = np.arange(0.1,5,0.05)
    
    nfwS    = Sigma_fit(p.Rp,p.Sigma*(1.e6**2),np.diag(CovS)*(1.e6**2),zmean,cosmo,True)
    nfw     = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean,cosmo,True)
    gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    gtr,gxr = GAMMA_components(rplot,zmean,ellip=er,M200 =nfw.M200,c200 = nfw.c200,cosmo=cosmo)
    
    mass = str(np.round(np.log10(nfw.M200),2))
    massS = str(np.round(np.log10(nfwS.M200),2))
    
    f, ax = plt.subplots()
    f0, ax0 = plt.subplots()
    f1, ax1 = plt.subplots()
    f2, ax2 = plt.subplots()
    f3, ax3 = plt.subplots()
    
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax0.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+massS)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

    ax1.legend()
    ax2.legend()

    ax0.plot(p.Rp,p.Sigma,'C1')
    ax0.plot(nfwS.xplot,nfwS.yplot,'C3',label='fited nfw')
    ax0.fill_between(p.Rp,p.Sigma+np.diag(CovS),p.DSigma_T-np.diag(CovS),color='C1',alpha=0.2)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylabel(r'$\Sigma$')
    ax0.set_xlabel('r [$h^{-1}$ Mpc]')
    ax0.set_ylim(2,500)
    ax0.set_xlim(0.1,10)
    ax0.xaxis.set_ticks([0.1,1,5,7])
    ax0.set_xticklabels([0.1,1,5,7])
    ax0.yaxis.set_ticks([5,10,100])
    ax0.set_yticklabels([5,10,100])
    ax0.legend()
        
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



def plt_profile_fitted(samp,RIN,ROUT,fittype=''):
    
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

    
    zmean = h['z_mean']
    q  = h['q2d_mean']
    qr = h['q2dr_mean']
    
    print('mean q standard',q)
    print('mean q reduced',qr)
    
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
    CovS   = cov.COV_S.reshape(len(GT),len(GT))
    
    CovGT  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGTr = cov.COV_GT_reduced.reshape(len(GT),len(GT))
    CovGTc = cov.COV_GT_control.reshape(len(GT),len(GT))
    
    CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
    CovGXr = cov.COV_GX_reduced.reshape(len(GT),len(GT))
    CovGXc = cov.COV_GX_control.reshape(len(GT),len(GT))

    rplot = np.arange(0.1,5,0.05)
        
    # MCMC results

    fitpar = fits.open(folder+'fitresults'+fittype+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fitpar_red = fits.open(folder+'fitresults'+fittype+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[0].header
  
    efit = (1. - fitpar['q']) / (1. + fitpar['q'])
    efit_red = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])
  
    DS = Delta_Sigma_NFW(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo=cosmo)
    DSr = Delta_Sigma_NFW(rplot,zmean,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo=cosmo)
    gt,gx   = GAMMA_components(rplot,zmean,ellip=efit,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo=cosmo)
    gtr,gxr   = GAMMA_components(rplot,zmean,ellip=efit_red,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo=cosmo)
    
    print('Results standard fit')
    print('log(M200) = ',fitpar['lM200'],' c200 = ',fitpar['c200'],' q ',fitpar['q'])
    print('Results reduced fit')
    print('log(M200) = ',fitpar_red['lM200'],' c200 = ',fitpar_red['c200'],' q ',fitpar_red['q'])
    
    ##############
    mass = str(np.round(fitpar['lM200'],2))
    
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
    ax.plot(rplot,DS,'C3',label='fited nfw')
    ax.plot(rplot,DSr,'C3--',label='fited nfw')
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
    ax.axvline(RIN/1000.)
    ax.axvline(ROUT/1000.)
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
    ax1.axvline(RIN/1000.)
    ax1.axvline(ROUT/1000.)
    
    # ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
    # ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')
    
    ax2.plot([0,5],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C2')
    ax2.plot(p.Rp,GXr,'C5--')
    ax2.plot(rplot,gx,'C3')
    ax2.plot(rplot,gxr,'C3--')
    ax2.axvline(RIN/1000.)
    ax2.axvline(ROUT/1000.)

    
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
