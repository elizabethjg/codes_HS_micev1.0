import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from astropy.constants import G,c,M_sun, pc
from fit_models_colossus import Delta_Sigma_fit

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)
import corner
folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'


def save_fitted(samp,RIN,ROUT,fittype='_2h_2q'):

    
    component=''
    terms='1h+2h'
    pname='NFW'
    
    p_name = 'profile_'+samp+'.fits'  
    profile = fits.open(folder+p_name)
    
    print(p_name)
      
    p   = profile[1].data
    h   = profile[0].header
    
    zmean = h['z_mean']
    
    rplot = p.Rp

        
    # MCMC results

    
    fitpar = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fitpar_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[0].header

    fitd = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[1].data
    fitd_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[1].data
  
    efit     = (1. - fitpar['q']) / (1. + fitpar['q'])
    efit_red = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])

    efit2h     = (1. - fitpar['q2h']) / (1. + fitpar['q2h'])
    efit_red2h = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])
    
    DS1h  = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h')
    DS2h  = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h')

    gt1h,gx1h   = GAMMA_components(rplot,zmean,ellip=efit,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h',pname=pname)
    gt1hr,gx1hr = GAMMA_components(rplot,zmean,ellip=efit_red,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo_params=params,terms='1h',pname=pname)

    gt2h,gx2h   = GAMMA_components(rplot,zmean,ellip=efit2h,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h',pname=pname)
    gt2hr,gx2hr = GAMMA_components(rplot,zmean,ellip=efit_red2h,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo_params=params,terms='2h',pname=pname)
    
    fitout = np.array([p.Rp,DS1h,DS2h,gt1h,gx1h,gt1hr,gx1hr,gt2h,gx2h,gt2hr,gx2hr])

    np.savetxt(folder+'fitprofile'+fittype+samp+'_'+str(int(RIN))+'_'+str(int(ROUT))+'.cat',fitout,fmt='%10.2f')


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
    
    
    ndots = p.shape[0]
    
    
    GT  = p.GAMMA_Tcos
    GT = p.GAMMA_Tcos_reduced
    GTc = p.GAMMA_Tcos_control
    
    GX  = p.GAMMA_Xsin
    GX = p.GAMMA_Xsin_reduced
    GXc = p.GAMMA_Xsin_control
    
    # '''
    CovDS  = cov.COV_ST.reshape(len(GT),len(GT))
    
    CovGT  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGT = cov.COV_GT_reduced.reshape(len(GT),len(GT))
    CovGTc = cov.COV_GT_control.reshape(len(GT),len(GT))
    
    CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
    CovGX = cov.COV_GX_reduced.reshape(len(GT),len(GT))
    CovGXc = cov.COV_GX_control.reshape(len(GT),len(GT))

    # FIT MONOPOLE
    rplot = np.arange(0.1,10,0.05)
    
    nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean)
    
    gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 = nfw.M200,c200=nfw.c200,cosmo_params=params)
    gtr,gxr   = GAMMA_components(rplot,zmean,ellip=er,M200 = nfw.M200,c200=nfw.c200,cosmo_params=params)
    
    mass = str(np.round(np.log10(nfw.M200),1))
    
    f, ax_all = plt.subplots(2,2, figsize=(12,8),sharex = True)
    f.subplots_adjust(hspace=0)
    ax,ax1,ax2,ax3 = ax_all[0,0],ax_all[0,1],ax_all[1,0],ax_all[1,1]

    '''
    f, ax = plt.subplots()
    f0, ax0 = plt.subplots()
    f1, ax1 = plt.subplots()
    f2, ax2 = plt.subplots()
    f3, ax3 = plt.subplots()
    '''
    
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

    
    ax.legend()
    
        
    ax.plot(nfw.xplot,nfw.yplot,'C0',label=samp1)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C0',alpha=0.4)
    
    ax1.plot(rplot,gt,'C0')
    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C0',alpha=0.4)
    
    ax2.plot(rplot,gx,'C0')    
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C0',alpha=0.4)        

    # ax3.plot(p.Rp,GTc,'C0', label = 'GT control')
    # ax3.plot(p.Rp,GXc,'C0--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C0',alpha=0.4)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C0',alpha=0.4)
    
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
    
    ndots = p.shape[0]
    
    
    GT  = p.GAMMA_Tcos
    GT = p.GAMMA_Tcos_reduced
    GTc = p.GAMMA_Tcos_control
    
    GX  = p.GAMMA_Xsin
    GX = p.GAMMA_Xsin_reduced
    GXc = p.GAMMA_Xsin_control
    
    # '''
    CovDS  = cov.COV_ST.reshape(len(GT),len(GT))
    
    CovGT  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGT = cov.COV_GT_reduced.reshape(len(GT),len(GT))
    CovGTc = cov.COV_GT_control.reshape(len(GT),len(GT))
    
    CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
    CovGX = cov.COV_GX_reduced.reshape(len(GT),len(GT))
    CovGXc = cov.COV_GX_control.reshape(len(GT),len(GT))

    # FIT MONOPOLE
    rplot = np.arange(0.1,5,0.05)
    
    nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean)
        
    gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 = nfw.M200,c200=nfw.c200,cosmo_params=params)
    gtr,gxr   = GAMMA_components(rplot,zmean,ellip=er,M200 = nfw.M200,c200=nfw.c200,cosmo_params=params)
    
    mass = str(np.round(np.log10(nfw.M200),1))
        
    ax.plot(1000,1000,'w.' ,label='$\log M_{200}=$'+mass)
    ax1.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax2.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)
    ax3.plot(1000,1000,'w.',label='$\log M_{200}=$'+mass)

            
    ax.plot(nfw.xplot,nfw.yplot,'C1',label=samp2)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Delta\Sigma$')
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(1,200)
    ax.set_xlim(0.1,10)
    ax.xaxis.set_ticks([0.1,1,5,7])
    ax.set_xticklabels([0.1,1,5,7])
    ax.yaxis.set_ticks([5,10,100])
    ax.set_yticklabels([5,10,100])
    ax.legend()
        
    ax1.plot(rplot,gt,'C1')
    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C1',alpha=0.4)

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
        
    ax2.plot([0,10],[0,0],'C7')
    ax2.plot(rplot,gx,'C1')
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C1',alpha=0.4)
    
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-20,20)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])
    
    
    ax3.plot([0,10],[0,0],'C7')
    # ax3.plot(p.Rp,GTc,'C1', label = 'GT control')
    # ax3.plot(p.Rp,GXc,'C1--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C1',alpha=0.4)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C1',alpha=0.4)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,10)
    ax3.set_ylim(-10,7)
    ax3.xaxis.set_ticks([0.1,1,5,7])
    ax3.set_xticklabels([0.1,1,5,7])
    

def plt_profile_wofit(samp):
    

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo_as = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
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
    rplot = p.Rp
    
    nfwS    = Sigma_fit(p.Rp,p.Sigma*(1.e6**2),np.diag(CovS)*(1.e6**2),zmean,cosmo_as,True)
    nfw    = Delta_Sigma_fit(p.Rp,p.DSigma_T,np.diag(CovDS),zmean)
        
    gt,gx   = GAMMA_components(rplot,zmean,ellip=e,M200 = nfw.M200,c200=nfw.c200,cosmo_params=params)
    gtr,gxr   = GAMMA_components(rplot,zmean,ellip=er,M200 = nfw.M200,c200=nfw.c200,cosmo_params=params)
    
    qfit = str(q)
    qfit_red = str(qr)
    
    mass = str(np.round(np.log10(nfw.M200),2))
    cfit = str(np.round((nfw.c200),2))
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
    ax0.plot(nfwS.xplot,nfwS.yplot,'C3')
    ax0.fill_between(p.Rp,p.Sigma+np.diag(CovS),p.Sigma-np.diag(CovS),color='C1',alpha=0.4)
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
    ax.plot(nfw.xplot,nfw.yplot,'C3',label='fitted nfw $\log M_{200}=$'+mass+' $c_{200} = $'+cfit)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Delta\Sigma [M_{\odot}pc^{-2} h ]$')
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
    ax1.plot(rplot,gt,'C3',label = '$q_{fit} = $'+qfit+', $q = $'+str(q))
    ax1.plot(rplot,gtr,'C3--',label = '$q_{fit} = $'+qfit_red+', $q = $'+str(qr))
    ax1.legend(loc=3,frameon=False)


    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C4',alpha=0.4)
    ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C0',alpha=0.4)
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
    
    ax2.plot([0,10],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C2')
    ax2.plot(p.Rp,GXr,'C5--')
    ax2.plot(rplot,gx,'C3')
    ax2.plot(rplot,gxr,'C3--')

    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C2',alpha=0.4)
    ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C5',alpha=0.4)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-20,20)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])
    
    
    ax3.plot([0,10],[0,0],'C7')
    ax3.plot(p.Rp,GTc,'k', label = 'GT control')
    ax3.plot(p.Rp,GXc,'C8--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C8',alpha=0.4)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C7',alpha=0.4)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,10)
    ax3.set_ylim(-10,7)
    ax3.xaxis.set_ticks([0.1,1,5,7])
    ax3.set_xticklabels([0.1,1,5,7])
    ax3.legend()


def plt_profile_fitted_2h(samp,RIN,ROUT,fittype='_onlyq',
                          substract = False,component='',
                          terms='1h',pname='NFW'):
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
    params = {'flat': True, 'H0': 70.0, 'Om0': 0.25, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.95}
    '''
    
    h = profile[1].header
    p = profile[1].data
    '''

    
    zmean = h['z_mean']
    q  = np.round(h['q2d_mean'],2)
    qr = np.round(h['q2dr_mean'],2)
    
    print('mean q standard',q)
    print('mean q reduced',qr)
    
    e = (1-q)/(1+q)
    er = (1-qr)/(1+qr)
    
    
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

    rplot = np.round(p.Rp,2)
        
    # MCMC results

    # fitpar = fits.open(folder+'fitresults_2h_250_2000_'+p_name)[0].header
    # fitpar_red = fits.open(folder+'fitresults_2h_250_5000_reduced_'+p_name)[0].header
    
    fitpar = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fitpar_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[0].header
  
    efit = (1. - fitpar['q']) / (1. + fitpar['q'])
    efit_red = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])
    DS = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms=terms)
    DSr = DS

    gt,gx   = GAMMA_components(rplot,zmean,ellip=efit,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms=terms,pname=pname)
    gtr,gxr   = GAMMA_components(rplot,zmean,ellip=efit_red,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo_params=params,terms=terms,pname=pname)
    
    print('Results standard fit')
    print('log(M200) = ',fitpar['lM200'],' c200 = ',fitpar['c200'],' q ',fitpar['q'])
    print('Results reduced fit')
    print('log(M200) = ',fitpar_red['lM200'],' c200 = ',fitpar_red['c200'],' q ',fitpar_red['q'])
    
    ##############
    mass = str(np.round(fitpar['lM200'],2))
    c200 = str(np.round(fitpar['c200'],2))
    qfit = str(np.round(fitpar['q'],2))

    mass_red = str(np.round(fitpar_red['lM200'],2))
    c200_red = str(np.round(fitpar_red['c200'],2))
    qfit_red = str(np.round(fitpar_red['q'],2))
    
    f, ax_all = plt.subplots(2,2, figsize=(12,8),sharex = True)
    f.subplots_adjust(hspace=0)
    ax,ax1,ax2,ax3 = ax_all[0,0],ax_all[0,1],ax_all[1,0],ax_all[1,1]
    
    if substract:
        GT  = GT - GTc
        GTr = GTr - GTc
        GX  = GX - GXc
        GXr = GXr - GXc
        p.DSigma_T = p.DSigma_T - p.DSigma_X

    ax.set_title(p_name+fittype+component)

    ax.plot(p.Rp,p.DSigma_T,'C1')
    ax.plot(rplot,DS,'C3',label='$\log M_{200}=$'+mass+', $c_{200} = $'+c200)
    ax.plot(rplot,DSr,'C3--',label='$\log M_{200}=$'+mass_red+', $c_{200} = $'+c200_red)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.4)
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
    ax.axvline(RIN/1000.,color='C7')
    ax.axvline(ROUT/1000.,color='C7')
    ax.legend(loc=3,frameon=False)
    
    # ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
    # ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')
    
    ax1.plot(p.Rp,GT,'C4')
    ax1.plot(p.Rp,GTr,'C0--')
    ax1.plot(rplot,gt,'C3',label = '$q_{fit} = $'+qfit+', $q = $'+str(q))
    ax1.plot(rplot,gtr,'C3--',label = '$q_{fit} = $'+qfit_red+', $q = $'+str(qr))
    ax1.legend(loc=3,frameon=False)

    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C4',alpha=0.4)
    ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C0',alpha=0.4)
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
    ax1.axvline(RIN/1000.,color='C7')
    ax1.axvline(ROUT/1000.,color='C7')
    
    # ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
    # ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')
    
    ax2.plot([0,10],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C2')
    ax2.plot(p.Rp,GXr,'C5--')
    ax2.plot(rplot,gx,'C3')
    ax2.plot(rplot,gxr,'C3--')
    ax2.axvline(RIN/1000.,color='C7')
    ax2.axvline(ROUT/1000.,color='C7')

    
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C2',alpha=0.4)
    ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C5',alpha=0.4)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-20,16)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])
    
    
    ax3.plot([0,10],[0,0],'C7')
    ax3.plot(p.Rp,GTc,'k', label = 'GT control')
    ax3.plot(p.Rp,GXc,'C8--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C8',alpha=0.4)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C7',alpha=0.4)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,10)
    ax3.set_ylim(-20,16)
    ax3.xaxis.set_ticks([0.1,1,5,7])
    ax3.set_xticklabels([0.1,1,5,7])
    ax3.legend(frameon=False)
    
    if substract == True:
        f.savefig(folder+'plots/profile_'+samp+fittype+component+'_substracted.png',bbox_inches='tight')
    else:
        f.savefig(folder+'plots/profile_'+samp+fittype+component+'.png',bbox_inches='tight')

def plt_profile_fitted_2h_2q(samp,RIN,ROUT,fittype='_2h_2q',
                          substract = False,component='',
                          terms='1h+2h',pname='NFW'):
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    params = {'flat': True, 'H0': 70.0, 'Om0': 0.25, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.95}
    
    zmean = h['z_mean']
    q  = np.round(h['q2d_mean'],2)
    qr = np.round(h['q2dr_mean'],2)
    
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
    
    CovDS  = cov.COV_ST.reshape(len(GT),len(GT))
    CovS   = cov.COV_S.reshape(len(GT),len(GT))
    
    CovGT  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGTr = cov.COV_GT_reduced.reshape(len(GT),len(GT))
    CovGTc = cov.COV_GT_control.reshape(len(GT),len(GT))
    
    CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
    CovGXr = cov.COV_GX_reduced.reshape(len(GT),len(GT))
    CovGXc = cov.COV_GX_control.reshape(len(GT),len(GT))

    rplot = np.round(p.Rp,2)
        
    # MCMC results
    
    fitpar = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fitpar_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[0].header

    fitd = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[1].data
    fitd_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[1].data
                
    rplot,DS1h,DS2h,gt1h,gx1h,gt1hr,gx1hr,gt2h,gx2h,gt2hr,gx2hr = np.loadtxt(folder+'fitprofile'+fittype+samp+'_'+str(int(RIN))+'_'+str(int(ROUT))+'.cat')

    DS  = DS1h  + DS2h
    gt  = gt1h  + gt2h
    gtr = gt1hr + gt2hr
    gx  = gx1h  + gx2h
    gxr = gx1hr + gx2hr
        
    ##############    
    
    print('Results standard fit')
    print('log(M200) = ',fitpar['lM200'],' c200 = ',fitpar['c200'],' q ',fitpar['q'])
    print('Results reduced fit')
    print('log(M200) = ',fitpar_red['lM200'],' c200 = ',fitpar_red['c200'],' q ',fitpar_red['q'])
    
    ##############
    mass = str(np.round(fitpar['lM200'],2))
    c200 = str(np.round(fitpar['c200'],2))
    qfit = str(np.round(fitpar['q'],2))
    qfit2h = str(np.round(fitpar['q2h'],2))

    mass_red = str(np.round(fitpar_red['lM200'],2))
    c200_red = str(np.round(fitpar_red['c200'],2))
    qfit_red = str(np.round(fitpar_red['q'],2))
    qfit2h_red = str(np.round(fitpar_red['q2h'],2))
    
    f, ax_all = plt.subplots(2,2, figsize=(12,8),sharex = True)
    f.subplots_adjust(hspace=0)
    ax,ax1,ax2,ax3 = ax_all[0,0],ax_all[0,1],ax_all[1,0],ax_all[1,1]
    
    if substract:
        GT  = GT - GTc
        GTr = GTr - GTc
        GX  = GX - GXc
        GXr = GXr - GXc
        p.DSigma_T = p.DSigma_T - p.DSigma_X

    ax.set_title(p_name+fittype+component)

    ax.plot(p.Rp,p.DSigma_T,'C1')
    ax.plot(rplot,DS,'C3',label='$\log M_{200}=$'+mass+', $c_{200} = $'+c200)
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C1',alpha=0.4)
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
    ax.axvline(RIN/1000.,color='C7')
    ax.axvline(ROUT/1000.,color='C7')
    ax.legend(loc=3,frameon=False)
    
    # ax1.plot(RMt[0]*0.7,RMt[1]/0.7,'k',label='redMaPPer')
    # ax1.errorbar(RMt[0]*0.7,RMt[1]/0.7,yerr=RMt[2]/0.7,fmt = 'none',ecolor='0.5')
    
    ax1.plot(p.Rp,GT,'C4')
    ax1.plot(p.Rp,GTr,'C0--')
    ax1.plot(rplot,gt,'C3',label = '$q_{fit} = $'+qfit+', $q_{2h} = $'+qfit2h+', $q = $'+str(q))
    ax1.plot(rplot,gtr,'C3--',label = '$q_{fit} = $'+qfit_red+', $q_{2h} = $'+qfit2h_red+', $q = $'+str(qr))
    # ax1.plot(rplot,gt1h,'C2',lw=2)
    # ax1.plot(rplot,gt2h,'C2',lw=1)
    # ax1.plot(rplot,gt1hr,'C2--',lw=2)
    # ax1.plot(rplot,gt2hr,'C2--',lw=1)

    ax1.legend(loc=3,frameon=False)

    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C4',alpha=0.4)
    ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C0',alpha=0.4)
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
    ax1.axvline(RIN/1000.,color='C7')
    ax1.axvline(ROUT/1000.,color='C7')
    
    # ax2.plot(RMt[0]*0.7,RMt[3]/0.7,'k',label='redMaPPer')
    # ax2.errorbar(RMt[0]*0.7,RMt[3]/0.7,yerr=RMt[4]/0.7,fmt = 'none',ecolor='0.5')
    
    ax2.plot([0,10],[0,0],'C7')
    ax2.plot(p.Rp,GX,'C2')
    ax2.plot(p.Rp,GXr,'C5--')
    ax2.plot(rplot,gx,'C3')
    ax2.plot(rplot,gxr,'C3--')
    # ax2.plot(rplot,gx1h,'C2',lw=2)
    # ax2.plot(rplot,gx2h,'C2',lw=1)
    # ax2.plot(rplot,gx1hr,'C2--',lw=2)
    # ax2.plot(rplot,gx2hr,'C2--',lw=1)
    ax2.axvline(RIN/1000.,color='C7')
    ax2.axvline(ROUT/1000.,color='C7')

    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C2',alpha=0.4)
    ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C5',alpha=0.4)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$')
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-20,16)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])
    
    
    ax3.plot([0,10],[0,0],'C7')
    ax3.plot(p.Rp,GTc,'k', label = 'GT control')
    ax3.plot(p.Rp,GXc,'C8--', label = 'GX control')
    ax3.fill_between(p.Rp,GXc+np.diag(CovGXc),GXc-np.diag(CovGXc),color='C8',alpha=0.4)
    ax3.fill_between(p.Rp,GTc+np.diag(CovGTc),GTc-np.diag(CovGTc),color='C7',alpha=0.4)
    ax3.set_xlabel('r [$h^{-1}$ Mpc]')
    ax3.set_xscale('log')
    ax3.set_xlim(0.1,10)
    ax3.set_ylim(-20,16)
    ax3.xaxis.set_ticks([0.1,1,5,7])
    ax3.set_xticklabels([0.1,1,5,7])
    ax3.legend(frameon=False)
    
    if substract == True:
        f.savefig(folder+'plots/profile_'+samp+'_2q_+'+fittype+component+'_substracted.png',bbox_inches='tight')
    else:
        f.savefig(folder+'plots/profile_'+samp+'_2q_+'+fittype+component+'.png',bbox_inches='tight')
        
    f, ax = plt.subplots(2,1, figsize=(12,4),sharex = True)
    f.subplots_adjust(hspace=0)
    ax[0].plot(fitd.q2h,'C7',alpha=0.7,label='2h')
    ax[0].plot(fitd.q,'C6',alpha=0.7,label='1h')
    ax[1].plot(fitd_red.q2h,'C7',alpha=0.7,label='2h')
    ax[1].plot(fitd_red.q,'C6',alpha=0.7,label='1h')
    ax[0].legend(frameon=False,loc=1)
    ax[0].set_ylabel('$q$')
    ax[1].set_ylabel('$q_r$')
    ax[1].set_xlabel('$N$')
    f.savefig(folder+'plots/fit_2q_'+samp+'.png',bbox_inches='tight')
    

def corner_plot(samp,RIN,ROUT,relax=True,
                fittype = '_2h_2q',component='',
                fname = 'pru.tab',model='NFW'):

    

    matplotlib.rcParams.update({'font.size': 14})
    
    if model == 'Einasto':
        fittype += '_Ein'

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data

    # MCMC results
          

    fitd = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[1].data
    fitd_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[1].data

    halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data
        
    Eratio = (2.*halos.K/abs(halos.U))
    
    mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
    mfit_NFW = (halos.cNFW_rho > 1.)*(halos.cNFW_S > 1.)*(halos.cNFW_rho < 10.)*(halos.cNFW_S < 10.)*(halos.lgMNFW_rho > 12)*(halos.lgMNFW_S > 12)
    mfit_Ein = (halos.cEin_rho > 1.)*(halos.cEin_S > 1.)*(halos.cEin_rho < 10.)*(halos.cEin_S < 10.)*(halos.lgMEin_rho > 12)*(halos.lgMEin_S > 12)*(halos.alpha_rho > 0.)*(halos.alpha_S > 0.)*(halos.alpha_rho < 0.7)*(halos.alpha_S < 0.7)
    mhalos = mhalos*mfit_Ein*mfit_NFW
    
    if relax:
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)
        mhalos = mhalos*mrelax
    
    halos = halos[mhalos]

    qh  = np.mean(halos.b2D/halos.a2D)
    qhr = np.mean(halos.b2Dr/halos.a2Dr)

    lMNFW = np.percentile(halos.lgMNFW_rho, [16,50,84])
    cNFW  = np.percentile(halos.cNFW_rho  , [16,50,84])
    lMEin = np.percentile(halos.lgMEin_rho, [16,50,84])
    cEin  = np.percentile(halos.cEin_rho , [16,50,84])
    alpha = np.percentile(halos.alpha_rho , [16,50,84])
    qdm   = np.percentile(halos.b2D/halos.a2D , [16,50,84])
    qdmr  = np.percentile(halos.b2Dr/halos.a2Dr , [16,50,84])

    lMNFW_mean = np.mean(halos.lgMNFW_rho)
    cNFW_mean  = np.mean(halos.cNFW_rho)
    lMEin_mean = np.mean(halos.lgMEin_rho)
    cEin_mean  = np.mean(halos.cEin_rho)
    alpha_mean = np.mean(halos.alpha_rho)
    qdm_mean   = np.mean(halos.b2D/halos.a2D)
    qdmr_mean  = np.mean(halos.b2Dr/halos.a2Dr)

    lMfit  = np.percentile(fitd.lM200[1500:], [16,50,84])
    cfit   = np.percentile(fitd.c200[1500:], [16,50,84])
    qfit   = np.percentile(fitd.q[1500:], [16,50,84])
    q2hfit = np.percentile(fitd.q2h[1500:], [16,50,84])
    qfit_red   = np.percentile(fitd_red.q[1500:], [16,50,84])
    q2hfit_red = np.percentile(fitd_red.q2h[1500:], [16,50,84])

    lM200 = np.median(fitd.lM200[1500:])
    mcmc     = np.array([fitd.q[1500:],fitd.q2h[1500:]]).T
    mcmc_red = np.array([fitd_red.q[1500:],fitd_red.q2h[1500:]]).T
    
    
    labels      = ['$q_{1h}$','$q_{2h}$']
        
    fs  = corner.corner(mcmc,labels=labels,smooth=1.,range=[(0.2,1.),(0.2,1.)],truths=[qh,qh],label_kwargs=({'fontsize':16}),truth_color='C2',quantiles=(0.16, 0.84))
    fr  = corner.corner(mcmc_red,labels=labels,smooth=1.,range=[(0.2,1.),(0.2,1.)],truths=[qhr,qhr],label_kwargs=({'fontsize':16}),truth_color='C2',quantiles=(0.16, 0.84))

    if model == 'Einasto':
        
        lMh = np.mean(halos.lgMEin_rho)
        ch  = np.mean(halos.cEin_rho)
        ah  = np.mean(halos.alpha_rho)
        
        labels_DS   = ['$\log M_{200}$','$c_{200}$',r'$\alpha$']
        
        mcmc_DS  = np.array([fitd.lM200[1500:],fitd.c200[1500:],fitd.alpha[1500:]]).T
        
        fds = corner.corner(mcmc_DS,labels=labels_DS,smooth=1.,range=[(lM200-0.07,lM200+0.07),(2.,4.5),(0.1,0.5)],truths=[lMh,ch,ah],label_kwargs=({'fontsize':16}),truth_color='C2',quantiles=(0.16, 0.84))
        
        fds.savefig(folder+'../final_plots/mcmc_'+samp+'_Ein_fds.pdf',bbox_inches='tight')
        fs.savefig(folder+'../final_plots/mcmc_'+samp+'_Ein_fs.pdf',bbox_inches='tight')
        fr.savefig(folder+'../final_plots/mcmc_'+samp+'_Ein_fr.pdf',bbox_inches='tight')
    
    
    else:
        
        lMh = np.mean(halos.lgMNFW_rho)
        ch  = np.mean(halos.cNFW_rho)
        labels_DS   = ['$\log M_{200}$','$c_{200}$']
        
        mcmc_DS  = np.array([fitd.lM200[1500:],fitd.c200[1500:]]).T
        
        fds = corner.corner(mcmc_DS,labels=labels_DS,smooth=1.,range=[(lM200-0.07,lM200+0.07),(2.,4.5)],truths=[lMh,ch],label_kwargs=({'fontsize':16}),truth_color='C2',quantiles=(0.16, 0.84))
        
        fds.savefig(folder+'../final_plots/mcmc_'+samp+'_fds.pdf',bbox_inches='tight')
        fs.savefig(folder+'../final_plots/mcmc_'+samp+'_fs.pdf',bbox_inches='tight')
        fr.savefig(folder+'../final_plots/mcmc_'+samp+'_fr.pdf',bbox_inches='tight')

    # Write tables
    
    mres = [lMNFW,cNFW,lMEin,cEin,alpha,lMfit,cfit]
    qres = [qdm,qfit,q2hfit,qdmr,qfit_red,q2hfit_red]
    mmean = [lMNFW_mean,cNFW_mean,lMEin_mean,cEin_mean,alpha_mean,qdm_mean,qdmr_mean]

    mres = mres + qres 
    
    fm=open(folder+'../'+'allres_'+fname,'a')
    fq=open(folder+'../'+'qres_'+fname,'a')

    # fm.write(samp+' & ')
    # fq.write(samp+' & ')

    for x in mres:
        # fm.write('$'+str('%.2f' % (x[1]))+'_{-'+str('%.2f' % (np.diff(x)[0]))+'}^{+'+str('%.2f' % (np.diff(x)[1]))+'}$ & ')
        fm.write(str('%.2f' % (x[1]))+'  '+str('%.2f' % (np.diff(x)[0]))+'  '+str('%.2f' % (np.diff(x)[1]))+'  ')
    # x = mres[-1]
    # fm.write(str('%.2f' % (x[1]))+'  '+str('%.2f' % (np.diff(x)[0]))+'  '+str('%.2f' % (np.diff(x)[1]))+'  \n')
    # fm.write('$'+str('%.2f' % (x[1]))+'_{-'+str('%.2f' % (np.diff(x)[0]))+'}^{+'+str('%.2f' % (np.diff(x)[1]))+r'}$ \\'+' \n')
    
    
    for x in mmean[:-1]:
        fm.write(str('%.2f' % x)+'  ')
    x = mmean[-1]
    fm.write(str('%.2f' % x)+'  \n')
    fm.close()
    
    for x in qres[:-1]:
        fq.write('$'+str('%.2f' % (x[1]))+'_{-'+str('%.2f' % (np.diff(x)[0]))+'}^{+'+str('%.2f' % (np.diff(x)[1]))+'}$ & ')
    x = qres[-1]
    fq.write('$'+str('%.2f' % (x[1]))+'_{-'+str('%.2f' % (np.diff(x)[0]))+'}^{+'+str('%.2f' % (np.diff(x)[1]))+r'}$ \\'+' \n')
    fq.close()

    

def plt_profile_fitted_final(samp,RIN,ROUT,axx3,fittype='_2h_2q'):


    matplotlib.rcParams.update({'font.size': 12})    
    ax,ax1,ax2 = axx3
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
        
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

    rplot = np.round(p.Rp,2)
        
    # MCMC results
          
    rplot,DS1h,DS2h,gt1h,gx1h,gt1hr,gx1hr,gt2h,gx2h,gt2hr,gx2hr = np.loadtxt(folder+'fitprofile'+fittype+samp+'_'+str(int(RIN))+'_'+str(int(ROUT))+'.cat')

    DS  = DS1h  + DS2h
    gt  = gt1h  + gt2h
    gtr = gt1hr + gt2hr
    gx  = gx1h  + gx2h
    gxr = gx1hr + gx2hr
        
    ##############    

    ax.plot(p.Rp,p.DSigma_T,'C7')
    ax.plot(rplot,DS1h,'C0',label='1h')
    ax.plot(rplot,DS2h,'C9',label='2h')
    ax.plot(rplot,DS,'C3',label='1h+2h')
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C7',alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Delta\Sigma$',labelpad=2)
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(0.5,200)
    ax.set_xlim(0.1,10)
    ax.xaxis.set_ticks([0.1,1,5,7])
    ax.set_xticklabels([0.1,1,5,7])
    ax.yaxis.set_ticks([1,10,100])
    ax.set_yticklabels([1,10,100])
    ax.axvline(RIN/1000.,color='k',ls=':')
    ax.axvline(ROUT/1000.,color='k',ls=':')
    # ax.legend(loc=3,frameon=False,ncol=2)
    
    
    ax1.plot(p.Rp,GT,'C7',label='standard')
    ax1.plot(p.Rp,GTr,'C6--',label='reduced')
    ax1.plot(rplot,gt,'C3')
    ax1.plot(rplot,gtr,'C3--')
    ax1.plot(rplot,gt1h,'C0')
    ax1.plot(rplot,gt2h,'C9')
    ax1.plot(rplot,gt1hr,'C0--')
    ax1.plot(rplot,gt2hr,'C9--')
    # ax1.legend(loc=3,frameon=False)
    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C7',alpha=0.4)
    ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C7',alpha=0.4)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [$h^{-1}$ Mpc]')
    ax1.set_ylabel(r'$\Gamma_T$',labelpad=2)
    ax1.set_ylim(0.5,100)
    ax1.set_xlim(0.1,10)
    ax1.xaxis.set_ticks([0.1,1,5,7])
    ax1.set_xticklabels([0.1,1,5,7])
    ax1.yaxis.set_ticks([1,10,100])
    ax1.set_yticklabels([1,10,100])
    ax1.axvline(RIN/1000.,color='k',ls=':')
    ax1.axvline(ROUT/1000.,color='k',ls=':')
        
    ax2.plot([0,10],[0,0],'k')
    ax2.plot(p.Rp,GX,'C7')
    ax2.plot(p.Rp,GXr,'C5--')
    ax2.plot(rplot,gx,'C3')
    ax2.plot(rplot,gxr,'C3--')
    ax2.plot(rplot,gx1h,'C0')
    ax2.plot(rplot,gx2h,'C9')
    ax2.plot(rplot,gx1hr,'C0--')
    ax2.plot(rplot,gx2hr,'C9--')
    ax2.axvline(RIN/1000.,color='k',ls=':')
    ax2.axvline(ROUT/1000.,color='k',ls=':')
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C7',alpha=0.4)
    ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C6',alpha=0.4)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times$',labelpad=2)
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-16,17)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])

def plt_profile_bias():
    

    samp = 'HM_Lz_relaxed'

    matplotlib.rcParams.update({'font.size': 12})    
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS  = p.DSigma_T
    GT  = p.GAMMA_Tcos    
    GX  = p.GAMMA_Xsin
    
    # '''
    CovDS  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT  = cov.COV_GT.reshape(len(GT),len(GT))    
    CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
        
    ##############    
    # misal
    ##############    
    
    p_name = 'profile_'+samp+'_mis30.fits'
    profile = fits.open(folder+p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS_al  = p.DSigma_T
    GT_al  = p.GAMMA_Tcos
    GX_al  = p.GAMMA_Xsin
    
    # '''
    CovDS_al  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT_al  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGX_al  = cov.COV_GX.reshape(len(GT),len(GT))

    ##############    
    # miscen
    ##############    
    
    p_name = 'profile_HM_Lz_soff_relaxed_miscen.fits'
    profile = fits.open(folder+p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS_c  = p.DSigma_T
    GT_c  = p.GAMMA_Tcos
    GX_c  = p.GAMMA_Xsin
    
    # '''
    CovDS_c  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT_c  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGX_c  = cov.COV_GX.reshape(len(GT),len(GT))

    ##############    
    # miscen + misal
    ##############    
    
    p_name = 'profile_'+samp+'_mis20_miscen.fits'
    profile = fits.open(folder+p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS_cal  = p.DSigma_T
    GT_cal  = p.GAMMA_Tcos
    GX_cal  = p.GAMMA_Xsin
    # '''
    CovDS_cal  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT_cal  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGX_cal  = cov.COV_GX.reshape(len(GT),len(GT))
    
    # '''    
    y  = [DS,GT,GX]
    Cy = [np.diag(CovDS),np.diag(CovGT),np.diag(CovGX)]

    y_c  = [DS_c,GT_c,GX_c]
    Cy_c = [np.diag(CovDS_c),np.diag(CovGT_c),np.diag(CovGX_c)]

    y_al  = [DS_al,GT_al,GX_al]
    Cy_al = [np.diag(CovDS_al),np.diag(CovGT_al),np.diag(CovGX_al)]

    y_cal  = [DS_cal,GT_cal,GX_cal]
    Cy_cal = [np.diag(CovDS_cal),np.diag(CovGT_cal),np.diag(CovGX_cal)]
    
    y2  = [y_al,y_c,y_cal]
    Cy2 = [Cy_al,Cy_c,Cy_cal]
    
    labels = [r'$\sigma(\Delta \theta) = 30^{\circ}$',r'$p_{cc} = 0.0$',r'$\sigma(\Delta \theta) = 20^{\circ}, p_{cc} = 0.75$']
    
    f, ax = plt.subplots(6,3, gridspec_kw={'height_ratios': [3, 1]*3},figsize=(14,10),sharex = True)
    f.subplots_adjust(hspace=0)
    
    for j in range(3):
        
        for i in range(3):
            ax[j*2,i].plot(p.Rp,y[i],'C7')        
            ax[j*2,i].fill_between(p.Rp,y[i]+Cy[i],y[i]-Cy[i],color='C7',alpha=0.4)
            ax[j*2,i].plot(p.Rp,y2[j][i],'C3',label = labels[j])        
            ax[j*2,i].fill_between(p.Rp,y2[j][i]+Cy2[j][i],y2[j][i]-Cy2[j][i],color='C3',alpha=0.4)
            ax[j*2+1,i].plot(p.Rp,(y[i]-y2[j][i])/y[i],'C7')        
            ax[j*2+1,i].plot([0,10],[0,0],'k')        
            ax[j*2+1,i].set_ylim(-0.5,1.)

        ax[j*2+1,i].set_ylim(-5,5)
        ax[-1,j].set_xlabel('r [$h^{-1}$ Mpc]')
        
        
        j *= 2
        
        ax[j,0].legend(frameon=False,loc=4)
        
        ax[j,0].set_xscale('log')
        ax[j,0].set_yscale('log')
        ax[j,0].set_ylabel(r'$\Delta\Sigma$',labelpad=2)
        ax[j,0].set_ylim(0.5,200)
        ax[j,0].set_xlim(0.1,10)
        ax[j,0].xaxis.set_ticks([0.1,1,5,7])
        ax[j,0].set_xticklabels([0.1,1,5,7])
        ax[j,0].yaxis.set_ticks([1,10,100])
        ax[j,0].set_yticklabels([1,10,100])
        
        ax[j,1].set_xscale('log')
        ax[j,1].set_yscale('log')
        ax[j,1].set_ylabel(r'$\Gamma_T$',labelpad=2)
        ax[j,1].set_ylim(0.5,100)
        ax[j,1].set_xlim(0.1,10)
        ax[j,1].xaxis.set_ticks([0.1,1,5,7])
        ax[j,1].set_xticklabels([0.1,1,5,7])
        ax[j,1].yaxis.set_ticks([1,10,100])
        ax[j,1].set_yticklabels([1,10,100])
            
        ax[j,2].plot([0,10],[0,0],'k')
        ax[j,2].set_ylabel(r'$\Gamma_\times$',labelpad=2)
        ax[j,2].set_xscale('log')
        ax[j,2].set_xlim(0.1,10)
        ax[j,2].set_ylim(-16,17)
        ax[j,2].xaxis.set_ticks([0.1,1,5,7])
        ax[j,2].set_xticklabels([0.1,1,5,7])
        
        f.savefig(folder+'../final_plots/profile_test_bias.pdf',bbox_inches='tight')



# save_fitted('HM_Lz_relaxed_miscen',250,5000)
# save_fitted('HM_Hz_relaxed_miscen',250,5000)
# save_fitted('LM_Lz_relaxed_miscen',250,5000)
# save_fitted('LM_Hz_relaxed_miscen',250,5000)

'''
f, ax_all = plt.subplots(4,3, figsize=(14,14),sharex = True)
f.subplots_adjust(hspace=0)

hsamples = ['HM_Lz','LM_Lz','HM_Hz','LM_Hz']
hsamples_relaxed = ['HM_Lz_relaxed','LM_Lz_relaxed','HM_Hz_relaxed','LM_Hz_relaxed']
hsamples_relaxed_misal = ['HM_Lz_relaxed_misalign','LM_Lz_relaxed_misalign','HM_Hz_relaxed_misalign','LM_Hz_relaxed_misalign']

for j in range(len(ax_all)):
    
    plt_profile_fitted_final(hsamples_relaxed[j],250,5000,ax_all[j],fittype='_Ein_2h_2q')
    ax_all[j,0].text(1,100,hsamples[j],fontsize=14)

ax_all[0,0].legend(loc=3,frameon=False)
ax_all[0,1].legend(loc=3,frameon=False)


    
f.savefig(folder+'../final_plots/profile_relaxed_Ein.pdf',bbox_inches='tight')
# '''

res = np.loadtxt(folder+'../allres_pru.tab').T
lMNFW    = res[0]
elMNFWm  = res[1]
elMNFWM  = res[2]
cNFW     = res[3]
ecNFWm   = res[4]
ecNFWM   = res[5]
lMEin    = res[6]
elMEinm  = res[7]
elMEinM  = res[8]
cEin     = res[9]
ecEinm   = res[10]
ecEinM   = res[11]
alpha    = res[12]
ealpham  = res[13]
ealphaM  = res[14]
lM       = res[15]
elMm     = res[16]
elMM     = res[17]
c        = res[18]
ecm      = res[19]
ecM      = res[20]
q        = res[21]
eqm      = res[22]
eqM      = res[23]
q1h      = res[24]
eq1hm    = res[25]
eq1hM    = res[26]
q2h      = res[27]
eq2hm    = res[28]
eq2hM    = res[29]
qr       = res[30]
eqrm     = res[31]
eqrM     = res[32]
q1hr     = res[33]
eq1hrm   = res[34]
eq1hrM   = res[35]
q2hr     = res[36]
eq2hrm   = res[37]
eq2hrM   = res[38]
lMNFW_mean = res[39]
cNFW_mean  = res[40]
lMEin_mean = res[41]
cEin_mean  = res[42]
alpha_mean = res[43]
qdm_mean   = res[44]
qdmr_mean  = res[45]
