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

folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'


def plt_profile_fitted_2h_2q(samp,RIN,ROUT,fittype='_2h_2q',
                          substract = False,component='',
                          terms='1h+2h',pname='NFW'):
    
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

    rplot = np.round(p.Rp,2)
        
    # MCMC results

    # fitpar = fits.open(folder+'fitresults_2h_250_2000_'+p_name)[0].header
    # fitpar_red = fits.open(folder+'fitresults_2h_250_5000_reduced_'+p_name)[0].header
    
    fitpar = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fitpar_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[0].header

    fitd = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[1].data
    fitd_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[1].data
  
    efit     = (1. - fitpar['q']) / (1. + fitpar['q'])
    efit_red = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])

    efit2h     = (1. - fitpar['q2h']) / (1. + fitpar['q2h'])
    efit_red2h = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])
    
    DS  = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms=terms)
    DSr = DS

    gt1h,gx1h   = GAMMA_components(rplot,zmean,ellip=efit,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h',pname=pname)
    gt1hr,gx1hr = GAMMA_components(rplot,zmean,ellip=efit_red,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo_params=params,terms='1h',pname=pname)

    gt2h,gx2h   = GAMMA_components(rplot,zmean,ellip=efit2h,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h',pname=pname)
    gt2hr,gx2hr = GAMMA_components(rplot,zmean,ellip=efit_red2h,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo_params=params,terms='2h',pname=pname)
    
    gt  = gt1h  + gt2h
    gtr = gt1hr + gt2hr
    gx  = gx1h  + gx2h
    gxr = gx1hr + gx2hr
    
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
    
    
def plt_map_fitted(samp,RIN,ROUT,fittype=''):
    
    m_name = '../maps/map_'+samp+'.fits'
    mapa = fits.open(folder+m_name)

    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    m   = mapa[1].data

    print(p_name)
    
    
    cosmo = LambdaCDM(H0=100*h['hcosmo'], Om0=0.25, Ode0=0.75)
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
    
    H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
    roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
    roc_mpc  = roc*((pc*1.0e6)**3.0)
    
    
    ndots = p.shape[0]
    
    
    S   = m.K
    Sr  = m.K_reduced
    Sc  = m.K_control
    
    GT  = m.GT
    GTr = m.GT_reduced
    GTc = m.GT_control
    
    GX  = m.GX
    GXr = m.GX_reduced
    GXc = m.GX_control
    
    x = m.xmpc
    y = m.ympc
    theta  = np.arctan2(m.ympc,m.xmpc)
    r = np.sqrt(m.xmpc**2 + m.ympc**2)
        
    # MCMC results

    fitpar = fits.open(folder+'fitresults'+fittype+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[0].header
    fitpar_red = fits.open(folder+'fitresults'+fittype+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[0].header
  
    efit = (1. - fitpar['q']) / (1. + fitpar['q'])
    efit_red = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])

    R = r*np.sqrt(fitpar['q']*(np.cos(theta))**2 + (np.sin(theta))**2 / fitpar['q'])
    Rr = r*np.sqrt(fitpar_red['q']*(np.cos(theta))**2 + (np.sin(theta))**2 / fitpar_red['q'])

    S0_fit = Sigma_NFW(r,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo=cosmo)/(1.e6**2)
    Se_fit  = Sigma_NFW(R,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo=cosmo)/(1.e6**2)
    Se_fit_red  = Sigma_NFW(Rr,zmean,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo=cosmo)/(1.e6**2)
    
    S2_fit  = quadrupole(r,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo=cosmo)
    
    gt0 = Delta_Sigma_NFW(r,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo=cosmo)

    gtc,gxs  = GAMMA_components(r,zmean,ellip=efit,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo=cosmo)
    gtc_red,gxs_red  = GAMMA_components(r,zmean,ellip=efit_red,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo=cosmo)
    
    GTfit = gt0 + gtc*np.cos(2*theta)
    GTfit_red = gt0 + gtc_red*np.cos(2*theta)

    GXfit = gxs*np.sin(2*theta)
    GXfit_red = gxs_red*np.sin(2*theta)
        
    Sfit = S0_fit+efit*S2_fit*np.cos(2*theta)
    Sfit_red = S0_fit+efit_red*S2_fit*np.cos(2*theta)
      
    print('Results standard fit')
    print('log(M200) = ',fitpar['lM200'],' c200 = ',fitpar['c200'],' q ',fitpar['q'])
    print('Results reduced fit')
    print('log(M200) = ',fitpar_red['lM200'],' c200 = ',fitpar_red['c200'],' q ',fitpar_red['q'])
    
    ##############
    # SURFACE DENSITY PLOT
    # JUST MONOPOLE
    f, ax = plt.subplots(1,2, figsize=(6.5,3), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0].scatter(x,y,c=Sc,vmin=-10,cmap='cividis',vmax=50.)
    ax[1].scatter(x,y,c=S0_fit,vmin=-10,cmap='cividis',vmax=50.)

    ax[0].set_title('S0_MICE')
    ax[1].set_title('S0_fit')

    ax[0].set_ylabel('y [Mpc]')
    ax[1].set_xlabel('x [Mpc]')
    ax[0].set_xlabel('x [Mpc]')

    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_S0_'+samp+fittype+'.png',bbox_inches='tight')    

    # Res
    f, ax = plt.subplots(1,1, figsize=(3.5,3), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax.scatter(x,y,c=Sc-S0_fit,vmin=-10,cmap='cividis',vmax=10.)

    ax.set_title('S0_MICE - S0_fit')

    ax.set_ylabel('y [Mpc]')
    ax.set_xlabel('x [Mpc]')

    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_S0_res_'+samp+fittype+'.png',bbox_inches='tight')    
    
    # S(R)
    
    f, ax = plt.subplots(2,2, figsize=(6.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0,0].scatter(x,y,c=S,vmin=-10,cmap='cividis',vmax=50.)
    ax[0,1].scatter(x,y,c=Se_fit,vmin=-10,cmap='cividis',vmax=50.)

    ax[1,0].scatter(x,y,c=Sr,vmin=-10,cmap='cividis',vmax=50.)
    ax[1,1].scatter(x,y,c=Se_fit_red,vmin=-10,cmap='cividis',vmax=50.)
    

    ax[0,0].set_ylabel('y [Mpc]')
    ax[1,0].set_ylabel('y [Mpc]')
    ax[1,1].set_xlabel('x [Mpc]')
    ax[1,0].set_xlabel('x [Mpc]')

    ax[0,0].set_title('S_MICE')
    ax[0,1].set_title('S_fit(R)')

    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_SR_'+samp+fittype+'.png',bbox_inches='tight')    
    
    # Res
    f, ax = plt.subplots(2, figsize=(3.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0].scatter(x,y,c=S - Se_fit,vmin=-10,cmap='cividis',vmax=10.)
    ax[1].scatter(x,y,c=Sr - Se_fit_red,vmin=-10,cmap='cividis',vmax=10.)

    ax[0].set_ylabel('y [Mpc]')
    ax[1].set_ylabel('y [Mpc]')
    ax[1].set_xlabel('x [Mpc]')

    ax[0].set_title('S_MICE - S_fit(R)')


    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_SR_res_'+samp+fittype+'.png',bbox_inches='tight')    

    # S(r) + S2

    f, ax = plt.subplots(2,2, figsize=(6.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0,0].scatter(x,y,c=S,vmin=-10,cmap='cividis',vmax=50.)
    ax[0,1].scatter(x,y,c=Sfit,vmin=-10,cmap='cividis',vmax=50.)

    ax[1,0].scatter(x,y,c=Sr,vmin=-10,cmap='cividis',vmax=50.)
    ax[1,1].scatter(x,y,c=Sfit_red,vmin=-10,cmap='cividis',vmax=50.)
    

    ax[0,0].set_ylabel('y [Mpc]')
    ax[1,0].set_ylabel('y [Mpc]')
    ax[1,1].set_xlabel('x [Mpc]')
    ax[1,0].set_xlabel('x [Mpc]')

    ax[0,0].set_title('S_MICE')
    ax[0,1].set_title('S0 + e*S2')

    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_S_'+samp+fittype+'.png',bbox_inches='tight')    
    
    # Res
    f, ax = plt.subplots(2, figsize=(3.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0].scatter(x,y,c=S - Se_fit,vmin=-10,cmap='cividis',vmax=10.)
    ax[1].scatter(x,y,c=Sr - Se_fit_red,vmin=-10,cmap='cividis',vmax=10.)

    ax[0].set_ylabel('y [Mpc]')
    ax[1].set_ylabel('y [Mpc]')
    ax[1].set_xlabel('x [Mpc]')

    ax[0].set_title('S_MICE - (S0 + e*S2)')


    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_S_res_'+samp+fittype+'.png',bbox_inches='tight')    
    
    # QUADRUPOLES
    # GTcos
    
    f, ax = plt.subplots(2,2, figsize=(6.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0,0].scatter(x,y,c=GT-gt0,vmin=-10,cmap='cividis',vmax=50.)
    ax[0,1].scatter(x,y,c=gtc*np.cos(2.*theta),vmin=-10,cmap='cividis',vmax=50.)

    ax[1,0].scatter(x,y,c=GTr-gt0,vmin=-10,cmap='cividis',vmax=50.)
    ax[1,1].scatter(x,y,c=gtc_red*np.cos(2.*theta),vmin=-10,cmap='cividis',vmax=50.)
    

    ax[0,0].set_ylabel('y [Mpc]')
    ax[1,0].set_ylabel('y [Mpc]')
    ax[1,1].set_xlabel('x [Mpc]')
    ax[1,0].set_xlabel('x [Mpc]')

    ax[0,0].set_title('GT_MICE')
    ax[0,1].set_title('GT_fit')

    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_GTcos_'+samp+fittype+'.png',bbox_inches='tight')    
    
    f, ax = plt.subplots(2,2, figsize=(6.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0,0].scatter(x,y,c=GT,vmin=-10,cmap='cividis',vmax=50.)
    ax[0,1].scatter(x,y,c=GTfit,vmin=-10,cmap='cividis',vmax=50.)

    ax[1,0].scatter(x,y,c=GTr,vmin=-10,cmap='cividis',vmax=50.)
    ax[1,1].scatter(x,y,c=GTfit_red,vmin=-10,cmap='cividis',vmax=50.)
    

    ax[0,0].set_ylabel('y [Mpc]')
    ax[1,0].set_ylabel('y [Mpc]')
    ax[1,1].set_xlabel('x [Mpc]')
    ax[1,0].set_xlabel('x [Mpc]')

    ax[0,0].set_title('GT_MICE')
    ax[0,1].set_title('GT_fit')

    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_GT_'+samp+fittype+'.png',bbox_inches='tight')    
    
    # Res
    f, ax = plt.subplots(2, figsize=(3.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0].scatter(x,y,c=GT - GTfit,vmin=-10,cmap='cividis',vmax=10.)
    ax[1].scatter(x,y,c=GTr - GTfit_red,vmin=-10,cmap='cividis',vmax=10.)

    ax[0].set_ylabel('y [Mpc]')
    ax[1].set_ylabel('y [Mpc]')
    ax[1].set_xlabel('x [Mpc]')

    ax[0].set_title('GT - (gt0 + e*gt2)')


    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_GT_res_'+samp+fittype+'.png',bbox_inches='tight')    



    # QUADRUPOLES
    # GX
    
    f, ax = plt.subplots(2,2, figsize=(6.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0,0].scatter(x,y,c=GX,vmin=-10,cmap='cividis',vmax=10.)
    ax[0,1].scatter(x,y,c=GXfit,vmin=-10,cmap='cividis',vmax=10.)

    ax[1,0].scatter(x,y,c=GXr,vmin=-10,cmap='cividis',vmax=10.)
    ax[1,1].scatter(x,y,c=GXfit_red,vmin=-10,cmap='cividis',vmax=10.)
    

    ax[0,0].set_ylabel('y [Mpc]')
    ax[1,0].set_ylabel('y [Mpc]')
    ax[1,1].set_xlabel('x [Mpc]')
    ax[1,0].set_xlabel('x [Mpc]')

    ax[0,0].set_title('GX_MICE')
    ax[0,1].set_title('GX_fit')

    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_GX_'+samp+fittype+'.png',bbox_inches='tight')    
    
    # Res
    f, ax = plt.subplots(2, figsize=(3.5,6), sharex=True, sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)

    im0 = ax[0].scatter(x,y,c=GX - GXfit,vmin=-10,cmap='cividis',vmax=10.)
    ax[1].scatter(x,y,c=GXr - GXfit_red,vmin=-10,cmap='cividis',vmax=10.)

    ax[0].set_ylabel('y [Mpc]')
    ax[1].set_ylabel('y [Mpc]')
    ax[1].set_xlabel('x [Mpc]')

    ax[0].set_title('GX - (e*gx2)')


    f.colorbar(im0, ax=ax, orientation='vertical', fraction=.05)
    f.savefig(folder+'../maps/plots/map_GX_res_'+samp+fittype+'.png',bbox_inches='tight')    

def plt_fitted_dist(samp,fittype,relax=True):
    
    halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    
    f     = fits.open(folder+'fitresults_'+fittype+'_'+p_name)[1].data
    f_red = fits.open(folder+'fitresults_'+fittype+'_reduced_'+p_name)[1].data

    hf     = fits.open(folder+'fitresults_'+fittype+'_'+p_name)[0].header
    hf_red = fits.open(folder+'fitresults_'+fittype+'_reduced_'+p_name)[0].header

    try:
        lM200 = np.percentile(f.lM200[1500:], [16, 50, 84])
        c200 = np.percentile(f.c200[1500:], [16, 50, 84])
    except:
        lM200 = np.array([hf['lM200'],hf['lM200'],hf['lM200']])
        c200 = np.array([hf['c200'],hf['c200'],hf['c200']])
        
    q = np.percentile(f.q[1500:], [16, 50, 84])
    qr = np.percentile(f_red.q[1500:], [16, 50, 84])
    
    Eratio = (2.*halos.K/abs(halos.U))
    
    mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
    
    if relax:
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)
        mhalos = mhalos*mrelax
    
    halos = halos[mhalos]

    qh = halos.b2D/halos.a2D
    qhr = halos.b2Dr/halos.a2Dr

    
    f, ax = plt.subplots(2,3, figsize=(10,6), sharey=True)
    
    ax[0,0].set_title('3D')
    ax[0,0].hist(halos.lgMNFW_rho,np.linspace(13.,14.5,50),histtype='step',label='NFW')
    ax[0,0].hist(halos.lgMEin_rho,np.linspace(13.,14.5,50),histtype='step',label='Einasto')
    ax[0,0].axvline(lM200[0],ls='--')
    ax[0,0].axvline(lM200[1],ls='-')
    ax[0,0].axvline(lM200[2],ls='--')
    ax[0,0].legend(frameon=False)

    ax[1,0].set_title('2D')
    ax[1,0].hist(halos.lgMNFW_S,np.linspace(13.,14.5,50),histtype='step')
    ax[1,0].hist(halos.lgMEin_S,np.linspace(13.,14.5,50),histtype='step')
    ax[1,0].axvline(lM200[0],ls='--')
    ax[1,0].axvline(lM200[1],ls='-')
    ax[1,0].axvline(lM200[2],ls='--')
    
    ax[0,1].set_title('3D')
    ax[0,1].hist(halos.cNFW_rho,np.linspace(1,10,50),histtype='step')
    ax[0,1].hist(halos.cEin_rho,np.linspace(1,10,50),histtype='step')
    ax[0,1].axvline(c200[0],ls='--')
    ax[0,1].axvline(c200[1],ls='-')
    ax[0,1].axvline(c200[2],ls='--')
    
    ax[1,1].set_title('2D')     
    ax[1,1].hist(halos.cNFW_S,np.linspace(1,10,50),histtype='step')
    ax[1,1].hist(halos.cEin_S,np.linspace(1,10,50),histtype='step')
    ax[1,1].axvline(c200[0],ls='--')
    ax[1,1].axvline(c200[1],ls='-')
    ax[1,1].axvline(c200[2],ls='--')
    
    ax[0,2].set_title('standard')
    ax[0,2].hist(qh,np.linspace(0,1,50),histtype='step')
    ax[0,2].hist(qh,np.linspace(0,1,50),histtype='step')
    ax[0,2].axvline(q[0],ls='--')
    ax[0,2].axvline(q[1],ls='-')
    ax[0,2].axvline(q[2],ls='--')
    
    ax[1,2].set_title('reduced')     
    ax[1,2].hist(qhr,np.linspace(0,1,50),histtype='step')
    ax[1,2].hist(qhr,np.linspace(0,1,50),histtype='step')
    ax[1,2].axvline(qr[0],ls='--')
    ax[1,2].axvline(qr[1],ls='-')
    ax[1,2].axvline(qr[2],ls='--')
    
    ax[1,0].set_xlabel('$\log M_{200}$')
    ax[1,1].set_xlabel('$c_{200}$')
    ax[1,2].set_xlabel('$q$')
    
    f.savefig(folder+'/../hist_'+samp+fittype+'.png',bbox_inches='tight')    
    

def try_einasto(samp,RIN,ROUT,fittype='_onlyq',
                          substract = False,component='',
                          terms='1h',pname='Einasto'):
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
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
        
    
    ndots = p.shape[0]

    maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))

    mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]

    CovDSfit  = (cov.COV_ST.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())
    CovDS  = cov.COV_ST.reshape(len(p.Rp),len(p.Rp))
    CovS   = cov.COV_S.reshape(len(p.Rp),len(p.Rp))
    
    CovGT  = cov.COV_GT.reshape(len(p.Rp),len(p.Rp))
    CovGTr = cov.COV_GT_reduced.reshape(len(p.Rp),len(p.Rp))
    CovGTc = cov.COV_GT_control.reshape(len(p.Rp),len(p.Rp))
    
    CovGX  = cov.COV_GX.reshape(len(p.Rp),len(p.Rp))
    CovGXr = cov.COV_GX_reduced.reshape(len(p.Rp),len(p.Rp))
    CovGXc = cov.COV_GX_control.reshape(len(p.Rp),len(p.Rp))

    
    GT  = p.GAMMA_Tcos
    GTr = p.GAMMA_Tcos_reduced
    GTc = p.GAMMA_Tcos_control
    
    GX  = p.GAMMA_Xsin
    GXr = p.GAMMA_Xsin_reduced
    GXc = p.GAMMA_Xsin_control

    rplot = p.Rp
    # '''
    fitEin = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSfit),zmean,'Einasto',1.e14,3.)
    
    qfit = 0.63
    # qfit2h = 0.34
    qfit2h = 0.54
    efit = (1-qfit)/(1+qfit)
    efit2h = (1-qfit2h)/(1+qfit2h)
    
    DS      = Delta_Sigma_Ein_2h(p.Rp,zmean,fitEin.M200,fitEin.c200,fitEin.alpha,cosmo_params=params,terms=terms)
    gt1h,gx1h   = GAMMA_components(p.Rp,zmean,ellip=efit,M200 = fitEin.M200,c200=fitEin.c200,cosmo_params=params,terms='1h',pname='Einasto',alpha=fitEin.alpha)
    gt2h,gx2h   = GAMMA_components(p.Rp,zmean,ellip=efit2h,M200 = fitEin.M200,c200=fitEin.c200,cosmo_params=params,terms='2h',pname='NFW',alpha=fitEin.alpha)
    
    gt = gt1h + gt2h
    gx = gx1h + gx2h
    
    gtr, gxr = gt, gx

    mass = str(np.round(np.log10(fitEin.M200),2))
    c200 = str(np.round(fitEin.c200,2))
    alpha = str(np.round(fitEin.alpha,2))
    qfit = str(np.round(qfit,2))
    qfit_red = qfit


    f, ax_all = plt.subplots(2,2, figsize=(12,8),sharex = True)
    f.subplots_adjust(hspace=0)
    ax,ax1,ax2,ax3 = ax_all[0,0],ax_all[0,1],ax_all[1,0],ax_all[1,1]
    
    ax.set_title(p_name+fittype+component)

    ax.plot(p.Rp,p.DSigma_T,'C1')
    ax.plot(rplot,DS,'C3',label='$\log M_{200}=$'+mass+', $c_{200} = $'+c200+r', $\alpha = $'+alpha)
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
    
    f.savefig(folder+'plots/profile_Ein_'+samp+fittype+component+'_'+terms+'.png',bbox_inches='tight')
