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
halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data        
Eratio = (2.*halos.K/abs(halos.U))


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

def extract_params(hsamples,RIN=250,relax=True):
    
    RIN = str(RIN)
        
    ROUT = ['2000','1000','2000','1000']
    
    qwh   = []
    NFW_h = []
    Ein_h = []
        
    NFW = []
    Ein = []
    o1h = []
    woc = []

    eNFW = []
    eEin = []
    eo1h = []
    ewoc = []
    
    
    for j in range(len(hsamples)):
        
        samp = hsamples[j]
        print(samp)

        # from individual halos

        h = fits.open(folder+'profile_'+samp+'.fits')[0].header
    
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mfit_NFW = (halos.cNFW_rho > 1.)*(halos.cNFW_S > 1.)*(halos.cNFW_rho < 10.)*(halos.cNFW_S < 10.)*(halos.lgMNFW_rho > 12)*(halos.lgMNFW_S > 12)
        mfit_Ein = (halos.cEin_rho > 1.)*(halos.cEin_S > 1.)*(halos.cEin_rho < 10.)*(halos.cEin_S < 10.)*(halos.lgMEin_rho > 12)*(halos.lgMEin_S > 12)*(halos.alpha_rho > 0.)*(halos.alpha_S > 0.)*(halos.alpha_rho < 0.7)*(halos.alpha_S < 0.7)
        mhalos = mhalos*mfit_Ein*mfit_NFW
    
        if relax:
            mrelax = (halos.offset < 0.1)*(Eratio < 1.35)
            mhalos = mhalos*mrelax
    
        halos_samp = halos[mhalos]
    
        qwh   += [h['q2d_mean']]
        # qwh   += [np.mean(halos_samp.q2d)]
        lMNFW = np.mean(halos_samp.lgMNFW_rho)
        cNFW  = np.mean(halos_samp.cNFW_rho)
        lMEin = np.mean(halos_samp.lgMEin_rho)
        cEin  = np.mean(halos_samp.cEin_rho)
        alpha = np.mean(halos_samp.alpha_rho)
        
        NFW_h += [[lMNFW,cNFW]]
        Ein_h += [[lMEin,cEin,alpha]]


        # 1 halo
        
        fstd = fits.open(folder+'fitresults_onlyq_'+RIN+'_'+ROUT[j]+'_profile_'+samp+'.fits')[1].data
        
        o1h  += [[np.median(fstd.lM200[1500:]),
                 np.median(fstd.q[1500:]),
                 np.median(fstd.c200[1500:])]]
                 
        eo1h += [[np.diff(np.percentile(fstd.lM200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.c200[1500:], [16,50,84]))]]

        # without_c

        fstd = fits.open(folder+'fitresults_2h_2q_woc_250_5000_profile_'+samp+'.fits')[1].data

        woc  += [[np.median(fstd.lM200[1500:]),
                 np.median(fstd.q[1500:]),
                 np.median(fstd.q2h[1500:])]]
                 
        ewoc += [[np.diff(np.percentile(fstd.lM200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q2h[1500:], [16,50,84]))]]
               
        # Einasto
        fstd = fits.open(folder+'fitresults_2h_2q_Ein_'+RIN+'_5000_profile_'+samp+'.fits')[1].data

        Ein  += [[np.median(fstd.lM200[1500:]),
                 np.median(fstd.q[1500:]),
                 np.median(fstd.c200[1500:]),
                 np.median(fstd.alpha[1500:]),
                 np.median(fstd.q2h[1500:])]]
                 
        eEin += [[np.diff(np.percentile(fstd.lM200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.c200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.alpha[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q2h[1500:], [16,50,84]))]]
        
        # NFW
        fstd = fits.open(folder+'fitresults_2h_2q_'+RIN+'_5000_profile_'+samp+'.fits')[1].data

        NFW  += [[np.median(fstd.lM200[1500:]),
                 np.median(fstd.q[1500:]),
                 np.median(fstd.c200[1500:]),
                 np.median(fstd.q2h[1500:])]]
                 
        eNFW += [[np.diff(np.percentile(fstd.lM200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.c200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q2h[1500:], [16,50,84]))]]
               
    return qwh,NFW_h,Ein_h,[NFW,eNFW],[Ein,eEin],[o1h,eo1h],[woc,ewoc]

def test_fitting(hsamples,RIN=250,relax=True):
    
    RIN = str(RIN)
        
    ROUT = ['2000','1000','2000','1000']
    
    
    for j in range(len(hsamples)):

        samp = hsamples[j]
        print(samp)

        f, ax = plt.subplots(2,1, figsize=(12,6),sharex = True)
        f.subplots_adjust(hspace=0)
        ax[0].set_title(samp)
        ax[0].axhspan(0.95,1.05,0,5,color='C7',alpha=0.5)
        ax[1].axhspan(0.95,1.05,0,5,color='C7',alpha=0.5)

        # alternative fit

        profile = fits.open(folder+'profile_'+samp+'.fits')
        h   = profile[0].header
        p   = profile[1].data
        cov = profile[2].data
        
        maskr   = (p.Rp > (float(RIN)/1000.))*(p.Rp < (float(ROUT[j])/1000.))
        mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]
        CovDSf = (cov.COV_ST.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())

        fitNFWDS = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSf),zmean)

        # from individual halos
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mfit_NFW = (halos.cNFW_rho > 1.)*(halos.cNFW_S > 1.)*(halos.cNFW_rho < 10.)*(halos.cNFW_S < 10.)*(halos.lgMNFW_rho > 12)*(halos.lgMNFW_S > 12)
        mfit_Ein = (halos.cEin_rho > 1.)*(halos.cEin_S > 1.)*(halos.cEin_rho < 10.)*(halos.cEin_S < 10.)*(halos.lgMEin_rho > 12)*(halos.lgMEin_S > 12)*(halos.alpha_rho > 0.)*(halos.alpha_S > 0.)*(halos.alpha_rho < 0.7)*(halos.alpha_S < 0.7)
        mhalos = mhalos*mfit_Ein*mfit_NFW
    
        if relax:
            mrelax = (halos.offset < 0.1)*(Eratio < 1.35)
            mhalos = mhalos*mrelax
    
        halos_samp = halos[mhalos]
    
        lMNFW = np.mean(halos_samp.lgMNFW_rho)
        cNFW  = np.mean(halos_samp.cNFW_rho)
        lMEin = np.mean(halos_samp.lgMEin_rho)
        cEin  = np.mean(halos_samp.cEin_rho)
        alpha = np.mean(halos_samp.alpha_rho)

        # 1 halo
        fstd = fits.open(folder+'fitresults_onlyq_'+RIN+'_'+ROUT[j]+'_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMNFW),label='1halo',alpha=0.5)
        ax[1].plot(fstd.c200[1500:]/cNFW,label='1halo',alpha=0.5)
                       
        # Einasto
        fstd = fits.open(folder+'fitresults_2h_2q_Ein_'+RIN+'_5000_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMEin),label='Einasto',alpha=0.5)
        ax[1].plot(fstd.c200[1500:]/cEin,label='Einasto',alpha=0.5)
        
        # NFW
        fstd = fits.open(folder+'fitresults_2h_2q_'+RIN+'_5000_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMNFW),label='NFW',alpha=0.5)
        ax[1].plot(fstd.c200[1500:]/cNFW,label='NFW',alpha=0.5)

        # without_c
        fstd = fits.open(folder+'fitresults_2h_2q_woc_250_5000_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMNFW),label='fix c',alpha=0.5)
        
        ax[0].axhline(fitNFWDS.M200/10**lMNFW,color='C3',lw=2)
        ax[1].axhline(fitNFWDS.c200/cNFW,color='C3',lw=2)
        ax[1].axhline(cEin/cNFW,color='C1',lw=2)
        
        ax[0].legend(frameon=False,ncol=4,loc=3)
        
        ax[0].set_ylim([0.6,1.15])
        ax[1].set_ylim([0.6,1.15])
        # f.savefig(folder+'../final_plots/test_fit_RIN'+RIN+'_'+samp+'.png',bbox_inches='tight')
        f.savefig(folder+'../test_plots/test_fit_RIN'+RIN+'_'+samp+'.png',bbox_inches='tight')


def plot_bias(hsamps,lhs,cstyle,nplot):
    qwh,NFW_h,Ein_h,NFW,Ein,o1h,woc = extract_params(hsamps)
    ###########
    #   q1h
    ###########
    fq = [NFW,Ein,woc,o1h]
    
    
    f,ax = plt.subplots(figsize=(14,3))
    plt.plot([0,5],[1,1],'C7--')
    
    
    xl = ['NFW','Ein','fix $c_{200}$','1h']
    
    param = 1
    
    ax.axhspan(0.95,1.05,0,5,color='C7',alpha=0.5)
    
    for hs in range(4):
        for fp in range(4):
            if fp == 0:
                plt.errorbar(fp+1+0.1*hs,fq[fp][0][hs][param]/qwh[hs],
                            yerr=np.array([fq[fp][1][hs][param]/qwh[hs]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
            else:
                plt.errorbar(fp+1+0.1*hs,fq[fp][0][hs][param]/qwh[hs],
                            yerr=np.array([fq[fp][1][hs][param]/qwh[hs]]).T,
                            fmt=cstyle[hs],markersize=10)
    
    plt.legend(frameon = False)
    plt.ylabel(r'$\tilde{q}_{1h}/\langle q \rangle$')
    ax.set_xticks(np.arange(4)+1)
    ax.set_xticklabels(xl)
    # f.savefig(folder+'../final_plots/model_q1h.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_q1h_'+nplot+'.png',bbox_inches='tight')
    
    
    ###########
    #   q2h
    ###########
    
    f,ax = plt.subplots(figsize=(14,3))
    plt.plot([0,5],[1,1],'C7--')
    
    param = -1
    
    ax.axhspan(0.95,1.05,0,5,color='C7',alpha=0.5)
    
    for hs in range(4):
        for fp in range(3):
            if fp == 0:
                plt.errorbar(fp+1+0.1*hs,fq[fp][0][hs][param]/qwh[hs],
                            yerr=np.array([fq[fp][1][hs][param]/qwh[hs]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
            else:
                plt.errorbar(fp+1+0.1*hs,fq[fp][0][hs][param]/qwh[hs],
                            yerr=np.array([fq[fp][1][hs][param]/qwh[hs]]).T,
                            fmt=cstyle[hs],markersize=10)
    
    plt.legend(frameon = False)
    plt.ylabel(r'$\tilde{q}_{2h}/\langle q \rangle$')
    ax.set_xticks(np.arange(3)+1)
    ax.set_xticklabels(xl[:-1])
    # f.savefig(folder+'../final_plots/model_q2h.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_q2h_'+nplot+'.png',bbox_inches='tight')
    
    ###########
    #   lM200
    ###########
    
    f,ax = plt.subplots(figsize=(14,3))
    plt.plot([0,5],[1,1],'C7--')
    
    param = 0
    
    ax.axhspan(0.95,1.05,0,5,color='C7',alpha=0.5)
    
    for hs in range(4):
        for fp in range(4):
            diff = 10**(fq[fp][0][hs][param] - NFW_h[hs][param])
            if fp == 0:
                plt.errorbar(fp+1+0.1*hs,(diff),
                            yerr=np.array([fq[fp][1][hs][param]*np.log(10)*diff]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
            else:
                plt.errorbar(fp+1+0.1*hs,(diff),
                            yerr=np.array([fq[fp][1][hs][param]*np.log(10)*diff]).T,
                            fmt=cstyle[hs],markersize=10)
    
    plt.legend(frameon = False)
    plt.ylabel(r'$\tilde{M_{200}}/M_{200}$')
    ax.set_xticks(np.arange(4)+1)
    ax.set_xticklabels(xl)
    # f.savefig(folder+'../final_plots/model_M200.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_M200_'+nplot+'.png',bbox_inches='tight')
    
    ##############
    #   lM200/q
    ##############
    
    
    f,ax = plt.subplots(figsize=(14,3))
    # plt.plot([0,1],[1,1],'C7--')
    
    param = 1
    
    ax.axhspan(0.95,1.05,0,1,color='C7',alpha=0.5)
    ax.axvspan(0.95,1.05,ymin=0.0,ymax=1.,color='C7',alpha=0.5)
    plt.axhline(1,color='C7',ls='--')
    plt.axvline(1,color='C7',ls='--')
    for hs in range(4):
        for fp in range(4):
            diff = 10**(fq[fp][0][hs][0] - NFW_h[hs][0])
            if fp == 0:
                plt.errorbar(diff,fq[fp][0][hs][param]/qwh[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qwh[hs]]).T,
                        fmt=cstyle[hs],markersize=10,label=lhs[hs])
            else:
                plt.errorbar(diff,fq[fp][0][hs][param]/qwh[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qwh[hs]]).T,
                        fmt=cstyle[hs],markersize=10)
    
    plt.legend(frameon = False)
    plt.xlabel(r'$\tilde{M_{200}}/M_{200}$')
    plt.ylabel(r'$\tilde{q}/\langle q \rangle$')
    # f.savefig(folder+'../final_plots/model_ratioq_M200.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_ratioq_M200_'+nplot+'.png',bbox_inches='tight')


# test_fitting(['HM_Lz_relaxed','LM_Lz_relaxed','HM_Hz_relaxed','LM_Hz_relaxed'])

# '''
lhs = ['HM-Lz','LM-Lz','HM-Hz','LM-Hz']
cstyle = ['C0^','C0v','C3^','C3v']
hsamps = ['HM_Lz_relaxed','LM_Lz_relaxed','HM_Hz_relaxed','LM_Hz_relaxed']

plot_bias(hsamps,lhs,cstyle,'original_samps')

