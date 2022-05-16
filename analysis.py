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
# folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles3/'
folder = '../profiles3/'
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
    fitd = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_'+p_name)[1].data

    efit     = (1. - fitpar['q']) / (1. + fitpar['q'])
    efit2h     = (1. - fitpar['q2h']) / (1. + fitpar['q2h'])
    
    DS1h  = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h')
    DS2h  = Delta_Sigma_NFW_2h(rplot,zmean,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h')

    gt1h,gx1h   = GAMMA_components(rplot,zmean,ellip=efit,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='1h',pname=pname)
    gt2h,gx2h   = GAMMA_components(rplot,zmean,ellip=efit2h,M200 = 10**fitpar['lM200'],c200=fitpar['c200'],cosmo_params=params,terms='2h',pname=pname)
    
    
    try:
        
        fitpar_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[0].header
        fitd_red = fits.open(folder+'fitresults'+fittype+component+'_'+str(int(RIN))+'_'+str(int(ROUT))+'_reduced_'+p_name)[1].data
  
        efit_red = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])
        efit_red2h = (1. - fitpar_red['q']) / (1. + fitpar_red['q'])

        gt1hr,gx1hr = GAMMA_components(rplot,zmean,ellip=efit_red,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo_params=params,terms='1h',pname=pname)
        gt2hr,gx2hr = GAMMA_components(rplot,zmean,ellip=efit_red2h,M200 = 10**fitpar_red['lM200'],c200=fitpar_red['c200'],cosmo_params=params,terms='2h',pname=pname)
        
    except:
        
        print('NON reduced')
        efit_red = efit
        efit_red2h = efit2h

        gt1hr,gx1hr = gt1h,gx1h
        gt2hr,gx2hr = gt2h,gx2h
    
    fitout = np.array([p.Rp,DS1h,DS2h,gt1h,gx1h,gt1hr,gx1hr,gt2h,gx2h,gt2hr,gx2hr])

    np.savetxt(folder+'fitprofile'+fittype+'_'+samp+'_'+str(int(RIN))+'_'+str(int(ROUT))+'.cat',fitout,fmt='%10.2f')

def plot_q_dist():
    
    import seaborn as sns
    
    hsamples = ['Lz','Mz','Hz']
    colors = ['C1','C3','C5']
    # hsamples = ['HM_Lz','HM_Mz','HM_Hz']
    # colors = ['C1','C3','C5']
    # hsamples = ['LM_Lz','LM_Mz','LM_Hz']
    # colors = ['C1','C3','C5']
    f, ax = plt.subplots(2,2, figsize=(6,5),sharex = True,sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)
    

    ax[0,0].set_xlim([0,1.15])
    ax[0,0].set_ylim([0,4.53])
    ax[0,0].text(1,4,'HM')
    ax[0,1].text(1,4,'LM')
    ax[0,0].text(0.1,4,'all halos')
    ax[1,0].text(0.1,4,'only relaxed')
    
    for j in range(len(hsamples)):

        ax[0,1].plot([-1]*2,[-1]*2,colors[j],label=hsamples[j])        
        samp = hsamples[j]
        # from individual halos

        # FIRST HM

        h = fits.open(folder+'profile_HM_'+samp+'.fits')[0].header
        print(h['N_LENSES'])
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)

        print(samp)
        print(h['LM_MIN'],h['LM_MAX'])
        print(h['z_min'],h['z_max'])
        
        sns.kdeplot(halos.q2d[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,0])
        ax[1,0].plot([np.mean(halos.q2d[mhalos*mrelax])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,0],ls='--')
        ax[1,0].plot([np.mean(halos.q2dr[mhalos*mrelax])]*2,[4.0,4.4],colors[j]+'--')
        sns.kdeplot(halos.q2d[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,0])
        ax[0,0].plot([np.mean(halos.q2d[mhalos])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,0],ls='--')
        ax[0,0].plot([np.mean(halos.q2dr[mhalos])]*2,[4.0,4.4],colors[j]+'--')

        samp = hsamples[j]
        # from individual halos

        # NOW LM

        h = fits.open(folder+'profile_LM_'+samp+'.fits')[0].header
        print(h['N_LENSES'])
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)

        print(samp)
        print(h['LM_MIN'],h['LM_MAX'])
        print(h['z_min'],h['z_max'])
        
        sns.kdeplot(halos.q2d[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,1])
        ax[1,1].plot([np.mean(halos.q2d[mhalos*mrelax])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,1],ls='--')
        ax[1,1].plot([np.mean(halos.q2dr[mhalos*mrelax])]*2,[4.0,4.4],colors[j]+'--')
        sns.kdeplot(halos.q2d[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,1])
        ax[0,1].plot([np.mean(halos.q2d[mhalos])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,1],ls='--')
        ax[0,1].plot([np.mean(halos.q2dr[mhalos])]*2,[4.0,4.4],colors[j]+'--')

    ax[0,1].plot([-1]*2,[-1]*2,'k--',label='standard')
    ax[0,1].plot([-1]*2,[-1]*2,'k-',label='reduced')
    
    ax[0,1].legend(loc=2,frameon=False)
    ax[1,1].set_xlabel('q')
    ax[1,0].set_xlabel('q')
    f.savefig(folder+'../final_plots/qdist.pdf',bbox_inches='tight')

def extract_params(hsamples,
                   RIN=250,
                   ROUToq = ['2000','1000','2000','1000'],
                   reduced=False,
                   cornplot=False):
    
    
    qh   = []
    qhr   = []
    NFW_h = []
    Ein_h = []
    NFW_hr = []
    Ein_hr = []
        
        
    NFW = []
    Ein = []
    o1h = []
    woc = []

    eNFW = []
    eNFWr = []
    eEin = []
    eo1h = []
    ewoc = []
    
    if reduced:
        ang = '_reduced'
        q2d = 'q2dr'
    else:
        ang = ''
        q2d = 'q2d'
    
    for j in range(len(hsamples)):
        
        samp = hsamples[j]
        # from individual halos

        h = fits.open(folder+'profile_'+samp+'.fits')[0].header
        
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)

        print(samp)
        print(h['LM_MIN'],h['LM_MAX'])
        print(h['z_min'],h['z_max'],h['z_mean'])
        print(h['N_LENSES'])
        print(mhalos.sum())
        print((mhalos*mrelax).sum())
        
        qhr   += [np.mean(halos[q2d][mhalos*mrelax])]
        qh    += [np.mean(halos[q2d][mhalos])]
               
    
        mfit_NFW = (halos.cNFW_rho > 1.)*(halos.cNFW_S > 1.)*(halos.cNFW_rho < 10.)*(halos.cNFW_S < 10.)*(halos.lgMNFW_rho > 12)*(halos.lgMNFW_S > 12)
        mfit_Ein = (halos.cEin_rho > 1.)*(halos.cEin_S > 1.)*(halos.cEin_rho < 10.)*(halos.cEin_S < 10.)*(halos.lgMEin_rho > 12)*(halos.lgMEin_S > 12)*(halos.alpha_rho > 0.)*(halos.alpha_S > 0.)*(halos.alpha_rho < 0.7)*(halos.alpha_S < 0.7)
        mhalos = mhalos*mfit_Ein*mfit_NFW
    
        halos_samp = halos[mhalos]
        
        # qwh   += [h['q2d_mean']]
        lMNFW = np.mean(halos[mhalos].lgMNFW_rho)
        cNFW  = np.mean(halos[mhalos].cNFW_rho)
        lMEin = np.mean(halos[mhalos].lgMEin_rho)
        cEin  = np.mean(halos[mhalos].cEin_rho)
        alpha = np.mean(halos[mhalos].alpha_rho)

        NFW_h += [[lMNFW,cNFW]]
        Ein_h += [[lMEin,cEin,alpha]]


        lMNFWr = np.mean(halos[mhalos*mrelax].lgMNFW_rho)
        cNFWr  = np.mean(halos[mhalos*mrelax].cNFW_rho)
        lMEinr = np.mean(halos[mhalos*mrelax].lgMEin_rho)
        cEinr  = np.mean(halos[mhalos*mrelax].cEin_rho)
        alphar = np.mean(halos[mhalos*mrelax].alpha_rho)
        
        NFW_hr += [[lMNFWr,cNFWr]]
        Ein_hr += [[lMEinr,cEinr,alphar]]


        # 1 halo
        
        fstd = fits.open(folder+'fitresults_onlyq_'+RIN[j]+'_'+ROUToq[j]+ang+'_profile_'+samp+'.fits')[1].data
        
        o1h  += [[np.median(fstd.lM200[1500:]),
                 np.median(fstd.q[1500:]),
                 np.median(fstd.c200[1500:])]]
                 
        eo1h += [[np.diff(np.percentile(fstd.lM200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.c200[1500:], [16,50,84]))]]
                 
        if cornplot:
            
            lM200 = np.median(fstd.lM200[1500:])
            labels_DS   = ['$\log M_{200}$','$c_{200}$']
        
            mcmc_DS  = np.array([fstd.lM200[1500:],fstd.c200[1500:]]).T
        
            fds = corner.corner(mcmc_DS,labels=labels_DS,smooth=1.,range=[(lM200-0.07,lM200+0.07),(2.,4.5)],truths=[lMNFWr,cNFWr],label_kwargs=({'fontsize':16}),truth_color='C2',quantiles=(0.16, 0.84))
            


        # without_c

        fstd = fits.open(folder+'fitresults_2h_2q_woc_'+RIN[j]+'_5000'+ang+'_profile_'+samp+'.fits')[1].data

        woc  += [[np.median(fstd.lM200[1500:]),
                 np.median(fstd.q[1500:]),
                 np.median(fstd.q2h[1500:])]]
                 
        ewoc += [[np.diff(np.percentile(fstd.lM200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q2h[1500:], [16,50,84]))]]
               
        # Einasto
        fstd = fits.open(folder+'fitresults_2h_2q_Ein_'+RIN[j]+'_5000'+ang+'_profile_'+samp+'.fits')[1].data

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
        fstd = fits.open(folder+'fitresults_2h_2q_'+RIN[j]+'_5000'+ang+'_profile_'+samp+'.fits')[1].data

        NFW  += [[np.median(fstd.lM200[1500:]),
                 np.median(fstd.q[1500:]),
                 np.median(fstd.c200[1500:]),
                 np.median(fstd.q2h[1500:])]]
                 
        eNFW += [[np.diff(np.percentile(fstd.lM200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.c200[1500:], [16,50,84])),
                 np.diff(np.percentile(fstd.q2h[1500:], [16,50,84]))]]

               
    return qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,[NFW,eNFW],[Ein,eEin],[o1h,eo1h],[woc,ewoc]

      
        

def plot_bias(hsamps,lhs,cstyle,nplot,RIN,ROUToq):
    
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc = extract_params(hsamps,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r = extract_params(hsamps,RIN,ROUToq,reduced=True)
    ###########
    #   q1h
    ###########
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_R]
    
    param = 1

    f, ax = plt.subplots(2,1, figsize=(16,14))

    for hs in range(len(hsamps)):
        ax[0].errorbar(qhr[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs],mfc='none')
        ax[0].errorbar(qhr_r[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs],mfc='none')
        ax[0].errorbar(qh[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
        ax[0].errorbar(qh_r[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])

        ax[1].errorbar(qhr[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFWr[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs],mfc='none')
        ax[1].errorbar(qhr_r[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFWr[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs],mfc='none')
        ax[1].errorbar(qh[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFWr[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
        ax[1].errorbar(qh_r[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFWr[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
    
    ax[0].plot([0.5,0.8],[0.5,0.8],'C7--')
    ax[1].plot([0.5,0.8],[0.5,0.8],'C7--')
    
    xl = ['NFW - 1h+2h','Ein - 1h+2h','NFW - 1h+2h - fix $c_{200}$','NFW 1h']
    # FOM
    
    f,ax = plt.subplots(figsize=(14,3))
    plt.plot([0,5],[0,0],'C7--')
    
    ax.axhspan(-0.05,0.05,0,5,color='C7',alpha=0.5)

    
    for hs in range(len(hsamps)):
        for fp in range(4):
            if fp == 0:
                plt.errorbar(fp+1+0.1*hs,(fq[fp][0][hs][param]-qhr[hs])/qhr[hs],
                            yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
            else:
                plt.errorbar(fp+1+0.1*hs,(fq[fp][0][hs][param]-qhr[hs])/qhr[hs],
                            yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                            fmt=cstyle[hs],markersize=10)
    
    plt.legend(frameon = False)
    plt.ylabel(r'$(\tilde{q}_{1h}-\langle q \rangle)/\langle q \rangle$')
    plt.axis([0,5,-0.12,0.12])
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
    
    for hs in range(len(hsamps)):
        for fp in range(3):
            if fp == 0:
                plt.errorbar(fp+1+0.1*hs,fq[fp][0][hs][param],
                            yerr=np.array([fq[fp][1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
            else:
                plt.errorbar(fp+1+0.1*hs,fq[fp][0][hs][param],
                            yerr=np.array([fq[fp][1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10)

    
    plt.legend(frameon = False)
    plt.ylabel(r'$\tilde{q}_{2h}/\langle q \rangle$')
    ax.set_xticks(np.arange(3)+1)
    ax.set_xticklabels(xl[:-1])
    # plt.axis([0,5,0.3,1.12])
    # f.savefig(folder+'../final_plots/model_q2h.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_q2h_'+nplot+'.png',bbox_inches='tight')
    
    ###########
    #   lM200
    ###########
    
    
    f,ax = plt.subplots(figsize=(14,3))
    plt.plot([0,5],[1,1],'C7--')
    
    param = 0
    
    ax.axhspan(0.95,1.05,0,5,color='C7',alpha=0.5)
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            diff  = 10**(fq[fp][0][hs][param] - NFW_h[hs][param])
            diff2 = 10**(fq[fp][0][hs][param] - NFW_hr[hs][param])
            if fp == 0:
                plt.errorbar(fp+1+0.1*hs,(diff),
                            yerr=np.array([fq[fp][1][hs][param]*np.log(10)*diff]).T,
                            fmt=cstyle[hs],markersize=10,label=lhs[hs])
                plt.errorbar(fp+1+0.1*hs,(diff2),
                            yerr=np.array([fq[fp][1][hs][param]*np.log(10)*diff2]).T,
                            fmt=cstyle[hs],markersize=5,label=lhs[hs])
            else:
                plt.errorbar(fp+1+0.1*hs,(diff),
                            yerr=np.array([fq[fp][1][hs][param]*np.log(10)*diff]).T,
                            fmt=cstyle[hs],markersize=10)
                plt.errorbar(fp+1+0.1*hs,(diff2),
                            yerr=np.array([fq[fp][1][hs][param]*np.log(10)*diff2]).T,
                            fmt=cstyle[hs],markersize=5)
                            

    
    plt.legend(frameon = False,loc=3)
    plt.ylabel(r'$\tilde{M_{200}}/M_{200}$')
    ax.set_xticks(np.arange(4)+1)
    ax.set_xticklabels(xl)
    # plt.axis([0,5,0.8,1.2])
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
    for hs in range(len(hsamps)):
        for fp in range(4):
            diff = 10**(fq[fp][0][hs][0] - NFW_h[hs][0])
            if fp == 0:
                plt.errorbar(diff,fq[fp][0][hs][param]/qshr[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qshr[hs]]).T,
                        fmt=cstyle[hs],markersize=10,label=lhs[hs])
            else:
                plt.errorbar(diff,fq[fp][0][hs][param]/qshr[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qshr[hs]]).T,
                        fmt=cstyle[hs],markersize=10)
                        
    # lMlC = 14.0544
    # elMlC = [0.00343227, 0.00317295]
                        
    # lMlsq = 14.0456
    # diffCl = 10**(lMlC - NFW_h[2][0])
    # difflsq = 10**(lMlsq - NFW_h[2][0])
    
    # plt.errorbar(diffCl,qlC/qwh[2],yerr=np.array([eqlC]).T,fmt='C3s',markersize=10)
    # plt.errorbar(difflsq,qlsq/qwh[2],yerr=np.array([eqlsq]).T,fmt='C3o',markersize=10)

    
    plt.legend(frameon = False)
    plt.xlabel(r'$\tilde{M_{200}}/M_{200}$')
    plt.ylabel(r'$\tilde{q}/\langle q \rangle$')
    plt.axis([0.8,1.2,0.8,1.12])
    # f.savefig(folder+'../final_plots/model_ratioq_M200.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_ratioq_M200_'+nplot+'.png',bbox_inches='tight')

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
          
    rplot,DS1h,DS2h,gt1h,gx1h,gt1hr,gx1hr,gt2h,gx2h,gt2hr,gx2hr = np.loadtxt(folder+'fitprofile'+fittype+'_'+samp+'_'+str(int(RIN))+'_'+str(int(ROUT))+'.cat')

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



# '''
lhs = ['HM-Lz','LM-Lz','HM-Hz','LM-Hz']
cstyle = ['C0^','C0v','C3^','C3v']
hsamps = ['HM_Lz_relaxed','LM_Lz_relaxed','HM_Hz_relaxed','LM_Hz_relaxed']
hsamps_nr = ['HM_Lz','LM_Lz','HM_Hz','LM_Hz']
ROUToq = ['2000','1000','2000','1000']

lhs_ext_rel = ['HM-Lz','LM-Lz','HM-Mz','LM-Mz','HM-Hz','LM-Hz','HM-HHz','LM-HHz']
cstyle_ext_rel = ['C0^','C0v','C1^','C1v','C3^','C3v','C4^','C4v']
ROUToq_ext_rel = ['2000','1000','2000','1000','2000','1000','2000','1000']

lhs_ext = ['HM-Lz','LM-Lz','HM-Mz','LM-Mz','HM-Hz','LM-Hz']
cstyle_ext = ['C1^','C1v','C3^','C3v','C5^','C5v']
ROUToq_ext = ['2000','1000','2000','1000','2000','1000']
RIN_mix00 = ['200','200','300','300','400','400']
RIN_mix50 = ['250','250','350','350','450','450']
RIN_mix = ['350','350','350','350','450','400']
RIN_mix350 = ['350','350','350','350','350','350']
RIN_mix400 = ['400']*6


hsamps_mix = ['HM_Lz','LM_Lz','HM_Mz','LM_Mz','HM_HHz_relaxed','LM_HHz_relaxed']
hsamps_ext = ['HM_Lz','LM_Lz','HM_Mz','LM_Mz','HM_HHz','LM_HHz']
hsamps = ['HM_Lz','LM_Lz','HM_Mz','LM_Mz','HM_Hz','LM_Hz']
hsamps_mis20 = ['HM_Lz_mis20','LM_Lz_mis20','HM_Mz_mis20','LM_Mz_mis20','HM_Hz_mis20','LM_Hz_mis20']
hsamps_miscen = ['HM_Lz_miscen','LM_Lz_miscen','HM_Mz_miscen','LM_Mz_miscen','HM_Hz_miscen','LM_Hz_miscen']
hsamps_misall = ['HM_Lz_mis20_miscen','LM_Lz_mis20_miscen','HM_Mz_mis20_miscen','LM_Mz_mis20_miscen','HM_Hz_mis20_miscen','LM_Hz_mis20_miscen']

hsamps_ext_rel = ['HM_Lz_relaxed','LM_Lz_relaxed',
                  'HM_Mz_relaxed','LM_Mz_relaxed',
                  'HM_Hz_relaxed','LM_Hz_relaxed',
                  'HM_HHz_relaxed','LM_HHz_relaxed']

# plot_bias(hsamps_mix2,lhs_ext,cstyle_ext,'mix00',RIN_mix00,ROUToq_ext,False)
# plot_bias(hsamps_mix2,lhs_ext,cstyle_ext,'mix50',RIN_mix50,ROUToq_ext,False)
# plot_bias(hsamps_mix2,lhs_ext,cstyle_ext,'mix350',RIN_mix350,ROUToq_ext,False)
# plot_bias(hsamps_mix2,lhs_ext,cstyle_ext,'mix',RIN_mix,ROUToq_ext,False)
plot_bias(hsamps,lhs_ext,cstyle_ext,'mix',RIN_mix,ROUToq_ext)
# plot_bias(hsamps_mis20,lhs_ext,cstyle_ext,'mix_mis20',RIN_mix,ROUToq_ext)
# plot_bias(hsamps_miscen,lhs_ext,cstyle_ext,'mix_miscen',RIN_mix,ROUToq_ext)
# plot_bias(hsamps_misall,lhs_ext,cstyle_ext,'mix_misall',RIN_mix,ROUToq_ext)
# '''
# plot_bias(hsamps_nr,lhs,cstyle,'nonrex_comprel_samps',250,ROUToq,True)
# plot_bias(hsamps_ext,lhs_ext,cstyle_ext,'final',350,ROUToq_ext,False)
# test_fitting(hsamps_mix2,RIN_mix2,ROUToq_ext,False)

# plot_bias(hsamps,lhs,cstyle,'relax_samps',250,ROUToq,True)
# plot_bias(hsamps_ext,lhs_ext,cstyle_ext,'extend_relax_samps',400,ROUToq_ext,True)

# test_fitting(hsamps_ext,400,ROUToq_ext,True)
# test_fitting(hsamps,250,ROUToq,True)
# test_fitting(hsamps_nr,250,ROUToq,False)
# '''

# for j in range(len(lhs_ext)):
        # samp = hsamps_mix2[j]
        # RIN = float(RIN_mix[j])
        # save_fitted(samp,RIN,5000,fittype='_2h_2q')
        # save_fitted(samp,RIN,5000,fittype='_2h_2q_Ein')
        # save_fitted(samp,RIN,5000,fittype='_2h_2q_woc')

'''
f, ax_all = plt.subplots(6,3, figsize=(16,14),sharex = True)
f.subplots_adjust(hspace=0)


for j in range(len(ax_all)):
    plt_profile_fitted_final(hsamps_ext[j],350,5000,ax_all[j],fittype='_2h_2q_Ein')
    ax_all[j,0].text(1,100,lhs_ext[j],fontsize=14)

ax_all[0,0].legend(loc=3,frameon=False)
ax_all[0,1].legend(loc=3,frameon=False)


    
f.savefig(folder+'../test_plots/profile_relaxed_Ein.png',bbox_inches='tight')


f, ax_all = plt.subplots(1,3, figsize=(16,4),sharex = True)
f.subplots_adjust(hspace=0)

plt_profile_fitted_final('HM_Lz',250,5000,ax_all,fittype='_2h_2q')


ax_all[0].legend(loc=3,frameon=False)

    
f.savefig(folder+'../test_plots/final/profile_HM_Lz.png',bbox_inches='tight')
'''
