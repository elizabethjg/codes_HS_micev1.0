import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from models_profiles import *
from astropy.constants import G,c,M_sun, pc
from fit_models_colossus import *

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)
import corner
folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'

plot_path = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/plot_paper_HSMice/'    

RIN  = 250.
ROUT = 2500.


def fit_profiles(samp,axes,relax=True):
    
    axDS,axM,axc,axa,axM2d,axc2d,axa2d = axes    
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
    
    params = {'flat': True, 'H0': 70.0, 'Om0': 0.25, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.95}

    
    zmean = h['z_mean']
            
    
    ndots = p.shape[0]

    maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))

    mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]

    CovDS  = cov.COV_ST.reshape(len(p),len(p))
    CovS   = cov.COV_S.reshape(len(p),len(p))

    CovDSf = (cov.COV_ST.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())
    CovSf = (cov.COV_S.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())

    fitNFWDS = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSf),zmean)
    fitEinDS = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSf),zmean,'Einasto',1.e14,3.)
    # fitEinDS  = Sigma_fit(p.Rp[maskr],p.Sigma[maskr],np.diag(CovSf),zmean,'Einasto',1.e14,3.)
    # fitNFWDS  = Sigma_fit(p.Rp[maskr],p.Sigma[maskr],np.diag(CovSf),zmean)

    out = np.array([fitEinDS.xplot,fitNFWDS.yplot,fitEinDS.yplot])
    
    np.savetxt(plot_path+'profile_'+samp+'.fit',out,fmt='%12.3f')

    # f,axDS = plt.subplots()
    axDS.plot(p.Rp,p.DSigma_T,'C7')
    axDS.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C7',alpha=0.4)
    # axDS.fill_between(p.Rp,p.Sigma+np.diag(CovS),p.Sigma-np.diag(CovS),color='C7',alpha=0.4)
    # axDS.plot(p.Rp,p.Sigma,'C7')
    axDS.plot(fitEinDS.xplot,fitEinDS.yplot,'seagreen',label='Einasto')
    axDS.plot(fitNFWDS.xplot,fitNFWDS.yplot,'orangered',label='NFW')
    axDS.set_xscale('log')
    axDS.set_yscale('log')
    axDS.set_ylabel(r'$\Delta\Sigma [M_\odot h/pc^2]$',labelpad=2)
    axDS.set_xlabel('r [$h^{-1}$ Mpc]')
    axDS.set_ylim(2,200)
    axDS.set_xlim(0.1,5)
    axDS.xaxis.set_ticks([0.1,1,5])
    axDS.set_xticklabels([0.1,1,5])
    axDS.yaxis.set_ticks([5,10,100])
    axDS.set_yticklabels([5,10,100])
    axDS.axvline(RIN/1000.,color='k',ls=':')
    axDS.axvline(ROUT/1000.,color='k',ls=':')
    # axDS.legend()
   
    # f,axS = plt.subplots()
    # axS.plot(p.Rp,p.Sigma,'C7')
    # axS.plot(fitEinS.xplot,fitEinS.yplot,'seagreen',label='Einasto')
    # axS.plot(fitNFWS.xplot,fitNFWS.yplot,'orangered',label='NFW')
    # axS.fill_between(p.Rp,p.Sigma+np.diag(CovS),p.Sigma-np.diag(CovS),color='C7',alpha=0.4)
    # axS.set_xscale('log')
    # axS.set_yscale('log')
    # axS.set_ylabel(r'$\Sigma$',labelpad=2)
    # axS.set_xlabel('r [$h^{-1}$ Mpc]')
    # axS.set_ylim(7,500)
    # axS.set_xlim(0.1,5)
    # axS.xaxis.set_ticks([0.1,1,5])
    # axS.set_xticklabels([0.1,1,5])
    # axS.yaxis.set_ticks([10,100])
    # axS.set_yticklabels([10,100])
    # axS.axvline(RIN/1000.,color='k',ls=':')
    # axS.axvline(ROUT/1000.,color='k',ls=':')
    
    out = np.array([np.log10(fitEinDS.M200),fitEinDS.c200,fitEinDS.alpha,fitEinDS.res,
                    np.log10(fitNFWDS.M200),fitNFWDS.c200,fitNFWDS.res])
    
    np.savetxt(folder+'../fitlens.res',out,fmt='%10.2f')
    
    
    halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data
                
    
    Eratio = (2.*halos.K/abs(halos.U))
    T       = (1. - halos.q**2)/(1. - halos.s**2)
    
    try:
        T_min = h['T_min']
        T_max = h['T_max']
    except:
        T_min = 0.
        T_max = 1.
    
    mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
    mfit_NFW = (halos.cNFW_rho > 1.)*(halos.cNFW_S > 1.)*(halos.cNFW_rho < 10.)*(halos.cNFW_S < 10.)*(halos.lgMNFW_rho > 12)*(halos.lgMNFW_S > 12)
    mfit_Ein = (halos.cEin_rho > 1.)*(halos.cEin_S > 1.)*(halos.cEin_rho < 10.)*(halos.cEin_S < 10.)*(halos.lgMEin_rho > 12)*(halos.lgMEin_S > 12)*(halos.alpha_rho > 0.)*(halos.alpha_S > 0.)*(halos.alpha_rho < 0.7)*(halos.alpha_S < 0.7)
    mshape   = (T >= T_min)*(T < T_max)*(halos.q2d >= h['q_min'])*(halos.q2d < h['q_max'])
    mhalos = mhalos*mfit_Ein*mfit_NFW*mshape
    
    if relax:
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)
        mhalos = mhalos*mrelax
    
    halos = halos[mhalos]

    # f,axM = plt.subplots()    
    axM.hist(halos.lgMNFW_rho,np.linspace(13.1,14.4,50),histtype='step',label='NFW',color='orangered')
    axM.hist(halos.lgMEin_rho,np.linspace(13.1,14.4,50),histtype='step',label='Einasto',color='seagreen')
    axM.axvline(np.mean(np.log10(10**halos.lgMEin_rho)),ls='-',color='seagreen')
    axM.axvline(np.mean(np.log10(10**halos.lgMNFW_rho)),ls='-',color='orangered')    
    axM.axvline(np.log10(fitEinDS.M200),ls='--',color='seagreen')
    axM.axvline(np.log10(fitNFWDS.M200),ls='--',color='orangered')
    
    # f,axc = plt.subplots()    
    axc.hist(halos.cNFW_rho,np.linspace(1.,8,50),histtype='step',label='NFW',color='orangered')
    axc.hist(halos.cEin_rho,np.linspace(1.,8,50),histtype='step',label='Einasto',color='seagreen')
    axc.axvline(np.mean(halos.cEin_rho),ls='-',color='seagreen')
    axc.axvline(np.mean(halos.cNFW_rho),ls='-',color='orangered')    
    axc.axvline(fitEinDS.c200,ls='--',color='seagreen')
    axc.axvline(fitNFWDS.c200,ls='--',color='orangered')

    # f,axa = plt.subplots()    
    axa.hist(halos.alpha_rho,np.linspace(0.05,0.6,50),histtype='step',label='Einasto',color='seagreen')
    axa.axvline(np.mean(halos.alpha_rho),ls='-',color='seagreen')
    axa.axvline(fitEinDS.alpha,ls='--',color='seagreen')

    
    '''
    ####### 2D
    print('#####2D#####')
    
    axM2d.hist(halos.lgMNFW_S,np.linspace(13.1,14.4,50),histtype='step',label='NFW',color='orangered')
    axM2d.hist(halos.lgMEin_S,np.linspace(13.1,14.4,50),histtype='step',label='Einasto',color='seagreen')
    axM2d.axvline(np.mean(np.log10(10**(halos.lgMEin_S))),ls='-',color='seagreen')
    axM2d.axvline(np.mean(np.log10(10**(halos.lgMNFW_S))),ls='-',color='orangered')    
    axM2d.axvline(np.log10(fitEinDS.M200),ls='--',color='seagreen')
    axM2d.axvline(np.log10(fitNFWDS.M200),ls='--',color='orangered')
    
    print('M ratio')
    print(10**(np.mean(np.log10(10**halos.lgMEin_S)) - np.log10(fitEinDS.M200)))
    print(10**(np.mean(np.log10(10**halos.lgMNFW_S)) - np.log10(fitNFWDS.M200)))

    f,axc2d = plt.subplots()    
    axc2d.hist(halos.cNFW_S,np.linspace(1.,8,50),histtype='step',label='NFW',color='orangered')
    axc2d.hist(halos.cEin_S,np.linspace(1.,8,50),histtype='step',label='Einasto',color='seagreen')
    axc2d.axvline(np.mean(halos.cEin_S),ls='-',color='seagreen')
    axc2d.axvline(np.mean(halos.cNFW_S),ls='-',color='orangered')    
    axc2d.axvline(fitEinDS.c200,ls='--',color='seagreen')
    axc2d.axvline(fitNFWDS.c200,ls='--',color='orangered')

    print('c ratio')
    print(np.mean(halos.cEin_S)/fitEinDS.c200)
    print(np.mean(halos.cNFW_S)/fitNFWDS.c200)

    # f,axa = plt.subplots()    
    axa2d.hist(halos.alpha_S,np.linspace(0.05,0.6,50),histtype='step',label='Einasto',color='seagreen')
    axa2d.axvline(np.mean(halos.alpha_S),ls='-',color='seagreen')
    axa2d.axvline(fitEinDS.alpha,ls='--',color='seagreen')

    print('alpha ratio')
    print(np.mean(halos.alpha_S)/fitEinDS.alpha)
    # '''    
    

    f=open(plot_path+'compare_lens.tab','a')
    f.write(samp+' & ')
    f.write('$'+str('%.2f' % (10**(np.log10(fitNFWDS.M200) - np.mean(np.log10(10**halos.lgMNFW_rho)))))+'$ & ')
    f.write('$'+str('%.2f' % (fitNFWDS.c200/np.mean(halos.cNFW_rho)))+'$ & ')
    f.write('$'+str('%.3f' % (fitNFWDS.res))+'$ & ')
    f.write('$'+str('%.2f' % (10**(np.log10(fitEinDS.M200) - np.mean(np.log10(10**halos.lgMEin_rho)))))+'$ & ')
    f.write('$'+str('%.2f' % (fitEinDS.c200/np.mean(halos.cEin_rho)))+'$ & ')
    f.write('$'+str('%.2f' % (fitEinDS.alpha/np.mean(halos.alpha_rho)))+'$ & ')
    f.write('$'+str('%.3f' % (fitEinDS.res))+'$ & ')
    f.write('$'+str('%.2f' % (fitNFWDS.c200/np.mean(halos.cEin_rho)))+r'$ \\ '+'\n') 
    f.close()

    # f=open(plot_path+'compare_lens_2d.tab','a')
    # f.write(samp+' & ')
    # f.write('$'+str('%.2f' % (10**(np.log10(fitNFWDS.M200) - np.mean(np.log10(10**halos.lgMNFW_S)))))+'$ & ')
    # f.write('$'+str('%.2f' % (fitNFWDS.c200/np.mean(halos.cNFW_S)))+'$ & ')
    # f.write('$'+str('%.3f' % (fitNFWDS.res))+'$ & ')
    # f.write('$'+str('%.2f' % (10**(np.log10(fitEinDS.M200) - np.mean(np.log10(10**halos.lgMEin_S)))))+'$ & ')
    # f.write('$'+str('%.2f' % (fitEinDS.c200/np.mean(halos.cEin_S)))+'$ & ')
    # f.write('$'+str('%.2f' % (fitEinDS.alpha/np.mean(halos.alpha_S)))+'$ & ')
    # f.write('$'+str('%.3f' % (fitEinDS.res))+r'$ \\ '+'\n') 
    # f.close()

def plot_old():

    fp1, axp1 = plt.subplots(4, 1, gridspec_kw={'height_ratios': [3, 1,1,1]},figsize=(5,8),sharex = True)
    fp2, axp2 = plt.subplots(4, 1, gridspec_kw={'height_ratios': [3, 1,1,1]},figsize=(5,8),sharex = True)
    
        
    axp = [axp1[0],axp1[1],axp1[2],axp1[3],axp2[0],axp2[1],axp2[2],axp2[3]]
    
    fp1.subplots_adjust(hspace=0)
    fp2.subplots_adjust(hspace=0)

    for j in range(8):
        axp[j].axvline(RIN/1000.,color='k',ls=':')
        axp[j].axvline(ROUT/1000.,color='k',ls=':')
        axp[j].plot([0,10],[1,1],'k')


    def make_plot(samp,cline,samp_name):

        pro = fits.open(folder+'profile_'+samp+'.fits')
        pro_old = fits.open(folder+'profile_'+samp+'_old.fits')
        x,ynfw,yein = np.loadtxt(plot_path+'profile_'+samp+'.fit')
        
        p   = pro[1].data
        cov = pro[2].data
        pold   = pro_old[1].data
        covold = pro_old[2].data
        
        CovDS  = cov.COV_ST.reshape(len(p),len(p))
        CovDSold  = covold.COV_ST.reshape(len(p),len(p))

        axp[0].fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color=cline,alpha=0.4)    
        axp[0].fill_between(pold.Rp,pold.DSigma_T+np.diag(CovDSold),pold.DSigma_T-np.diag(CovDSold),color=cline,alpha=0.4)    
        axp[0].plot(p.Rp,p.DSigma_T,cline,label=samp_name)
        axp[0].plot(pold.Rp,pold.DSigma_T,cline+'--')
        
        axp[1].plot(p.Rp,p.DSigma_T/pold.DSigma_T,cline)
        axp[2].plot(p.Rp,p.DSigma_T/ynfw,cline)
        axp[3].plot(p.Rp,p.DSigma_T/yein,cline)
        
        pro = fits.open(folder+'profile_'+samp+'_relaxed.fits')
        pro_old = fits.open(folder+'profile_'+samp+'_old_relaxed.fits')
        x,ynfw,yein = np.loadtxt(plot_path+'profile_'+samp+'_relaxed.fit')
        
        p   = pro[1].data
        cov = pro[2].data
        pold   = pro_old[1].data
        covold = pro_old[2].data
        
        CovDS  = cov.COV_ST.reshape(len(p),len(p))
        CovDSold  = covold.COV_ST.reshape(len(p),len(p))

        axp[4].fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color=cline,alpha=0.4)    
        axp[4].fill_between(pold.Rp,pold.DSigma_T+np.diag(CovDSold),pold.DSigma_T-np.diag(CovDSold),color=cline,alpha=0.4)    
        axp[4].plot(p.Rp,p.DSigma_T,cline)
        axp[4].plot(pold.Rp,pold.DSigma_T,cline+'--')
        
        axp[5].plot(p.Rp,p.DSigma_T/pold.DSigma_T,cline)
        axp[6].plot(p.Rp,p.DSigma_T/ynfw,cline)
        axp[7].plot(p.Rp,p.DSigma_T/yein,cline)


    hsamples = ['HM_Lz','LM_Lz','HM_Hz','LM_Hz']
    clines   = ['C0','C1','C2','C3']
    labels = ['HM-Lz','LM-Lz','HM-Hz','LM-Hz']

    for j in range(len(hsamples)):
        make_plot(hsamples[j],clines[j],labels[j])

    
    axp[0].plot([-10,-12],[-10,-12],'k',label='SS')
    axp[0].plot([-10,-12],[-10,-12],'k--',label='CM')
    axp[0].legend(frameon=False,loc=3,ncol=3)
    
    axp[0].set_xscale('log')
    axp[5].set_xscale('log')
    
    axp[0].set_yscale('log')
    axp[4].set_yscale('log')
    

    axp[0].set_ylabel(r'$\widetilde{\Delta\Sigma} [M_\odot h/pc^2]$',labelpad=2)
    axp[4].set_ylabel(r'$\widetilde{\Delta\Sigma} [M_\odot h/pc^2]$',labelpad=2)
    
    axp[1].set_ylabel(r'$\widetilde{\Delta\Sigma^{SS}} / \widetilde{\Delta\Sigma^{CM}}$',labelpad=2)
    axp[5].set_ylabel(r'$\widetilde{\Delta\Sigma^{SS}} / \widetilde{\Delta\Sigma^{CM}}$',labelpad=2)
    
    axp[2].set_ylabel(r'$\widetilde{\Delta\Sigma^{SS}} / \Delta\Sigma^{NFW}$',labelpad=2)
    axp[3].set_ylabel(r'$\widetilde{\Delta\Sigma^{SS}} / \Delta\Sigma^{Ein}$',labelpad=2)
    axp[6].set_ylabel(r'$\widetilde{\Delta\Sigma^{SS}} / \Delta\Sigma^{NFW}$',labelpad=2)
    axp[7].set_ylabel(r'$\widetilde{\Delta\Sigma^{SS}} / \Delta\Sigma^{Ein}$',labelpad=2)

    axp[3].set_xlabel('r [$h^{-1}$ Mpc]')
    axp[7].set_xlabel('r [$h^{-1}$ Mpc]')
    
    axp[0].set_ylim(2,200)
    axp[4].set_ylim(2,200)
    
    axp[1].set_ylim(0.6,2.5)
    axp[5].set_ylim(0.6,2.5)
    
    axp[2].set_ylim(0.1,2)
    axp[3].set_ylim(0.1,2)
    axp[6].set_ylim(0.1,2)
    axp[7].set_ylim(0.1,2)

    axp[0].set_xlim(0.1,5)
    axp[0].xaxis.set_ticks([0.1,1,5])
    axp[0].set_xticklabels([0.1,1,5])

    axp[4].set_xlim(0.1,5)
    axp[4].xaxis.set_ticks([0.1,1,5])
    axp[4].set_xticklabels([0.1,1,5])
    
    axp[0].text(1,100,'all halos')
    axp[4].text(1,100,'only relaxed')
    

    fp1.savefig(plot_path+'profile_lens_comparison_all.pdf',bbox_inches='tight')
    fp2.savefig(plot_path+'profile_lens_comparison_relax.pdf',bbox_inches='tight')
    
def plot_old_S():

    fp1, axp1 = plt.subplots(4, 1, gridspec_kw={'height_ratios': [3, 1,1,1]},figsize=(5,8),sharex = True)
    fp2, axp2 = plt.subplots(4, 1, gridspec_kw={'height_ratios': [3, 1,1,1]},figsize=(5,8),sharex = True)
    
        
    axp = [axp1[0],axp1[1],axp1[2],axp1[3],axp2[0],axp2[1],axp2[2],axp2[3]]
    
    fp1.subplots_adjust(hspace=0)
    fp2.subplots_adjust(hspace=0)

    for j in range(8):
        axp[j].axvline(RIN/1000.,color='k',ls=':')
        axp[j].axvline(ROUT/1000.,color='k',ls=':')
        axp[j].plot([0,10],[1,1],'k')


    def make_plot(samp,cline,samp_name):

        pro = fits.open(folder+'profile_'+samp+'.fits')
        pro_old = fits.open(folder+'profile_'+samp+'_old.fits')
        x,ynfw,yein = np.loadtxt(plot_path+'profile_'+samp+'_S.fit')
        
        p   = pro[1].data
        cov = pro[2].data
        pold   = pro_old[1].data
        covold = pro_old[2].data
        
        CovS  = cov.COV_S.reshape(len(p),len(p))
        CovSold  = covold.COV_S.reshape(len(p),len(p))

        axp[0].fill_between(p.Rp,p.Sigma+np.diag(CovS),p.Sigma-np.diag(CovS),color=cline,alpha=0.4)    
        axp[0].fill_between(pold.Rp,pold.Sigma+np.diag(CovSold),pold.Sigma-np.diag(CovSold),color=cline,alpha=0.4)    
        axp[0].plot(p.Rp,p.Sigma,cline,label=samp_name)
        axp[0].plot(pold.Rp,pold.Sigma,cline+'--')
        
        axp[1].plot(p.Rp,p.Sigma/pold.Sigma,cline)
        axp[2].plot(p.Rp,p.Sigma/ynfw,cline)
        axp[3].plot(p.Rp,p.Sigma/yein,cline)
        
        pro = fits.open(folder+'profile_'+samp+'_relaxed.fits')
        pro_old = fits.open(folder+'profile_'+samp+'_old_relaxed.fits')
        x,ynfw,yein = np.loadtxt(plot_path+'profile_'+samp+'_relaxed_S.fit')
        
        p   = pro[1].data
        cov = pro[2].data
        pold   = pro_old[1].data
        covold = pro_old[2].data
        
        CovDS  = cov.COV_S.reshape(len(p),len(p))
        CovDSold  = covold.COV_S.reshape(len(p),len(p))

        axp[4].fill_between(p.Rp,p.Sigma+np.diag(CovS),p.Sigma-np.diag(CovS),color=cline,alpha=0.4)    
        axp[4].fill_between(pold.Rp,pold.Sigma+np.diag(CovSold),pold.Sigma-np.diag(CovSold),color=cline,alpha=0.4)    
        axp[4].plot(p.Rp,p.Sigma,cline)
        axp[4].plot(pold.Rp,pold.Sigma,cline+'--')
        
        axp[5].plot(p.Rp,p.Sigma/pold.Sigma,cline)
        axp[6].plot(p.Rp,p.Sigma/ynfw,cline)
        axp[7].plot(p.Rp,p.Sigma/yein,cline)


    hsamples = ['HM_Lz','LM_Lz','HM_Hz','LM_Hz']
    clines   = ['C0','C1','C2','C3']
    labels = ['HM-Lz','LM-Lz','HM-Hz','LM-Hz']

    for j in range(len(hsamples)):
        make_plot(hsamples[j],clines[j],labels[j])

    
    axp[0].plot([-10,-12],[-10,-12],'k',label='SS')
    axp[0].plot([-10,-12],[-10,-12],'k--',label='CM')
    axp[0].legend(frameon=False,loc=3,ncol=3)
    
    axp[0].set_xscale('log')
    axp[5].set_xscale('log')
    
    axp[0].set_yscale('log')
    axp[4].set_yscale('log')
    

    axp[0].set_ylabel(r'$\widetilde{\Sigma} [M_\odot h/pc^2]$',labelpad=2)
    axp[4].set_ylabel(r'$\widetilde{\Sigma} [M_\odot h/pc^2]$',labelpad=2)
    
    axp[1].set_ylabel(r'$\widetilde{\Sigma^{SS}} / \widetilde{\Sigma^{CM}}$',labelpad=2)
    axp[5].set_ylabel(r'$\widetilde{\Sigma^{SS}} / \widetilde{\Sigma^{CM}}$',labelpad=2)
    
    axp[2].set_ylabel(r'$\widetilde{\Sigma^{SS}} / \Sigma^{NFW}$',labelpad=2)
    axp[3].set_ylabel(r'$\widetilde{\Sigma^{SS}} / \Sigma^{Ein}$',labelpad=2)
    axp[6].set_ylabel(r'$\widetilde{\Sigma^{SS}} / \Sigma^{NFW}$',labelpad=2)
    axp[7].set_ylabel(r'$\widetilde{\Sigma^{SS}} / \Sigma^{Ein}$',labelpad=2)

    axp[3].set_xlabel('r [$h^{-1}$ Mpc]')
    axp[7].set_xlabel('r [$h^{-1}$ Mpc]')
    
    axp[0].set_ylim(2,200)
    axp[4].set_ylim(2,200)
    
    axp[1].set_ylim(0.6,2.5)
    axp[5].set_ylim(0.6,2.5)
    
    axp[2].set_ylim(0.1,2)
    axp[3].set_ylim(0.1,2)
    axp[6].set_ylim(0.1,2)
    axp[7].set_ylim(0.1,2)

    axp[0].set_xlim(0.1,5)
    axp[0].xaxis.set_ticks([0.1,1,5])
    axp[0].set_xticklabels([0.1,1,5])

    axp[4].set_xlim(0.1,5)
    axp[4].xaxis.set_ticks([0.1,1,5])
    axp[4].set_xticklabels([0.1,1,5])
    
    axp[0].text(1,100,'all halos')
    axp[4].text(1,100,'only relaxed')
    

    fp1.savefig(plot_path+'profile_lens_comparison_S_all.pdf',bbox_inches='tight')
    fp2.savefig(plot_path+'profile_lens_comparison_S_relax.pdf',bbox_inches='tight')
    
def test_fitting(hsamples,
                   RIN,
                   ROUToq = ['2000','1000','2000','1000'],
                   relax=True):
    
        
    ROUT = ROUToq
        
    
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
        
        maskr   = (p.Rp > (float(RIN[j])/1000.))*(p.Rp < (float(ROUT[j])/1000.))
        mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]
        CovDSf = (cov.COV_ST.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())

        fitNFWDS = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSf),h['z_mean'])

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
        fstd = fits.open(folder+'fitresults_onlyq_'+RIN[j]+'_'+ROUT[j]+'_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMNFW),label='1halo',alpha=0.5)
        ax[1].plot(fstd.c200[1500:]/cNFW,label='1halo',alpha=0.5)
                       
        # Einasto
        fstd = fits.open(folder+'fitresults_2h_2q_Ein_'+RIN[j]+'_5000_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMEin),label='Einasto',alpha=0.5)
        ax[1].plot(fstd.c200[1500:]/cEin,label='Einasto',alpha=0.5)
        
        # NFW
        fstd = fits.open(folder+'fitresults_2h_2q_'+RIN[j]+'_5000_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMNFW),label='NFW',alpha=0.5)
        ax[1].plot(fstd.c200[1500:]/cNFW,label='NFW',alpha=0.5)

        # without_c
        fstd = fits.open(folder+'fitresults_2h_2q_woc_'+RIN[j]+'_5000_profile_'+samp+'.fits')[1].data
        ax[0].plot(10**(fstd.lM200[1500:]-lMNFW),label='fix c',alpha=0.5)
        
        ax[0].axhline(fitNFWDS.M200/10**lMNFW,color='C3',lw=2)
        ax[1].axhline(fitNFWDS.c200/cNFW,color='C3',lw=2)
        ax[1].axhline(cEin/cNFW,color='C1',lw=2)
        
        ax[0].legend(frameon=False,ncol=4,loc=3)
        
        ax[0].set_ylim([0.6,1.15])
        ax[1].set_ylim([0.6,1.15])
        
        ax[1].set_xlabel(r'Nit $\times$ Nchains')
        ax[1].set_ylabel(r'$c_{200}/\langle c_{200} \rangle$')
        ax[0].set_ylabel(r'$M_{200}/\langle M_{200} \rangle$')
        
        # f.savefig(folder+'../final_plots/test_fit_RIN'+RIN+'_'+samp+'.png',bbox_inches='tight')
        # f.savefig(folder+'../test_plots/test_fit_RIN'+RIN[j]+'_'+samp+'.png',bbox_inches='tight')
        f.savefig(folder+'../test_plots/test_fit_RIN'+RIN[j]+'_'+samp+'.png',bbox_inches='tight')

# '''    

f, ax = plt.subplots(4,3, figsize=(14,7))
f.subplots_adjust(hspace=0)

f2d, ax2d = plt.subplots(4,3, figsize=(14,7))
f2d.subplots_adjust(hspace=0)

fp, axp = plt.subplots(4,1, figsize=(4,9),sharey=True,sharex=True)
fp.subplots_adjust(hspace=0,wspace=0)

# hsamples_relaxed = ['HM_Lz_relaxed','LM_Lz_relaxed','HM_Hz_relaxed','LM_Hz_relaxed']
# hsamples_relaxed = ['HM_Lz_elong_relaxed','LM_Lz_elong_relaxed','HM_Hz_elong_relaxed','LM_Hz_elong_relaxed']
# hsamples_relaxed = ['HM_Lz_round_relaxed','LM_Lz_round_relaxed','HM_Hz_round_relaxed','LM_Hz_round_relaxed']
# hsamples_relaxed = ['HM_Lz_round_obl_relaxed','LM_Lz_round_obl_relaxed','HM_Hz_round_obl_relaxed','LM_Hz_round_obl_relaxed']
hsamples_relaxed = ['HM_Lz_round_pro_relaxed','LM_Lz_round_pro_relaxed','HM_Hz_round_pro_relaxed','LM_Hz_round_pro_relaxed']

ct = [25,100,35,130]
labels = ['HM-Lz','LM-Lz','HM-Hz','LM-Hz']

for j in range(len(hsamples_relaxed)):
    axes = [axp[j],ax[j,0],ax[j,1],ax[j,2],ax2d[j,0],ax2d[j,1],ax2d[j,2]]
    fit_profiles(hsamples_relaxed[j],axes,relax=True)
    ax[j,1].text(6,ct[j],labels[j])
    axp[j].text(0.8,100.,labels[j])
    for i in range(3):
        ax[j,i].set_ylabel('$N$')    

ax[j,0].set_xlabel(r'$M_{200}$')    
ax[j,1].set_xlabel(r'$c_{200}$')    
ax[j,2].set_xlabel(r'$\alpha$')    


# f.savefig(plot_path+'dist_lens_compare_relaxed_rescut.pdf',bbox_inches='tight')
# f2d.savefig(plot_path+'dist_lens_compare2d_relaxed.pdf',bbox_inches='tight')
# fp.savefig(plot_path+'profile_lens_relaxed.pdf',bbox_inches='tight')


'''
# all


f, ax = plt.subplots(4,3, figsize=(14,7))
f.subplots_adjust(hspace=0)

f2d, ax2d = plt.subplots(4,3, figsize=(14,7))
f2d.subplots_adjust(hspace=0)

fp, axp = plt.subplots(4,1, figsize=(4,9),sharey=True,sharex=True)
fp.subplots_adjust(hspace=0,wspace=0)

hsamples_relaxed = ['HM_Lz','LM_Lz','HM_Hz','LM_Hz']
ct = [30,150,40,200]
labels = ['HM-Lz','LM-Lz','HM-Hz','LM-Hz']

for j in range(len(hsamples_relaxed)):
    axes = [axp[j],ax[j,0],ax[j,1],ax[j,2],ax2d[j,0],ax2d[j,1],ax2d[j,2]]
    fit_profiles(hsamples_relaxed[j],axes,relax=False)
    ax[j,1].text(6,ct[j],labels[j])
    axp[j].text(0.8,100.,labels[j])
    for i in range(3):
        ax[j,i].set_ylabel('$N$')    

ax[j,0].set_xlabel(r'$M_{200}$')    
ax[j,1].set_xlabel(r'$c_{200}$')    
ax[j,2].set_xlabel(r'$\alpha$')    
    

f.savefig(plot_path+'dist_lens_compare.pdf',bbox_inches='tight')
f2d.savefig(plot_path+'dist_lens_compare2d.pdf',bbox_inches='tight')
fp.savefig(plot_path+'profile_lens.pdf',bbox_inches='tight')
# '''
