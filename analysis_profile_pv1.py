import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from astropy.constants import G,c,M_sun, pc
from fit_models_colossus import *

cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)
import corner
folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'

def fit_profiles(samp,RIN,ROUT,fittype='_onlyq',
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
            
    
    ndots = p.shape[0]

    maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))

    mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]

    CovDS  = cov.COV_ST.reshape(len(p),len(p))
    CovS   = cov.COV_S.reshape(len(p),len(p))

    CovDSf = (cov.COV_ST.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())
    CovSf = (cov.COV_S.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())

    fitEinDS = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSf),zmean,'Einasto',1.e14,3.)
    fitEinS  = Sigma_fit(p.Rp[maskr],p.Sigma[maskr],np.diag(CovSf),zmean,'Einasto',1.e14,3.)
    fitNFWDS = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSf),zmean)
    fitNFWS  = Sigma_fit(p.Rp[maskr],p.Sigma[maskr],np.diag(CovSf),zmean)

    f,axDS = plt.subplots()
    axDS.plot(p.Rp,p.DSigma_T,'C7')
    axDS.plot(fitEinDS.xplot,fitEinDS.yplot,'C0',label='Einasto')
    axDS.plot(fitNFWDS.xplot,fitNFWDS.yplot,'C9',label='NFW')
    axDS.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C7',alpha=0.4)
    axDS.set_xscale('log')
    axDS.set_yscale('log')
    axDS.set_ylabel(r'$\Delta\Sigma$',labelpad=2)
    axDS.set_xlabel('r [$h^{-1}$ Mpc]')
    axDS.set_ylim(7,200)
    axDS.set_xlim(0.1,5)
    axDS.xaxis.set_ticks([0.1,1,5])
    axDS.set_xticklabels([0.1,1,5])
    axDS.yaxis.set_ticks([10,100])
    axDS.set_yticklabels([10,100])
    axDS.axvline(RIN/1000.,color='k',ls=':')
    axDS.axvline(ROUT/1000.,color='k',ls=':')
    axDS.legend()
   
    f,axS = plt.subplots()
    axS.plot(p.Rp,p.Sigma,'C7')
    axS.plot(fitEinS.xplot,fitEinS.yplot,'C0',label='Einasto')
    axS.plot(fitNFWS.xplot,fitNFWS.yplot,'C9',label='NFW')
    axS.fill_between(p.Rp,p.Sigma+np.diag(CovS),p.Sigma-np.diag(CovS),color='C7',alpha=0.4)
    axS.set_xscale('log')
    axS.set_yscale('log')
    axS.set_ylabel(r'$\Sigma$',labelpad=2)
    axS.set_xlabel('r [$h^{-1}$ Mpc]')
    axS.set_ylim(7,500)
    axS.set_xlim(0.1,5)
    axS.xaxis.set_ticks([0.1,1,5])
    axS.set_xticklabels([0.1,1,5])
    axS.yaxis.set_ticks([10,100])
    axS.set_yticklabels([10,100])
    axS.axvline(RIN/1000.,color='k',ls=':')
    axS.axvline(ROUT/1000.,color='k',ls=':')
    
    out = np.array([np.log10(fitEinS.M200),fitEinS.c200,fitEinS.alpha,fitEinS.res,
                    np.log10(fitEinDS.M200),fitEinDS.c200,fitEinDS.alpha,fitEinDS.res,
                    np.log10(fitNFWS.M200),fitNFWS.c200,fitNFWS.alpha,fitNFWS.res,
                    np.log10(fitNFWDS.M200),fitNFWDS.c200,fitNFWDS.alpha,fitNFWDS.res])
    
    np.savetxt(folder+'../fitlens.res',out,fmt='%10.2f')
    
    
    halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data
                
    
    Eratio = (2.*halos.K/abs(halos.U))
    
    mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
    
    if relax:
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)
        mhalos = mhalos*mrelax
    
    halos = halos[mhalos]


    
    f, ax = plt.subplots(2,3, figsize=(10,6), sharey=True)
    
    ax[0,0].set_title('3D')
    ax[0,0].hist(halos.lgMNFW_rho,np.linspace(13.,14.5,50),histtype='step',label='NFW')
    ax[0,0].hist(halos.lgMEin_rho,np.linspace(13.,14.5,50),histtype='step',label='Einasto')
    ax[0,0].axvline(fitEinS.M200,ls='--')
    ax[0,0].axvline(fitNFWS.M200,ls='-')
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

