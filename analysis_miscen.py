import sys
import pylab
from astropy.io import fits
from astropy.cosmology import LambdaCDM
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
sys.path.append('/mnt/projects/lensing/lens_codes_v3.7')
from models_profiles import *
from fit_profiles_curvefit import *
from fit_models_colossus import Delta_Sigma_fit

folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'

def compute_profiles(samp,s_off): 
    
    folder = '/home/elizabeth/MICE/HS-lensing/profiles/'
    p_name = 'profile_'+samp+'.fits'  

    profile = fits.open(folder+p_name)
    
    print(p_name)
      
    p   = profile[1].data
    h   = profile[0].header
    cov = profile[2].data
    
    zmean = h['z_mean']

    RIN = 250
    ROUT = 2500

    # fitpar = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_profile_'+samp+'.fits')[0].header
    
    maskr   = (p.Rp > (RIN/1000.))*(p.Rp < (ROUT/1000.))
    mr = np.meshgrid(maskr,maskr)[1]*np.meshgrid(maskr,maskr)[0]
    
    CovDSf = (cov.COV_ST.reshape(len(p.Rp),len(p.Rp))[mr]).reshape(maskr.sum(),maskr.sum())
    
    fitpar = Delta_Sigma_fit(p.Rp[maskr],p.DSigma_T[maskr],np.diag(CovDSf),zmean)
    
    S0  = Sigma_NFW_2h(p.Rp,0.5,1.e14,3.)
    DS0 = Delta_Sigma_NFW_2h(p.Rp,0.5,1.e14,3.)
    
    Soff  = Sigma_NFW_miss(p.Rp,0.5,1.e14,c200=3.,s_off=s_off)
    DSoff  = Delta_Sigma_NFW_miss(p.Rp,0.5,1.e14,c200=3.,s_off=s_off)
    
    np.savetxt(folder+'missprof_pru.cat',np.array([S0,DS0,Soff,DSoff]),fmt='%10.2f')



def compute_profiles(samp): 
    

    p_name = 'profile_HM_Lz_soff1_relaxed_miscen.fits'  
    # p_name = 'profile_pru_relaxed_miscen.fits'  

    pmiss = fits.open(folder+p_name)[1].data
    p     = fits.open(folder+'profile_HM_Lz_relaxed.fits')[1].data      
    
    cov_miss = fits.open(folder+p_name)[2].data
    cov      = fits.open(folder+'profile_HM_Lz_relaxed.fits')[2].data      
    
    zmean = h['z_mean']

    RIN = 250
    ROUT = 5000

    fitpar = fits.open(folder+'fitresults_2h_2q_'+str(int(RIN))+'_'+str(int(ROUT))+'_profile_HM_Lz_relaxed.fits')[0].header
        
    S0,DS0,Soff,DSoff = np.loadtxt(folder+'missprof_HM_Lz_relaxed_1000.cat')
    # S0,DS0,Soff,DSoff = np.loadtxt(folder+'missprof_pru.cat')
    
    CovDS  = cov.COV_ST.reshape(len(p),len(p))
    CovS  = cov.COV_S.reshape(len(p),len(p))
    CovDSmiss  = cov_miss.COV_ST.reshape(len(p),len(p))
    CovSmiss  = cov_miss.COV_S.reshape(len(p),len(p))

    plt.figure()
    plt.fill_between(p.Rp,p.Sigma+np.diag(CovS),p.Sigma-np.diag(CovS),color='C3',alpha=0.7)
    plt.plot(p.Rp,p.Sigma,'C3')
    plt.plot(p.Rp,S0)
    plt.ylabel('$\Sigma$',fontsize=14)
    plt.xlabel('$R$',fontsize=14)
    plt.loglog()
    plt.savefig(folder+'plots/miss_S_1000.png')
    
    plt.figure()
    plt.fill_between(p.Rp,pmiss.Sigma+np.diag(CovSmiss),pmiss.Sigma-np.diag(CovSmiss),color='C3',alpha=0.7)
    plt.plot(p.Rp,pmiss.Sigma,'C3')
    plt.plot(p.Rp,Soff)    
    plt.ylabel('$\Sigma$',fontsize=14)
    plt.xlabel('$R$',fontsize=14)
    plt.loglog()
    plt.savefig(folder+'plots/miss_Soff_1000.png')
    
    plt.figure()
    plt.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C3',alpha=0.7)
    plt.plot(p.Rp,p.DSigma_T,'C3')
    plt.plot(p.Rp,DS0)
    plt.ylabel('$\Delta \Sigma$',fontsize=14)
    plt.xlabel('$R$',fontsize=14)
    plt.loglog()
    plt.savefig(folder+'plots/miss_DS_1000.png')

    
    plt.figure()
    plt.fill_between(p.Rp,pmiss.DSigma_T+np.diag(CovDSmiss),pmiss.DSigma_T-np.diag(CovDSmiss),color='C3',alpha=0.7)
    plt.plot(p.Rp,pmiss.DSigma_T,'C3')
    plt.plot(p.Rp,DSoff)        
    plt.ylabel('$\Delta \Sigma$',fontsize=14)
    plt.xlabel('$R$',fontsize=14)
    plt.loglog()
    plt.savefig(folder+'plots/miss_DSoff_1000.png')


def Delta_Sigma_NFW_miss(R,z,M200,s_off = None, tau = 0.2,
                         c200 = None, P_Roff = Rayleigh, cosmo_params=params):	
    
    '''
    Misscentred density contraste for NFW
    
    '''
  

    integral = []
    for r in R:
        argumento = lambda x: Sigma_NFW_2h(x,z,M200,c200)*x
        integral  += [integrate.quad(argumento, 0, r, epsabs=1.e-02, epsrel=1.e-02)[0]]

    DS    = (2./R**2)*integral - Sigma_NFW_2h(R,z,M200,c200)

    Delta_Sigma_NFW_2h(p.Rp,z,M200,c200)
