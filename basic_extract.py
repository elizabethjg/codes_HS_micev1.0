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
        lMNFW = np.log10(np.mean(10**(halos[mhalos].lgMNFW_rho)))
        cNFW  = np.mean(halos[mhalos].cNFW_rho)
        lMEin = np.log10(np.mean(10**(halos[mhalos].lgMEin_rho)))
        cEin  = np.mean(halos[mhalos].cEin_rho)
        alpha = np.mean(halos[mhalos].alpha_rho)

        NFW_h += [[lMNFW,cNFW]]
        Ein_h += [[lMEin,cEin,alpha]]


        lMNFWr = np.log10(np.mean(10**(halos[mhalos*mrelax].lgMNFW_rho)))
        cNFWr  = np.mean(halos[mhalos*mrelax].cNFW_rho)
        lMEinr = np.log10(np.mean(10**(halos[mhalos*mrelax].lgMNFW_rho)))
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

