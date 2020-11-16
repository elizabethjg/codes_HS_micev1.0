import sys
sys.path.append('/mnt/clemente/lensing')
sys.path.append('/mnt/clemente/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM
from maria_func import *
from profiles_fit import *
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from multiprocessing import Pool
from multiprocessing import Process
import argparse
from astropy.constants import G,c,M_sun,pc

#parameters
cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)

folder = '/mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/catalogs/'
S      = fits.open(folder+'MICE_sources.fits')[1].data


def partial_profile(RA0,DEC0,Z,angles,
                    RIN,ROUT,ndots,h,nboot=100):

        ndots = int(ndots)

        cosmo = LambdaCDM(H0=100*h, Om0=0.3, Ode0=0.7)
        dl  = cosmo.angular_diameter_distance(Z).value
        KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
        
        delta = ROUT/(3600*KPCSCALE)
        
        t0 = time.time()
        mask = (S.ra < (RA0+delta))&(S.ra > (RA0-delta))&(S.dec > (DEC0-delta))&(S.dec < (DEC0+delta))&(S.z_v > (Z+0.1))
                       
        catdata = S[mask]

        ds  = cosmo.angular_diameter_distance(catdata.z_v).value
        dls = cosmo.angular_diameter_distance_z1z2(Z, catdata.z_v).value
                
        
        BETA_array = dls/ds
        
        Dl = dl*1.e6*pc
        sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)


        rads, theta, test1,test2 = eq2p2(np.deg2rad(catdata.ra),
                                        np.deg2rad(catdata.dec),
                                        np.deg2rad(RA0),
                                        np.deg2rad(DEC0))


        theta2 = (2.*np.pi - theta) +np.pi/2.
        theta_ra = theta2
        theta_ra[theta2 > 2.*np.pi] = theta2[theta2 > 2.*np.pi] - 2.*np.pi
                       
        e1     = catdata.gamma1
        e2     = -1.*catdata.gamma2
        
        #Correct polar angle for e1, e2
        # theta = theta+np.pi/2.
        
        
        #get tangential ellipticities 
        et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
        #get cross ellipticities
        ex = (e1*np.sin(2*theta)-e2*np.cos(2*theta))*sigma_c

        
          
        del(e1)
        del(e2)
        
        r=np.rad2deg(rads)*3600*KPCSCALE
        del(rads)
        
        
        Ntot = len(catdata)
        del(catdata)    
        
        bines = np.logspace(np.log10(RIN),np.log10(ROUT),num=ndots+1)
        dig = np.digitize(r,bines)
        
        t = np.tile(theta_ra,(3,1)).T
        an = np.tile(angles,(len(theta_ra),1))
        at     = t - an
        
        DSIGMAwsum_T = []
        DSIGMAwsum_X = []
        N_inbin = []
        
        BOOTwsum_T   = np.zeros((nboot,ndots))
        BOOTwsum_X   = np.zeros((nboot,ndots))
        
        GAMMATcos_wsum = np.zeros((ndots,3))
        GAMMAXsin_wsum = np.zeros((ndots,3))
        
        BOOTwsum_Tcos  = np.zeros((nboot,ndots,3))
        BOOTwsum_Xsin  = np.zeros((nboot,ndots,3))
               
        for nbin in range(ndots):
                mbin = dig == nbin+1              
                
                DSIGMAwsum_T = np.append(DSIGMAwsum_T,et[mbin].sum())
                DSIGMAwsum_X = np.append(DSIGMAwsum_X,ex[mbin].sum())

                GAMMATcos_wsum[nbin,:] = np.sum((np.tile(et[mbin],(3,1))*np.cos(2.*at[mbin]).T),axis=1)
                GAMMAXsin_wsum[nbin,:] = np.sum((np.tile(ex[mbin],(3,1))*np.sin(2.*at[mbin]).T),axis=1)
                               
                N_inbin = np.append(N_inbin,len(et[mbin]))
                
                index = np.arange(mbin.sum())
                if mbin.sum() == 0:
                        continue
                else:
                        with NumpyRNGContext(1):
                                bootresult = bootstrap(index, nboot)
                        INDEX=bootresult.astype(int)
                        BOOTwsum_T[:,nbin] = np.sum(np.array(et[mbin])[INDEX],axis=1)
                        BOOTwsum_X[:,nbin] = np.sum(np.array(ex[mbin])[INDEX],axis=1)

                        
                        BOOTwsum_Tcos[:,nbin,:] = np.sum(((np.tile(et[mbin],(3,1))*np.cos(2.*at[mbin]).T))[:,INDEX],axis=2).T
                        BOOTwsum_Xsin[:,nbin,:] = np.sum(((np.tile(ex[mbin],(3,1))*np.sin(2.*at[mbin]).T))[:,INDEX],axis=2).T

        
        output = {'DSIGMAwsum_T':DSIGMAwsum_T,'DSIGMAwsum_X':DSIGMAwsum_X,
                   'BOOTwsum_T':BOOTwsum_T, 'BOOTwsum_X':BOOTwsum_X, 
                   'GAMMATcos_wsum': GAMMATcos_wsum, 'GAMMAXsin_wsum': GAMMAXsin_wsum,
                   'BOOTwsum_Tcos':BOOTwsum_Tcos, 'BOOTwsum_Xsin':BOOTwsum_Xsin, 
                   'N_inbin':N_inbin,'Ntot':Ntot}
        
        return output

def partial_profile_unpack(minput):
	return partial_profile(*minput)
        
'''
sample='pru'
lM_min=14.
lM_max=14.2
z_min = 0.1
z_max = 0.4
RIN = 100.
ROUT =5000.
ndots= 20
ncores=30
h=1.0
'''

def main(sample='pru',lM_min=14.,lM_max=14.2,
                z_min = 0.1, z_max = 0.4,
                RIN = 400., ROUT =5000.,
                ndots= 40, ncores=10, h=1.0):

        '''
        
        INPUT
        ---------------------------------------------------------
        sample         (str) sample name
        lM_min         (float) lower limit for log(Mass) - >=
        lM_max         (float) higher limit for log(Mass) - <
        z_min          (float) lower limit for z - >=
        z_max          (float) higher limit for z - <
        RIN            (float) Inner bin radius of profile
        ROUT           (float) Outer bin radius of profile
        ndots          (int) Number of bins of the profile
        ncores         (int) to run in parallel, number of cores
        h              (float) H0 = 100.*h
        '''

        cosmo = LambdaCDM(H0=100*h, Om0=0.3, Ode0=0.7)
        tini = time.time()
        
        print('Sample ',sample)
        print('Selecting groups with:')
        print(lM_min,' <= log(M) < ',lM_max)
        print(z_min,' <= z < ',z_max)
        print('Profile has ',ndots,'bins')
        print('from ',RIN,'kpc to ',ROUT,'kpc')
        print('h ',h)
              
        # Defining radial bins
        bines = np.logspace(np.log10(RIN),np.log10(ROUT),num=ndots+1)
        R = (bines[:-1] + np.diff(bines)*0.5)*1.e-3
        
        #reading cats
        
        L = fits.open(folder+'MICEv1.0_halo_cat_withshapes.fits')[1].data
        
        mregion = (L.ra < 80.)*(L.dec > 36.5)        
        mmass   = (L.lgm >= lM_min)*(L.lgm < lM_max)
        mz      = (L.z_v >= z_min)*(L.z_v < z_max)
        mlenses = mmass*mz*mregion
        Nlenses = mlenses.sum()

        if Nlenses < ncores:
                ncores = Nlenses
        
        print('Nlenses',Nlenses)
        print('CORRIENDO EN ',ncores,' CORES')
        
        L = L[mlenses]
        
        theta  = np.array([np.zeros(sum(mlenses)),np.arctan(L.a2dy/L.a2dx),np.arctan(L.a2dry/L.a2drx)]).T

        # Define K masks
        
        ra = L.ra
        dec = L.dec
        
        ramin  = np.min(ra)
        decmin = np.min(dec)
        dra  = ((np.max(ra)+1)  - ramin)/10.
        ddec = ((np.max(dec)+1) - decmin)/10.
        
        kmask = np.zeros((101,len(ra)))
        kmask[0] = np.ones(len(ra)).astype(bool)
        c    = 1
        
        for a in range(10): 
                for d in range(10): 
                        mra  = (ra  >= ramin + a*dra)*(ra < ramin + (a+1)*dra) 
                        mdec = (dec >= decmin + d*ddec)*(dec < decmin + (d+1)*ddec) 
        
                        kmask[c] = ~(mra*mdec)
                        c += 1

        
        # SPLIT LENSING CAT
        
        lbins = int(round(Nlenses/float(ncores), 0))
        slices = ((np.arange(lbins)+1)*ncores).astype(int)
        slices = slices[(slices < Nlenses)]
        Lsplit = np.split(L,slices)
        Tsplit = np.split(theta,slices)        
        Ksplit = np.split(kmask.T,slices)
        
        # WHERE THE SUMS ARE GOING TO BE SAVED
        
        DSIGMAwsum_T = np.zeros((101,ndots)) 
        DSIGMAwsum_X = np.zeros((101,ndots))
        
        BOOTwsum_T   = np.zeros((100,ndots))
        BOOTwsum_X   = np.zeros((100,ndots))
        
        GAMMATcos_wsum = np.zeros((101,ndots,3))
        GAMMAXsin_wsum = np.zeros((101,ndots,3))

        COS2_2theta_wsum = np.zeros((101,ndots,3))
        SIN2_2theta_wsum = np.zeros((101,ndots,3))
                                   
        BOOTwsum_Tcos  = np.zeros((100,ndots,3))
        BOOTwsum_Xsin  = np.zeros((100,ndots,3))
        
        Ninbin = np.zeros((101,ndots))
        
        Ntot         = np.array([])
        tslice       = np.array([])
        
        for l in range(len(Lsplit)):
                
                print('RUN ',l+1,' OF ',len(Lsplit))
                
                t1 = time.time()
                
                num = len(Lsplit[l])
                
                rin  = RIN*np.ones(num)
                rout = ROUT*np.ones(num)
                nd   = ndots*np.ones(num)
                h_array   = h*np.ones(num)
                
                if num == 1:
                        entrada = [Lsplit[l].ra[0], Lsplit[l].dec[0],
                                   Lsplit[l].z_v[0],Tsplit[l][0],
                                   RIN,ROUT,ndots,h]
                        
                        salida = [partial_profile_unpack(entrada)]
                else:          
                        entrada = np.array([Lsplit[l].ra,Lsplit[l].dec,
                                        Lsplit[l].z_v,Tsplit[l].tolist(),
                                        rin,rout,nd,h_array]).T
                        
                        pool = Pool(processes=(num))
                        salida = np.array(pool.map(partial_profile_unpack, entrada))
                        pool.terminate()
                                
                for j in range(len(salida)):
                        
                        profilesums = salida[j]
                        km          = np.tile(Ksplit[l][j],(3,ndots,1)).T
                                                
                        DSIGMAwsum_T += np.tile(profilesums['DSIGMAwsum_T'],(101,1))*km[:,:,0]
                        DSIGMAwsum_X += np.tile(profilesums['DSIGMAwsum_X'],(101,1))*km[:,:,0]
                        
                        Ninbin += np.tile(profilesums['N_inbin'],(101,1))*km[:,:,0]
                        
                        BOOTwsum_T   += profilesums['BOOTwsum_T']
                        BOOTwsum_X   += profilesums['BOOTwsum_X']
                        
                        GAMMATcos_wsum += np.tile(profilesums['GAMMATcos_wsum'],(101,1,1))*km
                        GAMMAXsin_wsum += np.tile(profilesums['GAMMAXsin_wsum'],(101,1,1))*km

                        COS2_2theta_wsum += np.tile((np.cos(2.*Tsplit[l][j]))**2,(101,1,1))*km
                        SIN2_2theta_wsum += np.tile((np.sin(2.*Tsplit[l][j]))**2,(101,1,1))*km
                                                
                        BOOTwsum_Tcos  += profilesums['BOOTwsum_Tcos']
                        BOOTwsum_Xsin  += profilesums['BOOTwsum_Xsin']
                        
                        
                        
                        Ntot         = np.append(Ntot,profilesums['Ntot'])
                
                t2 = time.time()
                ts = (t2-t1)/60.
                tslice = np.append(tslice,ts)
                print('TIME SLICE')
                print(ts)
                print('Estimated ramaining time')
                print(np.mean(tslice)*(len(Lsplit)-(l+1)))
        
        # COMPUTING PROFILE        
        Ninbin[DSIGMAwsum_T == 0] = 1.
                
        DSigma_T  = (DSIGMAwsum_T/Ninbin)
        DSigma_X  = (DSIGMAwsum_X/Ninbin)
        
        eDSigma_T =  np.std((BOOTwsum_T/np.tile(Ninbin[0],(100,1))),axis=0)
        eDSigma_X =  np.std((BOOTwsum_X/np.tile(Ninbin[0],(100,1))),axis=0)
        
        GAMMA_Tcos = (GAMMATcos_wsum/COS2_2theta_wsum).transpose(1,2,0)
        GAMMA_Xsin = (GAMMAXsin_wsum/SIN2_2theta_wsum).transpose(1,2,0)
        
        eGAMMA_Tcos =  np.std((BOOTwsum_Tcos/np.tile(Ninbin[0],(100,3,1)).transpose(0,2,1)),axis=0)
        eGAMMA_Xsin =  np.std((BOOTwsum_Xsin/np.tile(Ninbin[0],(100,3,1)).transpose(0,2,1)),axis=0)

        
        # AVERAGE LENS PARAMETERS
        
        zmean        = np.average(L.z_v,weights=Ntot)
        lM_mean      = np.average(L.lgm,weights=Ntot)
        
        q2d_mean     = np.average(L.q2d,weights=Ntot)
        q2dr_mean    = np.average(L.q2dr,weights=Ntot)
        q3d_mean     = np.average(L.q3d,weights=Ntot)
        q3dr_mean    = np.average(L.q3dr,weights=Ntot)
        s3d_mean     = np.average(L.s3d,weights=Ntot)
        s3dr_mean    = np.average(L.s3dr,weights=Ntot)
        
        # FITING AN NFW MODEL
        
        H        = cosmo.H(zmean).value/(1.0e3*pc) #H at z_pair s-1 
        roc      = (3.0*(H**2.0))/(8.0*np.pi*G) #critical density at z_pair (kg.m-3)
        roc_mpc  = roc*((pc*1.0e6)**3.0)
        
        try:
                nfw        = NFW_stack_fit(R,DSigma_T[0],eDSigma_T,zmean,roc)
        except:
                nfw          = [0.01,0.,100.,[0.,0.],[0.,0.],-999.,0.]

        M200_NFW   = (800.0*np.pi*roc_mpc*(nfw[0]**3))/(3.0*Msun)
        e_M200_NFW =((800.0*np.pi*roc_mpc*(nfw[0]**2))/(Msun))*nfw[1]
        le_M200    = (np.log(10.)/M200_NFW)*e_M200_NFW
 
        # WRITING OUTPUT FITS FILE
        
        table_hdu = [fits.Column(name='Rp', format='D', array=R),
                fits.Column(name='DSigma_T', format='D', array=DSigma_T[0]),
                fits.Column(name='error_DSigma_T', format='D', array=eDSigma_T),
                fits.Column(name='DSigma_X', format='D', array=DSigma_X[0]),
                fits.Column(name='error_DSigma_X', format='D', array=eDSigma_X),
                fits.Column(name='GAMMA_Tcos_control', format='D', array=GAMMA_Tcos[0,:,0]),
                fits.Column(name='GAMMA_Tcos', format='D', array=GAMMA_Tcos[0,:,1]),
                fits.Column(name='GAMMA_Tcos_reduced', format='D', array=GAMMA_Tcos[0,:,2]),
                fits.Column(name='error_GAMMA_Tcos_control', format='D', array=eGAMMA_Tcos[:,0]),
                fits.Column(name='error_GAMMA_Tcos', format='D', array=eGAMMA_Tcos[:,1]),
                fits.Column(name='error_GAMMA_Tcos_reduced', format='D', array=eGAMMA_Tcos[:,2]),
                fits.Column(name='GAMMA_Xsin_control', format='D', array=GAMMA_Xsin[0,:,0]),
                fits.Column(name='GAMMA_Xsin', format='D', array=GAMMA_Xsin[0,:,1]),
                fits.Column(name='GAMMA_Xsin_reduced', format='D', array=GAMMA_Xsin[0,:,2]),
                fits.Column(name='error_GAMMA_Xsin_control', format='D', array=eGAMMA_Xsin[:,0]),
                fits.Column(name='error_GAMMA_Xsin', format='D', array=eGAMMA_Xsin[:,1]),
                fits.Column(name='error_GAMMA_Xsin_reduced', format='D', array=eGAMMA_Xsin[:,2])]
                
        for j in np.arange(1,101):
                
                table_hdu += [fits.Column(name='DSigma_T_K'+str(j), format='D', array=DSigma_T[j])]
                
                table_hdu += [fits.Column(name='GAMMA_Tcos_control_K'+str(j), format='D', array=GAMMA_Tcos[:,0,j])]
                table_hdu += [fits.Column(name='GAMMA_Tcos_K'+str(j), format='D', array=GAMMA_Tcos[:,1,j])]
                table_hdu += [fits.Column(name='GAMMA_Tcos_reduced_K'+str(j), format='D', array=GAMMA_Tcos[:,2,j])]
                
                table_hdu += [fits.Column(name='GAMMA_Xsin_control_K'+str(j), format='D', array=GAMMA_Xsin[:,0,j])]
                table_hdu += [fits.Column(name='GAMMA_Xsin_K'+str(j), format='D', array=GAMMA_Xsin[:,1,j])]
                table_hdu += [fits.Column(name='GAMMA_Xsin_reduced_K'+str(j), format='D', array=GAMMA_Xsin[:,2,j])]
        
        
        tbhdu = fits.BinTableHDU.from_columns(table_hdu)
        
        h = tbhdu.header
        h.append(('N_LENSES',np.int(Nlenses)))
        h.append(('lM_min',np.int(lM_min)))
        h.append(('lM_max',np.int(lM_max)))
        h.append(('z_min',np.round(z_min,4)))
        h.append(('z_max',np.round(z_max,4)))
        h.append(('lM200_NFW',np.round(np.log10(M200_NFW),4)))
        h.append(('elM200_NFW',np.round(le_M200,4)))
        h.append(('CHI2_NFW',np.round(nfw[2],4)))
        h.append(('lM_mean',np.round(lM_mean,4)))
        h.append(('z_mean',np.round(zmean,4)))
        h.append(('q2d_mean',np.round(q2d_mean,4)))
        h.append(('q2dr_mean',np.round(q2dr_mean,4)))
        h.append(('q3d_mean',np.round(q3d_mean,4)))
        h.append(('q3dr_mean',np.round(q3dr_mean,4)))
        h.append(('s3d_mean',np.round(s3d_mean,4)))        
        h.append(('s3dr_mean',np.round(s3dr_mean,4))) 
        
        tbhdu.writeto(folder+'profile_'+sample+'.fits',overwrite=True)
                
        tfin = time.time()
        
        print('TOTAL TIME ',(tfin-tini)/60.)
        


if __name__ == '__main__':
        
        parser = argparse.ArgumentParser()
        parser.add_argument('-sample', action='store', dest='sample',default='pru')
        parser.add_argument('-lM_min', action='store', dest='lM_min', default=14.)
        parser.add_argument('-lM_max', action='store', dest='lM_max', default=14.2)
        parser.add_argument('-z_min', action='store', dest='z_min', default=0.1)
        parser.add_argument('-z_max', action='store', dest='z_max', default=0.4)
        parser.add_argument('-RIN', action='store', dest='RIN', default=400.)
        parser.add_argument('-ROUT', action='store', dest='ROUT', default=5000.)
        parser.add_argument('-nbins', action='store', dest='nbins', default=40)
        parser.add_argument('-ncores', action='store', dest='ncores', default=10)
        parser.add_argument('-h_cosmo', action='store', dest='h_cosmo', default=1.)
        args = parser.parse_args()
        
        sample     = args.sample
        lM_min     = float(args.lM_min) 
        lM_max     = float(args.lM_max) 
        z_min      = float(args.z_min) 
        z_max      = float(args.z_max) 
        RIN        = float(args.RIN)
        ROUT       = float(args.ROUT)
        nbins      = int(args.nbins)
        ncores     = int(args.ncores)
        h          = float(args.h_cosmo)
        
        main(sample,lM_min,lM_max,z_min,z_max,RIN,ROUT,nbins,ncores,h)
