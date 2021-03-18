import sys
sys.path.append('/mnt/projects/lensing')
sys.path.append('/mnt/projects/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM
from maria_func import *
from fit_profiles_curvefit import *
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from multiprocessing import Pool
from multiprocessing import Process
import argparse
from astropy.constants import G,c,M_sun,pc

###################
# ARREGLAR RPROX 14
##################

#parameters
cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)


parser = argparse.ArgumentParser()
parser.add_argument('-sample', action='store', dest='sample',default='pru')
parser.add_argument('-vmice', action='store', dest='vmice',default='2')
parser.add_argument('-rprox', action='store', dest='rprox',default='Rprox_lM14cut')
parser.add_argument('-rmin', action='store', dest='rmin',default=0.)
parser.add_argument('-rmax', action='store', dest='rmax',default=1000.)
parser.add_argument('-lM_min', action='store', dest='lM_min', default=14.)
parser.add_argument('-lM_max', action='store', dest='lM_max', default=15.5)
parser.add_argument('-z_min', action='store', dest='z_min', default=0.1)
parser.add_argument('-z_max', action='store', dest='z_max', default=0.4)
parser.add_argument('-q_min', action='store', dest='q_min', default=0.)
parser.add_argument('-q_max', action='store', dest='q_max', default=1.)
parser.add_argument('-RIN', action='store', dest='RIN', default=100.)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=10000.)
parser.add_argument('-nbins', action='store', dest='nbins', default=40)
parser.add_argument('-ncores', action='store', dest='ncores', default=10)
parser.add_argument('-h_cosmo', action='store', dest='h_cosmo', default=1.)
args = parser.parse_args()

sample     = args.sample
rprox      = args.rprox
rmin       = float(args.rmin)
rmax       = float(args.rmax) 
lM_min     = float(args.lM_min)
lM_max     = float(args.lM_max) 
z_min      = float(args.z_min) 
z_max      = float(args.z_max) 
q_min      = float(args.q_min) 
q_max      = float(args.q_max) 
RIN        = float(args.RIN)
ROUT       = float(args.ROUT)
ndots      = int(args.nbins)
ncores     = int(args.ncores)
vmice      = int(args.vmice)
hcosmo     = float(args.h_cosmo)

'''
sample='pru'
lM_min=14.
lM_max=14.4
z_min = 0.1
z_max = 0.12
q_min = 0.
q_max = 1.
RIN = 100.
ROUT = 10000.
ndots= 40
ncores = 40
hcosmo = 1.0 
vmice = 2
rmin = 0.
rmax = 1000.
rprox = 'Rprox_lM14cut'
'''


folder = '/mnt/projects/lensing/HALO_SHAPE/MICEv'+str(vmice)+'.0/'
S      = fits.open(folder+'catalogs/MICE_sources_HSN.fits')[1].data


def partial_profile(RA0,DEC0,Z,RIN,ROUT,ndots,h):

        ndots = int(ndots)

        cosmo = LambdaCDM(H0=100*h, Om0=0.25, Ode0=0.75)
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


        r=np.rad2deg(rads)*3600*KPCSCALE
                      
        e1     = catdata.gamma1
        e2     = -1.*catdata.gamma2
        
       
        #get tangential ellipticities 
        et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c

        k  = catdata.kappa*sigma_c
        
        r=np.rad2deg(rads)*3600*KPCSCALE      
                
        bines = np.logspace(np.log10(RIN),np.log10(ROUT),num=ndots+1)
        dig = np.digitize(r,bines)
                
        SIGMA = []
        DSIGMA = []
                       
        for nbin in range(ndots):
                mbin = dig == nbin+1              
                
                SIGMA = np.append(SIGMA,np.average(k[mbin]))
                DSIGMA = np.append(DSIGMA,np.average(k[mbin]))
                        
        return [SIGMA,DSIGMA]

def partial_profile_unpack(minput):
	return partial_profile(*minput)
            
        
        

def main(sample='pru', rprox = 'Rprox_lM14cut', 
                rmin = 0., rmax = 1000.,
                lM_min=14.,lM_max=14.2,
                z_min = 0.1, z_max = 0.4,
                q_min = 0., q_max = 1.0,
                RIN = 400., ROUT =5000.,
                ndots= 40, ncores=10, hcosmo=1.0):

        '''
        
        INPUT
        ---------------------------------------------------------
        sample         (str) sample name
        rprox          (str) radius to the neighbour criteria
        rmin           (float) lower limit for rprox - >=
        rmax           (float) higher limit for rprox - <
        lM_min         (float) lower limit for log(Mass) - >=
        lM_max         (float) higher limit for log(Mass) - <
        z_min          (float) lower limit for z - >=
        z_max          (float) higher limit for z - <
        q_min          (float) lower limit for q - >=
        q_max          (float) higher limit for q - <
        RIN            (float) Inner bin radius of profile
        ROUT           (float) Outer bin radius of profile
        ndots          (int) Number of bins of the profile
        ncores         (int) to run in parallel, number of cores
        h              (float) H0 = 100.*h
        '''

        cosmo = LambdaCDM(H0=100*hcosmo, Om0=0.25, Ode0=0.75)
        tini = time.time()
        
        print('Sample ',sample)
        print('Selecting groups with:')
        print(lM_min,' <= log(M) < ',lM_max)
        print(rmin,' <= '+rprox+' < ',rmax)
        print(z_min,' <= z < ',z_max)
        print(q_min,' <= q < ',q_max)
        print('Profile has ',ndots,'bins')
        print('from ',RIN,'kpc to ',ROUT,'kpc')
        print('h ',hcosmo)
        print('MICE version ',vmice)
              
        # Defining radial bins
        bines = np.logspace(np.log10(RIN),np.log10(ROUT),num=ndots+1)
        R = (bines[:-1] + np.diff(bines)*0.5)*1.e-3
        
        #reading cats
        
        L = fits.open(folder+'catalogs/MICE_halo_cat_withshapes.fits')[1].data
        
        try:
            mrcut   = (L[rprox] >= rmin)*(L[rprox] < rmax)
        except:
            print(rprox+' NOT FINDED')

            mrcut   = np.ones(len(L.ra)).astype(bool)
            
        mregion = (L.ra < 80.)*(L.dec < 50.)#*(L.dec > 36.5)        
        mmass   = (L.lgm >= lM_min)*(L.lgm < lM_max)
        mz      = (L.z_v >= z_min)*(L.z_v < z_max)
        mq      = (L.q2d >= q_min)*(L.q2d < q_max)
        mlenses = mmass*mz*mq*mrcut
        Nlenses = mlenses.sum()

        if Nlenses < ncores:
                ncores = Nlenses
        
        print('Nlenses',Nlenses)
        print('CORRIENDO EN ',ncores,' CORES')
        
        L = L[mlenses]
        
        ind = np.arange(len(L))        

        
        # SPLIT LENSING CAT
        
        lbins = int(round(Nlenses/float(ncores), 0))
        slices = ((np.arange(lbins)+1)*ncores).astype(int)
        slices = slices[(slices < Nlenses)]
        Lsplit = np.split(L,slices)
        Isplit = np.split(ind,slices)

                
        Ntot         = np.array([])

        table = [fits.Column(name='Rp', format='E', array=np.append(0,R))]
        
        for l in range(len(Lsplit)):
                
                print('RUN ',l+1,' OF ',len(Lsplit))
                
                t1 = time.time()
                
                num = len(Lsplit[l])
                
                rin  = RIN*np.ones(num)
                rout = ROUT*np.ones(num)
                nd   = ndots*np.ones(num)
                h_array   = hcosmo*np.ones(num)
                
                if num == 1:
                        entrada = [Lsplit[l].ra[0], Lsplit[l].dec[0],
                                   Lsplit[l].z_v[0],RIN,ROUT,ndots,hcosmo]
                        
                        salida = [partial_profile_unpack(entrada)]
                else:          
                        entrada = np.array([Lsplit[l].ra,Lsplit[l].dec,
                                        Lsplit[l].z_v,rin,rout,nd,h_array]).T
                        
                        pool = Pool(processes=(num))
                        salida = np.array(pool.map(partial_profile_unpack, entrada))
                        pool.terminate()
                                
                for j in range(len(salida)):
                        
                        S,DS = salida[j]
                        
                        S = np.append(Lsplit[l].unique_halo_id[j],S)
                        DS = np.append(Lsplit[l].unique_halo_id[j],DS)

                        table += [fits.Column(name='S'+str(Isplit[l][j]), format='E', array=S)]
                        table += [fits.Column(name='DS'+str(Isplit[l][j]), format='E', array=DS)]
                
                t2 = time.time()
                ts = (t2-t1)/60.
                tslice = np.append(tslice,ts)
                print('TIME SLICE')
                print(ts)
                print('Estimated ramaining time')
                print(np.mean(tslice)*(len(Lsplit)-(l+1)))
        
        
                
        
        
        # WRITING OUTPUT FITS FILE
                
        tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(table))
        
        h = fits.Header()
        h.append(('N_LENSES',np.int(Nlenses)))
        h.append(('MICE version',vmice))
        h.append((rprox+'_min',np.round(rmin,1)))
        h.append((rprox+'_max',np.round(rmax,1)))
        h.append(('lM_min',np.round(lM_min,2)))
        h.append(('lM_max',np.round(lM_max,2)))
        h.append(('z_min',np.round(z_min,2)))
        h.append(('z_max',np.round(z_max,2)))
        h.append(('q_min',np.round(q_min,2)))
        h.append(('q_max',np.round(q_max,2)))
        h.append(('hcosmo',np.round(hcosmo,4)))
        
        primary_hdu = fits.PrimaryHDU(header=h)
        
        hdul = fits.HDUList([primary_hdu, tbhdu])
        
        hdul.writeto(folder+'profiles/profile_'+sample+'_individual.fits',overwrite=True)
                
        tfin = time.time()
        
        print('TOTAL TIME ',(tfin-tini)/60.)
        


main(sample,rprox,rmin,rmax,lM_min,lM_max,z_min,z_max,q_min,q_max,RIN,ROUT,ndots,ncores,hcosmo)
