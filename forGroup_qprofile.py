import sys
sys.path.append('/mnt/projects/lensing')
sys.path.append('/mnt/projects/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import LambdaCDM
from astropy.wcs import WCS
from maria_func import *
from fit_profiles_curvefit import *
from astropy.stats import bootstrap
from astropy.utils import NumpyRNGContext
from multiprocessing import Pool
from multiprocessing import Process
import argparse
from astropy.constants import G,c,M_sun,pc
from scipy import stats
# For map
wcs = WCS(naxis=2)
wcs.wcs.crpix = [0., 0.]
wcs.wcs.cdelt = [1./3600., 1./3600.]
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]    

#parameters
cvel = c.value;   # Speed of light (m.s-1)
G    = G.value;   # Gravitational constant (m3.kg-1.s-2)
pc   = pc.value # 1 pc (m)
Msun = M_sun.value # Solar mass (kg)


parser = argparse.ArgumentParser()
parser.add_argument('-sample', action='store', dest='sample',default='pru')
parser.add_argument('-vmice', action='store', dest='vmice',default='2')
parser.add_argument('-lens_cat', action='store', dest='lcat',default='HALO_Props_MICE.fits')
parser.add_argument('-lM_min', action='store', dest='lM_min', default=14.)
parser.add_argument('-lM_max', action='store', dest='lM_max', default=15.5)
parser.add_argument('-z_min', action='store', dest='z_min', default=0.1)
parser.add_argument('-z_max', action='store', dest='z_max', default=0.4)
parser.add_argument('-q_min', action='store', dest='q_min', default=0.)
parser.add_argument('-q_max', action='store', dest='q_max', default=1.)
parser.add_argument('-rs_min', action='store', dest='rs_min', default=0.)
parser.add_argument('-rs_max', action='store', dest='rs_max', default=1.)
parser.add_argument('-relax', action='store', dest='relax', default='False')
parser.add_argument('-domap', action='store', dest='domap', default='False')
parser.add_argument('-RIN', action='store', dest='RIN', default=100.)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=10000.)
parser.add_argument('-nbins', action='store', dest='nbins', default=40)
parser.add_argument('-ncores', action='store', dest='ncores', default=10)
parser.add_argument('-h_cosmo', action='store', dest='h_cosmo', default=1.)
parser.add_argument('-ides_list', action='store', dest='idlist', default=None)
parser.add_argument('-resNFW_max', action='store', dest='resNFW_max', default=100.)
args = parser.parse_args()

sample     = args.sample
idlist     = args.idlist
lcat       = args.lcat
vmice      = args.vmice
lM_min     = float(args.lM_min)
lM_max     = float(args.lM_max) 
z_min      = float(args.z_min) 
z_max      = float(args.z_max) 
q_min      = float(args.q_min) 
q_max      = float(args.q_max) 
rs_min     = float(args.rs_min) 
rs_max     = float(args.rs_max) 
RIN        = float(args.RIN)
ROUT       = float(args.ROUT)
ndots      = int(args.nbins)
ncores     = int(args.ncores)
hcosmo     = float(args.h_cosmo)
resNFW_max = float(args.resNFW_max)
if args.relax == 'True':
    relax = True
elif args.relax == 'False':
    relax = False

if args.domap == 'True':
    domap = True
elif args.domap == 'False':
    domap = False

'''
lcat = 'HALO_Props_MICE.fits'
sample='pru'
lM_min=14.0
lM_max=14.5
z_min = 0.1
z_max = 0.2
q_min = 0.
q_max = 1.
rs_min = 0
rs_max = 1
RIN = 100.
ROUT = 10000.
ndots= 40
ncores = 32
hcosmo = 1.0 
vmice = '2'
rmin = 0.
rmax = 1000.
idlist = None
relax = False
domap = True
'''


folder = '/home/elizabeth/MICE/HS-lensing/'
S      = fits.open('/mnt/projects/lensing/HALO_SHAPE/MICEv'+vmice+'.0/catalogs/MICE_sources_HSN.fits')[1].data

def partial_map(RA0,DEC0,Z,angles,
                RIN,ROUT,ndots,h):

        
        lsize = int(np.sqrt(ndots))
        ndots = int(ndots)
        
        cosmo = LambdaCDM(H0=100*h, Om0=0.25, Ode0=0.75)
        dl  = cosmo.angular_diameter_distance(Z).value
        KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
        
        delta = (ROUT*1.5)/(3600*KPCSCALE)
        rarray = [[-1.*(ROUT*1.e-3),(ROUT*1.e-3)],[-1.*(ROUT*1.e-3),(ROUT*1.e-3)]]
        
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
        
       
        #get tangential ellipticities 
        et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
        #get cross ellipticities
        ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c
        # '''
        k  = catdata.kappa*sigma_c
                
        wcs.wcs.crval = [RA0,DEC0]
        dx, dy = wcs.wcs_world2pix(catdata.ra,catdata.dec, 0)

        dx = dx*KPCSCALE*1.e-3
        dy = dy*KPCSCALE*1.e-3
          
        del(e1)
        del(e2)
        
        r=np.rad2deg(rads)*3600*KPCSCALE
        del(rads)
        
        
        Ntot = len(catdata)
        del(catdata)    
        
        
        t = np.tile(theta_ra,(3,1)).T
        an = np.tile(angles,(len(theta_ra),1))
        at     = t - an
        
        GTsum = np.zeros((ndots,3))
        GXsum = np.zeros((ndots,3))
        Ksum  = np.zeros((ndots,3))
        Ninbin  = np.zeros((ndots,3))
        
        
        for j in range(3):
                x_t  = (dx*np.cos(an[:,j]) + dy*np.sin(an[:,j])) 
                y_t = (-1.*dx*np.sin(an[:,j]) + dy*np.cos(an[:,j])) 
                m   =  (abs(x_t) < (ROUT*1.e-3))*(abs(y_t) < (ROUT*1.e-3)) 
                out = stats.binned_statistic_2d(x_t[m], y_t[m], et[m], 
                                               'count',bins=lsize, range = rarray)
                Ninbin[:,j] = out.statistic.flatten()
                out = stats.binned_statistic_2d(x_t[m], y_t[m], et[m], 
                                               'sum',bins=lsize, range = rarray)
                GTsum[:,j] = out.statistic.flatten()
                out = stats.binned_statistic_2d(x_t[m], y_t[m], ex[m], 
                                               'sum',bins=lsize, range = rarray)
                GXsum[:,j] = out.statistic.flatten()
                out = stats.binned_statistic_2d(x_t[m], y_t[m], k[m], 
                                               'sum',bins=lsize, range = rarray)
                Ksum[:,j] = out.statistic.flatten()

        Ntot = m.sum()
   
        output = {'GTsum':GTsum,'GXsum':GXsum,
                   'Ksum': Ksum, 'N_inbin':Ninbin,
                   'Ntot': Ntot}
        
        return output

def partial_map_unpack(minput):
	return partial_map(*minput)

def partial_profile(RA0,DEC0,Z,angles,
                    RIN,ROUT,ndots,h):

        ndots = int(ndots)

        cosmo = LambdaCDM(H0=100*h, Om0=0.25, Ode0=0.75)
        dl  = cosmo.angular_diameter_distance(Z).value
        KPCSCALE   = dl*(((1.0/3600.0)*np.pi)/180.0)*1000.0
        
        delta = ROUT/(3600*KPCSCALE)
        
        t0 = time.time()
        mask = (S.ra < (RA0+delta))&(S.ra > (RA0-delta))&(S.dec > (DEC0-delta))&(S.dec < (DEC0+delta))&(S.z_v > (Z+0.1))
        t1 = time.time()  
        # print(t1-t0)             
        catdata = S[mask]

        ds  = cosmo.angular_diameter_distance(catdata.z_v).value
        dls = cosmo.angular_diameter_distance_z1z2(Z, catdata.z_v).value
                
        
        BETA_array = dls/ds
        
        Dl = dl*1.e6*pc
        sigma_c = (((cvel**2.0)/(4.0*np.pi*G*Dl))*(1./BETA_array))*(pc**2/Msun)

        t2 = time.time()
        # print(RA0,DEC0,t2-t1)

        rads, theta, test1,test2 = eq2p2(np.deg2rad(catdata.ra),
                                        np.deg2rad(catdata.dec),
                                        np.deg2rad(RA0),
                                        np.deg2rad(DEC0))
        
        theta2 = (2.*np.pi - theta) +np.pi/2.
        theta_ra = theta2
        theta_ra[theta2 > 2.*np.pi] = theta2[theta2 > 2.*np.pi] - 2.*np.pi
                       
        e1     = catdata.gamma1
        e2     = -1.*catdata.gamma2
        
       
        #get tangential ellipticities 
        et = (-e1*np.cos(2*theta)-e2*np.sin(2*theta))*sigma_c
        #get cross ellipticities
        ex = (-e1*np.sin(2*theta)+e2*np.cos(2*theta))*sigma_c
        # '''
        k  = catdata.kappa*sigma_c
          
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
        
        
        SIGMAwsum    = []
        DSIGMAwsum_T = []
        DSIGMAwsum_X = []
        N_inbin = []
              
        GAMMATcos_wsum = np.zeros((ndots,3))
        GAMMAXsin_wsum = np.zeros((ndots,3))
                
        COS2_2theta = np.zeros((ndots,3))
        SIN2_2theta = np.zeros((ndots,3))
                       
        for nbin in range(ndots):
                mbin = dig == nbin+1              
                
                SIGMAwsum    = np.append(SIGMAwsum,k[mbin].sum())
                DSIGMAwsum_T = np.append(DSIGMAwsum_T,et[mbin].sum())
                DSIGMAwsum_X = np.append(DSIGMAwsum_X,ex[mbin].sum())

                GAMMATcos_wsum[nbin,:] = np.sum((np.tile(et[mbin],(3,1))*np.cos(2.*at[mbin]).T),axis=1)
                GAMMAXsin_wsum[nbin,:] = np.sum((np.tile(ex[mbin],(3,1))*np.sin(2.*at[mbin]).T),axis=1)
                
                COS2_2theta[nbin,:] = np.sum((np.cos(2.*at[mbin]).T)**2,axis=1)
                SIN2_2theta[nbin,:] = np.sum((np.sin(2.*at[mbin]).T)**2,axis=1)
                               
                N_inbin = np.append(N_inbin,len(et[mbin]))
                
                index = np.arange(mbin.sum())
                '''
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

                '''
        
        output = {'SIGMAwsum':SIGMAwsum,'DSIGMAwsum_T':DSIGMAwsum_T,
                   'DSIGMAwsum_X':DSIGMAwsum_X,
                   'GAMMATcos_wsum': GAMMATcos_wsum, 'GAMMAXsin_wsum': GAMMAXsin_wsum,
                   'COS2_2theta_wsum':COS2_2theta,'SIN2_2theta_wsum':SIN2_2theta,
                   'N_inbin':N_inbin,'Ntot':Ntot}
        
        return output

def partial_profile_unpack(minput):
	return partial_profile(*minput)

def cov_matrix(array):
        
        K = len(array)
        Kmean = np.average(array,axis=0)
        bins = array.shape[1]
        
        COV = np.zeros((bins,bins))
        
        for k in range(K):
            dif = (array[k]- Kmean)
            COV += np.outer(dif,dif)        
        
        COV *= (K-1)/K
        return COV
        
def main(lcat, sample='pru',
         lM_min=14.,lM_max=14.2,
         z_min = 0.1, z_max = 0.4,
         q_min = 0., q_max = 1.0,
         rs_min = 0., rs_max = 1.0,
         resNFW_max = 100., relax=False,
         domap = False, RIN = 400., ROUT =5000.,
         ndots= 40, ncores=10, 
         idlist= None, hcosmo=1.0, vmice = '2'):

        '''
        
        INPUT
        ---------------------------------------------------------
        sample         (str) sample name
        lM_min         (float) lower limit for log(Mass) - >=
        lM_max         (float) higher limit for log(Mass) - <
        z_min          (float) lower limit for z - >=
        z_max          (float) higher limit for z - <
        q_min          (float) lower limit for q - >=
        q_max          (float) higher limit for q - <
        rs_min         (float) lower limit r_scale = r_c/r_max
        rs_max         (float) higher limit r_scale = r_c/r_max
        resNFW_max     (float) higher limit for resNFW_S
        relax          (bool)  Select only relaxed halos
        domap          (bool) Instead of computing a profile it 
                       will compute a map with 2D bins ndots lsize
        RIN            (float) Inner bin radius of profile
        ROUT           (float) Outer bin radius of profile
        ndots          (int) Number of bins of the profile
        ncores         (int) to run in parallel, number of cores
        h              (float) H0 = 100.*h
        '''

        cosmo = LambdaCDM(H0=100*hcosmo, Om0=0.25, Ode0=0.75)
        tini = time.time()
        
        print('Lens catalog ',lcat)
        print('Sample ',sample)
        print('Selecting groups with:')
        
        
        if idlist:
                print('From id list '+idlist)
        else:
                print(lM_min,' <= log(M) < ',lM_max)
                print(z_min,' <= z < ',z_max)
                print(q_min,' <= q < ',q_max)
                print(rs_min,' <= rs < ',rs_max)
                print('resNFW_S < ',resNFW_max)
                print('h ',hcosmo)
                        
        #reading cats
                
        L = fits.open(folder+lcat)[1].data               

        '''
        # To try all centre
        
        ra = np.rad2deg(np.arctan(L.xc/L.yc))
        ra[L.yc==0] = 90.
        dec = np.rad2deg(np.arcsin(L.zc/sqrt(L.xc**2 + L.yc**2 + L.zc**2)))

        L.ra_rc = ra
        L.dec_rc = dec
        '''
        ra = L.ra_rc
        dec = L.dec_rc
        
        Eratio = (2.*L.K/abs(L.U))
                               
        if idlist:
                ides = np.loadtxt(idlist).astype(int)
                mlenses = np.in1d(L.unique_halo_id,ides)
        else:
                
                rs      = L.offset
                mmass   = (L.lgM >= lM_min)*(L.lgM < lM_max)
                mz      = (L.z >= z_min)*(L.z < z_max)
                mq      = (L.q2d >= q_min)*(L.q2d < q_max)
                mrs     = (rs >= rs_min)*(rs < rs_max)
                mres    = L.resNFW_S < resNFW_max
                mlenses = mmass*mz*mq*mrs*mres
                
        if relax:
            print('Only relaxed halos are considered')
            sample = sample+'_relaxed'
            mlenses = mlenses*(L.offset < 0.1)*(Eratio < 1.35)
                
        Nlenses = mlenses.sum()

        if Nlenses < ncores:
                ncores = Nlenses
        
        print('Nlenses',Nlenses)
        print('CORRIENDO EN ',ncores,' CORES')
        
        L = L[mlenses]
                
        #Computing SMA axis
        theta  = np.array([np.zeros(sum(mlenses)),np.arctan(L.a2Dy/L.a2Dx),np.arctan(L.a2Dry/L.a2Drx)]).T                
        
        # Define K masks        
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
        
        if domap:

            print('Maps have ',ndots,'pixels')
            print('up to ',ROUT,'kpc')
            
            output_file = 'maps/map_'+sample+'.fits'

            # Defining 2D bins
            
            lsize = int(np.sqrt(ndots))
            
            ndots = int(lsize**2)
            
            Rmpc = 1.e-3*ROUT
            
            xb = np.linspace(-1.*Rmpc+(Rmpc/lsize),Rmpc-(Rmpc/lsize),lsize)
            
            ybin, xbin = np.meshgrid(xb,xb)
            
            ybin = ybin.flatten()
            xbin = xbin.flatten()

            # WHERE THE SUMS ARE GOING TO BE SAVED
            GTsum = np.zeros((101,ndots,3))
            GXsum = np.zeros((101,ndots,3))
            Ksum  = np.zeros((101,ndots,3))        
            
            Ninbin = np.zeros((101,ndots,3))
            
            # FUNCTION TO RUN IN PARALLEL
            partial = partial_map_unpack
            

        else:

            print('Profile has ',ndots,'bins')
            print('from ',RIN,'kpc to ',ROUT,'kpc')

            output_file = 'profiles/profile_'+sample+'.fits'

            # Defining radial bins
            bines = np.logspace(np.log10(RIN),np.log10(ROUT),num=ndots+1)
            R = (bines[:-1] + np.diff(bines)*0.5)*1.e-3

            # WHERE THE SUMS ARE GOING TO BE SAVED
            
            Ninbin = np.zeros((101,ndots))
            
            SIGMAwsum    = np.zeros((101,ndots)) 
            DSIGMAwsum_T = np.zeros((101,ndots)) 
            DSIGMAwsum_X = np.zeros((101,ndots))
                
            GAMMATcos_wsum = np.zeros((101,ndots,3))
            GAMMAXsin_wsum = np.zeros((101,ndots,3))
    
            COS2_2theta_wsum = np.zeros((101,ndots,3))
            SIN2_2theta_wsum = np.zeros((101,ndots,3))
            
            # FUNCTION TO RUN IN PARALLEL
            partial = partial_profile_unpack
            

        print('Saved in '+folder+output_file)
                                   
        
        
        Ntot         = np.array([])
        tslice       = np.array([])
        
        for l in range(len(Lsplit)):
                
                print('RUN ',l+1,' OF ',len(Lsplit))
                
                t1 = time.time()
                
                num = len(Lsplit[l])
                
                rin  = RIN*np.ones(num)
                rout = ROUT*np.ones(num)
                nd   = ndots*np.ones(num)
                h_array   = hcosmo*np.ones(num)
                
                if num == 1:
                        entrada = [Lsplit[l].ra_rc[0], Lsplit[l].dec_rc[0],
                                   Lsplit[l].z[0],Tsplit[l][0],
                                   RIN,ROUT,ndots,hcosmo]
                        
                        salida = [partial(entrada)]
                else:          
                        entrada = np.array([Lsplit[l].ra_rc,Lsplit[l].dec_rc,
                                        Lsplit[l].z,Tsplit[l].tolist(),
                                        rin,rout,nd,h_array]).T
                        
                        pool = Pool(processes=(num))
                        salida = np.array(pool.map(partial, entrada))
                        pool.terminate()
                                
                for j in range(len(salida)):
                        
                        profilesums = salida[j]
                        Ntot         = np.append(Ntot,profilesums['Ntot'])
                        
                        if domap:
                            
                            km      = np.tile(Ksplit[l][j],(3,ndots,1)).T
                            Ninbin += np.tile(profilesums['N_inbin'],(101,1,1))*km

                            GTsum += np.tile(profilesums['GTsum'],(101,1,1))*km
                            GXsum += np.tile(profilesums['GXsum'],(101,1,1))*km
                            Ksum  += np.tile(profilesums['Ksum'],(101,1,1))*km
                            
                        else:

                            km      = np.tile(Ksplit[l][j],(3,ndots,1)).T

                            Ninbin += np.tile(profilesums['N_inbin'],(101,1))*km[:,:,0]
                                                
                            SIGMAwsum    += np.tile(profilesums['SIGMAwsum'],(101,1))*km[:,:,0]
                            DSIGMAwsum_T += np.tile(profilesums['DSIGMAwsum_T'],(101,1))*km[:,:,0]
                            DSIGMAwsum_X += np.tile(profilesums['DSIGMAwsum_X'],(101,1))*km[:,:,0]
                            
                            GAMMATcos_wsum += np.tile(profilesums['GAMMATcos_wsum'],(101,1,1))*km
                            GAMMAXsin_wsum += np.tile(profilesums['GAMMAXsin_wsum'],(101,1,1))*km
    
                            COS2_2theta_wsum += np.tile(profilesums['COS2_2theta_wsum'],(101,1,1))*km
                            SIN2_2theta_wsum += np.tile(profilesums['SIN2_2theta_wsum'],(101,1,1))*km
                        
                
                t2 = time.time()
                ts = (t2-t1)/60.
                tslice = np.append(tslice,ts)
                print('TIME SLICE')
                print(ts)
                print('Estimated ramaining time')
                print(np.mean(tslice)*(len(Lsplit)-(l+1)))

        # AVERAGE LENS PARAMETERS AND SAVE IT IN HEADER
        
        zmean        = np.average(L.z,weights=Ntot)
        lM_mean      = np.log10(np.average(10**L.lgM,weights=Ntot))
        c200_mean    = np.average(L.cNFW_S,weights=Ntot)
        lM200_mean   = np.log10(np.average(10**L.lgMNFW_S,weights=Ntot))
        
        q2d_mean     = np.average(L.q2d,weights=Ntot)
        q2dr_mean    = np.average(L.q2dr,weights=Ntot)
        q3d_mean     = np.average(L.q,weights=Ntot)
        q3dr_mean    = np.average(L.qr,weights=Ntot)
        s3d_mean     = np.average(L.s,weights=Ntot)
        s3dr_mean    = np.average(L.sr,weights=Ntot)

        h = fits.Header()
        h.append(('N_LENSES',np.int(Nlenses)))
        h.append(('Lens cat',lcat))
        h.append(('MICE version sources',vmice))
        h.append(('rs_min',np.round(rs_min,1)))
        h.append(('rs_max',np.round(rs_max,1)))
        h.append(('lM_min',np.round(lM_min,2)))
        h.append(('lM_max',np.round(lM_max,2)))
        h.append(('z_min',np.round(z_min,2)))
        h.append(('z_max',np.round(z_max,2)))
        h.append(('q_min',np.round(q_min,2)))
        h.append(('q_max',np.round(q_max,2)))
        h.append(('resNFW_max',np.round(resNFW_max,2)))
        h.append(('hcosmo',np.round(hcosmo,4)))
        h.append(('lM_mean',np.round(lM_mean,4)))
        h.append(('lM200_mean',np.round(lM200_mean,4)))
        h.append(('c200_mean',np.round(c200_mean,4)))
        h.append(('z_mean',np.round(zmean,4)))
        h.append(('q2d_mean',np.round(q2d_mean,4)))
        h.append(('q2dr_mean',np.round(q2dr_mean,4)))
        h.append(('q3d_mean',np.round(q3d_mean,4)))
        h.append(('q3dr_mean',np.round(q3dr_mean,4)))
        h.append(('s3d_mean',np.round(s3d_mean,4)))        
        h.append(('s3dr_mean',np.round(s3dr_mean,4))) 
        
        if domap:
            
            # COMPUTING PROFILE        
                
            GT  = (GTsum/Ninbin)
            GX  = (GXsum/Ninbin)
            K   = (Ksum/Ninbin)
                    
            # COMPUTE COVARIANCE
            
            COV_Gtc = cov_matrix(GT[1:,:,0])
            COV_Gt  = cov_matrix(GT[1:,:,1])
            COV_Gtr = cov_matrix(GT[1:,:,2])
            
            COV_Gxc = cov_matrix(GX[1:,:,0])
            COV_Gx  = cov_matrix(GX[1:,:,1])
            COV_Gxr = cov_matrix(GX[1:,:,2])
    
            COV_Kc = cov_matrix(GX[1:,:,0])
            COV_K  = cov_matrix(GX[1:,:,1])
            COV_Kr = cov_matrix(GX[1:,:,2])
            
            table_pro = [fits.Column(name='xmpc', format='E', array=xbin),
                    fits.Column(name='ympc', format='E', array=ybin),
                    fits.Column(name='GT_control', format='E', array=GT[0,:,0]),
                    fits.Column(name='GT', format='E', array=GT[0,:,1]),
                    fits.Column(name='GT_reduced', format='E', array=GT[0,:,2]),
                    fits.Column(name='GX_control', format='E', array=GX[0,:,0]),
                    fits.Column(name='GX', format='E', array=GX[0,:,1]),
                    fits.Column(name='GX_reduced', format='E', array=GX[0,:,2]),
                    fits.Column(name='K_control', format='E', array=K[0,:,0]),
                    fits.Column(name='K', format='E', array=K[0,:,1]),
                    fits.Column(name='K_reduced', format='E', array=K[0,:,2])]
                    
                        
            table_cov = [fits.Column(name='COV_GT_control', format='E', array=COV_Gtc.flatten()),
                        fits.Column(name='COV_GT', format='E', array=COV_Gt.flatten()),
                        fits.Column(name='COV_GT_reduced', format='E', array=COV_Gtr.flatten()),
                        fits.Column(name='COV_GX_control', format='E', array=COV_Gxc.flatten()),
                        fits.Column(name='COV_GX', format='E', array=COV_Gx.flatten()),
                        fits.Column(name='COV_GX_reduced', format='E', array=COV_Gxr.flatten()),
                        fits.Column(name='COV_K_control', format='E', array=COV_Kc.flatten()),
                        fits.Column(name='COV_K', format='E', array=COV_K.flatten()),
                        fits.Column(name='COV_K_reduced', format='E', array=COV_Kr.flatten())]
    
    
        else:
            # COMPUTING PROFILE        
            Ninbin[DSIGMAwsum_T == 0] = 1.
                    
            Sigma     = (SIGMAwsum/Ninbin)
            DSigma_T  = (DSIGMAwsum_T/Ninbin)
            DSigma_X  = (DSIGMAwsum_X/Ninbin)
            
            GAMMA_Tcos = (GAMMATcos_wsum/COS2_2theta_wsum).transpose(1,2,0)
            GAMMA_Xsin = (GAMMAXsin_wsum/SIN2_2theta_wsum).transpose(1,2,0)
            
            # COMPUTE COVARIANCE
    
            COV_S   = cov_matrix(Sigma[1:,:])
            COV_St  = cov_matrix(DSigma_T[1:,:])
            COV_Sx  = cov_matrix(DSigma_X[1:,:])
            
            COV_Gtc = cov_matrix(GAMMA_Tcos[:,0,1:].T)
            COV_Gt  = cov_matrix(GAMMA_Tcos[:,1,1:].T)
            COV_Gtr = cov_matrix(GAMMA_Tcos[:,2,1:].T)
            
            COV_Gxc = cov_matrix(GAMMA_Xsin[:,0,1:].T)
            COV_Gx  = cov_matrix(GAMMA_Xsin[:,1,1:].T)
            COV_Gxr = cov_matrix(GAMMA_Xsin[:,2,1:].T)
            
            # FITTING NFW PROFILE AND SAVE IT IN HEADER
            
            nfw    = Delta_Sigma_fit(R,DSigma_T[0],np.diag(COV_St),zmean,cosmo,True)
    
            M200_NFW   = nfw.M200
            c200_NFW   = nfw.c200
            e_c200_NFW = nfw.error_c200
            e_M200_NFW = nfw.error_M200
            le_M200    = (np.log(10.)/M200_NFW)*e_M200_NFW
        
            h.append(('lM200_NFW',np.round(np.log10(M200_NFW),4)))
            h.append(('elM200_NFW',np.round(le_M200,4)))
            h.append(('c200_NFW',np.round(c200_NFW,4)))
            h.append(('ec200_NFW',np.round(e_c200_NFW,4)))
            h.append(('CHI2_NFW',np.round(nfw.chi2,4)))

            # WRITING OUTPUT FITS FILE
            
            table_pro = [fits.Column(name='Rp', format='E', array=R),
                    fits.Column(name='Sigma', format='E', array=Sigma[0]),
                    fits.Column(name='DSigma_T', format='E', array=DSigma_T[0]),
                    fits.Column(name='DSigma_X', format='E', array=DSigma_X[0]),
                    fits.Column(name='GAMMA_Tcos_control', format='E', array=GAMMA_Tcos[:,0,0]),
                    fits.Column(name='GAMMA_Tcos', format='E', array=GAMMA_Tcos[:,1,0]),
                    fits.Column(name='GAMMA_Tcos_reduced', format='E', array=GAMMA_Tcos[:,2,0]),
                    fits.Column(name='GAMMA_Xsin_control', format='E', array=GAMMA_Xsin[:,0,0]),
                    fits.Column(name='GAMMA_Xsin', format='E', array=GAMMA_Xsin[:,1,0]),
                    fits.Column(name='GAMMA_Xsin_reduced', format='E', array=GAMMA_Xsin[:,2,0])]
                    
                        
            table_cov = [fits.Column(name='COV_ST', format='E', array=COV_St.flatten()),
                        fits.Column(name='COV_S', format='E', array=COV_S.flatten()),
                        fits.Column(name='COV_SX', format='E', array=COV_Sx.flatten()),
                        fits.Column(name='COV_GT_control', format='E', array=COV_Gtc.flatten()),
                        fits.Column(name='COV_GT', format='E', array=COV_Gt.flatten()),
                        fits.Column(name='COV_GT_reduced', format='E', array=COV_Gtr.flatten()),
                        fits.Column(name='COV_GX_control', format='E', array=COV_Gxc.flatten()),
                        fits.Column(name='COV_GX', format='E', array=COV_Gx.flatten()),
                        fits.Column(name='COV_GX_reduced', format='E', array=COV_Gxr.flatten())]
        
        tbhdu_pro = fits.BinTableHDU.from_columns(fits.ColDefs(table_pro))
        tbhdu_cov = fits.BinTableHDU.from_columns(fits.ColDefs(table_cov))
        
        
        primary_hdu = fits.PrimaryHDU(header=h)
        
        hdul = fits.HDUList([primary_hdu, tbhdu_pro, tbhdu_cov])
        
        hdul.writeto(folder+output_file,overwrite=True)
                
        tfin = time.time()
        
        print('TOTAL TIME ',(tfin-tini)/60.)
        


main(lcat,sample,lM_min,lM_max,z_min,z_max,q_min,q_max,rs_min,rs_max,resNFW_max,relax,domap,RIN,ROUT,ndots,ncores,idlist,hcosmo,vmice)
