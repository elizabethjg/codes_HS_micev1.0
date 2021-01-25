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
from fit_profiles_curvefit import *
from multiprocessing import Pool
from multiprocessing import Process
import argparse
from astropy.constants import G,c,M_sun,pc
from astropy.wcs import WCS
from scipy import stats
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
parser.add_argument('-rprox', action='store', dest='rprox',default='Rprox_lM14cut')
parser.add_argument('-rmin', action='store', dest='rmin',default=0.)
parser.add_argument('-rmax', action='store', dest='rmax',default=5000.)
parser.add_argument('-lM_min', action='store', dest='lM_min', default=14.4)
parser.add_argument('-lM_max', action='store', dest='lM_max', default=14.6)
parser.add_argument('-z_min', action='store', dest='z_min', default=0.1)
parser.add_argument('-z_max', action='store', dest='z_max', default=0.4)
parser.add_argument('-q_min', action='store', dest='q_min', default=0.)
parser.add_argument('-q_max', action='store', dest='q_max', default=1.)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=5000.)
parser.add_argument('-nbins', action='store', dest='nbins', default=1000)
parser.add_argument('-ncores', action='store', dest='ncores', default=32)
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
ROUT       = float(args.ROUT)
ndots      = int(args.nbins)
ncores     = int(args.ncores)
vmice      = int(args.vmice)
hcosmo     = float(args.h_cosmo)

'''
sample='pru'
lM_min=14.
lM_max=14.2
z_min = 0.1
z_max = 0.4
q_min = 0.
q_max = 0.6
RIN = 400.
ROUT = 5000.
ndots= 40
ncores = 40
hcosmo = 1.0 
vmice = 2
'''


folder = '/mnt/clemente/lensing/HALO_SHAPE/MICE_v'+str(vmice)+'.0/catalogs/'
S      = fits.open(folder+'MICE_sources.fits')[1].data


def partial_map(RA0,DEC0,Z,angles,ROUT,ndots,h):

        lsize = int(np.sqrt(ndots))

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
        k  = catdata.kappa
                
        wcs.wcs.crval = [RA0,DEC0]
        dx, dy = wcs.wcs_world2pix(catdata.ra,catdata.dec, 0)

        dx = dx*KPCSCALE*1.e-3
        dy = dy*KPCSCALE*1.e-3
          
        # del(e1)
        # del(e2)
        
        r=np.rad2deg(rads)*3600*KPCSCALE
        del(rads)
        
        
        Ntot = len(catdata)
        # del(catdata)    
        
        
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
                   'Ksum': Ksum, 'Ninbin':Ninbin,
                   'Ntot': Ntot}
        
        return output

def partial_map_unpack(minput):
	return partial_map(*minput)

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

            
        
        

def main(sample='pru', rprox = 'Rprox_lM14cut', 
                rmin = 0., rmax = 1000.,
                lM_min=14.,lM_max=14.2,
                z_min = 0.1, z_max = 0.4,
                q_min = 0., q_max = 1.0,
                ROUT =5000., ndots= 40, 
                ncores=10, hcosmo=1.0):

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
        print('Map pixel size ',ndots,'bins')
        print('max',ROUT,'kpc')
        print('h ',hcosmo)
        print('MICE version ',vmice)
              
        # Defining 2D bins
        
        lsize = int(np.sqrt(ndots))
        
        ndots = int(lsize**2)
        
        Rmpc = 1.e-3*ROUT
        
        xb = np.linspace(-1.*Rmpc+(Rmpc/lsize),Rmpc-(Rmpc/lsize),lsize)
        
        ybin, xbin = np.meshgrid(xb,xb)
        
        ybin = ybin.flatten()
        xbin = xbin.flatten()
                
        #reading cats
        
        L = fits.open(folder+'MICE_halo_cat_withshapes.fits')[1].data
        
        try:
            mrcut   = (L[rprox] >= rmin)*(L[rprox] < rmax)
        except:
            print(rprox+' NOT FINDED')

            mrcut   = np.ones(len(L.ra)).astype(bool)
            
        mregion = (L.ra < 80.)*(L.dec < 50.)#*(L.dec > 36.5)        
        mmass   = (L.lgm >= lM_min)*(L.lgm < lM_max)
        mz      = (L.z_v >= z_min)*(L.z_v < z_max)
        mq      = (L.q2d >= q_min)*(L.q2d < q_max)
        mlenses = mmass*mz*mregion*mq*mrcut
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
                      
        GTsum = np.zeros((101,ndots,3))
        GXsum = np.zeros((101,ndots,3))
        Ksum  = np.zeros((101,ndots,3))
                                   
        Ninbin = np.zeros((101,ndots,3))
        
        tslice       = np.array([])
        Ntot         = np.array([])
        
        for l in range(len(Lsplit)):
                
                print('RUN ',l+1,' OF ',len(Lsplit))
                
                t1 = time.time()
                
                num = len(Lsplit[l])
                
                rout = ROUT*np.ones(num)
                nd   = (ndots*np.ones(num)).astype(int)
                h_array   = hcosmo*np.ones(num)
                
                if num == 1:
                        entrada = [Lsplit[l].ra[0], Lsplit[l].dec[0],
                                   Lsplit[l].z_v[0],Tsplit[l][0],
                                   ROUT,ndots,hcosmo]
                        
                        salida = [partial_map_unpack(entrada)]
                else:          
                        entrada = np.array([Lsplit[l].ra,Lsplit[l].dec,
                                        Lsplit[l].z_v,Tsplit[l].tolist(),
                                        rout,nd,h_array]).T
                        
                        pool = Pool(processes=(num))
                        salida = np.array(pool.map(partial_map_unpack, entrada))
                        pool.terminate()
                                
                for j in range(len(salida)):
                        
                        profilesums = salida[j]
                        km          = np.tile(Ksplit[l][j],(3,ndots,1)).T
                                                
                        GTsum += np.tile(profilesums['GTsum'],(101,1,1))*km
                        GXsum += np.tile(profilesums['GXsum'],(101,1,1))*km
                        Ksum  += np.tile(profilesums['Ksum'],(101,1,1))*km
                        
                        Ninbin += np.tile(profilesums['Ninbin'],(101,1,1))*km
                        Ntot         = np.append(Ntot,profilesums['Ntot'])
                        
                
                t2 = time.time()
                ts = (t2-t1)/60.
                tslice = np.append(tslice,ts)
                print('TIME SLICE')
                print(ts)
                print('Estimated ramaining time')
                print(np.mean(tslice)*(len(Lsplit)-(l+1)))
        
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
        
        
        # AVERAGE LENS PARAMETERS
        
        zmean        = np.average(L.z_v,weights=Ntot)
        lM_mean      = np.average(L.lgm,weights=Ntot)
        
        q2d_mean     = np.average(L.q2d,weights=Ntot)
        q2dr_mean    = np.average(L.q2dr,weights=Ntot)
        q3d_mean     = np.average(L.q3d,weights=Ntot)
        q3dr_mean    = np.average(L.q3dr,weights=Ntot)
        s3d_mean     = np.average(L.s3d,weights=Ntot)
        s3dr_mean    = np.average(L.s3dr,weights=Ntot)
        
        if vmice == 2:
            lM_v2_mean = np.average(L.lmhalo,weights=Ntot)
        else:
            lM_v2_mean = lM_mean
            
        
        
        # WRITING OUTPUT FITS FILE
        
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
        
        tbhdu_pro = fits.BinTableHDU.from_columns(fits.ColDefs(table_pro))
        tbhdu_cov = fits.BinTableHDU.from_columns(fits.ColDefs(table_cov))
        
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
        h.append(('lM_mean',np.round(lM_mean,4)))
        h.append(('lM_v2_mean',np.round(lM_v2_mean,4)))
        h.append(('z_mean',np.round(zmean,4)))
        h.append(('q2d_mean',np.round(q2d_mean,4)))
        h.append(('q2dr_mean',np.round(q2dr_mean,4)))
        h.append(('q3d_mean',np.round(q3d_mean,4)))
        h.append(('q3dr_mean',np.round(q3dr_mean,4)))
        h.append(('s3d_mean',np.round(s3d_mean,4)))        
        h.append(('s3dr_mean',np.round(s3dr_mean,4))) 
        
        primary_hdu = fits.PrimaryHDU(header=h)
        
        hdul = fits.HDUList([primary_hdu, tbhdu_pro, tbhdu_cov])
        
        hdul.writeto(folder+'profile_'+sample+'.fits',overwrite=True)
                
        tfin = time.time()
        
        print('TOTAL TIME ',(tfin-tini)/60.)
        


main(sample,rprox,rmin,rmax,lM_min,lM_max,z_min,z_max,q_min,q_max,ROUT,ndots,ncores,hcosmo)
