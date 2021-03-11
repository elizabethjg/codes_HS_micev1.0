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
parser.add_argument('-rmax', action='store', dest='rmax',default=1000.)
parser.add_argument('-lM_min', action='store', dest='lM_min', default=14.4)
parser.add_argument('-lM_max', action='store', dest='lM_max', default=14.6)
parser.add_argument('-z_min', action='store', dest='z_min', default=0.1)
parser.add_argument('-z_max', action='store', dest='z_max', default=0.4)
parser.add_argument('-q_min', action='store', dest='q_min', default=0.)
parser.add_argument('-q_max', action='store', dest='q_max', default=1.)
parser.add_argument('-ROUT', action='store', dest='ROUT', default=5000.)
parser.add_argument('-nbins', action='store', dest='nbins', default=5000)
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
lM_min=13.6
lM_max=13.65
z_min = 0.1
z_max = 0.11
q_min = 0.
q_max = 1.0
RIN = 400.
ROUT = 5000.
ndots= 5000
ncores = 40
hcosmo = 1.0 
vmice = 2
'''


folder = '/mnt/projects/lensing/HALO_SHAPE/MICEv'+str(vmice)+'.0/catalogs/'
S      = fits.open(folder+'catalogs/MICE_sources.fits')[1].data


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
        mlenses = mmass*mz*mregion*mq*mrcut
        Nlenses = mlenses.sum()
        
        print('Nlenses',Nlenses)
        
        L = L[mlenses]
        
        for j in range(mlenses.sum()):
                
                idhalo = L.unique_halo_id_raw[j]       
                Nlenses = 1
        
                if Nlenses < ncores:
                        ncores = Nlenses
                
                print('Nlenses',Nlenses)
                print('LENS ',j)
                
                
                theta  = np.array([0.,np.arctan(L.a2dy[j]/L.a2dx[j]),np.arctan(L.a2dry[j]/L.a2drx[j])]).T
        
                
                entrada = [L.ra[j],L.dec[j],L.z_v[j],theta,ROUT,ndots,hcosmo]
                                               
                
                profilesums = partial_map_unpack(entrada)
                                        
                GTsum = profilesums['GTsum']
                GXsum = profilesums['GXsum']
                Ksum  = profilesums['Ksum']
                
                Ninbin = profilesums['Ninbin']
                                
                                
                # COMPUTING PROFILE        
                        
                GT  = (GTsum/Ninbin)
                GX  = (GXsum/Ninbin)
                K   = (Ksum/Ninbin)
                        
                
                
                # AVERAGE LENS PARAMETERS
                
                zmean        = L.z_v[j]
                lM_mean      = L.lgm[j]
                
                q2d_mean     = L.q2d[j]
                q2dr_mean    = L.q2dr[j]
                q3d_mean     = L.q3d[j]
                q3dr_mean    = L.q3dr[j]
                s3d_mean     = L.s3d[j]
                s3dr_mean    = L.s3dr[j]
                
                if vmice == 2:
                        lM_v2_mean = L.lmhalo[j]
                else:
                        lM_v2_mean = lM_mean
                
                
                
                # WRITING OUTPUT FITS FILE
                
                table_pro = [fits.Column(name='xmpc', format='E', array=xbin),
                        fits.Column(name='ympc', format='E', array=ybin),
                        fits.Column(name='GT_control', format='E', array=GT[:,0]),
                        fits.Column(name='GT', format='E', array=GT[:,1]),
                        fits.Column(name='GT_reduced', format='E', array=GT[:,2]),
                        fits.Column(name='GX_control', format='E', array=GX[:,0]),
                        fits.Column(name='GX', format='E', array=GX[:,1]),
                        fits.Column(name='GX_reduced', format='E', array=GX[:,2]),
                        fits.Column(name='K_control', format='E', array=K[:,0]),
                        fits.Column(name='K', format='E', array=K[:,1]),
                        fits.Column(name='K_reduced', format='E', array=K[:,2])]
                        
                                        
                tbhdu_pro = fits.BinTableHDU.from_columns(fits.ColDefs(table_pro))

                
                h = fits.Header()
                h.append(('LENS',np.int(idhalo)))
                h.append(('MICE version',vmice))
                h.append(('ndots',ndots))
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
                
                hdul = fits.HDUList([primary_hdu, tbhdu_pro])
                
                hdul.writeto(folder+'mapas/mapa_ind_'+sample+'_'+str(idhalo)+'.fits',overwrite=True)
                
        tfin = time.time()
        
        print('TOTAL TIME ',(tfin-tini)/60.)
        


main(sample,rprox,rmin,rmax,lM_min,lM_max,z_min,z_max,q_min,q_max,ROUT,ndots,ncores,hcosmo)
