import sys
import pylab
from astropy.io import fits
from scipy import stats
sys.path.append('/home/eli/lens_codes_v3.7')
sys.path.append('/home/elizabeth/lens_codes_v3.7')
from colossus.cosmology import cosmology  
from colossus.lss import peaks
params = {'flat': True, 'H0': 70.0, 'Om0': 0.25, 'Ob0': 0.044, 'sigma8': 0.8, 'ns': 0.95}
cosmology.addCosmology('MICE', params)
cosmo = cosmology.setCosmology('MICE')


def q_75(y):
    return np.quantile(y, 0.75)

def q_25(y):
    return np.quantile(y, 0.25)


def binned(x,y,nbins=10):
    
    bined = stats.binned_statistic(x,y,statistic='mean', bins=nbins)
    x_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    ymean     = bined.statistic

    bined = stats.binned_statistic(x,y,statistic='median', bins=nbins)
    x_b = 0.5*(bined.bin_edges[:-1] + bined.bin_edges[1:])
    q50     = bined.statistic
    
    bined = stats.binned_statistic(x,y,statistic=q_25, bins=nbins)
    q25     = bined.statistic

    bined = stats.binned_statistic(x,y,statistic=q_75, bins=nbins)
    q75     = bined.statistic

    bined = stats.binned_statistic(x,y,statistic='count', bins=nbins)
    N     = bined.statistic

    bined = stats.binned_statistic(x,y,statistic='std', bins=nbins)
    sigma = bined.statistic
    
    dig   = np.digitize(x,bined.bin_edges)
    mz    = np.ones(len(x))
    for j in range(nbins):
        mbin = dig == (j+1)
        mz[mbin] = y[mbin] >= q50[j]   
    mz = mz.astype(bool)
    return x_b,q50,q25,q75,mz,ymean,sigma/np.sqrt(N)
            

def make_plot(X,Y,Z,zlim=0.3,nbins=20,plt=plt,error = False):
    plt.figure()
    x,q50,q25,q75,nada,ymean,ers = binned(X,Y,nbins)
    plt.scatter(X,Y, c=Z, alpha=0.3,s=1,vmax=zlim)
    if error:
        plt.plot(x,ymean,'C3')
        plt.plot(x,ymean+ers,'C3--')
        plt.plot(x,ymean-ers,'C3--')    
    else:
        plt.plot(x,q50,'C3')
        plt.plot(x,q25,'C3--')
        plt.plot(x,q75,'C3--')
    plt.colorbar()

def make_plot2(X,Y,color='C0',nbins=20,plt=plt,label='',error = False,lw=1,lt='-'):
    x,q50,q25,q75,nada,ymean,ers = binned(X,Y,nbins)
    if error:
        plt.plot(x,ymean,lt,color=color,label=label,lw=lw)
        plt.fill_between(x,ymean+ers,ymean-ers,color=color,alpha=0.2)
    else:
        plt.plot(x,q50,lt,color=color,label=label,lw=lw)
        plt.fill_between(x,q75,q25,color=color,alpha=0.2)


folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'
halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data        
Eratio = (2.*halos.K/abs(halos.U))

mrelax = (halos.offset < 0.1)*(Eratio < 1.35)
mfit_NFW = (halos.cNFW_rho > 1.)*(halos.cNFW_S > 1.)*(halos.cNFW_rho < 10.)*(halos.cNFW_S < 10.)*(halos.lgMNFW_rho > 12)*(halos.lgMNFW_S > 12)
mfit_Ein = (halos.cEin_rho > 1.)*(halos.cEin_S > 1.)*(halos.cEin_rho < 10.)*(halos.cEin_S < 10.)*(halos.lgMEin_rho > 12)*(halos.lgMEin_S > 12)*(halos.alpha_rho > 0.)*(halos.alpha_S > 0.)*(halos.alpha_rho < 0.7)*(halos.alpha_S < 0.7)
mhalos = mrelax*mfit_Ein*mfit_NFW*(halos.z <= 0.5)

halos = halos[mhalos]

nu200 = peaks.peakHeight(10**halos.lgMNFW_rho, halos.z)
nufof = peaks.peakHeight(10**halos.lgM, halos.z)

nu = nu200

mzL = (halos.z > 0.1)*(halos.z < 0.2)
mzH = (halos.z > 0.2)*(halos.z < 0.4)

plt.figure()
make_plot2(nu,halos.s,'C2',15)
make_plot2(nu[mzL],halos.s[mzL],'C0',15)
make_plot2(nu[mzH],halos.s[mzH],'C3',15)
plt.xlabel(r'$\nu(M_{200},z)$')
plt.ylabel(r'$s = c/a$')

plt.figure()
make_plot2(nu,halos.q2d,'C2',15)
make_plot2(nu[mzL],halos.q2d[mzL],'C0',15)
make_plot2(nu[mzH],halos.q2d[mzH],'C3',15)
plt.xlabel(r'$\nu(M_{200},z)$')
plt.ylabel(r'$q_{2D} = b/a$')
