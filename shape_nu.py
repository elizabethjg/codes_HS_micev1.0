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

def make_plot2(X,Y,color='C0',nbins=20,plt=plt,label='',error = False,lw=1,lt='-',log=False):
    if log:
        Y = 10**Y
    
    x,q50,q25,q75,nada,ymean,ers = binned(X,Y,nbins)
    if error:
        if log:
            ers = ers*np.log10(np.exp(1))/ymean
            plt.plot(x,np.log10(ymean),lt,color=color,label=label,lw=lw)
            plt.fill_between(np.log10(x),np.log10(ymean)+ers,ymean-ers,color=color,alpha=0.2)
        else:
            plt.plot(x,ymean,lt,color=color,label=label,lw=lw)
            plt.fill_between(x,ymean+ers,ymean-ers,color=color,alpha=0.2)
    else:
        plt.plot(x,q50,lt,color=color,label=label,lw=lw)
        plt.fill_between(x,q75,q25,color=color,alpha=0.2)


# folder = '/home/eli/Documentos/Astronomia/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles2/'
folder = '/home/elizabeth/Documentos/proyectos/HALO-SHAPE/MICE/HS-lensing/profiles3/'
halos = fits.open('../HALO_Props_MICE.fits')[1].data        
mfit_NFW = (halos.cNFW_rho > 1.)*(halos.cNFW_S > 1.)*(halos.cNFW_S_E > 1.)*(halos.cNFW_rho_E > 1.)*(halos.lgMNFW_rho > 12)*(halos.lgMNFW_S > 12)*(halos.lgMNFW_S_E > 12)*(halos.lgMNFW_rho_E > 12)*(halos.cNFW_rho < 15)*(halos.cNFW_S < 15)*(halos.cNFW_S_E < 15)*(halos.cNFW_rho_E < 15)
mfit_Ein = (halos.cEin_rho > 1.)*(halos.cEin_S > 1.)*(halos.cEin_S_E > 1.)*(halos.cEin_rho_E > 1.)*(halos.lgMEin_rho > 12)*(halos.lgMEin_S > 12)*(halos.lgMEin_S_E > 12)*(halos.lgMEin_rho_E > 12)*(halos.cEin_rho < 15)*(halos.cEin_S < 15)*(halos.cEin_S_E < 15)*(halos.cEin_rho_E < 15)*(halos.alpha_rho > 0.)*(halos.alpha_rho_E > 0.)

mfit = mfit_NFW*mfit_Ein#*(~mraro)
mhalos = mfit#*(halos.z <= 0.5)

halos = halos[mhalos]

Eratio = (2.*halos.K/abs(halos.U))

mrelax = (halos.offset < 0.1)*(Eratio < 1.35)


nu200 = peaks.peakHeight(10**halos.lgMNFW_rho, halos.z)
nu200_E = peaks.peakHeight(10**halos.lgMNFW_rho_E, halos.z)
nu200e = peaks.peakHeight(10**halos.lgMEin_rho, halos.z)
nu200e_E = peaks.peakHeight(10**halos.lgMEin_rho_E, halos.z)
nufof = peaks.peakHeight(10**halos.lgM, halos.z)

nu = nufof
nu = nu200
nu = nu200_E

mzL = (halos.z > 0.1)*(halos.z < 0.2)
mzH = (halos.z > 0.2)*(halos.z < 0.4)

a = -0.257
b = -0.219

f, ax = plt.subplots(2,1, figsize=(5,7),sharey=True,sharex=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[1].plot([0,0],[0,0],'gold',lw=3,label='Bonamigo et al. (2015)')
ax[1].plot([0,0],[0,0],'orangered',lw=3,label='NFW')
ax[1].plot([0,0],[0,0],'seagreen',lw=3,label='Einasto')
ax[1].plot([0,0],[0,0],'C9',lw=3,label='FOF')
ax[1].plot([0,0],[0,0],'k',lw=3,label='spherical')
ax[1].plot([0,0],[0,0],'k--',lw=3,label='ellipsoidal')
ax[1].legend(loc=3,frameon=False,ncol=2)

make_plot2(np.log10(nufof),np.log10(halos.s),'C9',15,label=r'$M_{FOF}$',lw=3,error=True,plt=ax[0],log=True)
make_plot2(np.log10(nu200),np.log10(halos.s),'orangered',15,label=r'$M^S_{200}$',lw=3,error=True,plt=ax[0],log=True)
make_plot2(np.log10(nu200_E),np.log10(halos.s),'orangered',15,label=r'$M^E_{200}$',lw=3,lt='--',error=True,plt=ax[0],log=True)
make_plot2(np.log10(nu200e),np.log10(halos.s),'seagreen',15,label=r'$M^S_{200}$',lw=3,error=True,plt=ax[0],log=True)
make_plot2(np.log10(nu200e_E),np.log10(halos.s),'seagreen',15,label=r'$M^E_{200}$',lw=3,lt='--',error=True,plt=ax[0],log=True)

make_plot2(np.log10(nufof[mrelax]),np.log10(halos.s[mrelax]),'C9',15,label=r'$M_{FOF}$',lw=3,error=True,plt=ax[1],log=True)
make_plot2(np.log10(nu200[mrelax]),np.log10(halos.s[mrelax]),'orangered',15,label=r'$M^S_{200}$',lw=3,error=True,plt=ax[1],log=True)
make_plot2(np.log10(nu200_E[mrelax]),np.log10(halos.s[mrelax]),'orangered',15,label=r'$M^E_{200}$',lw=3,lt='--',error=True,plt=ax[1],log=True)
make_plot2(np.log10(nu200e[mrelax]),np.log10(halos.s[mrelax]),'seagreen',15,label=r'$M^S_{200}$',lw=3,error=True,plt=ax[1],log=True)
make_plot2(np.log10(nu200e_E[mrelax]),np.log10(halos.s[mrelax]),'seagreen',15,label=r'$M^E_{200}$',lw=3,lt='--',error=True,plt=ax[1],log=True)

ax[0].plot(np.log10(nu200),a*np.log10(nu200)+b,'gold',lw=3,label='Bonamigo et al. (2015)')
ax[1].plot(np.log10(nu200),a*np.log10(nu200)+b,'gold',lw=3,label='Bonamigo et al. (2015)')
ax[1].set_xlabel(r'$\log \nu(M,z)$')
ax[0].set_ylabel(r'$\log(S)$')
ax[1].set_ylabel(r'$\log(S)$')

ax[0].text(0.45,-0.25,'all halos')
ax[1].text(0.45,-0.25,'only relaxed halos')
plt.axis([0.15,0.7,-0.55,-0.22])
f.savefig('../nu_shape.pdf',bbox_inches='tight')

plt.figure()
make_plot2(np.log10(nufof),np.log10(halos.q2d),'C3',15,label=r'$M_{FOF}$')
make_plot2(np.log10(nu200),np.log10(halos.q2d),'C2',15,label=r'$M^S_{200}$')
make_plot2(np.log10(nu200_E),np.log10(halos.q2d),'C1',15,label=r'$M^E_{200}$')
plt.xlabel(r'$\log \nu(M_{200},z)$')
plt.ylabel(r'$\log(q = b/a)$')
plt.axis([0.15,0.65,-0.4,-0.1])
plt.legend(loc=3)

plt.figure()
make_plot2(halos.lgM+3.*np.log10(1+halos.z),halos.s,'C3',15,label=r'$M_{FOF}$')
make_plot2(halos.lgMNFW_rho+3.*np.log10(1+halos.z),halos.s,'C2',15,label=r'$M^S_{200}$')
make_plot2(halos.lgMNFW_rho_E+3.*np.log10(1+halos.z),halos.s,'C1',15,label=r'$M^E_{200}$')
plt.xlabel(r'$\log (M_{200},z)$')
plt.ylabel(r'$(s = b/a)$')
# plt.axis([0.15,0.65,-0.4,-0.1])
plt.legend(loc=3)



'''
plt.figure()
plt.plot(np.log10(nu),a*np.log10(nu)+b,'C4',lw=3,label='Bonamigo et al. 2015')
make_plot2(np.log10(nu),np.log10(halos.s),'C2',15)
# make_plot2(np.log10(nu[mzL]),np.log10(halos.s[mzL]),'C0',15)
# make_plot2(np.log10(nu[mzH]),np.log10(halos.s[mzH]),'C3',15)
plt.xlabel(r'$\log \nu(M_{200},z)$')
plt.ylabel(r'$\log(s = c/a)$')

plt.figure()
make_plot2(np.log10(nu),np.log10(halos.q2d),'C2',15)
# make_plot2(np.log10(nu[mzL]),np.log10(halos.q2d[mzL]),'C0',15)
# make_plot2(np.log10(nu[mzH]),np.log10(halos.q2d[mzH]),'C3',15)
plt.xlabel(r'$\log \nu(M_{200},z)$')
plt.ylabel(r'$\log(q = b/a)$')
'''
