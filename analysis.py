import sys
sys.path.append('/home/elizabeth/lens_codes_v3.7')
import pylab
from astropy.io import fits
from basic_extract import extract_params
import corner
from matplotlib import cm
from binned_plots import make_plot2

folder = '../profiles3/'

def plot_q_dist():

    halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data        
    Eratio = (2.*halos.K/abs(halos.U))
    
    import seaborn as sns
    
    hsamples = ['Lz','Mz','Hz']
    colors = ['C1','C3','C5']

    fi  = np.arctan(halos.a2Dy/halos.a2Dx)
    fir = np.arctan(halos.a2Dry/halos.a2Drx)

    Dfi = np.rad2deg(fi - fir)

    f, ax = plt.subplots(2,2, figsize=(6,5),sharex = True,sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)
    

    ax[0,0].set_xlim([0,1.15])
    ax[0,0].set_ylim([0,4.53])
    ax[0,0].text(1,4,'HM')
    ax[0,1].text(1,4,'LM')
    ax[0,0].text(0.1,4,'all halos')
    ax[1,0].text(0.1,4,'only relaxed')
    
    for j in range(len(hsamples)):

        ax[0,1].plot([-1]*2,[-1]*2,colors[j],label=hsamples[j])        
        samp = hsamples[j]
        # from individual halos

        # FIRST HM
        
        print('HM_'+samp)

        h = fits.open(folder+'profile_HM_'+samp+'.fits')[0].header
        print(h['N_LENSES'])
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)

        
        print(h['LM_MIN'],h['LM_MAX'])
        print(h['z_min'],h['z_max'])
        print(np.mean(Dfi[mhalos]),np.std(Dfi[mhalos]))
        print(np.mean(Dfi[mhalos*mrelax]),np.std(Dfi[mhalos*mrelax]))
        
        sns.kdeplot(halos.q2d[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,0])
        ax[1,0].plot([np.mean(halos.q2d[mhalos*mrelax])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,0],ls='--')
        ax[1,0].plot([np.mean(halos.q2dr[mhalos*mrelax])]*2,[4.0,4.4],colors[j]+'--')
        sns.kdeplot(halos.q2d[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,0])
        ax[0,0].plot([np.mean(halos.q2d[mhalos])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,0],ls='--')
        ax[0,0].plot([np.mean(halos.q2dr[mhalos])]*2,[4.0,4.4],colors[j]+'--')

        samp = hsamples[j]
        # from individual halos

        # NOW LM
        
        print('LM_'+samp)

        h = fits.open(folder+'profile_LM_'+samp+'.fits')[0].header
        print(h['N_LENSES'])
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)

        print(samp)
        print(h['LM_MIN'],h['LM_MAX'])
        print(h['z_min'],h['z_max'])
        print(np.mean(Dfi[mhalos]),np.std(Dfi[mhalos]))
        print(np.mean(Dfi[mhalos*mrelax]),np.std(Dfi[mhalos*mrelax]))
        
        sns.kdeplot(halos.q2d[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,1])
        ax[1,1].plot([np.mean(halos.q2d[mhalos*mrelax])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1,1],ls='--')
        ax[1,1].plot([np.mean(halos.q2dr[mhalos*mrelax])]*2,[4.0,4.4],colors[j]+'--')
        sns.kdeplot(halos.q2d[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,1])
        ax[0,1].plot([np.mean(halos.q2d[mhalos])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0,1],ls='--')
        ax[0,1].plot([np.mean(halos.q2dr[mhalos])]*2,[4.0,4.4],colors[j]+'--')

    ax[0,1].plot([-1]*2,[-1]*2,'k--',label='standard')
    ax[0,1].plot([-1]*2,[-1]*2,'k-',label='reduced')
    
    ax[0,1].legend(loc=2,frameon=False)
    ax[1,1].set_xlabel('$q$')
    ax[1,0].set_xlabel('$q$')
    ax[0,0].set_ylabel('$P(q)$')
    ax[1,0].set_ylabel('$P(q)$')
    f.savefig(folder+'../final_plots/qdist.pdf',bbox_inches='tight')  
    
    plt.figure()
    mhalos = (halos.lgM >= 13.5)*(halos.lgM < 14.5)*(halos.z > 0.1)*(halos.z < 0.4)
    plt.scatter(halos.q2d[mhalos],halos.q2dr[mhalos],c=halos.offset[mhalos],vmax=0.3,s=1,cmap='inferno')
    cbar = plt.colorbar()
    plt.plot([0.1,1],[0.1,1],'C7--')
    plt.xlabel('$q$')
    plt.ylabel('$q_r$')
    cbar.set_label(r'$r_c/r_{max}$')
    
    plt.savefig(folder+'../final_plots/qrelax.png',bbox_inches='tight')
    
    plt.figure(figsize=(6,2))
    
    make_plot2(halos.q2d[mhalos],halos.offset[mhalos],color='C5',nbins=10,plt=plt,label='standard',error = False,lw=1,lt='-')
    make_plot2(halos.q2dr[mhalos],halos.offset[mhalos],color='C1',nbins=10,plt=plt,label='reduced',error = False,lw=1,lt='--')
    plt.xlabel('$q$')
    plt.ylabel(r'$r_c/r_{max}$')
    plt.legend()
    
    plt.savefig(folder+'../final_plots/qrelax2.pdf',bbox_inches='tight')
 


def plot_q_dist_testq():

    halos = fits.open(folder+'../HALO_Props_MICE.fits')[1].data        
    Eratio = (2.*halos.K/abs(halos.U))
    
    hsamps = ['all', 'all_qmin', 'all_qmed', 'all_qmax']
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps,['350']*4,['1500']*4)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps,['350']*4,['1500']*4,reduced=True)

    
    import seaborn as sns
    
    hsamples = ['all']
    colors = ['k']

    f, ax = plt.subplots(2,1, figsize=(6,5),sharex = True,sharey=True)
    f.subplots_adjust(hspace=0,wspace=0)
    
    ax[0].set_xlim([0,1.15])
    ax[0].set_ylim([0,4.53])
    ax[0].text(0.1,4,'all halos')
    ax[1].text(0.1,4,'only relaxed')
    
    for j in range(len(hsamples)):

        ax[0].plot([-1]*2,[-1]*2,colors[j],label=hsamples[j])        
        samp = hsamples[j]
        # from individual halos

        # FIRST HM

        h = fits.open(folder+'profile_'+samp+'.fits')[0].header
        print(h['N_LENSES'])
        mhalos = (halos.lgM >= h['LM_MIN'])*(halos.lgM < h['LM_MAX'])*(halos.z >= h['z_min'])*(halos.z < h['z_max'])
        mrelax = (halos.offset < 0.1)*(Eratio < 1.35)

        print(samp)
        print(h['LM_MIN'],h['LM_MAX'])
        print(h['z_min'],h['z_max'])
        
        sns.kdeplot(halos.q2d[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1])
        ax[1].plot([np.mean(halos.q2d[mhalos*mrelax])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos*mrelax].byteswap().newbyteorder(),color=colors[j],ax=ax[1],ls='--')
        ax[1].plot([np.mean(halos.q2dr[mhalos*mrelax])]*2,[4.0,4.4],colors[j]+'--')
        sns.kdeplot(halos.q2d[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0])
        ax[0].plot([np.mean(halos.q2d[mhalos])]*2,[4.0,4.4],colors[j])
        sns.kdeplot(halos.q2dr[mhalos].byteswap().newbyteorder(),color=colors[j],ax=ax[0],ls='--')
        ax[0].plot([np.mean(halos.q2dr[mhalos])]*2,[4.0,4.4],colors[j]+'--')

        # for i in range(len(hsamps)):
            

    ax[0].plot([-1]*2,[-1]*2,'k--',label='standard')
    ax[0].plot([-1]*2,[-1]*2,'k-',label='reduced')
    
    ax[0].legend(loc=2,frameon=False)
    ax[1].set_xlabel('q')
    ax[1].set_xlabel('q')
    f.savefig(folder+'../test_plots/qdist_all_q.png',bbox_inches='tight')   



def plot_bias(hsamps,lhs,cstyle,nplot,RIN,ROUToq,D = 1):
    
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps,RIN,ROUToq,reduced=True)
    ###########
    #   q1h
    ###########
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_r]
    
    param = 1

    f, ax = plt.subplots(2,1, figsize=(6,8),sharex=True)
    f.subplots_adjust(hspace=0)
    
    for hs in range(len(hsamps)):
        
        ax[0].plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        
        ax[0].errorbar(qhr[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,mfc='none')
        ax[0].errorbar(qhr_r[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10,mfc='none')
        ax[0].errorbar(qh[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10)
        ax[0].errorbar(qh_r[hs],NFW[0][hs][param],
                            yerr=np.array([NFW[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10)
                            

        ax[1].errorbar(qhr[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFW_r[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,mfc='none')
        ax[1].errorbar(qhr_r[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFW_r[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10,mfc='none')
        ax[1].errorbar(qh[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFW_r[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10)
        ax[1].errorbar(qh_r[hs],NFW_r[0][hs][param],
                            yerr=np.array([NFW_r[1][hs][param]]).T,alpha=0.5,
                            fmt=cstyle[hs],markersize=10)
    
    ax[0].plot([0.5,0.8],[0.5,0.8],'C7--')
    ax[1].plot([0.5,0.8],[0.5,0.8],'C7--')
    ax[0].set_xlim([0.5,0.82])
    ax[0].set_ylim([0.57,0.67])
    ax[1].set_ylim([0.52,0.65])
    # ax[0].set_xlim([0.1,0.9])
    # ax[0].set_ylim([0.1,0.9])
    # ax[1].set_ylim([0.1,0.9])
    ax[1].plot([-1,-1],[-1,-1],'k^',label='standard - all')
    ax[1].plot([-1,-1],[-1,-1],'k^',label='standard - relaxed',mfc='none')
    ax[1].plot([-1,-1],[-1,-1],'k^',label='reduced - all',alpha=0.5)
    ax[1].plot([-1,-1],[-1,-1],'k^',label='reduced - relaxed',mfc='none',alpha = 0.5)

    ax[0].set_ylabel(r'$\tilde{q_{1h}}(\hat{\phi})$')
    ax[1].set_ylabel(r'$\tilde{q_{1h}}(\hat{\phi}_r)$')
    ax[0].set_xlabel(r'$\langle q \rangle$')
    ax[1].set_xlabel(r'$\langle q \rangle$')

    ax[0].legend(frameon=False,loc=4,ncol=3)
    ax[1].legend(loc=1,ncol=2)
    # f.savefig(folder+'../final_plots/qcomparison.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/qcomparison_'+nplot+'.png',bbox_inches='tight')
    # f.savefig(folder+'../final_plots/qcomparison_rel.png',bbox_inches='tight')
    
    
    f, ax = plt.subplots(1,1, figsize=(6,4))

    for hs in range(len(hsamps)):
        
        ax.plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        param = -1
        ax.errorbar(NFW[0][hs][param],NFW_r[0][hs][param],
                            xerr=np.array([NFW[1][hs][param]]).T,
                            yerr=np.array([NFW_r[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10,mfc='none')
        param = 1
        ax.errorbar(NFW[0][hs][param],NFW_r[0][hs][param],
                            xerr=np.array([NFW[1][hs][param]]).T,
                            yerr=np.array([NFW_r[1][hs][param]]).T,
                            fmt=cstyle[hs],markersize=10)
    
    ax.plot([0.1,1.2],[0.1,1.2],'C7--')
    ax.set_xlim([0.15,1.1])
    ax.set_ylim([0.15,1.1])
    # ax.set_xlim([0.25,0.7])
    # ax.set_ylim([0.25,0.72])
    ax.plot([-1,-1],[-1,-1],'k^',label=r'$\tilde{q}_{1h}$')
    ax.plot([-1,-1],[-1,-1],'k^',label=r'$\tilde{q}_{2h}$',mfc='none')
    ax.legend(loc=2,ncol = 4)
    ax.set_xlabel(r'$\tilde{q}(\hat{\phi})$')
    ax.set_ylabel(r'$\tilde{q}(\hat{\phi}_r)$')
    # f.savefig(folder+'../final_plots/q2h_rel.png',bbox_inches='tight')
    # f.savefig(folder+'../final_plots/q2h.png',bbox_inches='tight')
    # f.savefig(folder+'../final_plots/q2h.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/q2h_'+nplot+'.png',bbox_inches='tight')


    xl = ['NFW - 1h+2h','Ein - 1h+2h','NFW - 1h+2h - fix $c_{200}$','NFW 1h']
    
    # FOM
    f, ax = plt.subplots(2,1, figsize=(14,6),sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    ax = [ax[1],ax[0]]
    
    ax[0].plot([0,5],[0,0],'C7--')
    
    ax[0].axhspan(-0.05,0.05,0,5,color='C7',alpha=0.5)
    ax[0].axhspan(-0.025,0.025,0,5,color='C7',alpha=0.5)

    
    for hs in range(len(hsamps)):
        ax[0].plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        for fp in range(4):
            ax[0].errorbar(fp+1+0.1*hs,(fq[fp][0][hs][param]-qhr[hs])/qhr[hs],
                         yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                         fmt=cstyle[hs],markersize=10,mfc='none')
            # ax[0].errorbar(fp+1+0.1*hs,(fq[fp][0][hs][param]-qh[hs])/qh[hs],
                         # yerr=np.array([fq[fp][1][hs][param]/qh[hs]]).T,
                         # fmt=cstyle[hs],markersize=10)
    
    ax[0].legend(frameon = False,loc=2,ncol=7)
    ax[0].set_ylabel(r'$(\tilde{q}_{1h}(\hat{\phi})-\langle q \rangle)/\langle q \rangle$')
    ax[0].axis([0,5,-0.12,0.12])
    ax[0].set_xticks(np.arange(4)+1)
    ax[0].set_xticklabels(xl)
    # f.savefig(folder+'../final_plots/model_q1h.pdf',bbox_inches='tight')
    # f.savefig(folder+'../test_plots/model_q1hr_'+nplot+'.png',bbox_inches='tight')
    
    # FOM
    param = 1
    ax[1].plot([0,5],[0,0],'C7--')
    
    ax[1].axhspan(-0.05,0.05,0,5,color='C7',alpha=0.5)
    ax[1].axhspan(-0.025,0.025,0,5,color='C7',alpha=0.5)
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            ax[1].errorbar(fp+1+0.1*hs,(fqr[fp][0][hs][param]-qh[hs])/qh[hs],
                          yerr=np.array([fqr[fp][1][hs][param]/qh[hs]]).T,
                          fmt=cstyle[hs],markersize=10)

    ax[1].plot([-1,-1],[-1,-1],'k^',label=r'all halos')
    ax[1].plot([-1,-1],[-1,-1],'k^',label=r'relaxed',mfc='none')
    
    ax[1].legend(frameon = False,loc=1,ncol=3)
    ax[1].set_ylabel(r'$(\tilde{q}_{1h}(\hat{\phi}_r)-\langle q \rangle)/\langle q \rangle$')
    ax[1].axis([0,5,-0.14,0.14])
    ax[1].set_xticks(np.arange(4)+1)
    ax[1].set_xticklabels(xl)
    # f.savefig(folder+'../final_plots/model_q1h.pdf',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_q1h_'+nplot+'.png',bbox_inches='tight')
    
    # FOM - mass
    f, ax = plt.subplots(2,1, figsize=(14,6),sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    ax = [ax[1],ax[0]]
    
    ax[0].plot([0,5],[0,0],'C7--')
    
    ax[0].axhspan(-0.05,0.05,0,5,color='C7',alpha=0.5)
    ax[0].axhspan(-0.025,0.025,0,5,color='C7',alpha=0.5)

    
    for hs in range(len(hsamps)):
        ax[0].plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        for fp in range(4):
            diff = (10**fq[fp][0][hs][0] - 10**NFW_h[hs][0])/10**NFW_h[hs][0]
            ax[0].errorbar(fp+1+0.1*hs,diff,
                         yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                         fmt=cstyle[hs],markersize=10,mfc='none')
    
    ax[0].legend(frameon = False,loc=3,ncol=3)
    ax[0].set_ylabel(r'$(\tilde{M}_{200}-\langle M_{200} \rangle)/\langle M_{200} \rangle$')
    ax[0].set_xticks(np.arange(4)+1)
    ax[0].set_xticklabels(xl)
    # f.savefig(folder+'../final_plots/model_q1h.pdf',bbox_inches='tight')
    # f.savefig(folder+'../test_plots/model_q1hr_'+nplot+'.png',bbox_inches='tight')
    
    # FOM
    param = 1
    ax[1].plot([0,5],[0,0],'C7--')
    
    ax[1].axhspan(-0.05,0.05,0,5,color='C7',alpha=0.5)
    ax[1].axhspan(-0.025,0.025,0,5,color='C7',alpha=0.5)
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            diff = (10**fq[fp][0][hs][0] - 10**NFW_hr[hs][0])/10**NFW_hr[hs][0]
            ax[1].errorbar(fp+1+0.1*hs,diff,
                          yerr=np.array([fqr[fp][1][hs][param]/qh[hs]]).T,
                          fmt=cstyle[hs],markersize=10)

    ax[1].plot([-1,-1],[-1,-1],'k^',label=r'all halos')
    ax[1].plot([-1,-1],[-1,-1],'k^',label=r'relaxed',mfc='none')
    
    ax[1].legend(frameon = False,loc=3,ncol=3)
    ax[1].set_ylabel(r'$(\tilde{M}_{200}-\langle M_{200} \rangle)/\langle M_{200} \rangle$')
    ax[1].axis([0,5,-0.23,0.23])
    ax[1].set_xticks(np.arange(4)+1)
    ax[1].set_xticklabels(xl)
    
    # f.savefig(folder+'../final_plots/model_M200.png',bbox_inches='tight')
    f.savefig(folder+'../test_plots/model_M200_'+nplot+'.png',bbox_inches='tight')
    # f.savefig(folder+'../final_plots/model_q1h_rel.png',bbox_inches='tight')
    # f.savefig(folder+'../test_plots/model_q1h_'+nplot+'.png',bbox_inches='tight')
    
    

def plot_M200q(hsamps,lhs,cstyle,RIN,ROUToq):
    
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps,RIN,ROUToq,reduced=True)
    ###########
    #   q1h
    ###########
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_r]
    
    param = 1
    D = 1.
    
    ##############
    #   lM200/q
    ##############
    mec = ['k','C2','C6','C9']
    xl = ['NFW - 1h+2h','Ein - 1h+2h','NFW - 1h+2h - fix $c_{200}$','NFW 1h']
    
    f, ax = plt.subplots(1,1, figsize=(10,6),sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    param = 1
    ax = [ax]
    
    ax[0].axhspan(-0.05,0.05,0,1,color='C7',alpha=0.3)
    ax[0].axvspan(-0.05,0.05,ymin=0.0,ymax=1.,color='C7',alpha=0.2)
    ax[0].axhline(0,color='C7',ls='--')
    ax[0].axvline(0,color='C7',ls='--')
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            ql = fqr[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fqr[fp][0][hs][0] - 10**Ein_h[hs][0])/10**Ein_h[hs][0]
            else:
                diff = (10**fqr[fp][0][hs][0] - 10**NFW_h[hs][0])/10**NFW_h[hs][0]
            ax[0].errorbar(diff,(ql-qh[hs])/qh[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                        fmt=cstyle[hs],markersize=15,mec=mec[fp])
                        
    
    

    # ---------------------------------
    D = 0.78
    
    hsamps_misall = ['HM_Lz_mis20_miscen',
                     'LM_Lz_mis20_miscen',
                     'HM_Mz_mis20_miscen',
                     'LM_Mz_mis20_miscen',
                     'HM_Hz_mis20_miscen',
                     'LM_Hz_mis20_miscen']
                     
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps_misall,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps_misall,RIN,ROUToq,reduced=True)
    
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_r]

    for hs in range(len(hsamps)):
        ax[0].plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        for fp in range(4):
            ql = fqr[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fqr[fp][0][hs][0] - 10**Ein_h[hs][0])/10**Ein_h[hs][0]
            else:
                diff = (10**fqr[fp][0][hs][0] - 10**NFW_h[hs][0])/10**NFW_h[hs][0]
            ax[0].errorbar(diff,(ql-qh[hs])/qh[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                        fmt=cstyle[hs],markersize=10,mec=mec[fp],alpha = 0.5)
            if hs == len(hsamps)-1:
                ax[0].plot([-1,-1],[-1,-1],'^',label=xl[fp],mfc='none',mec=mec[fp])
        
                            
    
    ax[0].set_ylabel(r'$(\tilde{q}_{1h}(\hat{\phi}_r)-\langle q \rangle)/\langle q \rangle$')
    ax[0].set_xlim([-0.2,0.2])
    ax[0].set_ylim([-0.2,0.2])
                    
    ax[0].plot([-1,-1],[-1,-1],'k^',label='$p_{cc} = 1, \Delta \phi = 0$',markersize=15)
    ax[0].plot([-1,-1],[-1,-1],'k^',label='$p_{cc} = 0.75, \Delta \phi = 20^\circ$',markersize=10,alpha=0.5)


    ax[0].legend(frameon = False,ncol = 2,loc=3)
    ax[0].set_xlabel(r'$(\tilde{M}_{200}-\langle M_{200} \rangle)/\langle M_{200} \rangle$')
    
    
    
    f.savefig(folder+'../final_plots/model_ratioq_M200.pdf',bbox_inches='tight')
    # f.savefig(folder+'../final_plots/model_ratioq_M200.png',bbox_inches='tight')

def plot_M200q_qcut():

    hsamps = ['HM_Lz','LM_Lz','HM_Mz','LM_Mz','HM_Hz','LM_Hz']

    
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps,RIN,ROUToq,reduced=True)
    ###########
    #   q1h
    ###########
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_r]
    
    param = 1
    D = 1.
    
    ##############
    #   lM200/q
    ##############
    mec = ['k','C2','C6','C9']
    xl = ['NFW - 1h+2h','Ein - 1h+2h','NFW - 1h+2h - fix $c_{200}$','NFW 1h']
    
    f, ax = plt.subplots(1,1, figsize=(10,6),sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    param = 1
    ax = [ax]
    
    ax[0].axhspan(-0.05,0.05,0,1,color='C7',alpha=0.3)
    ax[0].axvspan(-0.05,0.05,ymin=0.0,ymax=1.,color='C7',alpha=0.2)
    ax[0].axhline(0,color='C7',ls='--')
    ax[0].axvline(0,color='C7',ls='--')
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            ql = fq[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fq[fp][0][hs][0] - 10**Ein_h[hs][0])/10**Ein_h[hs][0]
            else:
                diff = (10**fq[fp][0][hs][0] - 10**NFW_h[hs][0])/10**NFW_h[hs][0]
            ax[0].errorbar(diff,(ql-qh[hs])/qh[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                        fmt=cstyle[hs],markersize=15,mec=mec[fp],alpha=0.7)
                        
    
    

    # ---------------------------------
    
    hsamps_misall = ['HM_Lz_qcut',
                     'LM_Lz_qcut',
                     'HM_Mz_qcut',
                     'LM_Mz_qcut',
                     'HM_Hz_qcut',
                     'LM_Hz_qcut']
                     
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps_misall,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps_misall,RIN,ROUToq,reduced=True)
    
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_r]

    for hs in range(len(hsamps)):
        ax[0].plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        for fp in range(4):
            ql = fq[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fqr[fp][0][hs][0] - 10**Ein_h[hs][0])/10**Ein_h[hs][0]
            else:
                diff = (10**fqr[fp][0][hs][0] - 10**NFW_h[hs][0])/10**NFW_h[hs][0]
            ax[0].errorbar(diff,(ql-qh[hs])/qh[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qh[hs]]).T,
                        fmt=cstyle[hs],markersize=10,mec=mec[fp])
            if hs == len(hsamps)-1:
                ax[0].plot([-1,-1],[-1,-1],'^',label=xl[fp],mfc='none',mec=mec[fp])
        
                            
    
    ax[0].set_ylabel(r'$(\tilde{q}_{1h}(\hat{\phi})-\langle q \rangle)/\langle q \rangle$')
    ax[0].set_xlim([-0.2,0.2])
    ax[0].set_ylim([-0.2,0.2])
                    
    ax[0].plot([-1,-1],[-1,-1],'k^',label='$q > 0.5$',markersize=10)
    ax[0].plot([-1,-1],[-1,-1],'k^',label='all halos',markersize=15,alpha=0.7)


    ax[0].legend(frameon = False,ncol = 2,loc=3)
    ax[0].set_xlabel(r'$(\tilde{M}_{200}-\langle M_{200} \rangle)/\langle M_{200} \rangle$')
    
    
    
    f.savefig(folder+'../final_plots/model_ratioq_M200.pdf',bbox_inches='tight')


def plot_M200q_prev(hsamps,lhs,cstyle,RIN,ROUToq):
    
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps,RIN,ROUToq,reduced=True)
    ###########
    #   q1h
    ###########
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_r]
    
    param = 1
    D = 1.
    
    ##############
    #   lM200/q
    ##############
    mec = ['k','C2','C6','C9']
    xl = ['NFW - 1h+2h','Ein - 1h+2h','NFW - 1h+2h - fix $c_{200}$','NFW 1h']
    
    f, ax = plt.subplots(2,1, figsize=(14,6),sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    param = 1
    
    ax[0].axhspan(-0.05,0.05,0,1,color='C7',alpha=0.3)
    ax[0].axvspan(-0.05,0.05,ymin=0.0,ymax=1.,color='C7',alpha=0.2)
    ax[0].axhline(0,color='C7',ls='--')
    ax[0].axvline(0,color='C7',ls='--')
    ax[1].axhspan(-0.05,0.05,0,1,color='C7',alpha=0.3)
    ax[1].axvspan(-0.05,0.05,ymin=0.0,ymax=1.,color='C7',alpha=0.2)
    ax[1].axhline(0,color='C7',ls='--')
    ax[1].axvline(0,color='C7',ls='--')
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            ql = fq[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fq[fp][0][hs][0] - 10**Ein_hr[hs][0])/10**Ein_hr[hs][0]
            else:
                diff = (10**fq[fp][0][hs][0] - 10**NFW_hr[hs][0])/10**NFW_hr[hs][0]
            ax[0].errorbar(diff,(ql-qhr[hs])/qhr[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                        fmt=cstyle[hs],markersize=15,mec=mec[fp])
                        
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            ql = fqr[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fqr[fp][0][hs][0] - 10**Ein_h[hs][0])/10**Ein_h[hs][0]
            else:
                diff = (10**fqr[fp][0][hs][0] - 10**NFW_h[hs][0])/10**NFW_h[hs][0]
            ax[1].errorbar(diff,(ql-qh[hs])/qh[hs],
                        yerr=np.array([fqr[fp][1][hs][param]/qh[hs]]).T,
                        fmt=cstyle[hs],markersize=15,mec=mec[fp])

    

    # ---------------------------------
    D = 0.78
    hsamps_misall = ['HM_Lz_mis20_miscen','LM_Lz_mis20_miscen','HM_Mz_mis20_miscen','LM_Mz_mis20_miscen','HM_Hz_mis20_miscen','LM_Hz_mis20_miscen']
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps_misall,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps_misall,RIN,ROUToq,reduced=True)
    
    fq = [NFW,Ein,woc,o1h]
    fqr = [NFW_r,Ein_r,woc_r,o1h_r]

    for hs in range(len(hsamps)):
        ax[0].plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        for fp in range(4):
            ql = fq[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fq[fp][0][hs][0] - 10**Ein_hr[hs][0])/10**Ein_hr[hs][0]
            else:
                diff = (10**fq[fp][0][hs][0] - 10**NFW_hr[hs][0])/10**NFW_hr[hs][0]
            ax[0].errorbar(diff,(ql-qhr[hs])/qhr[hs],
                        yerr=np.array([fq[fp][1][hs][param]/qhr[hs]]).T,
                        fmt=cstyle[hs],markersize=10,mec=mec[fp],alpha = 0.5)
                        
    
    ax[0].plot([-1,-1],[-1,-1],'k^',label='$p_{cc} = 1, \Delta \phi = 0$',markersize=15)
    ax[0].plot([-1,-1],[-1,-1],'k^',label='$p_{cc} = 0.75, \Delta \phi = 20^\circ$',markersize=10,alpha=0.5)
    
    ax[0].legend(frameon = False,ncol = 4,loc=3)
    
    ax[0].set_ylabel(r'$(\tilde{q}_{1h}(\hat{\phi})-\langle q \rangle)/\langle q \rangle$')
    ax[0].set_xlim([-0.15,0.1])
    ax[0].set_ylim([-0.2,0.2])
    
    for hs in range(len(hsamps)):
        for fp in range(4):
            ql = fqr[fp][0][hs][param]
            el = ((1. - ql)/(1. + ql))/D
            ql = ((1. - el)/(1. + el))
            if fp == 1:
                diff = (10**fqr[fp][0][hs][0] - 10**Ein_h[hs][0])/10**Ein_h[hs][0]
            else:
                diff = (10**fqr[fp][0][hs][0] - 10**NFW_h[hs][0])/10**NFW_h[hs][0]
            ax[1].errorbar(diff,(ql-qh[hs])/qh[hs],
                        yerr=np.array([fqr[fp][1][hs][param]/qh[hs]]).T,
                        fmt=cstyle[hs],markersize=10,mec=mec[fp],alpha = 0.5)
            if hs == 0:
                ax[1].plot([-1,-1],[-1,-1],'^',label=xl[fp],mfc='none',mec=mec[fp])

    
    ax[1].legend(frameon = False,ncol=2,loc=3)
    ax[1].set_ylabel(r'$(\tilde{q}_{1h}(\hat{\phi}_r)-\langle q \rangle)/\langle q \rangle$')
    ax[1].set_xlabel(r'$(\tilde{M}_{200}-\langle M_{200} \rangle)/\langle M_{200} \rangle$')
    
    
    f.savefig(folder+'../final_plots/model_ratioq_M200.pdf',bbox_inches='tight')

def plt_profile_fitted_final(samp,RIN,ROUT,axx3,fittype='_2h_2q'):


    matplotlib.rcParams.update({'font.size': 11})    
    ax,ax1,ax2 = axx3
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
        
    ndots = p.shape[0]
    
    GT  = p.GAMMA_Tcos
    GTr = p.GAMMA_Tcos_reduced
    GTc = p.GAMMA_Tcos_control
    
    GX  = p.GAMMA_Xsin
    GXr = p.GAMMA_Xsin_reduced
    GXc = p.GAMMA_Xsin_control
    
    # '''
    CovDS  = np.sqrt(cov.COV_ST.reshape(len(GT),len(GT)))
    CovS   = np.sqrt(cov.COV_S.reshape(len(GT),len(GT)))
             
    CovGT  = np.sqrt(cov.COV_GT.reshape(len(GT),len(GT)))
    CovGTr = np.sqrt(cov.COV_GT_reduced.reshape(len(GT),len(GT)))
    CovGTc = np.sqrt(cov.COV_GT_control.reshape(len(GT),len(GT)))
             
    CovGX  = np.sqrt(cov.COV_GX.reshape(len(GT),len(GT)))
    CovGXr = np.sqrt(cov.COV_GX_reduced.reshape(len(GT),len(GT)))
    CovGXc = np.sqrt(cov.COV_GX_control.reshape(len(GT),len(GT)))

    rplot = np.round(p.Rp,2)
        
    # MCMC results
          
    rplot,DS1h,DS2h,gt1h,gx1h,gt1hr,gx1hr,gt2h,gx2h,gt2hr,gx2hr = np.loadtxt(folder+'fitprofile'+fittype+'_'+samp+'_'+str(int(RIN))+'_'+str(int(ROUT))+'.cat')

    DS  = DS1h  + DS2h
    gt  = gt1h  + gt2h
    gtr = gt1hr + gt2hr
    gx  = gx1h  + gx2h
    gxr = gx1hr + gx2hr
        
    ##############    

    ax.plot(p.Rp,p.DSigma_T,'C7')
    ax.plot(rplot,DS1h,'C1',label='1h')
    ax.plot(rplot,DS2h,'C8',label='2h')
    ax.plot(rplot,DS,'C3',label='1h+2h')
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C7',alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Delta\Sigma [h M_\odot/pc^2]$',labelpad=2)
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(0.5,200)
    ax.set_xlim(0.1,10)
    ax.xaxis.set_ticks([0.1,1,5,7])
    ax.set_xticklabels([0.1,1,5,7])
    ax.yaxis.set_ticks([1,10,100])
    ax.set_yticklabels([1,10,100])
    ax.axvline(RIN/1000.,color='k',ls=':')
    ax.axvline(ROUT/1000.,color='k',ls=':')
    # ax.legend(loc=3,frameon=False,ncol=2)
    
    
    ax1.plot(p.Rp,GT,'C7',label='standard')
    ax1.plot(p.Rp,GTr,'C6--',label='reduced')
    ax1.plot(rplot,gt,'C3')
    ax1.plot(rplot,gtr,'C3--')
    ax1.plot(rplot,gt1h,'C1')
    ax1.plot(rplot,gt2h,'C8')
    ax1.plot(rplot,gt1hr,'C1--')
    ax1.plot(rplot,gt2hr,'C8--')
    # ax1.legend(loc=3,frameon=False)
    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C7',alpha=0.4)
    ax1.fill_between(p.Rp,GTr+np.diag(CovGTr),GTr-np.diag(CovGTr),color='C7',alpha=0.4)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [$h^{-1}$ Mpc]')
    ax1.set_ylabel(r'$\Gamma_T [h M_\odot/pc^2]$',labelpad=1.2)
    ax1.set_ylim(0.5,100)
    ax1.set_xlim(0.1,10)
    ax1.xaxis.set_ticks([0.1,1,5,7])
    ax1.set_xticklabels([0.1,1,5,7])
    ax1.yaxis.set_ticks([1,10,100])
    ax1.set_yticklabels([1,10,100])
    ax1.axvline(RIN/1000.,color='k',ls=':')
    ax1.axvline(ROUT/1000.,color='k',ls=':')
        
    ax2.plot([0,10],[0,0],'k')
    ax2.plot(p.Rp,GX,'C7')
    ax2.plot(p.Rp,GXr,'C5--')
    ax2.plot(rplot,gx,'C3')
    ax2.plot(rplot,gxr,'C3--')
    ax2.plot(rplot,gx1h,'C1')
    ax2.plot(rplot,gx2h,'C8')
    ax2.plot(rplot,gx1hr,'C1--')
    ax2.plot(rplot,gx2hr,'C8--')
    ax2.axvline(RIN/1000.,color='k',ls=':')
    ax2.axvline(ROUT/1000.,color='k',ls=':')
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C7',alpha=0.4)
    ax2.fill_between(p.Rp,GXr+np.diag(CovGXr),GXr-np.diag(CovGXr),color='C6',alpha=0.4)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times [h M_\odot/pc^2]$',labelpad=1.2)
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-16,17)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])


def plt_profile_fitted_withnoise(samp,RIN,ROUT,axx3):


    matplotlib.rcParams.update({'font.size': 11})    
    ax,ax1,ax2 = axx3
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    p   = profile[1].data
    cov = profile[2].data
        
    ndots = p.shape[0]
    
    GT  = p.GAMMA_Tcos    
    GX  = p.GAMMA_Xsin
    
    # '''
    CovDS  = np.sqrt(cov.COV_ST.reshape(len(GT),len(GT)))
    CovGT  = np.sqrt(cov.COV_GT.reshape(len(GT),len(GT)))             
    CovGX  = np.sqrt(cov.COV_GX.reshape(len(GT),len(GT)))


    p_name = 'profile_'+samp+'_des.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    # '''
    h   = profile[0].header
    pn   = profile[1].data
    cov = profile[2].data
        
    ndots = pn.shape[0]
    
    GTn  = pn.GAMMA_Tcos    
    GXn  = pn.GAMMA_Xsin
    
    # '''
    CovDSn  = np.sqrt(cov.COV_ST.reshape(len(GTn),len(GTn)))
    CovGTn  = np.sqrt(cov.COV_GT.reshape(len(GTn),len(GTn)))             
    CovGXn  = np.sqrt(cov.COV_GX.reshape(len(GTn),len(GTn)))

                
    ##############    

    ax.plot(p.Rp,p.DSigma_T,'C7')
    ax.fill_between(p.Rp,p.DSigma_T+np.diag(CovDS),p.DSigma_T-np.diag(CovDS),color='C7',alpha=0.4)
    ax.errorbar(pn.Rp,pn.DSigma_T,yerr=np.diag(CovDSn),fmt='C7o',markersize=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'$\Delta\Sigma [h M_\odot/pc^2]$',labelpad=2)
    ax.set_xlabel('r [$h^{-1}$ Mpc]')
    ax.set_ylim(0.5,200)
    ax.set_xlim(0.1,10)
    ax.xaxis.set_ticks([0.1,1,5,7])
    ax.set_xticklabels([0.1,1,5,7])
    ax.yaxis.set_ticks([1,10,100])
    ax.set_yticklabels([1,10,100])
    ax.axvline(RIN/1000.,color='k',ls=':')
    ax.axvline(ROUT/1000.,color='k',ls=':')
    # ax.legend(loc=3,frameon=False,ncol=2)
    
    
    ax1.plot(p.Rp,GT,'C7',label='standard')
    ax1.fill_between(p.Rp,GT+np.diag(CovGT),GT-np.diag(CovGT),color='C7',alpha=0.4)
    ax1.errorbar(pn.Rp,GTn,yerr=np.diag(CovGTn),fmt='C7o',markersize=2)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('r [$h^{-1}$ Mpc]')
    ax1.set_ylabel(r'$\Gamma_T [h M_\odot/pc^2]$',labelpad=1.2)
    ax1.set_ylim(0.5,100)
    ax1.set_xlim(0.1,10)
    ax1.xaxis.set_ticks([0.1,1,5,7])
    ax1.set_xticklabels([0.1,1,5,7])
    ax1.yaxis.set_ticks([1,10,100])
    ax1.set_yticklabels([1,10,100])
    ax1.axvline(RIN/1000.,color='k',ls=':')
    ax1.axvline(ROUT/1000.,color='k',ls=':')
        
    ax2.plot([0,10],[0,0],'k')
    ax2.plot(p.Rp,GX,'C7')
    ax2.errorbar(pn.Rp,GXn,yerr=np.diag(CovGXn),fmt='C7o',markersize=2)
    ax2.axvline(RIN/1000.,color='k',ls=':')
    ax2.axvline(ROUT/1000.,color='k',ls=':')
    ax2.fill_between(p.Rp,GX+np.diag(CovGX),GX-np.diag(CovGX),color='C7',alpha=0.4)
    ax2.set_xlabel('r [$h^{-1}$ Mpc]')
    ax2.set_ylabel(r'$\Gamma_\times [h M_\odot/pc^2]$',labelpad=1.2)
    ax2.set_xscale('log')
    ax2.set_xlim(0.1,10)
    ax2.set_ylim(-16,17)
    ax2.xaxis.set_ticks([0.1,1,5,7])
    ax2.set_xticklabels([0.1,1,5,7])


def misal():
    
    def g(x,disp): 
        return (1./(np.sqrt(2*np.pi)*disp))*np.exp(-0.5*(x/disp)**2)



    from scipy import integrate

    lhs = ['HM-Lz','LM-Lz','HM-Mz','LM-Mz','HM-Hz','LM-Hz']
    hsamps = ['HM_Lz','LM_Lz','HM_Mz','LM_Mz','HM_Hz','LM_Hz']
    ROUToq = ['2000','1000']*3
    RIN = ['350','350','350','350','450','400']

    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps,RIN,ROUToq)
    qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFW_r,Ein_r,o1h_r,woc_r,zh = extract_params(hsamps,RIN,ROUToq,reduced=True)

    e1h = (1. - np.array(NFW[0]).T[1])/(1. + np.array(NFW[0]).T[1])
    e2h = (1. - np.array(NFW[0]).T[-1])/(1. + np.array(NFW[0]).T[-1])
    e1h_r = (1. - np.array(NFW_r[0]).T[1])/(1. + np.array(NFW_r[0]).T[1])
    e2h_r = (1. - np.array(NFW_r[0]).T[-1])/(1. + np.array(NFW_r[0]).T[-1])
    
    edm    = (1. - np.array(qh))/(1. + np.array(qh))
    edm_rel = (1. - np.array(qhr))/(1. + np.array(qhr))


    hsamps_mis10  = ['HM_Lz_mis10','LM_Lz_mis10','HM_Mz_mis10','LM_Mz_mis10','HM_Hz_mis10','LM_Hz_mis10']
    hsamps_mis20  = ['HM_Lz_mis20','LM_Lz_mis20','HM_Mz_mis20','LM_Mz_mis20','HM_Hz_mis20','LM_Hz_mis20']
    hsamps_mis30  = ['HM_Lz_mis30','LM_Lz_mis30','HM_Mz_mis30','LM_Mz_mis30','HM_Hz_mis30','LM_Hz_mis30']
    
    hsamps_miscen = ['HM_Lz_miscen',
                     'LM_Lz_miscen',
                     'HM_Mz_miscen',
                     'LM_Mz_miscen',
                     'HM_Hz_miscen',
                     'LM_Hz_miscen']
    
    hsamps_misall  = ['HM_Lz_mis20_miscen',
                      'LM_Lz_mis20_miscen',
                      'HM_Mz_mis20_miscen',
                      'LM_Mz_mis20_miscen',
                      'HM_Hz_mis20_miscen',
                      'LM_Hz_mis20_miscen']

    sampall = [hsamps_mis10,
               hsamps_mis20,
               hsamps_mis30,
               hsamps_miscen,
               hsamps_misall]

    mod  = ['NFW','Ein','NFW_fixc','NFW_1h']
    sav  = [[],[],[],[]]
    Dil  = []

    for j in range(5):
        
        argumento = lambda x: g(x,(j+1)*10.)*np.cos(2*np.deg2rad(x))
        D         = integrate.quad(argumento, -90., 90.)[0]
        print('EXPECTED DILUTION ',D)
        Dil += [D]

        qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFWmis,Einmis,o1hmis,wocmis,zh = extract_params(sampall[j],RIN,ROUToq)
        qh_r,NFW_h,Ein_h,qhr_r,NFW_hr,Ein_hr,NFWmis_r,Einmis_r,o1hmis_r,wocmis_r,zh = extract_params(sampall[j],RIN,ROUToq,reduced=True)
        
        fm   = [NFWmis,Einmis,wocmis,o1hmis]
        fm_r = [NFWmis_r,Einmis_r,wocmis_r,o1hmis_r]
        
        
        for m in range(len(mod)-1):
        
            e1hmis   = (1. - np.array(fm[m][0]).T[1])/(1. + np.array(fm[m][0]).T[1])
            e2hmis   = (1. - np.array(fm[m][0]).T[-1])/(1. + np.array(fm[m][0]).T[-1])
            e1hmis_r = (1. - np.array(fm_r[m][0]).T[1])/(1. + np.array(fm_r[m][0]).T[1])
            e2hmis_r = (1. - np.array(fm_r[m][0]).T[-1])/(1. + np.array(fm_r[m][0]).T[-1])
            
            sav[m] += [e1hmis/e1h,e2hmis/e2h,e1hmis/edm_rel,
                      e1hmis_r/e1h_r,e2hmis_r/e2h_r,e1hmis_r/edm]
            
            print('EXPECTED DILUTION ',D)
            print(mod[m])
            print('dm ratio')
            print(e1hmis/edm_rel)
            print(e1hmis_r/edm)
    
            f=open(folder+'../misal_res_'+mod[m]+'.tab','a')
            f.write(' ----- mis '+str((j+1)*10)+'\n') 
            
            for i in range(len(lhs)):
                f.write(lhs[i]+' & ')
                f.write('$'+str('%.2f' % (e1hmis/e1h)[i])+'$ & ')
                f.write('$'+str('%.2f' % (e2hmis/e2h)[i])+'$ & ')
                f.write('$'+str('%.2f' % (e1hmis/edm_rel)[i])+'$ & ')
                f.write('$'+str('%.2f' % (e1hmis_r/e1h_r)[i])+'$ & ')
                f.write('$'+str('%.2f' % (e2hmis_r/e2h_r)[i])+'$ & ')
                f.write('$'+str('%.2f' % (e1hmis_r/edm)[i])+r'$ \\ '+'\n') 

            f.close()
            
        m = m + 1
        
        e1hmis   = (1. - np.array(fm[m][0]).T[1])/(1. + np.array(fm[m][0]).T[1])
        e1hmis_r = (1. - np.array(fm_r[m][0]).T[1])/(1. + np.array(fm_r[m][0]).T[1])
        
        print('EXPECTED DILUTION ',D)
        print('dm ratio')
        print(e1hmis/edm_rel)
        print(e1hmis_r/edm)
        
        f=open(folder+'../misal_res_'+mod[m]+'.tab','a')
        
        for i in range(len(lhs)):
            f.write(lhs[i]+' & ')
            f.write('$'+str('%.2f' % (e1hmis/e1h)[i])+'$ & ')
            f.write(' -  & ')
            f.write('$'+str('%.2f' % (e1hmis/edm_rel)[i])+'$ & ')
            f.write('$'+str('%.2f' % (e1hmis_r/e1h_r)[i])+'$ & ')
            f.write(' - ')
            f.write('$'+str('%.2f' % (e1hmis_r/edm)[i])+r'$ \\ '+'\n') 
            
        f.close()
        
        sav[m] += [e1hmis/e1h,np.ones(len(e1h))*-1.,e1hmis/edm_rel,
                      e1hmis_r/e1h_r,np.ones(len(e1h))*-1.,e1hmis_r/edm]

    Dil[-1] = Dil[1]
    Dil[-2] = 1.

    for m in range(len(mod)):
        print(mod[m])
        for j in range(5):
            print((j+1)*10)
            for i in range(6):
                # ratio = Dil[j]/sav[m][i+j*6]
                ratio = sav[m][i+j*6][:-1]/Dil[j]
                print(mean(ratio),np.std(ratio))
            
        
def plt_profile_bias():
    

    samp = 'HM_Lz'

    matplotlib.rcParams.update({'font.size': 12})    
    
    p_name = 'profile_'+samp+'.fits'
    profile = fits.open(folder+p_name)

    print(p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS  = p.DSigma_T
    GT  = p.GAMMA_Tcos    
    GX  = p.GAMMA_Xsin
    
    # '''
    CovDS  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT  = cov.COV_GT.reshape(len(GT),len(GT))    
    CovGX  = cov.COV_GX.reshape(len(GT),len(GT))
        
    ##############    
    # misal
    ##############    
    
    p_name = 'profile_'+samp+'_mis30.fits'
    profile = fits.open(folder+p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS_al  = p.DSigma_T
    GT_al  = p.GAMMA_Tcos
    GX_al  = p.GAMMA_Xsin
    
    # '''
    CovDS_al  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT_al  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGX_al  = cov.COV_GX.reshape(len(GT),len(GT))

    ##############    
    # miscen
    ##############    
    
    p_name = 'profile_HM_Lz_Nmis_miscen.fits'
    profile = fits.open(folder+p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS_c  = p.DSigma_T
    GT_c  = p.GAMMA_Tcos
    GX_c  = p.GAMMA_Xsin
    
    # '''
    CovDS_c  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT_c  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGX_c  = cov.COV_GX.reshape(len(GT),len(GT))

    ##############    
    # miscen + misal
    ##############    
    
    p_name = 'profile_'+samp+'_mis20_miscen.fits'
    profile = fits.open(folder+p_name)
    
    p   = profile[1].data
    cov = profile[2].data
    
    DS_cal  = p.DSigma_T
    GT_cal  = p.GAMMA_Tcos
    GX_cal  = p.GAMMA_Xsin
    # '''
    CovDS_cal  = cov.COV_ST.reshape(len(GT),len(GT))
    CovGT_cal  = cov.COV_GT.reshape(len(GT),len(GT))
    CovGX_cal  = cov.COV_GX.reshape(len(GT),len(GT))
    
    # '''    
    y  = [DS,GT,GX]
    Cy = [np.diag(CovDS),np.diag(CovGT),np.diag(CovGX)]

    y_c  = [DS_c,GT_c,GX_c]
    Cy_c = [np.diag(CovDS_c),np.diag(CovGT_c),np.diag(CovGX_c)]

    y_al  = [DS_al,GT_al,GX_al]
    Cy_al = [np.diag(CovDS_al),np.diag(CovGT_al),np.diag(CovGX_al)]

    y_cal  = [DS_cal,GT_cal,GX_cal]
    Cy_cal = [np.diag(CovDS_cal),np.diag(CovGT_cal),np.diag(CovGX_cal)]
    
    y2  = [y_al,y_c,y_cal]
    Cy2 = [Cy_al,Cy_c,Cy_cal]
    
    labels = [r'$\sigma(\Delta \phi) = 30^{\circ}, p_{cc} = 1.0$',r'$\sigma(\Delta \phi) = 0^{\circ}, p_{cc} = 0.0$',r'$\sigma(\Delta \phi) = 20^{\circ}, p_{cc} = 0.75$']
    
    # f, ax = plt.subplots(6,3, gridspec_kw={'height_ratios': [3, 1]*3},figsize=(14,10),sharex = True)
    f, ax = plt.subplots(3,3,figsize=(14,8),sharex = True)
    f.subplots_adjust(hspace=0)
    
    for j in range(3):
        
        for i in range(3):
            ax[j,i].plot(p.Rp,y[i],'C7')        
            ax[j,i].fill_between(p.Rp,y[i]+Cy[i],y[i]-Cy[i],color='C7',alpha=0.4)
            ax[j,i].plot(p.Rp,y2[j][i],'C3',label = labels[j])        
            ax[j,i].fill_between(p.Rp,y2[j][i]+Cy2[j][i],y2[j][i]-Cy2[j][i],color='C3',alpha=0.4)
            # ax[j*2,i].plot(p.Rp,y[i],'C7')        
            # ax[j*2,i].fill_between(p.Rp,y[i]+Cy[i],y[i]-Cy[i],color='C7',alpha=0.4)
            # ax[j*2,i].plot(p.Rp,y2[j][i],'C3',label = labels[j])        
            # ax[j*2,i].fill_between(p.Rp,y2[j][i]+Cy2[j][i],y2[j][i]-Cy2[j][i],color='C3',alpha=0.4)
            # ax[j*2+1,i].plot(p.Rp,(y[i]-y2[j][i]),'C7')        
            # ax[j*2+1,i].plot([0,10],[0,0],'k')        
            # ax[j*2+1,i].set_ylim(-0.5,1.)


        ax[-1,j].set_xlabel('r [$h^{-1}$ Mpc]')
        
        
        # j *= 2
        
        ax[j,0].legend(frameon=False,loc=4)
        
        ax[j,0].set_xscale('log')
        ax[j,0].set_yscale('log')
        ax[j,0].set_ylabel(r'$\Delta\Sigma [h M_\odot/pc^2]$',labelpad=0.5)
        ax[j,0].set_ylim(0.5,200)
        ax[j,0].set_xlim(0.1,10)
        ax[j,0].xaxis.set_ticks([0.1,1,5,7])
        ax[j,0].set_xticklabels([0.1,1,5,7])
        ax[j,0].yaxis.set_ticks([1,10,100])
        ax[j,0].set_yticklabels([1,10,100])
        
        ax[j,1].set_xscale('log')
        ax[j,1].set_yscale('log')
        ax[j,1].set_ylabel(r'$\Gamma_T [h M_\odot/pc^2]$',labelpad=0.5)
        ax[j,1].set_ylim(0.5,100)
        ax[j,1].set_xlim(0.1,10)
        ax[j,1].xaxis.set_ticks([0.1,1,5,7])
        ax[j,1].set_xticklabels([0.1,1,5,7])
        ax[j,1].yaxis.set_ticks([1,10,100])
        ax[j,1].set_yticklabels([1,10,100])
            
        ax[j,2].plot([0,10],[0,0],'k')
        ax[j,2].set_ylabel(r'$\Gamma_\times [h M_\odot/pc^2]$',labelpad=0.5)
        ax[j,2].set_xscale('log')
        ax[j,2].set_xlim(0.1,10)
        ax[j,2].set_ylim(-16,17)
        ax[j,2].xaxis.set_ticks([0.1,1,5,7])
        ax[j,2].set_xticklabels([0.1,1,5,7])
        
        f.savefig(folder+'../final_plots/profile_test_bias.pdf',bbox_inches='tight')


def plot_withnoise():


    hsamps  = ['HM_Lz_qcut',
               'LM_Lz_qcut',
               'HM_Mz_qcut',
               'LM_Mz_qcut',
               'HM_Hz_qcut',
               'LM_Hz_qcut']

    RIN = ['350','350','350','350','450','400']
    ROUToq = ['2000','1000']*3
    
    qh,NFW_h,Ein_h,qhr,NFW_hr,Ein_hr,NFW,Ein,o1h,woc,zh = extract_params(hsamps,RIN,ROUToq)
    
    qfit_des = []
    qm_des   = []
    eq_des   = []

    qfit_sn = []
    qm_sn   = []
    eq_sn   = []

    for j in range(len(hsamps)):
        
        qs_sn = fits.open(folder+'fitresults_2h_2q_'+RIN[j]+'_5000_profile_'+hsamps[j]+'.fits')[1].data.q
        qfit_sn += [qs_sn]
        qm_sn += [np.median(qs_sn)]
        eq_sn += [np.diff(np.percentile(qs_sn, [16,50,84]))]

        qs_des = fits.open(folder+'fitresults_2h_2q_'+RIN[j]+'_5000_profile_'+hsamps[j]+'_des.fits')[1].data.q
        qfit_des += [qs_des]
        qm_des += [np.median(qs_des)]
        eq_des += [np.diff(np.percentile(qs_des, [16,50,84]))]
        
    qfit = [qfit_sn,qfit_des,qfit_des]
    qm   = [qm_sn,qm_des,qm_des]
    eq   = [eq_sn,eq_des,eq_des]


    xl = ['without shape noise', '6.3 gal/arcmin2','27 gal/arcmin2']

    
    # FOM
    f, ax = plt.subplots(1,1, figsize=(14,6),sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    
    ax = [ax]
    
    ax[0].plot([0,5],[0,0],'C7--')
    
    ax[0].axhspan(-0.05,0.05,0,5,color='C7',alpha=0.5)
    ax[0].axhspan(-0.025,0.025,0,5,color='C7',alpha=0.5)

    
    for hs in range(len(hsamps)):
        ax[0].plot([-1,-1],[-1,-1],cstyle[hs],label=lhs[hs])
        for fp in range(3):
            ax[0].errorbar(fp+1+0.1*hs,(qm[fp][hs]-qh[hs])/qh[hs],
                         yerr=np.array([eq[fp][hs]/qh[hs]]).T,
                         fmt=cstyle[hs],markersize=15)
                         
            ax[0].violinplot((qfit[fp][hs]-qh[hs])/qh[hs],positions=[fp+1+0.1*hs],showextrema=False, showmedians=True)
    
    ax[0].legend(frameon = False,loc=2,ncol=7)
    ax[0].set_ylabel(r'$(\tilde{q}_{1h}(\hat{\phi})-\langle q \rangle)/\langle q \rangle$')
    ax[0].axis([0.5,4,-0.6,0.6])
    ax[0].set_xticks(np.arange(3)+1)
    ax[0].set_xticklabels(xl)
    # f.savefig(folder+'../final_plots/model_q1h.pdf',bbox_inches='tight')
    # f.savefig(folder+'../test_plots/model_q1hr_'+nplot+'.png',bbox_inches='tight')


# '''

cstyle = ['C9v','C8v','C2v','ko','C1^','C3^','C5^']
hsamps = ['all_q1','all_q2','all_q3','all','all_q4','all_q5','all_q6']
hsamps = ['all_q1_relaxed','all_q2_relaxed','all_q3_relaxed','all_relaxed','all_q4_relaxed','all_q5_relaxed','all_q6_relaxed']

lhs = ['HM-Lz','LM-Lz','HM-Mz','LM-Mz','HM-Hz','LM-Hz']
cstyle = ['C1^','C1v','C3^','C3v','C5^','C5v']
ROUToq = ['2000','1000']*3
RIN = ['350','350','350','350','450','400']

lhsq = ['all','q < 0.5','0.5 < q < 0.7','q > 0.7']
cstyleq = ['ko','C1v','C3o','C5^']
ROUToqq = ['1500']*4
RINq = ['350']*4
hsampsq = ['all','all_qmin','all_qmed','all_qmax']
hsampsq_rel = ['all_relaxed','all_qmin_relaxed','all_qmed_relaxed','all_qmax_relaxed']

hsamps = ['HM_Lz_qcut','LM_Lz_qcut','HM_Mz_qcut','LM_Mz_qcut','HM_Hz_qcut','LM_Hz_qcut']
hsamps = ['HM_Lz','LM_Lz','HM_Mz','LM_Mz','HM_Hz','LM_Hz']
hsamps_rel = ['HM_Lz_relaxed','LM_Lz_relaxed','HM_Mz_relaxed','LM_Mz_relaxed','HM_Hz_relaxed','LM_Hz_relaxed']
hsamps_mis20 = ['HM_Lz_mis20','LM_Lz_mis20','HM_Mz_mis20','LM_Mz_mis20','HM_Hz_mis20','LM_Hz_mis20']
hsamps_miscen = ['HM_Lz_miscen','LM_Lz_miscen','HM_Mz_miscen','LM_Mz_miscen','HM_Hz_miscen','LM_Hz_miscen']
hsamps_misall = ['HM_Lz_mis20_miscen','LM_Lz_mis20_miscen','HM_Mz_mis20_miscen','LM_Mz_mis20_miscen','HM_Hz_mis20_miscen','LM_Hz_mis20_miscen']

# plot_bias(hsampsq,lhsq,cstyleq,'all_q',RINq,ROUToq)
# plot_bias(hsamps,lhs,cstyle,'all',RIN,ROUToq)
# plot_bias(hsamps_rel,lhs,cstyle,'all',RIN,ROUToq)
# plot_bias(hsamps_misall,lhs,cstyle,'bias',RIN,ROUToq,0.78)
# plot_M200q(hsamps,lhs,cstyle,RIN,ROUToq)

# plot_M200q(hsamps_rel,lhs,cstyle,RIN,ROUToq)
# '''
# from basic_extract import save_fitted

# for j in range(len(hsamps_rel)):
        # samp = hsamps_rel[j]
        # rin = float(RIN[j])
        # save_fitted(samp,rin,5000,fittype='_2h_2q')
        # save_fitted(samp,rin,5000,fittype='_2h_2q_Ein')
        # save_fitted(samp,rin,5000,fittype='_2h_2q_woc')
        
'''

f, ax_all = plt.subplots(6,3, figsize=(14,16),sharex = True)
f.subplots_adjust(hspace=0)


for j in range(len(ax_all)):
    rin = float(RIN[j])
    plt_profile_fitted_final(hsamps_rel[j],rin,5000,ax_all[j],fittype='_2h_2q')
    ax_all[j,0].text(1,100,lhs[j],fontsize=14)

ax_all[0,0].legend(loc=3,frameon=False)
ax_all[0,1].legend(loc=3,frameon=False)


    
f.savefig(folder+'../final_plots/profile_rel.png',bbox_inches='tight')



f, ax_all = plt.subplots(4,3, figsize=(14,16),sharex = True)
f.subplots_adjust(hspace=0)


for j in range(len(ax_all)):
    plt_profile_fitted_final(hsampsq[j],350,5000,ax_all[j],fittype='_2h_2q_Ein')
    ax_all[j,0].text(1,100,lhsq[j],fontsize=14)

ax_all[0,0].legend(loc=3,frameon=False)
ax_all[0,1].legend(loc=3,frameon=False)


    
f.savefig(folder+'../test_plots/profile_all_q_Ein.png',bbox_inches='tight')
'''


hsampsq  = ['HM_Lz_qcut',
               'LM_Lz_qcut',
               'HM_Mz_qcut',
               'LM_Mz_qcut',
               'HM_Hz_qcut',
               'LM_Hz_qcut']


f, ax_all = plt.subplots(4,3, figsize=(14,16),sharex = True)
f.subplots_adjust(hspace=0)


for j in range(len(ax_all)):
    plt_profile_fitted_withnoise(hsampsq[j],350,5000,ax_all[j])
    ax_all[j,0].text(1,100,lhsq[j],fontsize=14)

ax_all[0,0].legend(loc=3,frameon=False)
ax_all[0,1].legend(loc=3,frameon=False)


    
f.savefig(folder+'../test_plots/profile_all_q_Ein.png',bbox_inches='tight')
