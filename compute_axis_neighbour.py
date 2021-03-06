import sys
sys.path.append('/mnt/projects/lensing/lens_codes_v3.7')
sys.path.append('/home/eli/lens_codes_v3.7')
import time
import numpy as np
from multiprocessing import Pool
from multiprocessing import Process
from astropy.io import fits
from member_distribution import *
from scipy import spatial

t1 = time.time() 

folder = '../../MICEv2.0/catalogs/'
cat = fits.open(folder+'MICE_halo_cat_withshapes.fits')[1].data
pesos = True

t2 = time.time()

ncores = 32

xc = cat.xc
yc = cat.yc
zc = cat.zc
lgm = cat.lgm

tree=spatial.cKDTree(np.array([xc,yc,zc]).T)

def axis_neigh(indices,nfile,pesos = False):

        if pesos:

                f = open('out_w_'+str(nfile),'w')
        else:
                f = open('out_'+str(nfile),'w')

        times = []

        show = 5000

        for j in indices:
                t1 = time.time()
                dist,ind=tree.query(np.array([xc[j],yc[j],zc[j]]).T,k=6)
                
                if pesos:
                        v3d,w3d,v2d,w2d = compute_axis(xc[ind],yc[ind],zc[ind],10**lgm[ind])
                else:
                        v3d,w3d,v2d,w2d = compute_axis(xc[ind],yc[ind],zc[ind])
                
                abc = str(np.sqrt(w3d[0]))+'   '+str(np.sqrt(w3d[1]))+'     '+str(np.sqrt(w3d[2]))+'   '
                av  = str(v3d[0,0])+'   '+str(v3d[1,0])+'     '+str(v3d[2,0])+'   '
                bv  = str(v3d[0,1])+'   '+str(v3d[1,1])+'     '+str(v3d[2,1])+'   '
                cv  = str(v3d[0,2])+'   '+str(v3d[1,2])+'     '+str(v3d[2,2])+'   '
                
                ab  = str(np.sqrt(w2d[0]))+'   '+str(np.sqrt(w2d[1]))+'    '
                avp  = str(v2d[0,0])+'   '+str(v2d[1,0])+'     '
                bvp  = str(v2d[0,1])+'   '+str(v2d[1,1])+'     '


                
                f.write(str(j)+'  '+abc+av+bv+cv+ab+avp+bvp+'\n')
                
                t2 = time.time()
                
                times += [(t2-t1)*60]
                
                if j == show:
                
                        print('Estimated ramaining time')
                        print((np.mean(times)*(len(indices)-show))/60.)
                        
                        show += 5000
        
        f.close()


def axis_neigh_unpack(minput):
	return axis_neigh(*minput)

index = np.arange(len(cat))

slicer = int(round(len(index)/float(ncores), 0))
slices = ((np.arange(ncores-1)+1)*slicer).astype(int)
slices = slices[(slices <= len(index))]

index_splitted = np.split(index,slices)

files = np.arange(ncores)
if pesos:
        pesos = np.ones(ncores).astype(bool)
else:
        pesos = np.zeros(ncores).astype(bool)

entrada = np.array([index_splitted,files,pesos]).T

pool = Pool(processes=(ncores))
pool.map(axis_neigh_unpack, entrada)
pool.terminate()


t3 = time.time()

print('Tiempo de lectura')
print(t2 - t1)
print('Tiempo de ejecución')
print(t3 - t2)
print('Tiempo total')
print(t3 - t1)
