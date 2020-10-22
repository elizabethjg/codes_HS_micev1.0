import os
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
sys.path.append('/mnt/clemente/lensing')

from lensing37.gentools import classonly
from lensing37.LensCat.main import Survey, read_columns, cat_paths


class MICE(Survey):
	
	
	def __init__(self):
		pass


	@classonly
	def load(cls):
		'''Method to load the catalogue
		'''

		path = '/mnt/clemente/lensing/HALO_SHAPE/MICE_v1.0/catalogs/MICE_sources.fits'
		columns = fits.open(path)[1].columns

		# Somehow we load the data and save it in cls.data
		f      = fits.open(path)
		dl = Table(f[1].data).to_pandas()  
		catid = np.arange(dl.shape[0]).astype(np.int32)
		catid = pd.DataFrame({'CATID': catid})
		cls.data = pd.concat([catid,dl], axis=1)


	@classonly
	def drop(cls):
		'''Method to drop the Survey and GriSPy grid from memory
		'''
		try:
			del	cls.data
		except AttributeError:
			pass
		
		try:
			del cls.gsp
		except AttributeError:
			pass
