import sys
import numpy as np
import os.path
from scipy.io import loadmat
from readGeotiff import readGeotiff
import pickle

def mad(arr, axis=None, flag=1, keepdims=False):
	if flag==0:
		mean = np.nanmean(arr, axis=axis, keepdims=True)
		mad = np.nanmean(np.abs(arr-mean), axis=axis, keepdims=keepdims)
	elif flag==1:
		median = np.nanmedian(arr, axis=axis, keepdims=True)
		mad = np.nanmedian(np.abs(arr-median), axis=axis, keepdims=keepdims)
	return mad

def getNeighborOffsets(fileNames, nB):
	# initialize output
	dz = np.full((nB.shape[0], 1), np.nan)
	dz_mad = np.full((nB.shape[0], 1), np.nan)

	fileNames = sorted(fileNames) #TMP
	ext = [os.path.splitext(e)[1].lower() for e in fileNames]

	# Get offsets between tile above
	for n in np.arange(nB.shape[0]):
		print(f'subtile {n} of {nB.shape[0]}')

		# read this tile coordinates
		if ext[nB[n,0]] == '.mat':
			m0 = loadmat(fileNames[nB[n,0]])
		elif ext[nB[n,0]] == '.tif':
			m0 = readGeotiff(fileNames[nB[n,0]], mapInfoOnlyFlag=True)

		# read coordinates for tile above
		if ext[nB[n,1]] == '.mat':
			m1 = loadmat(fileNames[nB[n,1]])
		elif ext[nB[n,1]] == '.tif':
			m1 = readGeotiff(fileNames[nB[n,1]], mapInfoOnlyFlag=True)

		# find indices of overlap
		_, cols0, cols1 = np.intersect1d(m0.x, m1.x, return_indices=True)
		_, rows0, rows1 = np.intersect1d(m0.y, m1.y, return_indices=True)

		# check to make sure buffers are equal size between tiles
		if len(cols0) != len(cols1):
			raise Exception('buffer cols not same length')
		elif len(rows0) != len(rows1):
			raise Exception('buffer rows not same length')

		# read buffer for this tile
		if ext[nB[n,0]] == '.mat':
			x0 = m0.x
			y0 = m0.y
			x0 = x0[cols0]
			y0 = y0[rows0]
			z0 = m0.z[rows0, cols0]
		elif ext[nB[n,0]] == '.tif':
			#print(f'{cols0=}, {rows0=}')
			#print(f'{[cols0[0], cols0[-1]+1, rows0[-1], rows0[0]+1]=}')
			z0 = readGeotiff(fileNames[nB[n,0]], pixel_subset=[cols0[0], cols0[-1], rows0[-1], rows0[0]])
			#print(f'{z0.x=}, {z0.y=}, {z0.z=}')
			#print(f'{z0.y.shape=}, {z0.x.shape=}')
			x0 = z0.x
			y0 = z0.y
			z0 = z0.z
		# read buffer for neighbor tile
		if ext[nB[n,1]] == '.mat':
			x1 = m1.x
			y1 = m1.y
			x1 = x1[cols1]
			y1 = y1[rows1]
			z1 = m1.z[rows1,cols1]
		elif ext[nB[n,1]] == '.tif':
			#print(f'{[cols1[0], cols1[-1]+1, rows1[-1], rows1[0]+1]=}')
			z1 = readGeotiff(fileNames[nB[n,1]], pixel_subset=[cols1[0], cols1[-1], rows1[-1], rows1[0]])
			#print(f'{z1.y.shape=}, {z1.x.shape=}')
			x1 = z1.x
			y1 = z1.y
			z1 = z1.z
		z0[z0 == -9999] = np.nan
		z1[z1 == -9999] = np.nan
		#print(f'{z0=}, {z1=}')
		#sys.exit(1)
		#print(f'{y0=}, {y1=}, {x0=}, {x1=}')
		#with open('neighbors.pcl', 'wb') as f:
		#	pickle.dump({'y0': y0, 'y1': y1, 'x0': x0, 'x1': x1, 'z1': z1, 'z0': z0}, f)
		#sys.exit(1)
		if not np.any(y0!=y1) and not np.any(x0!=x1):
			dzn = z0.flatten() - z1.flatten()
			if np.sum(~np.isnan(dzn))/dzn.size > 0.1:
				dz[n] = np.nanmedian(dzn)
				dz_mad[n] = mad(dzn, flag=1)
		else:
			print('inconsistent buffers, offset not computed')

	return dz, dz_mad
