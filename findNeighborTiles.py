import sys
from scipy.io import loadmat
import os.path
import numpy as np
from readGeotiff import readGeotiff
import pickle

def findNeighborTiles(fileNames):
	
	# get filename extensions for determining import function
	ext = [os.path.splitext(e)[1].lower() for e in fileNames]

	# loop through each tile and get extents
	x0 = np.full((len(fileNames), 1), np.nan)
	x1 = np.full((len(fileNames), 1), np.nan)
	y0 = np.full((len(fileNames), 1), np.nan)
	y1 = np.full((len(fileNames), 1), np.nan)
	print('looping through tile files to get extents')
	for i in np.arange(len(fileNames)):
		if ext[i] == '.mat':
			m = loadmat(fileNames[i])
		elif ext[i] == '.tif':
			m = readGeotiff(fileNames[i], mapInfoOnlyFlag=True)

		x0[i] = np.min(m.x)
		x1[i] = np.max(m.x)
		y0[i] = np.min(m.y)
		y1[i] = np.max(m.y)

	# using ranges, find neighbors on each side
	ntop = np.full((len(fileNames),2), np.nan)
	nright = np.full((len(fileNames),2), np.nan)
	ntopRight = np.full((len(fileNames),2), np.nan)
	nbottomRight = np.full((len(fileNames),2), np.nan)

	# mean ranges
	mnx = (x0+x1) / 2.0
	mny = (y0+y1) / 2.0

	print('looping through ranges to find neighbors')
	for i in np.arange(len(fileNames)):
		#if fileNames[i] != '/fs/ess/PZS0720/skhuvis/mosaicking/Iturralde/tiles/Iturralde_19_0_20220615_619000_8582000.tif':
		#	continue
		# top
		n = np.where(np.logical_and(y0[i] < y0, np.logical_and(y1[i] >= y0,  np.logical_and(mnx[i] > x0, mnx[i] < x1))))[0]
		if n.size > 0:
			ntop[i,:] = [i,n]

		# right
		n = np.where(np.logical_and(x0[i] < x0, np.logical_and(x1[i] >= x0,  np.logical_and(mny[i] > y0, mny[i] < y1))))[0]
		if n.size > 0:
			nright[i,:] = [i,n]

		# top-right
		n = np.where(np.logical_and(y0[i] < y0, np.logical_and(y1[i] >= y0,  np.logical_and(x0[i] < x0, np.logical_and(x1[i] >= x0, np.logical_and(mnx > x1[i], mny > y1[i]))))))[0]
		if n.size > 0:
			ntopRight[i,:] = [i,n]

		# bottom-right
		n = np.where(np.logical_and(y0[i] > y0, np.logical_and(y0[i] <= y1, np.logical_and(x0[i] < x0, np.logical_and(x1[i] >= x0, np.logical_and(mnx > x1[i], mny < y0[i]))))))[0]
		if n.size > 0:
			nbottomRight[i,:] = [i,n]
	ntop = np.delete(ntop, np.isnan(ntop[:,0]), axis=0).astype(int)
	nright = np.delete(nright, np.isnan(nright[:,0]), axis=0).astype(int)
	ntopRight = np.delete(ntopRight, np.isnan(ntopRight[:,0]), axis=0).astype(int)
	nbottomRight = np.delete(nbottomRight, np.isnan(nbottomRight[:,0]), axis=0).astype(int)

	X = np.concatenate((x0, x1, y0, y1), axis=1)

	return ntop, nright, ntopRight, nbottomRight, X
