import numpy as np
import pickle
import sys

def adjustOffsets(NTiles, n1, n2, dz, dze):
	# remove nan pairs or with large offsets and/or mad values
	input_array = np.concatenate((n1, n2, dz, dze), axis=1)
	n = np.isnan(input_array).any(axis=1)
	n1 = np.delete(n1, n)
	n2 = np.delete(n2, n)
	dz = np.delete(dz, n)
	dze = np.delete(dze, n)

	Npairs = len(n1)

	# Build design and weight matrices
	A = np.zeros((Npairs, NTiles))
	A[np.arange(Npairs), n1] = 1
	A[np.arange(Npairs), n2] = -1
	
	# locate missing tiles
	n_missing = np.logical_not(np.any(A, axis=0))
	A = np.delete(A, n_missing, axis=1)

	# add delta=0 lines
	A = np.concatenate((A,np.diag(np.ones(A.shape[1]))))
	dz = np.concatenate((dz, np.zeros((A.shape[1])))).reshape(-1,1)
	dze = np.concatenate((np.ones(dze.shape), np.ones((A.shape[1])) * 100))

	# calculate weights
	wz = 1.0 / dze.reshape(-1,1)

	# build offset vector
	dZ = np.full((NTiles, 1), np.nan)

	a=wz*A
	b=wz*dz
	x,_,_,_=np.linalg.lstsq(a,b,rcond=None)
	dZ[~n_missing] = x

	return dZ
