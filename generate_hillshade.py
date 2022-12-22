#!/usr/bin/env python3

from scipy.io import loadmat
import pickle
import numpy as np
import matplotlib.pyplot as plt
import os.path
import sys

def cart2pol(y, x):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(x, y)
    return phi, rho

def setborder(grid, bs, bvalue):
    grid[:bs+1, :] = bvalue # top rows
    grid[-1-bs:, :] = bvalue # bottom rows
    grid[:,:bs+1] = bvalue # left cols
    grid[:,-1-bs:] = bvalue # right cols
    return grid

def hillshade(dem, X, Y, azimuth=315, altitude=45, zf=1, plotit=False):	

    if np.unique(np.diff(X)).size > 1 or np.unique(np.diff(Y)).size > 1:
        raise Exception('hillshade function assumes constant spacing in x and y directions')

    # get cell spacing in x and y direction
    dx = np.abs(X[1]-X[0])
    dy = np.abs(Y[1]-Y[0])

    # Determine extents from X and Y vectors
    extent =[np.min(X), np.max(X), np.min(Y), np.max(Y)]

    # lighting azimuth
    azimuth = 360.0 - azimuth + 90 # convert to mathematic unit
    if azimuth >= 360.0:
        azimuth = azimuth - 360.0
    azimuth = azimuth * (np.pi/180.0) # convert to radians

    # lighting altitude
    altitude = (90.0 - altitude) * (np.pi / 180.0) # convert zenith angle to radians

    # calc slope and aspect (radians)
    fy, fx = np.gradient(dem,dx,dy, edge_order=2) # uses simple, unweighted gradient of immediate neighbors #TMP
    asp, grad = cart2pol(fy, fx) # convert to cartesian coordinates #TMP
    grad = np.arctan(zf*grad) # steepest slope
    # convert asp
    asp[asp<np.pi] = asp[asp<np.pi] + (np.pi/2.0)
    asp[asp<0] = asp[asp<0] + (2.0*np.pi)

    # hillshade calculation
    h = 255.0 * ( (np.cos(altitude) * np.cos(grad)) + (np.sin(altitude) * np.sin(grad) * np.cos(azimuth-asp)) ) # ESRIs algorithm
    h[h<0] = 0 # set hillshade values to min of 0

    h = setborder(h, 0, np.nan) # set border cells to NaN

    if plotit:
        fig = plt.figure(0)
        plt.imshow(h, cmap='gray', vmin=0, vmax=255)
        plt.locator_params(axis='both', nbins=5)
        fig.savefig('hillshade.png')
        plt.close('all')
    return h

def main(datfile):
    if datfile.endswith('.mat'):
        m = loadmat(datfile)
        suff = '_hillshade_matlab.png'
    elif datfile.endswith('.pcl'):
        with open(datfile, 'rb') as f:
            m = pickle.load(f)
        suff = '_hillshade_python.png'
    else:
        raise Exception(f'File {datfile} is not a valid format format. Accepted formats are: .mat and .pcl.')
    outname = os.path.basename(datfile).split('.')[0] + suff
    if os.path.exists(outname):
        print(f'{outname} already exists, skipping')
        return
    hillshade(m['z'], m['x'].flatten(), m['y'].flatten(), zf=4, plotit=True)
    os.rename('hillshade.png', outname)
    print(f'Wrote to {outname}')

if __name__ == "__main__":
    datfile = sys.argv[1]
    main(datfile)
