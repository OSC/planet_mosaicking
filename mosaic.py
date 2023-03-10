#!/usr/bin/env python3

import glob
import pickle
from mosaicDEMTiles import mosaicDEMTiles

tileDir = '/fs/ess/PZS0720/skhuvis/mosaicking/Iturralde/tiles'
fileNames = sorted(glob.glob(tileDir + '/*0.tif'))
[x,y,z] = mosaicDEMTiles(fileNames)
with open('mosaic.pcl', 'wb') as f:
	pickle.dump({'x': x, 'y': y, 'z': z}, f)
