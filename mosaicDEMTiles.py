import sys
import pickle
from osgeo import ogr
from findNeighborTiles import findNeighborTiles
from readGeotiff import readGeotiff
from getNeighborOffsets import getNeighborOffsets
from adjustOffsets import adjustOffsets
import numpy as np
from inpaint_nans import inpaint_nans

if sys.version_info[0] < 3:
	import raster_array_tools as rat
else:
	from lib import raster_array_tools as rat

def blendTile(z, zsub):
	# fill missing data in this subtile with current data in mosaic subset
	z[np.logical_and(np.isnan(z), np.logical_not(np.isnan(zsub)))] = zsub[np.logical_and(np.isnan(z), np.logical_not(np.isnan(zsub)))]
	# Create the blending array by setting zeros at the far edge of
	# subtile/tile overlap and ones at the other edge, linearly
	# interpolating between the two.
	# find where data is missing in both subtile and tile
	buffA = np.logical_not(np.logical_and(~np.isnan(z), ~np.isnan(zsub))).astype(float)
	# set pixels with data in both subtile and tiles to NaN as an
	# interpolation flag for inpaint_nans
	buffA[~(buffA.astype(bool))] = np.nan
	
	# set boundaries of blend array to zero
	buffA[0,np.isnan(buffA[0,:])] = 0
	buffA[-1, np.isnan(buffA[-1,:])] = 0
	
	buffA[np.isnan(buffA[:,0]), 0] = 0
	buffA[np.isnan(buffA[:,-1]), -1] = 0
	
	# interpolate linearly across NaNs
	buffA = inpaint_nans(buffA, 2)
	
	# find where there is data in the tile subset
	notMissing = np.logical_not(np.isnan(zsub))
	
	# blend the subtile and mosaic subset, applying the edge-distance
	# weighting
	z[notMissing] = z[notMissing] * buffA[notMissing] + zsub[notMissing]*(1-buffA[notMissing])
	#with open('blending.pcl', 'wb') as f:
	#	pickle.dump({'z': z, 'buffA': buffA, 'notMissing': notMissing, 'zsub': zsub}, f)
	#sys.exit(1)

	return z

def mosaicDEMTiles(fileNames, dx=None, extent=[]):
	
	fileNames = sorted(fileNames) #TMP
	# find neighbors
	ntop, nright, _, _, X = findNeighborTiles(fileNames)

	# Get offsets between tile above
	dzup, dzup_mad = getNeighborOffsets(fileNames, ntop)

	# Get offsets between tile on right
	dzrt, dzrt_mad = getNeighborOffsets(fileNames, nright)

	# Adjustment
	dZ = adjustOffsets(len(fileNames), np.concatenate((ntop[:,0], nright[:,0])).reshape(-1,1), np.concatenate((ntop[:,1], nright[:,1])).reshape(-1,1), np.concatenate((dzup, dzrt)).reshape(-1,1), np.concatenate((dzup_mad, dzrt_mad)).reshape(-1,1))

	# Define Mosaic
	if extent == []:
		# define output grid based on file extents if not provided
		x0 = np.min(X[:,0])
		x1 = np.max(X[:,1])
		y0 = np.min(X[:,2])
		y1 = np.max(X[:,3])
	else:
		x0 = extent[0]
		x1 = extent[1]
		y0 = extent[2]
		y1 = extent[3]

	# define output resolution based on first DEM file if not specified
	if not dx:
		m = readGeotiff(fileNames[0], mapInfoOnlyFlag=True)
		dx = m.x[1] - m.x[0]

	# make a polyshape out of boundary for checking subtile overlap
	tile_wkt = f"POLYGON (({x0} {y0}, {x0} {y1}, {x1} {y1}, {x1} {y0}, {x0} {y0}))"
	tile_poly = ogr.CreateGeometryFromWkt(tile_wkt)

	# build tile coordinate vectors
	x = np.arange(x0, x1+float(dx), float(dx))
	y = np.arange(y1, y0-float(dx), -float(dx))

	# build mosaic output array
	z = np.full((len(y), len(x)), np.nan)

	## Add DEMs with adjustments to mosaic

	# initialize subtile count for use in n-weighted alignment
	subtile_n = 1

	# initialize count of pixels with data in mosaic for error checking
	NTiles = len(fileNames)
	Nn = 0
	for filen in np.arange(NTiles):
		# first only add subtiles with adjustments
		if np.isnan(dZ[filen]):
			continue

		print(f'adding subtile {filen} of {NTiles}: {fileNames[filen]}')

		# read DEM geotiff map coordinates to get extent
		m = readGeotiff(fileNames[filen], mapInfoOnlyFlag=True)

		# make a polyshape out of this subtile boundary
		minx = np.min(m.x)
		maxx = np.max(m.x)
		miny = np.min(m.y)
		maxy = np.max(m.y)
		filen_wkt = f"POLYGON (({minx} {miny}, {minx} {maxy}, {maxx} {maxy}, {maxx} {miny}, {minx} {miny}))"
		filen_poly = ogr.CreateGeometryFromWkt(filen_wkt)

		# test overlap
		if not tile_poly.Intersects(filen_poly):
			print('subtile out of bounds, skipping')
			continue

		# read portion of DEM covering mosaic
		m = readGeotiff(fileNames[filen], map_subset=[np.min(x), np.max(x), np.min(y), np.max(y)])
		#sys.exit(1)

		# set nodata to nan (would be better to read from tif if given)
		m.z[m.z==-9999] = np.nan

		# add adjustment offset
		m.z = m.z - dZ[filen]

		# find indices of overlap between tile and mosaic
		ncols = np.where(np.logical_and(x>=minx, x<=maxx))[0]
		nrows = np.where(np.logical_and(y>=miny, y<=maxy))[0]

		# check if tile grid is same as mosaic grid
		if m.x[1] - m.x[0] != dx or m.y[1]-m.y[0] != -dx or m.x[0] != x[ncols[0]] or m.y[0] != y[nrows[0]] or len(m.x) != len(ncols) or len(m.y) != len(nrows):
			# if not, interpolate to mosaic grid
			m.z = rat.interp2_gdal(m.x, m.y, m.z, x[ncols], y[nrows], 'bilinear')

		# subset area of tile from mosaic
		zsub = z[np.meshgrid(nrows, ncols)].T

		# find overlapping pixels with values between DEM and mosaic subset
		n_overlap = np.logical_and(~np.isnan(zsub.flatten()), ~np.isnan(m.z.flatten()))

		if np.any(n_overlap):
			# if there's overlapping pixels, blend tile into mosaic with edge
			# feathering
			m.z = blendTile(m.z, zsub)
		else:
			# if no pixels overlap, just add the subset data into the subtile,
			# replacing the NaNs
			m.z[np.logical_not(np.isnan(zsub))] = zsub[np.logical_not(np.isnan(zsub))]

		# place the blended subtile into the tile
		z[np.meshgrid(nrows, ncols)] = m.z.T
		#with open('interim.pcl', 'wb') as f:
		#	pickle.dump({'z': z}, f)

		# The below was error checking during initial testing.
		# Probably not needed now
		# count the number of pixels with data after this merge
		Nn1 = np.sum(np.logical_not(np.isnan(z)))

		# if the new number of pixels with data is now less than before,
		# then data was overwritten with NaNs, which is an error
		if Nn1 < Nn:
			raise Exception('more nan in mosaic with this iteration')

		# reset count of pixels with data for next iteration:
		subtile_n = subtile_n + 1

	## Add DEMs with nan dZ adjustment values by subtracting buffer difference

	# index vector of offsets with nans. We will remove these as they are added
	# to the mosaic (or are found to not contain usable data).
	nf = np.where(np.isnan(dZ))[0].astype(int)

	length_nf = len(nf)
	# count of which dem in the nf vector to attempt to add to the mosaic. If
	# there is no overlapping data currently in the mosaic, the count will
	# increase to the next dem and so forth, rotating back to 1 once it cycles
	# through. If the count goes through a full cycle of 1:length(nf), the
	# remaining nf dems will be added without registration.
	count=0

	# flag to indicate if that addition of all dems have failed due to lack of
	# overlap and to add without registration.
	noRegFlag = False

	while nf.size > 0:
		if count >= nf.size: # cycle complete, reset counter
			count = 0
			# check if no more subtiles added this cycle
			if nf.size == length_nf:
				noRegFlag = True # turn off registration if no overlap
			# reset length nf of this cycle to commpare with next
			length_nf = nf.size
		filen = nf[count]
		print(f'{nf.size} remaining subtiles, attempting to add: {fileNames[nf[count]]}')

		# THIS PART IS REPEATED FROM ABOVE UNLESS OTHERWISE INDICATED
		m = readGeotiff(fileNames[filen], mapInfoOnlyFlag=True)

		# Check if subtile overlaps mosaic boundary
		minx = np.min(m.x)
		maxx = np.max(m.x)
		miny = np.min(m.y)
		maxy = np.max(m.y)
		filen_wkt = f"POLYGON (({minx} {miny}, {minx} {maxy}, {maxx} {maxy}, {maxx} {miny}, {minx} {miny}))"
		filen_poly = ogr.CreateGeometryFromWkt(filen_wkt)

		# test overlap
		if not tile_poly.Intersects(filen_poly):
			nf = np.delete(nf, count)
			print('subtile out of bounds, skipping')
			continue

		m = readGeotiff(fileNames[filen], map_subset=[minx, maxx, miny, maxy])
		m.z[m.z==-9999] = np.nan

		# find indices of overlap between tile and mosaic
		ncols = np.where(np.logical_and(x>=minx, x<=maxx))[0]
		nrows = np.where(np.logical_and(y>=miny, y<=maxy))[0]

		# check if tile grid is same as mosaic grid
		if m.x[1] - m.x[0] != dx or m.y[1]-m.y[0] != -dx or m.x[0] != x[ncols[0]] or m.y[0] != y[nrows[0]] or len(m.x) != len(ncols) or len(m.y) != len(nrows):
			# if not, interpolate mosaic grid
			m.z = rat.interp2_gdal(m.x, m.y, m.z, x[ncols], y[nrows], 'bilinear')

		zsub = z[np.meshgrid(nrows, ncols)].T

		# find overlapping non-nan and non-water pixels
		n_overlap = np.logical_and(np.logical_and(~np.isnan(zsub.flatten()), ~np.isnan(m.z.flatten())))

		# check if overlapping pixels exist to determine if blending is needed
		if np.any(n_overlap):
			# THIS PART IS DIFFERENT
			# get median difference between this subtile and mosaic subset
			dz_med = np.nanmedian(m.z[n_overlap] - zsub[n_overlap])
			m.z = m.z - dz_med

			# blend tile into mosaic with edge feathering
			m.z = blendTile(m.z, zsub)
		else:
			if noRegFlag:
				m.z[np.logical_not(np.isnan(zsub))] = zsub[np.logical_not(np.isnan(zsub))]
			else:
				# skip this subtile for this cycle
				count = count + 1
				continue

		# place the blended subtile into the tile
		z[np.meshgrid(nrows, ncols)] = m.z.T

		# count the number of pixels with data after this merge
		Nn1 = np.sum(~np.isnan(z))

		# if the new number of pixels with data is now less than before,
		# then data was overwritten with NaNs, which is an error
		if Nn1 < Nn:
			raise Exception('more nans in mosaic with this iteration')

		# reset count of pixels with data for next iteration
		Nn = Nn1

		# remove this subtile from index
		nf = np.delete(nf, count)

	return x, y, z
