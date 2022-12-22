import rasterio
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)
# create return variable classes
class Mapinf:
	def __init__(self):
		self.dx = None
		self.dy = None
		self.mapx = None
		self.mapy = None
class Info:
	def __init__(self):
		self.cols = None
		self.rows = None
		self.isize = None
		self.bands = None
		self.map_info = Mapinf()
		
class I:
	def __init__(self):
		self.z = None
		self.x = None
		self.y = None
		self.info = None
		self.Tinfo = None

# read a geotiff file by name. return a populated class above
def readGeotiff(name, reportMapSubsetError=False, pixel_subset=None, map_subset=None, mapInfoOnlyFlag=False):
	if name is None: 
		return None

	# example variables
	# dataset.bounds = BoundingBox(left=-2017348.0, bottom=708478.0, right=-1996638.0, top=731298.0)
	# dataset.res = (10.0,10.0)
	dataset = rasterio.open(name)
	if dataset is None:
		raise Exception("tiff file opened improperly. readGeotiff.py")

	info = Info()
	info.cols = dataset.width
	info.rows = dataset.height
	info.imsize = dataset.offsets #potentially dataset.scales
	info.bands = dataset.count

	info.map_info.dx = dataset.res[0]
	info.map_info.dy = dataset.res[1]
	info.map_info.mapx = dataset.bounds[0]
	info.map_info.mapy = dataset.bounds[3]

	subrows = np.array([0, info.rows-1])
	subcols = np.array([0, info.cols-1])

	minx = info.map_info.mapx
	maxy = info.map_info.mapy

	x = minx + (np.arange(info.cols) * info.map_info.dx)
	y = maxy - (np.arange(info.rows) * info.map_info.dy)

	if pixel_subset is not None:
		subrows = np.array(pixel_subset[2:4])
		subcols = np.array(pixel_subset[0:2])
	if map_subset is not None:
		map_subset = np.array(map_subset)
		subcols = (map_subset[0:2]-info.map_info.mapx)/info.map_info.dx+[0,1]
		subrows = (info.map_info.mapy - map_subset[-1:1:-1])/info.map_info.dy+[0,1]
		subcols = np.around(subcols)
		subcols = subcols.astype(int)
		subrows = np.around(subrows)
		subrows = subrows.astype(int)

		subcols[np.where(subcols < 0)] = 0
		subrows[np.where(subrows < 0)] = 0
		subcols[np.where(subcols >= info.cols)] = info.cols-1
		subrows[np.where(subrows >= info.rows)] = info.rows-1

	returnI = I()
	returnI.x = x[int(subcols[0]):int(subcols[1])+1]
	returnI.y = y[int(subrows[0]):int(subrows[1])+1]
	returnI.z = []	

	if not mapInfoOnlyFlag:
		dataset2 = rasterio.open(name)
		if dataset2 is None:
			raise Exception("tiff file opened improperly. readGeotiff.py, dataset2")
		returnI.z = dataset2.read(1)
		returnI.z = np.array(returnI.z)
		returnI.z = returnI.z[subrows[0]:subrows[1]+1, subcols[0]:subcols[1]+1]

	returnI.info = info
	returnI.Tinfo = dataset	
	returnI.z = np.array(returnI.z)
	return returnI
