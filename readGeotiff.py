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

	subrows = np.array([1, info.rows])
	subcols = np.array([1, info.cols])

	minx = info.map_info.mapx
	maxy = info.map_info.mapy

	x = minx + (np.arange(0,info.cols) * info.map_info.dx)
	y = maxy - (np.arange(0,info.rows) * info.map_info.dy)

	if pixel_subset is not None:
		subrows = np.array(pixel_subset[2:4])
		subcols = np.array(pixel_subset[0:2])
	if map_subset is not None:
		map_subset = np.array(map_subset)
		subcols = (map_subset[0:2:]-info.map_info.mapx)/info.map_info.dy+1
		subrows = (info.map_info.mapy - map_subset[2:4:])/info.map_info.dy+1
		subcols = np.around(subcols)
		subcols = subcols.astype(int)
		subrows = np.around(subrows)
		subrows = subrows.astype(int)
		if reportMapSubsetError:
			if any(subcols < 1) or any(subrows < 1) or any(subcols > info.cols) or any(subrows > info.rows):
				raise Exception('readGeotiff map_subset extends beyond raster extend')
		subcols[np.where(subcols < 1)] = 1
		subrows[np.where(subrows < 1)] = 1
		subcols[np.where(subcols > info.cols)] = info.cols
		subrows[np.where(subrows > info.rows)] = info.rows

	returnI = I()
	#print(f'{x.shape=}, {y.shape=}, {subcols=}, {subrows=}')
	returnI.x = x[int(subcols[0]):int(subcols[1])+1]
	returnI.y = y[int(subrows[0]):int(subrows[1])+1]
	returnI.z = []	

	if not mapInfoOnlyFlag:
		dataset2 = rasterio.open(name)
		if dataset2 is None:
			raise Exception("tiff file opened improperly. readGeotiff.py, dataset2")
		returnI.z = dataset2.read(1)
		#print(f'b{returnI.z.shape=}')
		returnI.z = np.array(returnI.z)
		#print(f'c{returnI.z.shape=}')
		#print(f'{subrows=}, {subcols=}')
		returnI.z = returnI.z[subrows[0]:subrows[1]+1, subcols[0]:subcols[1]+1]
		#print(f'd{returnI.z.shape=}')

	returnI.info = info
	returnI.Tinfo = dataset	
	returnI.z = np.array(returnI.z)
	return returnI
