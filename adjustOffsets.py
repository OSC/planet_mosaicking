import numpy as np
import sys

def adjustOffsets(NTiles, n1, n2, dz, dze):
	# remove nan pairs or with large offsets and/or mad values
	input_array = np.concatenate((n1, n2, dz, dze), axis=1)
	print(f'{input_array=}')
	n = np.isnan(input_array).any(axis=1)
	print(f'{n=}')
	sys.exit(1)
	return dZ
