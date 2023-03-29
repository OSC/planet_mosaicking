# Python mosaicking for Planet data

This code is largely ported from [MATLAB code written by Ian Howat](https://github.com/ihowat/setsm_postprocessing "setsm_postprocessing on GitHub"), glaciologist and professor of Earth Sciences at The Ohio State University (OSU). Some helper functions were copied from [Python code for postprocessing setsm by PGC](https://github.com/PolarGeospatialCenter/setsm_postprocessing_python).

## Python requirements
Scripts in this repo are intended to work with Python 3.9+. 

### Python package dependencies
* NumPy
* SciPy
* scikit-image
* OpenCV
* GDAL (including OSGeo, OGR, OSR)
* Shapely
* pickle
* rosterio

These dependencies can be handled with [conda](https://conda.io/docs/index.html "conda landing page"). Once conda is installed, you can create a new environment with
```
conda create --name s2s -c conda-forge python=3 gdal=3 numpy scipy scikit-image opencv shapely rasterio --yes
```

Note: You can install Miniconda [here](https://repo.anaconda.com/miniconda) if your machine does not have a conda version of Python available.

## Running
The file `run.slurm` contains a sample SLURM submission script. You will likely need to modify the path defined by `tileDir` and the file pattern used to define the list of tif files for `fileNames`.
