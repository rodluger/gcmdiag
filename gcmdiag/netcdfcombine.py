#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`netcdfcombine.py`
--------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
from scipy.io import netcdf
import os
import glob

def NetCDFCombine(path = '.', outfile = 'combined.nc'):
  '''
  Concatenates the arrays in a set of NetCDF files along the
  `lat` and `latb` dimensions. The input files should be in
  the `path` directory and will be combined into `outfile`.
  
  '''
  
  # Create a new NetCDF file
  newfile = netcdf.netcdf_file(outfile, mode = 'w')

  # Get all files
  files = glob.glob(os.path.join(path, '*.nc.*'))

  # Copy dimension info from the first file to the new file
  # The parallelization takes place over the latitude dimension
  f = netcdf.netcdf_file(files[0])
  # The `time` dimension is `None` for some reason; let's fix this
  f.dimensions['time'] = f.variables['time'].shape[0]
  newdims = dict(f.dimensions)
  newdims['lat'] = newdims['lat'] * len(files)
  newdims['latb'] = newdims['latb'] * len(files) + 1
  for name, length in newdims.items(): 
    newfile.createDimension(name, length)

  # Create a dict with the new variables
  tmpdims = dict(f.dimensions)
  tmpdims['lat'] = 0
  tmpdims['latb'] = 0
  newvars = {}
  for v in f.variables.keys():
    newvars.update({v: np.empty(tuple([tmpdims[d] for d in f.variables[v].dimensions]))})
  f.close()

  # Now populate them with each of the runs
  for file in files:
    f = netcdf.netcdf_file(file)
  
    # Stitch the arrays together
    for v in newvars.keys():
      dims = np.array(f.variables[v].dimensions)
      axis = np.where((dims == 'lat') | (dims == 'latb'))[0]
      if len(axis):
        # This is an array we should concatenate
        newvars[v] = np.concatenate([newvars[v], f.variables[v][:]], axis = axis[0])
      else:
        # This is an array that is the same for all files
        newvars[v] = np.array(f.variables[v][:])
    f.close()

  # Finally, create the new variables
  f = netcdf.netcdf_file(files[0])
  for v in newvars.keys():
    var = newfile.createVariable(v, f.variables[v].data.dtype, f.variables[v].dimensions)
    var[:] = newvars[v]
  
    # Copy over attributes
    for atr, val in f.variables[name]._attributes.items():
      setattr(var, atr, val)
      
  # Close the files    
  f.close()
  newfile.close()