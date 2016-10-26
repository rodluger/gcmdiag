#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`rect.py` Rectify the pressure grid
-------------------------------------------

Converts a NetCDF file from sigma coordinates to pressure coordinates,
then interpolates all arrays to a rectified pressure grid. The user
can choose between linear interpolation (slow) and cubic interpolation
(really, really slow).

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .constants import PREF
import matplotlib.pyplot as pl
import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d
from scipy.io import netcdf
import sys

__all__ = ['Rectify']

def Rectify(file, interpolation = 'linear'):
  '''
  
  '''
  
  f = netcdf.netcdf_file(file)
  varnames = [k for k in f.variables.keys() if f.variables[k].dimensions == ('time', 'pfull', 'lat', 'lon')]
  
  # Create the true irregular pressure grid from the sigma coordinate
  # This is a 4D array (time, pressure, lat, lon). Units are mb.
  print("Computing new pressure grid...")
  pfull_ = np.stack([f.variables['ps'][:] * p / PREF for p in f.variables['pfull'][:]], axis = 1)

  # Rectify pfull
  pfull = np.nanmedian(pfull_, axis = (0,2,3))

  # Rectify each of the 4D variables
  newvars = {}
  for n, name in enumerate(varnames):
    var = np.array(f.variables[name][:])
    for i, _ in enumerate(f.variables['time'][:]):
      sys.stdout.write('\rVariable %d/%d: Processing time %d/%d...' % (n + 1, len(varnames), i + 1, len(f.variables['time'][:])))
      sys.stdout.flush()
      for k, _ in enumerate(f.variables['lat'][:]):
        for l, _ in enumerate(f.variables['lon'][:]):
          var[i,:,k,l] = interp1d(pfull_[i,:,k,l], var[i,:,k,l], kind=interpolation, bounds_error=False)(pfull)  
    newvars.update({name: var})
    
  # Create a new NetCDF file with just the rectified grids
  outfile = '%s.rect.nc' % file[:-3]
  print("\nSaving to `%s`..." % outfile)
  fnew = netcdf.netcdf_file(outfile, mode = 'w')
  
  # Create the dimensions and the variables
  for name, length in dict(f.dimensions).items(): 
    fnew.createDimension(name, length)
  for name in f.variables.keys():
    var = fnew.createVariable(name, f.variables[name].data.dtype, f.variables[name].dimensions)
    if name in varnames:
      var[:] = newvars[name]
    else:
      var[:] = f.variables[name][:]
    for atr, val in f.variables[name]._attributes.items():
      setattr(var, atr, val)
  
  # Update the pressure grid
  fnew.variables['pfull'][:] = pfull
  fnew.variables['pfull'].long_name = 'True pressure'
  
  # Close
  fnew.close()
  f.close()
  print("Done!")