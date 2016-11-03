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
from .utils import array
from .constants import PREF
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import netcdf
import sys

__all__ = ['Rectify']

def Rectify(file, interpolation = 'linear', rect_file = None):
  '''
  
  '''
  
  # Open the file
  f = netcdf.netcdf_file(file)
  
  # Get the variables we need to modify
  pfvars = []
  phvars = []
  newvars = {}
  for name in list(f.variables.keys()):
    if f.variables[name].dimensions == ('time', 'pfull', 'lat', 'lon'):
      pfvars.append(name)
    elif f.variables[name].dimensions == ('time', 'phalf', 'lat', 'lon'):
      phvars.append(name)
    
  # Create the true irregular pressure grid from the sigma coordinate
  # This is a 4D array (time, pressure, lat, lon). Units are mb.
  print("Computing new pressure grid...")
  pfull_ = np.stack([f.variables['ps'][:] * p / PREF for p in f.variables['pfull'][:]], axis = 1)
  phalf_ = np.stack([f.variables['ps'][:] * p / PREF for p in f.variables['phalf'][:]], axis = 1)

  # Rectify pfull
  pfull = np.nanmedian(pfull_, axis = (0,2,3))
  phalf = np.nanmedian(phalf_, axis = (0,2,3))
    
  # Get the dimensions
  ntime = f.variables['time'].shape[0]
  nlat = f.variables['lat'].shape[0]
  nlon = f.variables['lon'].shape[0]
    
  # Rectify each of the 4D phalf variables and interpolate to the pfull grid
  newvars.update({'phalf': phalf, 'phalf_unit': np.array('mb'), 'phalf_dims': np.array(['phalf']), 'phalf_desc': np.array('True half-pressure grid')})
  for n, name in enumerate(phvars):
    var_old = np.array(f.variables[name][:])
    var_new = np.zeros((ntime, len(pfull), nlat, nlon))
    for i in range(ntime):
      sys.stdout.write('\rVariable %d/%d: Processing time %d/%d...' % (n + 1, len(pfvars) + len(phvars), i + 1, ntime))
      sys.stdout.flush()
      for k in range(nlat):
        for l in range(nlon):
          var_new[i,:,k,l] = interp1d(phalf_[i,:,k,l], var_old[i,:,k,l], kind=interpolation, 
                                      bounds_error=False, fill_value="extrapolate")(pfull)  
    newvars.update({name: var_new})
    newvars.update({name + '_unit': np.array(f.variables[name].units.decode('utf-8'))})
    newvars.update({name + '_dims': np.array(['time', 'pfull', 'lat', 'lon'])})
    newvars.update({name + '_desc': np.array(f.variables[name].long_name.decode('utf-8'))})  
    
  # Rectify each of the 4D pfull variables
  newvars.update({'pfull': pfull, 'pfull_unit': np.array('mb'), 'pfull_dims': np.array(['pfull']), 'pfull_desc': np.array('True pressure grid')})
  for n, name in enumerate(pfvars):
    var = np.array(f.variables[name][:])
    for i in range(ntime):
      sys.stdout.write('\rVariable %d/%d: Processing time %d/%d...' % (len(phvars) + n + 1, len(pfvars) + len(phvars), i + 1, ntime))
      sys.stdout.flush()
      for k in range(nlat):
        for l in range(nlon):
          var[i,:,k,l] = interp1d(pfull_[i,:,k,l], var[i,:,k,l], kind=interpolation, 
                                  bounds_error=False, fill_value="extrapolate")(pfull)  
    newvars.update({name: var})
    newvars.update({name + '_unit': np.array(f.variables[name].units.decode('utf-8'))})
    newvars.update({name + '_dims': np.array(['time', 'pfull', 'lat', 'lon'])})
    newvars.update({name + '_desc': np.array(f.variables[name].long_name.decode('utf-8'))})
    
  # Close original file
  f.close()
  
  # Create a new NetCDF file with just the rectified grids
  if rect_file is None:
    rect_file = '%s.rect.npz' % file[:-3]
  print("\nSaving to `%s`..." % rect_file)
  np.savez(rect_file, **newvars)
  
  f.close()
  print("Done!")