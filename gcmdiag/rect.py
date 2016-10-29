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
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import netcdf
import sys

__all__ = ['Rectify']

def Rectify(file, interpolation = 'linear', debug = False, outfile = None):
  '''
  
  '''
  
  # Open the file
  f = netcdf.netcdf_file(file)
  
  # Get the variables we need to modify
  pres_variables = {}
  for var in list(f.variables.keys()):
    if f.variables[var].dimensions == ('time', 'pfull', 'lat', 'lon'):
      pres_variables.update({var: [f.variables[var][:].dtype, dict(f.variables[var]._attributes)]})
   
  if not debug: 
    # Create the true irregular pressure grid from the sigma coordinate
    # This is a 4D array (time, pressure, lat, lon). Units are mb.
    print("Computing new pressure grid...")
    pfull_ = np.stack([f.variables['ps'][:] * p / PREF for p in f.variables['pfull'][:]], axis = 1)

    # Rectify pfull
    pfull = np.nanmedian(pfull_, axis = (0,2,3))
  
  # Get the other variables
  other_variables = {}
  for var in list(f.variables.keys()):
    if var not in pres_variables.keys():
      other_variables.update({var: [f.variables[var][:].dtype, list(f.variables[var].dimensions), 
                                    dict(f.variables[var]._attributes), np.array(f.variables[var][:])]})
  
  # Get the dimensions
  dimensions = dict(f.dimensions)
  ntime = f.variables['time'].shape[0]
  dimensions['time'] = ntime
  nlat = f.variables['lat'].shape[0]
  nlon = f.variables['lon'].shape[0]
    
  # Rectify each of the 4D variables
  newvars = {}
  for n, name in enumerate(pres_variables):
    var = np.array(f.variables[name][:])
    if not debug:
      for i in range(ntime):
        sys.stdout.write('\rVariable %d/%d: Processing time %d/%d...' % (n + 1, len(pres_variables), i + 1, ntime))
        sys.stdout.flush()
        for k in range(nlat):
          for l in range(nlon):
            var[i,:,k,l] = interp1d(pfull_[i,:,k,l], var[i,:,k,l], kind=interpolation, 
                                    bounds_error=False, fill_value="extrapolate")(pfull)  
    newvars.update({name: var})
  
  # Close original file
  f.close()
  
  # Create a new NetCDF file with just the rectified grids
  if outfile is None:
    outfile = '%s.rect.nc' % file[:-3]
  print("\nSaving to `%s`..." % outfile)
  fnew = netcdf.netcdf_file(outfile, mode = 'w')
  
  # Create the dimensions and the variables
  for name, length in dimensions.items(): 
    fnew.createDimension(name, length)
  
  for name in pres_variables.keys():
    var = fnew.createVariable(name, pres_variables[name][0], ('time', 'pfull', 'lat', 'lon'))
    var[:] = np.array(newvars[name])
    for atr, val in pres_variables[name][1].items():
      setattr(var, atr, val)
  for name in other_variables.keys():
    var = fnew.createVariable(name, other_variables[name][0], other_variables[name][1])
    var[:] = other_variables[name][3]
    for atr, val in other_variables[name][2].items():
      setattr(var, atr, val)
  
  if not debug:
    # Update the pressure grid
    fnew.variables['pfull'][:] = pfull
    fnew.variables['pfull'].long_name = 'True pressure'
  
  # Close
  fnew.close()
  f.close()
  print("Done!")