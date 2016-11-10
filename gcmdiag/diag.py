#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`diag.py` Diagnostic Functions
--------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .constants import *
from .utils import array
from .rect import Rectify
from scipy.io import netcdf
import numpy as np
import os

__all__ = ['NetCDF']

class GCMOutput(object):
  '''
  
  '''
  
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    raise NotImplementedError('Not a user-facing class!')
  
  @property
  def streamfunc(self):
    '''
    Returns the streamfunction Psi averaged over time and longitude
    
    '''
    
    vbar = self.vcomp.avg('time', 'lon')
    cosTheta = np.cos(self.lat * np.pi / 180.).reshape(1,-1)
    integrand = vbar * cosTheta
    psi = np.zeros((len(self.pfull), len(self.lat)))
    for i in range(len(self.pfull)):
      psi[i] = -np.trapz(integrand[:i], self.pfull[:i], axis = 0)
    return array(psi, name = 'streamfunc', desc = 'stream function', unit = '')
  
  @property
  def eddyangmom(self):
    '''
    Returns the time and zonal mean eddy angular momentum flux
  
    '''
  
    uvp = self.vcomp.prime('time', 'lon') * self.ucomp.prime('time', 'lon')
    F = uvp.avg('time', 'lon') * REARTH * np.cos(self.lat * np.pi / 180.).reshape(1, -1)
    F.name = 'eddyangmom'
    F.desc = 'eddy angular momentum flux'
    F.unit = ''
    return F
    
  @property
  def eddyheat(self):
    '''
    Returns the time and zonal mean meridional eddy heat flux
  
    '''
  
    F = (CPAIR * self.vcomp.prime('time', 'lon') * self.temp.prime('time', 'lon')).avg('time', 'lon')
    F.name = 'eddyheat'
    F.desc = 'meridional eddy heat flux'
    F.unit = 'K m/s'
    return F
  
  @property
  def relangmom(self):
    '''
    Returns the relative angular momentum of the atmosphere

    '''
      
    # Relative angular momentum
    M = self.ucomp.avg('time', 'lon') * REARTH * np.cos(self.lat * np.pi / 180.).reshape(1,-1)
    M.name = 'relangmom'
    M.desc = 'relative angular momentum'
    M.unit = ''
    return M

  @property
  def totalangmom(self):
    '''
    Returns the total angular momentum of the atmosphere
    
    '''

    M = self.relangmom + OMEGA * REARTH ** 2 * np.cos(self.lat * np.pi / 180.).reshape(1,-1) ** 2
    M.name = 'totalangmom'
    M.desc = 'total angular momentum'
    M.unit = ''
    return M
  
  @property
  def toaimbalance(self):
    '''
    Returns the top-of-atmosphere radiative imbalance (DSW - ULW)
    
    '''
    swdn = self.swdn_toa[0]
    olr = self.olr.avg('time')
    TOA = swdn - olr
    TOA.dims = olr.dims
    TOA.name = 'toaimbalance'
    TOA.desc = 'TOA imbalance'
    return TOA
  
  @property
  def eddydrystaticenergy(self):
    '''
    Returns the instantaneous eddy dry static energy flux
    
    '''
    
    vprime = self.vcomp.prime('time', 'lon')
    tprime = self.temp.prime('time', 'lon')
    zprime = self.hght.prime('time', 'lon')
    
    DSE = vprime * (CPAIR * tprime + GRAV * zprime)
    DSE.name = 'eddydrystaticenergy'
    DSE.desc = 'eddy dry static energy flux'
    DSE.unit = ''
    return DSE

  @property
  def eddylatentheat(self):
    '''
    Returns the instantaneous eddy latent heat flux
    
    '''
    
    vprime = self.vcomp.prime('time', 'lon')
    qprime = self.sphum.prime('time', 'lon')
    
    LH = HLV * vprime * qprime
    LH.name = 'eddylatentheat'
    LH.desc = 'eddy latent heat flux'
    LH.unit = ''
    return LH

  @property
  def eddymoiststaticenergy(self):
    '''
    Returns the instantaneous eddy moist static energy flux
    
    '''

    MSE = self.eddylatentheat + self.eddydrystaticenergy
    MSE.name = 'eddymoiststaticenergy'
    MSE.desc = 'eddy moist static energy flux'
    MSE.unit = ''
    return MSE

  @property
  def drystaticenergy(self):
    '''
    Returns the instantaneous dry static energy flux
    
    '''
    
    DSE = self.vcomp * (CPAIR * self.temp + GRAV * self.hght)
    DSE.name = 'drystaticenergy'
    DSE.desc = 'dry static energy flux'
    DSE.unit = ''
    return DSE

  @property
  def latentheat(self):
    '''
    Returns the instantaneous latent heat flux
    
    '''
    
    LH = HLV * self.vcomp * self.sphum
    LH.name = 'latentheat'
    LH.desc = 'latent heat flux'
    LH.unit = ''
    return LH

  @property
  def moiststaticenergy(self):
    '''
    Returns the instantaneous moist static energy flux
    
    '''

    MSE = self.latentheat + self.drystaticenergy
    MSE.name = 'moiststaticenergy'
    MSE.desc = 'moist static energy flux'
    MSE.unit = ''
    return MSE

  @property
  def convintmeridflux(self):
    '''
    Returns the convergence of the vertically-integrated time and zonal average meridional 
    flux of moist static energy
    
    '''
    
    # Vertically-integrated meridional flux (time and zonal average)
    z = (1. / GRAV) * self.moiststaticenergy.avg('time', 'lon').integral(self.pfull * 100.)
    cimf = z.grad(REARTH * self.lat * np.pi / 180.)
    cimf.name = 'convintmeridflux'
    cimf.desc = 'time-mean, zonal-mean convergence of the vertically-integrated meridional flux of moist static energy'
    cimf.unit = ''
    return cimf

class NetCDF(GCMOutput):
  '''
  A smart netcdf data container
  
  '''

  def __init__(self, file, rectify = False, interpolation = 'linear', 
               avg_file = None, rect_file = None, burnin = 0.25):
    '''
    
    '''
    
    # Open file to get variable names
    self.file = file
    f = netcdf.netcdf_file(self.file)
    self._vars = {}
    for var in f.variables.keys():
      self._vars.update({var: None})
    f.close()
    self.burnin = burnin
    
    # If there's an average file, peek at it
    if avg_file is None:
      if os.path.exists(self.file.replace('daily', 'average')):
        avg_file = self.file.replace('daily', 'average')
    self.avg_file = avg_file
    self._avg_vars = {}
    if self.avg_file:
      f = netcdf.netcdf_file(self.avg_file)
      for var in f.variables.keys():
        self._avg_vars.update({var: None})
      f.close()
    
    # Rectify the pressure grid?
    if rectify:
      self.rect_file = '%s.rect.npz' % self.file[:-3]
      if not os.path.exists(self.rect_file):
        Rectify(self.file, interpolation = interpolation, rect_file = rect_file)
      rect = np.load(self.rect_file)
      for var in rect.keys():
        if not (var.endswith('_unit') or var.endswith('_dims') or var.endswith('_desc')):
          setattr(self, var, array(rect[var], unit = rect[var + '_unit'][()], 
                                   name = var, 
                                   dims = tuple(rect[var + '_dims']), 
                                   desc = rect[var + '_desc'][()],
                                   burnin = self.burnin))

  def __getattr__(self, var):  
    '''
    
    '''

    if var in self._vars:
      if self._vars[var] is None:
        f = netcdf.netcdf_file(self.file)
        val = np.array(f.variables[var][:])
        unit = f.variables[var].units.decode('utf-8')
        desc = f.variables[var].long_name.decode('utf-8')
        dims = f.variables[var].dimensions
        self._vars.update({var: array(val, unit = unit, name = var, dims = dims, desc = desc, burnin = self.burnin)})
        f.close()
      return self._vars[var]
    elif var in self._avg_vars:
      if self._avg_vars[var] is None:
        f = netcdf.netcdf_file(self.avg_file)
        val = np.array(f.variables[var][:])
        unit = f.variables[var].units.decode('utf-8')
        desc = f.variables[var].long_name.decode('utf-8')
        dims = f.variables[var].dimensions
        self._avg_vars.update({var: array(val, unit = unit, name = var, dims = dims, desc = desc, burnin = self.burnin)})
        f.close()
      return self._avg_vars[var]
    else:
      return None