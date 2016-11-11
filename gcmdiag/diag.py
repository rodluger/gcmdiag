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
  def streamfunction(self):
    '''
    Returns the streamfunction Psi averaged over time and longitude
    
    '''
    
    vbar = self.vcomp.avg('time', 'lon')
    cosTheta = np.cos(self.lat * np.pi / 180.).reshape(1,-1)
    integrand = vbar * cosTheta
    psi = np.zeros((len(self.pfull), len(self.lat)))
    for i in range(len(self.pfull)):
      psi[i] = -np.trapz(integrand[:i], self.pfull[:i] * 100., axis = 0)
    psi = array(psi)
    psi.name = 'streamfunction'
    psi.desc = 'streamfunction'
    psi.unit = 'kg / s^3'
    return psi
  
  @property
  def eddy_angular_momentum_flux(self):
    '''
    Returns the time and zonal mean eddy angular momentum flux
  
    '''
  
    uvp = self.vcomp.prime('time', 'lon') * self.ucomp.prime('time', 'lon')
    F = uvp.avg('time', 'lon') * REARTH * np.cos(self.lat * np.pi / 180.).reshape(1, -1)
    F.name = 'eddy_angular_momentum_flux'
    F.desc = 'eddy angular momentum flux'
    F.unit = 'm^3 / s^2'
    return F
    
  @property
  def eddy_heat_flux(self):
    '''
    Returns the time and zonal mean meridional eddy heat flux
  
    '''
  
    F = (CPAIR * self.vcomp.prime('time', 'lon') * self.temp.prime('time', 'lon')).avg('time', 'lon')
    F.name = 'eddy_heat_flux'
    F.desc = 'meridional eddy heat flux'
    F.unit = 'K m/s'
    return F
  
  @property
  def relative_angular_momentum(self):
    '''
    Returns the relative angular momentum of the atmosphere

    '''
      
    # Relative angular momentum
    M = self.ucomp.avg('time', 'lon') * REARTH * np.cos(self.lat * np.pi / 180.).reshape(1,-1)
    M.name = 'relative_angular_momentum'
    M.desc = 'relative angular momentum'
    M.unit = 'm^2 / s'
    return M

  @property
  def total_angular_momentum(self):
    '''
    Returns the total angular momentum of the atmosphere
    
    '''

    M = self.relangmom + OMEGA * REARTH ** 2 * np.cos(self.lat * np.pi / 180.).reshape(1,-1) ** 2
    M.name = 'total_angular_momentum'
    M.desc = 'total angular momentum'
    M.unit = 'm^2 / s'
    return M
  
  @property
  def toa_imbalance(self):
    '''
    Returns the top-of-atmosphere radiative imbalance (DSW - ULW)
    
    '''
    swdn = self.swdn_toa[0]
    olr = self.olr.avg('time')
    TOA = swdn - olr
    TOA.dims = olr.dims
    TOA.name = 'toa_imbalance'
    TOA.desc = 'TOA imbalance'
    return TOA
  
  @property
  def eddy_dry_static_energy_flux(self):
    '''
    Returns the instantaneous eddy dry static energy flux
    
    '''
    
    vprime = self.vcomp.prime('time', 'lon')
    tprime = self.temp.prime('time', 'lon')
    zprime = self.hght.prime('time', 'lon')
    
    DSE = vprime * (CPAIR * tprime + GRAV * zprime)
    DSE.name = 'eddy_dry_static_energy_flux'
    DSE.desc = 'eddy dry static energy flux'
    DSE.unit = 'm^3 / s^3'
    return DSE

  @property
  def eddy_latent_heat_flux(self):
    '''
    Returns the instantaneous eddy latent heat flux
    
    '''
    
    vprime = self.vcomp.prime('time', 'lon')
    qprime = self.sphum.prime('time', 'lon')
    
    LH = HLV * vprime * qprime
    LH.name = 'eddy_latent_heat_flux'
    LH.desc = 'eddy latent heat flux'
    LH.unit = 'm^3 / s^3'
    return LH

  @property
  def eddy_moist_static_energy_flux(self):
    '''
    Returns the instantaneous eddy moist static energy flux
    
    '''

    MSE = self.eddy_latent_heat_flux + self.eddy_dry_static_energy_flux
    MSE.name = 'eddy_moist_static_energy_flux'
    MSE.desc = 'eddy moist static energy flux'
    MSE.unit = 'm^3 / s^3'
    return MSE

  @property
  def dry_static_energy_flux(self):
    '''
    Returns the instantaneous dry static energy flux
    
    '''
    
    DSE = self.vcomp * (CPAIR * self.temp + GRAV * self.hght)
    DSE.name = 'dry_static_energy_flux'
    DSE.desc = 'dry static energy flux'
    DSE.unit = 'm^3 / s^3'
    return DSE

  @property
  def latent_heat_flux(self):
    '''
    Returns the instantaneous latent heat flux
    
    '''
    
    LH = HLV * self.vcomp * self.sphum
    LH.name = 'latent_heat_flux'
    LH.desc = 'latent heat flux'
    LH.unit = 'm^3 / s^3'
    return LH

  @property
  def moist_static_energy_flux(self):
    '''
    Returns the instantaneous moist static energy flux
    
    '''

    MSE = self.latent_heat_flux + self.dry_static_energy_flux
    MSE.name = 'moist_static_energy_flux'
    MSE.desc = 'moist static energy flux'
    MSE.unit = 'm^3 / s^3'
    return MSE
  
  @property
  def atmospheric_energy_flux(self):
    '''
    Returns the atmospheric energy flux, defined in equation (9.8) of
    Pierrehumbert's "Principles of Planetary Climate."
    
    '''
    
    # This is the vertically-integrated meridional moist static energy flux (time and zonal average)
    phi = (1. / REARTH) * self.moist_static_energy_flux.integral(self.pfull * 100. / GRAV).avg('time', 'lon')    
    phi.name = 'atmospheric_energy_flux'
    phi.desc = 'atmospheric energy flux'
    phi.unit = 'W / m^2'
    return phi

  @property
  def eddy_atmospheric_energy_flux(self):
    '''
    Returns the eddy atmospheric energy flux, similar to equation (9.8) of
    Pierrehumbert's "Principles of Planetary Climate."
    
    '''
    
    # This is the vertically-integrated meridional moist static energy flux (time and zonal average)
    phi = (1. / REARTH) * self.eddy_moist_static_energy_flux.integral(self.pfull * 100. / GRAV).avg('time', 'lon')
    phi.name = 'eddy_atmospheric_energy_flux'
    phi.desc = 'eddy atmospheric energy flux'
    phi.unit = 'W / m^2'
    return phi
    
  @property
  def energy_flux_gradient(self):
    '''
    Returns the latitudinal gradient of the atmospheric energy flux. More specifically, this is the convergence of the 
    vertically-integrated time and zonal average meridional flux of moist static energy. 
    Check out equation (21) in Trenberth (1997) and equation (9.8) in Pierrehumbert's "Principles of Planetary Climate."
    
    '''
    
    # Vertically-integrated meridional flux (time and zonal average)
    lat = self.lat * np.pi / 180.
    z = self.atmospheric_energy_flux * np.cos(lat)
    cimf = (1. / np.cos(lat)) * z.grad(lat)
    
    # DEBUG
    cimf -= self.vcomp.integral(self.pfull * 100.).avg('time', 'lon')
    # /DEBUG
    
    
    cimf.name = 'energy_flux_gradient'
    cimf.desc = 'time-mean, zonal-mean convergence of the vertically-integrated meridional flux of moist static energy'
    cimf.unit = 'W / m^2'
    return cimf

  @property
  def eddy_energy_flux_gradient(self):
    '''
    Returns the latitudinal gradient of the eddy atmospheric energy flux. More specifically, this is the convergence of the 
    vertically-integrated time and zonal average meridional eddy flux of moist static energy. 
    Check out equation (21) in Trenberth (1997) and equation (9.8) in Pierrehumbert's "Principles of Planetary Climate."
    
    '''
    
    # Vertically-integrated meridional flux (time and zonal average)
    lat = self.lat * np.pi / 180.
    z = self.eddy_atmospheric_energy_flux * np.cos(lat)
    cimf = (1. / np.cos(lat)) * z.grad(lat)
    cimf.name = 'energy_flux_gradient'
    cimf.desc = 'time-mean, zonal-mean convergence of the vertically-integrated meridional eddy flux of moist static energy'
    cimf.unit = 'W / m^2'
    return cimf

  @property
  def integral_vbar_zbar(self):
    '''
    The vertically-integrated ``vbar * g * zbar``
    
    '''
    
    vgz = self.vcomp.avg('time', 'lon') * GRAV * self.hght.avg('time', 'lon')
    vgz = vgz.integral(self.pfull * 100.)
    vgz.name = 'integral_vbar_zbar'
    vgz.desc = 'vertically-integrated vbar * g * zbar'
    vgz.unit = 'm^2 kg / s^5'
    return vgz

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