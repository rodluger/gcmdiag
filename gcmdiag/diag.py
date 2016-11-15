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

__all__ = ['NetCDF', 'GCMOutput']

class GCMOutput(object):
  '''
  
  '''
  
  def __init__(self, *args, **kwargs):
    '''
    
    '''
    
    raise NotImplementedError('Not a user-facing class!')
  
  def streamfunction(self):
    '''
    Returns the time and zonal mean streamfunction Psi
    
    :returns ndarray Z: A :py:obj:`(pfull, lat)`-shaped array.
    
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
  
  def angular_momentum_flux(self, eddy = False):
    '''
    Returns the instantaneous angular momentum flux (eddy or total)
    
    :returns ndarray Z: A :py:obj:`(time, pfull, lat, lon)`-shaped array.
    
    '''
    
    if eddy:
      uvp = self.vcomp.prime('time', 'lon') * self.ucomp.prime('time', 'lon')
    else:
      uvp = self.vcomp * self.ucomp
    
    # Get the cosine of the latitude in the correct shape
    shape = np.ones(len(uvp.shape), dtype = int)
    shape[np.argmax(np.array(uvp.dims) == 'lat')] = -1
    shape = tuple(shape)
    coslat = np.cos(self.lat * np.pi / 180.).reshape(*shape)
    
    F = uvp * REARTH * coslat
    F.name = '%sangular_momentum_flux' % ('eddy_' if eddy else '')
    F.desc = '%sangular momentum flux' % ('eddy ' if eddy else '')
    F.unit = 'm^3 / s^2'
    
    return F
    
  def heat_flux(self, eddy = False, zero_mass_flux = True):
    '''
    Returns the instantaneous heat flux (eddy or total). If :py:obj:`eddy` 
    is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, subtracts off a
    pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    :returns ndarray Z: A :py:obj:`(time, pfull, lat, lon)`-shaped array.
    
    '''
  
    if eddy:
      vprime = self.vcomp.prime('time', 'lon')
      tprime = self.temp.prime('time', 'lon')
      HF = vprime * CPAIR * tprime
    else:
      if zero_mass_flux:
        # First we subtract off a pressure-integrated mean velocity. This makes sure
        # the net mass flux across each latitude circle is zero, as it should be for 
        # the time-mean flow on a season-less planet.
        vmu = self.vcomp.integral(self.pfull / self.pfull[-1]).avg('time', 'lon').reshape(1,1,-1,1)
      else:
        vmu = 0.
      HF = (self.vcomp - vmu) * CPAIR * self.temp
    
    HF.name = '%sheat_flux' % ('eddy_' if eddy else '')
    HF.desc = '%sheat flux' % ('eddy ' if eddy else '')
    HF.unit = 'm^3 / s^3'
    return HF
  
  def geopotential_energy_flux(self, eddy = False, zero_mass_flux = True):
    '''
    Returns the instantaneous geopotential energy flux (eddy or total). If :py:obj:`eddy` 
    is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, subtracts off a
    pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    :returns ndarray Z: A :py:obj:`(time, pfull, lat, lon)`-shaped array.
    
    '''
    
    if eddy:
      vprime = self.vcomp.prime('time', 'lon')
      zprime = self.hght.prime('time', 'lon')
      GEF = vprime * GRAV * zprime
    else:
      if zero_mass_flux:
        # First we subtract off a pressure-integrated mean velocity. This makes sure
        # the net mass flux across each latitude circle is zero, as it should be for 
        # the time-mean flow on a season-less planet.
        vmu = self.vcomp.integral(self.pfull / self.pfull[-1]).avg('time', 'lon')
        # Now we need to reshape it
        shape = np.ones(len(self.vcomp.shape), dtype = int)
        shape[np.argmax(np.array(self.vcomp.dims) == 'lat')] = -1
        shape = tuple(shape)
        vmu = vmu.reshape(*shape)
      else:
        vmu = 0.
      GEF = (self.vcomp - vmu) * GRAV * self.hght

    GEF.name = '%sgeopotential_energy_flux' % ('eddy_' if eddy else '')
    GEF.desc = '%sgeopotential energy flux' % ('eddy ' if eddy else '')    
    GEF.unit = 'm^3 / s^3'
    return GEF
    
  def angular_momentum(self, relative = True):
    '''
    Returns the total or relative instantaneous angular momentum of the atmosphere
    
    :returns ndarray Z: A :py:obj:`(time, pfull, lat, lon)`-shaped array.
    
    '''
      
    # Get the cosine of the latitude in the correct shape
    shape = np.ones(len(self.ucomp.shape), dtype = int)
    shape[np.argmax(np.array(self.ucomp.dims) == 'lat')] = -1
    shape = tuple(shape)
    coslat = np.cos(self.lat * np.pi / 180.).reshape(*shape)
    
    # Relative angular momentum
    M = self.ucomp * coslat
    
    if not relative:
      # Add the rotational component to get the total
      M += OMEGA * REARTH ** 2 * coslat ** 2
    
    M.name = '%s_angular_momentum' % ('relative' if relative else 'total')
    M.desc = '%s angular momentum' % ('relative' if relative else 'total')
    M.unit = 'm^2 / s'
    return M
  
  def toa_imbalance(self):
    '''
    Returns the instantaneous top-of-atmosphere radiative imbalance (DSW - ULW)
    
    :returns ndarray Z: A :py:obj:`(time, lat, lon)`-shaped array.
    
    '''
    
    # We need to reshape the downwards sw flux
    shape = list(self.olr.shape)
    shape[np.argmax(np.array(self.olr.dims) == 'time')] = 1
    shape = tuple(shape)
    swdn = self.swdn_toa[0].reshape(*shape)
    
    # Compute the difference 
    TOA = swdn - self.olr
    TOA.dims = self.olr.dims
    TOA.name = 'toa_imbalance'
    TOA.desc = 'TOA imbalance'
    return TOA
  
  def dry_static_energy_flux(self, eddy = False, zero_mass_flux = True):
    '''
    Returns the instantaneous dry static energy flux (eddy or total). If :py:obj:`eddy` 
    is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, subtracts off a
    pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    :returns ndarray Z: A :py:obj:`(time, pfull, lat, lon)`-shaped array.
    
    '''
    
    if eddy:
      vprime = self.vcomp.prime('time', 'lon')
      tprime = self.temp.prime('time', 'lon')
      zprime = self.hght.prime('time', 'lon')
      DSE = vprime * (CPAIR * tprime + GRAV * zprime)
    else:
      if zero_mass_flux:
        # First we subtract off a pressure-integrated mean velocity. This makes sure
        # the net mass flux across each latitude circle is zero, as it should be for 
        # the time-mean flow on a season-less planet.
        vmu = self.vcomp.integral(self.pfull / self.pfull[-1]).avg('time', 'lon')
        # Now we need to reshape it
        shape = np.ones(len(self.vcomp.shape), dtype = int)
        shape[np.argmax(np.array(self.vcomp.dims) == 'lat')] = -1
        shape = tuple(shape)
        vmu = vmu.reshape(*shape)
      else:
        vmu = 0.
      DSE = (self.vcomp - vmu) * (CPAIR * self.temp + GRAV * self.hght)

    DSE.name = '%sdry_static_energy_flux' % ('eddy_' if eddy else '')
    DSE.desc = '%sdry static energy flux' % ('eddy ' if eddy else '')    
    DSE.unit = 'm^3 / s^3'
    return DSE

  def latent_heat_flux(self, eddy = False, zero_mass_flux = True):
    '''
    Returns the instantaneous latent heat flux (eddy or total). If :py:obj:`eddy` 
    is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, subtracts off a
    pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    :returns ndarray Z: A :py:obj:`(time, pfull, lat, lon)`-shaped array.
    
    '''
    
    if eddy:
      vprime = self.vcomp.prime('time', 'lon')
      qprime = self.sphum.prime('time', 'lon')
      LH = HLV * vprime * qprime
    else:    
      if zero_mass_flux:
        # First we subtract off a pressure-integrated mean velocity. This makes sure
        # the net mass flux across each latitude circle is zero, as it should be for 
        # the time-mean flow on a season-less planet.
        vmu = self.vcomp.integral(self.pfull / self.pfull[-1]).avg('time', 'lon')
        # Now we need to reshape it
        shape = np.ones(len(self.vcomp.shape), dtype = int)
        shape[np.argmax(np.array(self.vcomp.dims) == 'lat')] = -1
        shape = tuple(shape)
        vmu = vmu.reshape(*shape)
      else:
        vmu = 0.
      LH = HLV * (self.vcomp - vmu) * self.sphum
    
    LH.name = '%slatent_heat_flux' % ('eddy_' if eddy else '')
    LH.desc = '%slatent heat flux' % ('eddy ' if eddy else '')
    LH.unit = 'm^3 / s^3'
    return LH

  def moist_static_energy_flux(self, eddy = False, zero_mass_flux = True):
    '''
    Returns the instantaneous moist static energy flux (eddy or total). If :py:obj:`eddy` 
    is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, subtracts off a
    pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    :returns ndarray Z: A :py:obj:`(time, pfull, lat, lon)`-shaped array.
    
    '''

    MSE = self.latent_heat_flux(eddy = eddy, zero_mass_flux = zero_mass_flux) + \
          self.dry_static_energy_flux(eddy = eddy, zero_mass_flux = zero_mass_flux)
    MSE.name = '%smoist_static_energy_flux' % ('eddy_' if eddy else '')
    MSE.desc = '%smoist static energy flux' % ('eddy ' if eddy else '')
    MSE.unit = 'm^3 / s^3'
    return MSE
    
  def atmospheric_energy_flux(self, eddy = False, zero_mass_flux = True, dry = True, latent = True):
    '''
    Returns the instantaneous atmospheric energy flux, defined in equation (9.8) of
    Pierrehumbert's "Principles of Planetary Climate." If :py:obj:`eddy` 
    is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, subtracts off a
    pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    :returns ndarray Z: A :py:obj:`(time, lat, lon)`-shaped array.
    
    '''
    
    # This is the instantaneous vertically-integrated meridional moist static energy flux
    if dry and latent:
      Z = self.moist_static_energy_flux(eddy = eddy, zero_mass_flux = zero_mass_flux)
    elif dry:
      Z = self.dry_static_energy_flux(eddy = eddy, zero_mass_flux = zero_mass_flux)
    elif latent:
      Z = self.latent_heat_flux(eddy = eddy, zero_mass_flux = zero_mass_flux)
    else:
      raise ValueError('At least one of `dry`, `latent` must be `True`.')
    phi = (1. / REARTH) * Z.integral(self.pfull * 100. / GRAV)   
    phi.name = '%satmospheric_energy_flux' % ('eddy_' if eddy else '')
    phi.desc = '%satmospheric energy flux' % ('eddy ' if eddy else '')
    phi.unit = 'W / m^2'
    return phi
    
  def energy_flux_gradient(self, eddy = False, zero_mass_flux = True, dry = True, latent = True):
    '''
    Returns the latitudinal gradient of the atmospheric energy flux (eddy or total). More 
    specifically, this is the convergence of the vertically-integrated time and zonal 
    average meridional flux of moist static energy. Check out equation (21) in Trenberth 
    (1997) and equation (9.8) in Pierrehumbert's "Principles of Planetary Climate."
    If :py:obj:`eddy` is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, 
    subtracts off a pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    :returns ndarray Z: A :py:obj:`(lat)`-shaped array.
    
    '''
    
    # Vertically-integrated meridional flux (time and zonal average)
    lat = self.lat * np.pi / 180.
    z = self.atmospheric_energy_flux(eddy = eddy, zero_mass_flux = zero_mass_flux, 
                                     dry = dry, latent = latent).avg('time', 'lon') * np.cos(lat)
    cimf = (1. / np.cos(lat)) * z.grad(lat)
    cimf.name = '%senergy_flux_gradient' % ('eddy_' if eddy else '')
    cimf.desc = 'time-mean, zonal-mean convergence of the vertically-integrated ' + \
                '%smeridional flux of moist static energy' % ('eddy ' if eddy else '')
    cimf.unit = 'W / m^2'
    return cimf

  def net_vgz(self, eddy = False, zero_mass_flux = True):
    '''
    The vertically-integrated geopotential energy flux.
    If :py:obj:`eddy` is :py:obj:`False` and :py:obj:`zero_mass_flux` is :py:obj:`True`, 
    subtracts off a pressure-integrated mean velocity to zero out the latitudinal mass flux.
    
    '''
    
    vgz = self.geopotential_energy_flux(eddy = eddy, zero_mass_flux = zero_mass_flux).integral(self.pfull * 100.)
    vgz.name = 'integral_%svbar_zbar' % ('eddy_' if eddy else '')
    vgz.desc = 'vertically-integrated %svbar * g * zbar' % ('eddy ' if eddy else '')
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