#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`utils.py` General utilities
------------------------------------

'''

import numpy as np

__all__ = ['array']

class array(np.ndarray):
  '''
  A custom subclass of numpy ndarray with some extra
  attributes.
  
  '''
  
  def __new__(cls, input_array, unit = None, name = None, dims = None, desc = None):
    # Input array is an already formed ndarray instance
    # We first cast to be our class type
    obj = np.asarray(input_array).view(cls)
    # add the new attribute to the created instance
    obj.unit = unit
    obj.name = name
    obj.desc = desc
    obj.dims = dims
    # Finally, we must return the newly created object:
    return obj

  def __array_finalize__(self, obj):
    # see InfoArray.__array_finalize__ for comments
    if obj is None: return
    self.unit = getattr(obj, 'unit', None)
    self.name = getattr(obj, 'name', None)
    self.desc = getattr(obj, 'desc', None)
    self.dims = getattr(obj, 'dims', None)
  
  def __array_wrap__(self, out_arr, context = None):
    # Call the parent
    return np.ndarray.__array_wrap__(self, out_arr, context)
  
  def avg(self, *axes):
    '''
    Returns the mean along one or more of the axes, which
    are input as strings.
    
    '''
    
    axis = tuple(np.where([d in axes for d in self.dims])[0])
    newdims = list(self.dims)
    [newdims.remove(a) for a in axes]
    if len(axis):
      return array(np.nanmean(self, axis = axis), unit = self.unit, 
                   name = self.name, dims = newdims, desc = self.desc)
    else:
      return self
  
  def prime(self, *axes):
    '''
    Returns the deviations from the mean along the axes `axes`
    
    '''
    
    axis = tuple(np.where([d in axes for d in self.dims])[0])
    shape = list(self.shape)
    for a in axis: 
      shape[a] = 1
    return self - array(np.nanmean(self, axis = axis).reshape(*shape),
                        unit = self.unit, name = self.name, dims = self.dims,
                        desc = self.desc)
  
  def divergence(self):
    '''
    Computes the divergence of the array
  
    '''
    
    return array(np.nansum(np.gradient(self), axis = 0), unit = '', name = 'div(%s)' % self.name, dims = self.dims, desc = 'div(%s)' % self.name)
  
  def integral(self, x):
    '''
    Returns the integral along axis axis `x`.
    
    '''
    
    newdims = list(self.dims)
    newdims.remove(x.name)
    res = np.trapz(self, x, axis = np.where([d == x.name for d in self.dims])[0])
    return array(res, name = 'int(%s)' % self.name, unit = '', dims = newdims, desc = 'int(%s)' % self.name)
  
  @property
  def label(self):
    '''
    
    '''
    
    if len(self.unit):
      return '%s (%s)' % (self.desc, self.unit)
    else:
      return self.desc
