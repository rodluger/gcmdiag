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
  
  def __new__(cls, input_array, unit = None, name = None, dims = None, desc = None, burnin = 0.):
    # Input array is an already formed ndarray instance
    # We first cast to be our class type
    obj = np.asarray(input_array).view(cls)
    # add the new attribute to the created instance
    obj.unit = unit
    obj.name = name
    obj.desc = desc
    obj.dims = dims
    obj.burnin = burnin
    # Finally, we must return the newly created object:
    return obj

  def __array_finalize__(self, obj):
    # see InfoArray.__array_finalize__ for comments
    if obj is None: return
    self.unit = getattr(obj, 'unit', None)
    self.name = getattr(obj, 'name', None)
    self.desc = getattr(obj, 'desc', None)
    self.dims = getattr(obj, 'dims', None)
    self.burnin = getattr(obj, 'burnin', 0.)
  
  def __array_wrap__(self, out_arr, context = None):
    # Call the parent
    return np.ndarray.__array_wrap__(self, out_arr, context)
  
  def avg(self, *axes):
    '''
    Returns the mean along one or more of the axes, which
    are input as strings. If taking the mean along the time
    axis, discards the burn-in.
    
    '''
    
    axis = tuple(np.where([d in axes for d in self.dims])[0])
    newdims = list(self.dims)
    [newdims.remove(a) for a in axes]
    if len(axis):
      if 'time' in axes:
        n = np.argmax(axes == 'time')
        arr = np.delete(self, int(self.burnin * self.shape[n]), axis = n)
      else:
        arr = self
      return array(np.nanmean(arr, axis = axis), unit = self.unit, 
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
    if 'time' in axes:
      n = np.argmax(axes == 'time')
      arr = np.delete(self, int(self.burnin * self.shape[n]), axis = n)
    else:
      arr = self
    return self - array(np.nanmean(arr, axis = axis).reshape(*shape),
                        unit = self.unit, name = self.name, dims = self.dims,
                        desc = self.desc)
  
  def grad(self, x):
    '''
    Computes the gradient along the axis `x` (a 1D `array` instance).
    Based on http://stackoverflow.com/a/25877826.
    
    '''
    
    # Get the axis index
    try:
      axis = np.where([d == x.name for d in self.dims])[0][0]
    except IndexError:
      raise ValueError('No axis named `%s`.' % x.name)
    
    # Reshape x
    shape = [1 for i in range(len(self.shape))]
    shape[axis] = -1
    shape = tuple(shape)
    x = x.reshape(*shape)   
    sz = self.shape[axis]  
     
    # Compute a central difference.
    x0 = np.delete(x, [sz - 1, sz - 2], axis = axis)
    x1 = np.delete(x, [0, sz - 1], axis = axis)
    x2 = np.delete(x, [0, 1], axis = axis)
    y0 = np.delete(self, [sz - 1, sz - 2], axis = axis)
    y1 = np.delete(self, [0, sz - 1], axis = axis)
    y2 = np.delete(self, [0, 1], axis = axis)
    f = (x2 - x1) / (x2 - x0)
    result = (1 - f) * (y2 - y1) / (x2 - x1) + f * (y1 - y0) / (x1 - x0)
    
    # Add zeros to the two edges for simplicity
    eshape = list(result.shape)
    eshape[axis] = 1
    eshape = tuple(eshape)
    edge = np.zeros(eshape)
    result = np.concatenate([edge, result, edge], axis = axis)
    
    return array(result, unit = '', name = 'grad(%s)' % self.name, dims = self.dims, desc = 'grad(%s)' % self.name)

  def integral(self, x):
    '''
    Returns the integral along the axis `x` (an `array` instance).
    
    '''
    
    newdims = list(self.dims)
    newdims.remove(x.name)
    try:
      axis = np.where([d == x.name for d in self.dims])[0][0]
    except IndexError:
      raise ValueError('No axis named `%s`.' % x.name)
    res = np.trapz(self, x, axis = axis)
    return array(res, name = 'int(%s)' % self.name, unit = '', dims = newdims, desc = 'int(%s)' % self.name)
  
  @property
  def label(self):
    '''
    
    '''
    
    if len(self.unit):
      return '%s (%s)' % (self.desc, self.unit)
    else:
      return self.desc
