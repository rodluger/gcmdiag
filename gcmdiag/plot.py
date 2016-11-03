#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`plot.py` Plotting Routines
-----------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import matplotlib.pyplot as pl
import numpy as np

__all__ = ['ColorMap', 'Compare']

def ColorMap(x, y, z, c = None, ax = None, invert_y = True, cax = None, **kwargs):
  '''
  
  '''
  
  if ax is None:
    ax = pl.gca()
  
  # Plot the colormap
  xgrid, ygrid = np.meshgrid(x, y)
  im = ax.pcolormesh(xgrid, ygrid, z, shading='gouraud', **kwargs)
  ax.set_title(z.label, fontsize = 14)
  if cax:
    cb = pl.colorbar(im, cax = cax)
  else:
    cb = pl.colorbar(im, cax = ax)
  
  # Plot a contour of something
  if c is not None:
    cont = ax.contour(xgrid, ygrid, c, 5, colors = 'k', **kwargs)
    pl.clabel(cont, fontsize=9, inline=1, fmt="%d")
  
  # Set limits and stuff
  if invert_y:
    ax.invert_yaxis()
    ax.set_ylim(y.max(), y.min())
  else:
    ax.set_ylim(y.min(), y.max())
  ax.set_xlim(x.min(), x.max())
  ax.set_xlabel(x.label, fontsize = 12)
  ax.set_ylabel(y.label, fontsize = 12)
  
  return ax, cb

def Compare(x, y, z, titles = None, c = None, **kwargs):
  '''
  
  '''
  
  if c is None:
    c = [None for i in z]
  vmin = kwargs.pop('vmin', np.min([zn.min() for zn in z]))
  vmax = kwargs.pop('vmax', np.max([zn.max() for zn in z]))
  fig = pl.figure(figsize = (14, 10))
  fig.subplots_adjust(top = 0.875)
  ax = [pl.subplot2grid((1, 11 * len(z) + 1), (0, 11 * n), colspan = 10) for n in range(len(z))]
  cax = pl.subplot2grid((1, 11 * len(z) + 1), (0, 11 * len(z)))
  for n in range(len(z)):
    ColorMap(x[n], y[n], z[n], ax = ax[n], c = c[n], vmin = vmin, vmax = vmax, cax = cax, **kwargs)
    if titles is not None:
      ax[n].set_title(titles[n])
    if n != 0:
      ax[n].set_ylabel('')
      ax[n].set_yticklabels([])
  pl.suptitle(z[0].label, fontsize = 18)
  return fig, ax