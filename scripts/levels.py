import gcmdiag as gcm
from gcmdiag import GCMDIAG_OUT
from gcmdiag.constants import *
import os
import matplotlib.pyplot as pl
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

for run in ['1bar', 'p5bar', 'p25bar']:

  # Load
  print("Running %s..." % run)
  file = os.path.join(GCMDIAG_OUT, 'N2-%s/day1000h00/day1000h00.atmos_daily.nc' % run)
  data = gcm.NetCDF(file, rectify = True)

  # Compute variables
  lon = np.array(data.lon)
  lat = np.array(data.lat)
  pfull = np.array(data.pfull)
  v = data.vcomp.avg('time')
  u = data.ucomp.avg('time')
  t = data.temp.avg('time')
  mr = u * REARTH * np.cos(data.lat * np.pi / 180.).reshape(1,-1,1)
  mt = mr + OMEGA * REARTH ** 2 * np.cos(data.lat * np.pi / 180.).reshape(1,-1,1) ** 2
  uvp = data.vcomp.prime('time', 'lon') * data.ucomp.prime('time', 'lon')
  em = uvp.avg('time') * REARTH * np.cos(data.lat * np.pi / 180.).reshape(1,-1,1)
  eh = (data.vcomp.prime('time', 'lon') * data.temp.prime('time', 'lon')).avg('time')
  dse = data.drystaticenergy.avg('time')
  lhf = data.latentheat.avg('time')
  del data
  
  sig0p5 = int(np.argmin(np.abs(pfull - pfull[-1] / 2.)))
  
  for level, ind in zip(['surface', 'sig0p5'], [-1, sig0p5]):
  
    params = [[lon, lat, v[ind], None, False, 'Meridional wind'],
              [lon, lat, u[ind], None, False, 'Zonal wind'],
              [lon, lat, t[ind], None, False, 'Temperature'],
              [lon, lat, mt[ind], None, False, 'Total angular momentum'],
              [lon, lat, mr[ind], None, False, 'Relative angular momentum'],
              [lon, lat, em[ind], None, False, 'Eddy angular momentum flux'],
              [lon, lat, eh[ind], None, False, 'Eddy heat flux'],
              [lon, lat, dse[ind], None, False, 'Dry static energy flux'],
              [lon, lat, lhf[ind], None, False, 'Latent heat flux']]
    
    # Set up plot
    fig, ax = pl.subplots(3,3, figsize = (14, 10))
    fig.subplots_adjust(wspace = 0.33, hspace = 0.33)
    ax = ax.flatten()
    div = [make_axes_locatable(axis) for axis in ax]
    cax = [d.append_axes("right", size="7%", pad=0.05) for d in div]
    
    for n, _ in enumerate(params):
      _x, _y, _z, _c, _i, _t = params[n]
      _, cbar = gcm.ColorMap(_x, _y, _z, c = _c, invert_y = _i, ax = ax[n], cax = cax[n])
      cbar.ax.tick_params(labelsize=8)
      ax[n].set_title(_t, fontsize = 12)
      [tick.label.set_fontsize(8) for tick in ax[n].yaxis.get_major_ticks() + ax[n].xaxis.get_major_ticks()]
      ax[n].yaxis.label.set_fontsize(8)
      ax[n].xaxis.label.set_fontsize(8)
    pl.suptitle(level, fontsize = 18)
    
    fig.savefig('images/levels_%s_%s.png' % (level, run))
    pl.close()