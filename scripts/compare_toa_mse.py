import gcmdiag as gcm
from gcmdiag import GCMDIAG_OUT
import os
import matplotlib.pyplot as pl
import numpy as np

# Do all three original runs
for run in ['1bar', 'p5bar', 'p25bar']:

  # Load
  print("Running %s..." % run)
  file = os.path.join(GCMDIAG_OUT, 'N2-%s/day1000h00/day1000h00.atmos_daily.nc' % run)
  data = gcm.NetCDF(file, rectify = True)
  
  # Get arrays
  x = data.lon
  y = data.lat
  toa = data.toaimbalance
  mse = data.convintmeridflux
  
  # Plot
  vmin = np.min([toa.min(), mse.min()])
  vmax = np.max([toa.max(), mse.max()])
  fig = pl.figure(figsize = (14, 10))
  fig.subplots_adjust(top = 0.875)
  ax = [pl.subplot2grid((1, 11 * 2 + 1), (0, 11 * n), colspan = 10) for n in range(2)]
  cax = pl.subplot2grid((1, 11 * 2 + 1), (0, 11 * 2))

  ColorMap(x, y, toa, ax = ax[0], cax = cax, vmin = vmin, vmax = vmax)
  ax[0].set_title('TOA imbalance')
  
  ColorMap(x, y, mse, ax = ax[1], cax = cax, vmin = vmin, vmax = vmax)
  ax[1].set_title('ConvIntMeridFlux')
  ax[1].set_ylabel('')
  ax[1].set_yticklabels([])
  
  pl.suptitle('Energy conservation', fontsize = 18)
  fig.savefig('images/compare_toa_mse_%s.png' % run)