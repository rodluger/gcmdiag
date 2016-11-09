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
  fig = pl.figure(figsize = (14, 10))
  ax = [pl.subplot2grid((1, 22), (0, 0), colspan = 9),
        pl.subplot2grid((1, 22), (0, 12), colspan = 9)]
  cax = [pl.subplot2grid((1, 22), (0, 9)),
         pl.subplot2grid((1, 22), (0, 21))]
  
  fig.subplots_adjust(top = 0.875)

  gcm.ColorMap(x, y, toa, ax = ax[0], cax = cax[0])
  ax[0].set_title('TOA imbalance')
  
  gcm.ColorMap(x, y, mse, ax = ax[1], cax = cax[1])
  ax[1].set_title('ConvIntMeridFlux')
  
  pl.suptitle('Energy conservation', fontsize = 18)
  fig.savefig('images/compare_toa_mse_%s.png' % run)