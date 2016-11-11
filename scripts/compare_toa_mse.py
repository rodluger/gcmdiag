import gcmdiag as gcm
from gcmdiag import GCMDIAG_OUT
from gcmdiag.constants import *
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
  toa = data.toa_imbalance
  efg = data.energy_flux_gradient
  efg_eddy = data.eddy_energy_flux_gradient
  
  # Plot
  fig = pl.figure(figsize = (14, 10))
  ax = [pl.subplot2grid((1, 22), (0, 0), colspan = 9),
        pl.subplot2grid((1, 22), (0, 12), colspan = 9)]
  cax = pl.subplot2grid((1, 22), (0, 9))
  
  fig.subplots_adjust(top = 0.875)
  
  # Colormap
  gcm.ColorMap(x, y, toa, ax = ax[0], cax = cax)
  ax[0].set_title('TOA imbalance')
  
  # 1D comparison
  ax[1].plot(efg, y, 'b-', lw = 2, label = r'$\nabla\Phi$')
  ax[1].plot(efg_eddy, y, 'b--', lw = 2, alpha = 0.5, label = r'$\nabla\Phi_\mathrm{eddy}$')
  ax[1].plot(toa.avg('lon'), y, 'r-', lw = 2, label = 'TOA')
  ax[1].set_ylabel(ax[0].get_ylabel())
  ax[1].set_ylim(ax[0].get_ylim())
  ax[1].set_xlabel('Energy flux (W/m^2)')
  ax[1].yaxis.tick_right()
  ax[1].yaxis.set_label_position("right")
  ax[1].legend(loc = 'upper right')
  ax[1].set_title(r'Comparison to $\nabla\Phi$')
  
  pl.suptitle('Energy conservation', fontsize = 18)
  fig.savefig('images/compare_toa_mse_%s.png' % run)
  