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
  toa = data.toa_imbalance().avg('time', 'lon')
  efg = data.energy_flux_gradient()
  efg_eddy = data.energy_flux_gradient(eddy = True)
  efg_mean = efg - efg_eddy
  efg_dry = data.energy_flux_gradient(dry = True, latent = False)
  efg_latent = data.energy_flux_gradient(dry = False, latent = True)
  
  # Plot
  fig = pl.figure(figsize = (14, 10))
  ax = [pl.subplot2grid((1, 21), (0, 0), colspan = 10),
        pl.subplot2grid((1, 21), (0, 11), colspan = 10)]
  fig.subplots_adjust(top = 0.875)
  
  # Mean and eddy parts
  ax[0].plot(toa, y, 'r-', lw = 2, label = 'TOA')
  ax[0].plot(efg, y, 'b-', lw = 2, label = r'$\nabla\Phi$')
  ax[0].plot(efg_eddy, y, 'b--', lw = 2, alpha = 0.5, label = r'$\nabla\Phi_\mathrm{eddy}$')
  ax[0].plot(efg_mean, y, 'b-.', lw = 2, alpha = 0.5, label = r'$\nabla\Phi_\mathrm{mean}$')
  ax[0].set_ylabel(ax[0].get_ylabel())
  ax[0].set_ylim(ax[0].get_ylim())
  ax[0].set_xlabel('Energy flux (W/m^2)')
  ax[0].yaxis.tick_right()
  ax[0].yaxis.set_label_position("right")
  ax[0].legend(loc = 'upper right')
  ax[0].set_title('Mean and eddy')
  
  # Dry and latent
  ax[0].plot(toa, y, 'r-', lw = 2, label = 'TOA')
  ax[1].plot(efg, y, 'b-', lw = 2, label = r'$\nabla\Phi$')
  ax[1].plot(efg_dry, y, 'b--', lw = 2, alpha = 0.5, label = r'$\nabla\Phi_\mathrm{dry}$')
  ax[1].plot(efg_latent, y, 'b--', lw = 2, alpha = 0.5, label = r'$\nabla\Phi_\mathrm{latent}$')
  ax[1].set_ylabel(ax[0].get_ylabel())
  ax[1].set_ylim(ax[0].get_ylim())
  ax[1].set_xlabel('Energy flux (W/m^2)')
  ax[1].yaxis.tick_right()
  ax[1].yaxis.set_label_position("right")
  ax[1].legend(loc = 'upper right')
  ax[1].set_title('Latent and dry')
  
  pl.suptitle('Energy conservation', fontsize = 18)
  fig.savefig('images/energy_balance_%s.png' % run)