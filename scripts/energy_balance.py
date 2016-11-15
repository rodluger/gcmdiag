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
  fig, ax = pl.subplots(1, 3, figsize = (16, 6))
  fig.subplots_adjust(wspace = 0.05)
  ax = ax.flatten()
  
  # Energy balance
  ax[0].plot(y, toa, 'r-', lw = 2, label = 'TOA')
  ax[0].plot(y, efg, 'b-', lw = 2, label = r'$\nabla\Phi$')
  ax[0].set_xlabel('Latitude (deg)')
  ax[0].set_xlim(-90, 90)
  ax[0].set_ylabel('Energy flux (W/m^2)')
  ax[0].legend(loc = 'upper left')
  ax[0].set_title('Mean and eddy')
  
  # Mean and eddy
  ax[1].plot(y, efg_eddy, 'r-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{eddy}$')
  ax[1].plot(y, efg_mean, 'b-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{mean}$')
  ax[1].set_xlabel('Latitude (deg)')
  ax[1].set_xlim(-90, 90)
  ax[1].set_ylabel('')
  ax[1].set_yticklabels([])
  ax[1].legend(loc = 'upper left')
  ax[1].set_title('Mean and eddy')
  
  # Dry and latent
  ax[2].plot(y, efg_dry, 'r-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{dry}$')
  ax[2].plot(y, efg_latent, 'b-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{latent}$')
  ax[2].set_xlabel('Latitude (deg)')
  ax[2].set_xlim(-90, 90)
  ax[2].set_ylabel('')
  ax[2].set_yticklabels([])
  ax[2].legend(loc = 'upper left')
  ax[2].set_title('Latent and dry')
  
  # Scale all plots the same
  ymin = np.min([axis.get_ylim()[0] for axis in ax])
  ymax = np.max([axis.get_ylim()[1] for axis in ax])
  [axis.set_ylim(ymin, ymax) for axis in ax]
  
  pl.suptitle('Energy conservation', fontsize = 18)
  fig.savefig('images/energy_balance_%s.png' % run)