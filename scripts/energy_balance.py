import gcmdiag as gcm
from gcmdiag import GCMDIAG_OUT
from gcmdiag.constants import *
import os
import matplotlib.pyplot as pl
import numpy as np

# The runs we'll plot
runs = ['1bar', 'p5bar', 'p25bar']

# Plot
fig, ax = pl.subplots(3, 3, figsize = (12, 12))
fig.subplots_adjust(wspace = 0.05, hspace = 0.075)

# Do all three original runs
for i, run in enumerate(runs):
  
  '''
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
  '''
  
  # Energy balance
  #ax[i,0].plot(y, toa, 'r-', lw = 2, label = 'TOA')
  #ax[i,0].plot(y, efg, 'b-', lw = 2, label = r'$\nabla\Phi$')
  ax[i,0].set_xlim(-90, 90)
  if i < 2: ax[i,0].set_xticklabels([])
  ax[i,0].set_ylabel(r'Energy flux (W/m$^2$)')
  ax[i,0].legend(loc = 'upper left')
  
  # Mean and eddy
  #ax[i,1].plot(y, efg_eddy, 'r-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{eddy}$')
  #ax[i,1].plot(y, efg_mean, 'b-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{mean}$')
  ax[i,1].set_xlim(-90, 90)
  ax[i,1].set_yticklabels([])
  if i < 2: ax[i,1].set_xticklabels([])
  ax[i,1].legend(loc = 'upper left')
  
  # Dry and latent
  #ax[i,2].plot(y, efg_dry, 'r-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{dry}$')
  #ax[i,2].plot(y, efg_latent, 'b-', lw = 2, alpha = 0.75, label = r'$\nabla\Phi_\mathrm{latent}$')
  ax[i,2].set_xlim(-90, 90)
  ax[i,2].set_yticklabels([])
  if i < 2: ax[i,2].set_xticklabels([])
  ax[i,2].legend(loc = 'upper left')
  
# Appearance
ax[0,0].set_title('Energy balance')
ax[0,1].set_title('Mean and eddy')
ax[0,2].set_title('Latent and dry')
[ax[n,0].text(-0.4, 0.5, runs[n], rotation = 90, fontsize = 16, transform = ax[n,0].transAxes, ha = 'center', va = 'center') for n in range(3)]
[ax[-1,n].set_xlabel('Latitude (deg)') for n in range(3)]
[axis.axhline(0, color = 'k', alpha = 0.5, zorder = -99, ls = '--') for axis in ax.flatten()]
[axis.set_ylim(-200, 250) for axis in ax.flatten()]

# Save
fig.savefig('images/energy_balance.png')