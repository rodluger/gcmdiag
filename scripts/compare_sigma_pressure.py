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

  # Rectify and compare
  x = []
  y = []
  z = []
  c = []
  for rectify in [False, True]:
  
    data = gcm.NetCDF(file, rectify = rectify)
    x.append(data.lat)
    y.append(data.pfull)
    z.append(data.total_angular_momentum)
    c.append(data.streamfunction)
    del data
  
  titles = ['Sigma coordinates', 'Pressure coordinates']
  fig, ax = gcm.Compare(x, y, z, titles = titles, c = c, invert_y = True)
  pl.suptitle('Hadley Cell', fontsize = 18)
  fig.savefig('images/compare_sigma_pressure_%s.png' % run)