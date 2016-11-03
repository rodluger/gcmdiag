import gcmdiag as gcm
from gcmdiag import GCMDIAG_SRC
import os
import matplotlib.pyplot as pl
import numpy as np

# Load
file = os.path.join(GCMDIAG_SRC, 'N2-1bar/day1000h00/day1000h00.atmos_daily.nc')

# Rectify and compare
x = []
y = []
z = []
c = []
for rectify in [False, True]:
  
  data = gcm.NetCDF(file, rectify = rectify)
  x.append(data.lat)
  y.append(data.pfull)
  z.append(data.totalangmom)
  c.append(data.streamfunc)
  del data
  
titles = ['Sigma coordinates', 'Pressure coordinates']
fig, ax = gcm.Compare(x, y, z, titles = titles, c = c, invert_y = True)
pl.show()