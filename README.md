# gcmdiag

## Quick guide

To merge multiple NetCDF files that were generated during a parallel run:

```python
import gcmdiag as gcm
gcm.NetCDFCombine(path = 'path/to/files/', outfile = 'path/to/output.nc')
```

To convert a NetCDF file from sigma coordinates to pressure coordinates and interpolate all arrays:

```python
import gcmdiag as gcm
gcm.Rectify(path = 'path/to/file.nc', interpolation = 'linear')
```
