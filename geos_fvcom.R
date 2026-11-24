#!/usr/bin/Rscript
library(sf) 
library(ncdf4)
library(terra)

# read in netCDF
my_res = 5e-3 # resolution [degrees]
ncid <- nc_open('data/eri_2026.nc')
eta = ncvar_get(ncid, 'zeta')[,1] # should be a vector (assign in a time loop if need be)

nds = st_read('data/spatial/eri_nodes.kml', quiet=T)
nds['eta'] = eta
v = vect(nds['eta'])

r_template = rast(v, resolution=my_res)

# (rasterize) put points into cells (missing = NA)
#r = rasterize(v, r_template, field='eta', fun=mean)

# interpolate onto rast
r = interpIDW(r_template, v, field='eta', radius=my_res*4, power=2, maxPoints=10)


# now find gradient
M = nrow(r)
N = ncol(r)


## Option 2: Use `terrain()` (handles lon/lat properly)
#
#With a geographic (lon/lat) CRS, the resolution is in degrees and the east-west cell width changes with latitude, so Option 1 needs extra work. `terrain()` already computes distances in metres for lon/lat rasters, so you can split its output into components:
#
#```r
slp <- terrain(r, v = "slope",  unit = "radians", neighbors = 8)
asp <- terrain(r, v = "aspect", unit = "radians", neighbors = 8)
#
## aspect = direction the slope faces (downhill), clockwise from north
dzdx <- -tan(slp) * sin(asp)   # east component
dzdy <- -tan(slp) * cos(asp)   # north component
