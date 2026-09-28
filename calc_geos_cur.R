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
#```
#
#Check: if values increase to the east, the slope faces west (aspect = 270°). Then sin(aspect) = −1, so dz/dx comes out positive, as it should.
#
#- The units are z-units per metre. For example, a field in metres gives m/m.
#- Flat cells have an undefined aspect (`NA`), so set their gradient to 0 if you need to.
#
#If you're working with SSH (sea surface height) to get geostrophic velocities, Option 2 is the easiest correct route on a lon/lat grid. After that, u = −(g/f)·∂η/∂y and v = (g/f)·∂η/∂x.
#






### Defining f from latitude
#
#The Coriolis parameter is
#
#**f = 2 Ω sin(φ)**
#
#where Ω = 7.2921×10⁻⁵ rad/s is Earth's rotation rate and φ is latitude. Your raster is in lon/lat, so you can build an f raster directly from the cell-centre y coordinates:
#
#```r
#lat   <- init(r, "y")                 # raster of cell-centre latitudes (degrees)
#Omega <- 7.2921e-5                    # rad/s
#f     <- 2 * Omega * sin(lat * pi/180)  # s^-1
#
#g <- 9.81
#u <- -(g / f) * dzdy   # eastward
#v <-  (g / f) * dzdx   # northward
#```
#
#Around Lake Erie (~42°N), f ≈ 9.8×10⁻⁵ s⁻¹. The lake only spans about 41.4–42.9°N, so f changes by roughly 3% across it. A single constant, `f <- 2*Omega*sin(42*pi/180)`, would work nearly as well. The per-cell version costs nothing extra, though.
#
### Why multiply by g
#
#Geostrophic balance means the Coriolis force balances the pressure-gradient force:
#
#f **k** × **u** = −(1/ρ) ∇p
#
#The surface slope ∂η/∂x is just a ratio (m/m), not a force. Hydrostatic balance turns it into one. Under a surface that's higher by η, the pressure at any level is higher by ρgη, so
#
#−(1/ρ) ∇p = −g ∇η
#
#g converts the tilt of the surface into a pressure-gradient acceleration. You can see it in the units:
#
#g [m s⁻²] × slope [m/m] ÷ f [s⁻¹] = **m/s**
#
#Without g, you'd get a slope divided by a frequency, which isn't a velocity.
#
### Caveats for this setup
#- **Geostrophy may be weak in Erie.** The lake is shallow and mostly driven by wind, so friction and Ekman effects matter. Geostrophic currents from η will only be part of the flow. Comparing against FVCOM's own `u`/`v` output would show how much.
#- **IDW interpolation and gradients.** `interpIDW` can leave small bumps around each node, and taking the gradient amplifies them into noisy velocities. You may want to smooth first (e.g. `focal(r, w=3, fun=mean)`), or compute gradients on the native FVCOM mesh.
#- **Flat cells**, where aspect is `NA`, should have their gradient set to 0, as the note in the file already says.







