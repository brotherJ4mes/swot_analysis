#!/bin/Rscript
library(terra)
library(ncdf4)

#varname='ssha_unfiltered'
#fin = 'SWOT_L2_LR_SSH_Basic_032_132_20250502T171047_20250502T180215_PIC2_01.nc'
var_name='ssh_karin_2' # name for L2
#fin = 'SWOT_L2_LR_SSH_Basic_032_160_20250503T171118_20250503T180246_PIC2_01.nc'
#fin = 'SWOT_L2_LR_SSH_Expert_032_119_20250502T060159_20250502T065327_PIC2_01.nc'
#fin = '../data_swot/calval/out.nc'
fin = 'ont_L2_LR_SSH_basic.nc'

print('reading SWOT...')
ncid = nc_open(fin)
lons = as.vector(ncvar_get(ncid,'longitude'))
lats = as.vector(ncvar_get(ncid,'latitude'))
ssh  = as.vector(ncvar_get(ncid, var_name))
#nc_close(ncid)

lons[lons>180] = lons[lons>180] - 360
df = data.frame(lon=lons, lat=lats, eta=ssh)
#df = df[!is.na(df$eta),]

print('reading lakes...')
lks = st_read('homes_shoreline.kml')
print('fixing lakes (write this out)...')
lks = st_make_valid(lks)

print('vectorizing...')
pts = vect(df, crs='+proj=longlat +datum=WGS84')

print('cropping...')
pts = mask(pts, vect(lks))
#pts = na.omit(pts, field='eta')

#plot(pts, 'eta', type='continuous')
#print('skipping intersecting!')
##print('intersecting...')
#sel_water = is.related(pts, vect(lks), 'intersects')
#pts = pts[sel_water,]

egm = rast('egm2008_GL.tif')
egm_pts = extract(egm, pts)
pts$eta = pts$eta - egm_pts$geoid_undulation


print('rasterizing...')
r_eta = rast(ext(pts), resolution = 2.5e-2, crs = crs(pts)) # res in degrees (guessed to match LR resolution)
r_eta = interpIDW(r_eta, pts, field = "eta", radius = 50, power = 2)
print('done')


# Unscaled Horn kernels (numerator only: units of eta)
wgts_x <- matrix(c(-1, 0, 1,
                   -2, 0, 2,
                   -1, 0, 1)/8, nrow = 3, byrow = TRUE)

wgts_y <- matrix(c( 1,  2,  1,
                    0,  0,  0,
                   -1, -2, -1)/8, nrow = 3, byrow = TRUE)

# Cell-center latitude in radians
lat_rad = init(r_eta, "y") * (pi / 180)

# Physical cell widths in meters (WGS84 equatorial radius R = 6,378,137 m)
R_e = 6378137
dx_m = (res(r_eta)[1] * pi / 180) * R_e * cos(lat_rad)
dy_m = (res(r_eta)[2] * pi / 180) * R_e

# True physical partial derivatives in [units of eta / meter]
d_eta_dx <- focal(r_eta, w = wgts_x, fun = "sum") / dx_m
d_eta_dy <- focal(r_eta, w = wgts_y, fun = "sum") / dy_m

grad_xy <- c(d_eta_dx, d_eta_dy)
names(grad_xy) <- c("d_eta_dx", "d_eta_dy")

omega = 0.00007292 # earths rotation (for f)
lats = crds(grad_xy)[,2]*pi/180
f = 2*omega*sin(lats)

#dn/dx = -fv; dn/dy = fu

# plot eta
plot(r_eta)
plot(lks, add=T)
plot(pts, add=T, pch=20, cex=.1)

# TODO: fix dimension mismatch 
u = grad_xy$d_eta_dy/-mean(f)
v = grad_xy$d_eta_dx/mean(f)
speed = sqrt(u^2+v^2)

print('masking the raster...')
r_eta = mask(r_eta, pts)

vectorplot(raster::stack(c(u,v)), isField = "dXY", par.settings = viridisTheme())


# Smooth interpolation via Thin Plate Spline (requires fields)
#tps <- Tps(crds(pts), pts$etapat
#r <- interpolate(r_tmpl, tps)

# Option A: Component-wise gradient via central differences
#dx <- res(r)[1]
#dy <- res(r)[2]
#
#kx <- matrix(c(0,  0, 0,
#              -1,  0, 1,
#               0,  0, 0) / (2 * dx), nrow = 3, byrow = TRUE)
#
#ky <- matrix(c(0,  1, 0,
#               0,  0, 0,
#               0, -1, 0) / (2 * dy), nrow = 3, byrow = TRUE)
#
#dz_dx <- focal(r, w = kx)
#dz_dy <- focal(r, w = ky)
#grad_mag <- sqrt(dz_dx^2 + dz_dy^2)

# Option B: Magnitude and direction directly via terra::terrain()
#slp <- terrain(r_surf, v = c("slope", "aspect"), unit = "radians")
#grad_mag_terrain <- tan(slp$slope) # |grad z| = tan(slope)


#plot(lks)
#plot(pts, 'eta', add=T)


#print('plotting...')
#plot(pts_lk, 'adj', type='continuous')
#plot(lks, add=T)




