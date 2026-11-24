#!/bin/Rscript
library(terra)
library(sf)
library(ncdf4)

#varname='ssha_unfiltered'
#fin = 'SWOT_L2_LR_SSH_Basic_032_132_20250502T171047_20250502T180215_PIC2_01.nc'
var_name='ssh_karin_2' # name for L2
#fin = 'SWOT_L2_LR_SSH_Basic_032_160_20250503T171118_20250503T180246_PIC2_01.nc'
#fin = 'SWOT_L2_LR_SSH_Expert_032_119_20250502T060159_20250502T065327_PIC2_01.nc'
#fin = '../data_swot/calval/out.nc'
#fin = 'ont_L2_LR_SSH_basic.nc'
fin = commandArgs(trail=T)
out_dir = 'png'

print('reading SWOT...')
ncid = nc_open(fin)
lons = as.vector(ncvar_get(ncid,'longitude'))
lats = as.vector(ncvar_get(ncid,'latitude'))
ssh  = as.vector(ncvar_get(ncid, var_name))
valid_dt = mean(as.POSIXct(ncvar_get(ncid, 'time'), origin='2000-01-01', tz='UTC'), na.rm=T)
date_str = format(valid_dt, '%Y-%b-%d %HZ')
print(date_str)

#nc_close(ncid)

lons[lons>180] = lons[lons>180] - 360
df = data.frame(lon=lons, lat=lats, eta=ssh)
#df = df[!is.na(df$eta),]

print('reading lakes...')
lks = st_read('homes_shoreline.kml', quiet=T)
#print('fixing lakes (write this out)...')
#lks = st_make_valid(lks)

print('vectorizing...')
pts = vect(df, crs='+proj=longlat +datum=WGS84')

print('removing points over land ...')
pts = mask(pts, vect(lks))
pts = na.omit(pts, field='eta')

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

grad_xy = mask(grad_xy, pts)

omega = 0.00007292 # earths rotation (for f)
lats = crds(grad_xy, na.rm=T)[,2]*pi/180
f = 2*omega*sin(lats)

#dn/dx = -fv; dn/dy = fu
# TODO: fix dimension mismatch 

u = grad_xy$d_eta_dy/-mean(f)
v = grad_xy$d_eta_dx/mean(f)
speed = sqrt(u^2+v^2)



# vectorize for plotting
r_eta = mask(r_eta, pts)
uvec =  u[!is.na(u)]
vvec =  v[!is.na(v)]
svec =  speed[!is.na(speed)]

x0 = crds(u, na.rm=T)[,1]
y0 = crds(u, na.rm=T)[,2]

scale=1000 # set big enough to eliminate magnitude

fout = paste(out_dir, format(valid_dt, '%Y%m%dT%H.png'), sep='/')
png(fout, width=1400, height=700, pointsize=24)
layout(matrix(c(1,1,2,3),2,2), widths=c(.7,.4))
par(mar=c(0,0,0,0))
plot(r_eta, main=date_str)
plot(lks, col=NA, add=T, border='lightblue')
arrows(x0=x0, y0=y0, x1=x0+uvec/scale, y1=y0+vvec/scale, length=.1, angle=20)


par(mar=c(2,0,2,0))
xlim = c(-2,2)*100
ubar = round(global(abs(u), 'mean', na.rm=T)*100, digits=1)
vbar = round(global(abs(v), 'mean', na.rm=T)*100, digits=1)
hist(u*100, main='u (cm/s)', yaxt='n', xlab=NA, xlim=xlim, ylab=NA)
mtext(sprintf('mean: %4.1f', ubar), adj=.95, line=-2)
hist(v*100, main='v (cm/s)', yaxt='n', xlab=NA, xlim=xlim, ylab=NA)
mtext(sprintf('mean: %4.1f', vbar), adj=.95, line=-2)
dev.off()

#arrows(x0=x0, y0=y0, x1=x0+uvec/svec/scale, y1=y0+vvec/svec/scale, length=.1, angle=20)

#print('masking the raster...')
#r_eta = mask(r_eta, pts)



#vectorplot(raster::stack(c(u,v)), isField = "dXY", par.settings = viridisTheme())

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




