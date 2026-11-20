#!/bin/Rscript
library(terra)
library(ncdf4)

#varname='ssha_unfiltered'
#fin = 'SWOT_L2_LR_SSH_Basic_032_132_20250502T171047_20250502T180215_PIC2_01.nc'
var_name='ssh_karin_2' # name for L2
#fin = 'SWOT_L2_LR_SSH_Basic_032_160_20250503T171118_20250503T180246_PIC2_01.nc'
#fin = 'SWOT_L2_LR_SSH_Expert_032_119_20250502T060159_20250502T065327_PIC2_01.nc'
fin = '../data_swot/calval/out.nc'

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
lks = st_geometry(st_read('homes_shoreline.kml'))

print('vectorizing...')
pts = vect(df, crs='+proj=longlat +datum=WGS84')
#print('cropping...')
#pts = crop(pts, vect(lks), ext=T)
pts = na.omit(pts, field='eta')

#plot(pts, 'eta', type='continuous')
#print('intersecting...')
#sel_water = is.related(pts, vect(lks), 'intersects')
#pts_lk = pts[sel_water,]

print('skipping intersecting!')
#sel_water = is.related(pts, vect(lks), 'intersects')
pts_lk = pts

egm = rast('egm2008_GL.tif')
egm_pts = extract(egm, pts_lk)
pts_lk$adj = pts_lk$eta - egm_pts$geoid_undulation




r_tmpl <- rast(ext(pts), resolution = 100, crs = crs(pts))

# Smooth interpolation via Thin Plate Spline (requires fields)
#tps <- Tps(crds(pts pts$eta)
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
#slp <- terrain(r, v = c("slope", "aspect"), unit = "radians")
#grad_mag_terrain <- tan(slp$slope) # |grad z| = tan(slope)


#plot(lks)
#plot(pts, 'eta', add=T)



#print('plotting...')
#plot(pts_lk, 'adj', type='continuous')
#plot(lks, add=T)




