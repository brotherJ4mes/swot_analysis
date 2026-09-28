#!/bin/Rscript
library(terra)
library(ncdf4)

#varname='ssha_unfiltered'
#fin = 'SWOT_L2_LR_SSH_Basic_032_132_20250502T171047_20250502T180215_PIC2_01.nc'
varname='ssh_karin' # name for L2
#fin = 'SWOT_L2_LR_SSH_Basic_032_160_20250503T171118_20250503T180246_PIC2_01.nc'
fin = 'SWOT_L2_LR_SSH_Expert_032_119_20250502T060159_20250502T065327_PIC2_01.nc'

print('reading SWOT...')
ncid = nc_open(fin)
lons = as.vector(ncvar_get(ncid,'longitude'))
lats = as.vector(ncvar_get(ncid,'latitude'))
ssh  = as.vector(ncvar_get(ncid, varname))
nc_close(ncid)

lons[lons>180] = lons[lons>180] - 360

df = data.frame(lon=lons, lat=lats, eta=ssh)
#df = df[!is.na(df$eta),]

print('reading lakes...')
lks = st_geometry(st_read('~/work/homes_shoreline.shp'))

print('vectorizing...')
pts = vect(df, crs='+proj=longlat +datum=WGS84')
print('cropping...')
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


plot(lks)
plot(pts, 'eta', add=T)



print('plotting...')
#plot(pts_lk, 'adj', type='continuous')
#plot(lks, add=T)




