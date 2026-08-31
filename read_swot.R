#!/bin/Rscript
library(terra)
library(ncdf4)


#varname='ssha_unfiltered'
fin = 'SWOT_L2_LR_SSH_Basic_032_132_20250502T171047_20250502T180215_PIC2_01.nc'
varname='ssh_karin' # name for L2
#fin = 'SWOT_L3_LR_SSH_Basic_032_160_20250503T171119_20250503T180245_v3.0.nc'

ncid = nc_open(fin)
lons = as.vector(ncvar_get(ncid,'longitude'))
lats = as.vector(ncvar_get(ncid,'latitude'))
ssh  = as.vector(ncvar_get(ncid, varname))
nc_close(ncid)


lons[lons>180] = lons[lons>180] - 360

df = data.frame(lon=lons, lat=lats, eta=ssh)
#df = df[!is.na(df$eta),]

lks = st_geometry(st_read('~/work/geospatial_data/homes_shoreline.shp'))

pts = vect(df, crs='+proj=longlat +datum=WGS84')
pts = crop(pts, vect(lks), ext=T)
pts = na.omit(pts, field='eta')

#plot(pts, 'eta', type='continuous')
sel_water = is.related(pts, vect(lks), 'intersects')
pts_lk = pts[sel_water,]


plot(pts_lk, 'eta', type='continuous')
plot(lk, add=T)








