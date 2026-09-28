#!/usr/bin/Rscript
library(sf) 
library(ncdf4)

# read in netCDF
ncid <- nc_open('data/eri_2026.nc')
eta = ncvar_get(ncid, 'zeta')[,1] # should be a vector (assign in a time loop if need be)


### =========== OPTION 1 (simple) ===============
### plot variable on Voronoi polygon
nds = st_read('data/spatial/eri_nodes.kml', quiet=T)
nds['eta'] = eta
plot(nds['eta'], axes=T, lwd=.1)


#### =========== OPTION 2  (medium) ===============
#### plot variable on elements (after straight avging over nodes)
#els = st_read('spatial_files/deb_elems.shp', quiet=T)
#
## verts is an num_elems x 4 matrix with the indices for each vertex (NA in 4th column for tris)
#verts <- matrix(as.numeric(unlist(strsplit(els$nodes, split=','))), ncol=4, byrow=T)
#verts[verts == -99] <- NA
#els['varname'] <- apply(matrix(plot_var[c(verts)], ncol=4, byrow=F), 1, mean, na.rm=T)
#plot(els['varname'], axes=T, lwd=.1)
#
