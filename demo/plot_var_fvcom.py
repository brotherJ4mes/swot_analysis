#!/usr/bin/python3.9
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
plt.ion()

# read in netCDF
nc=netCDF4.Dataset('model_output/eri_24.nc')
plot_var = nc['zeta'][:] # select var of interest (vector; num_nodes x 1)
#plot_var = nc['elev'][:] # substitute for variable of interest


### =========== OPTION 1 (simple) ===============
### plot variable on Voronoi node
vor = gpd.read_file('spatial_files/eri_nodes.shp')
vor['zeta'] = plot_var # replace varname with "depth", "zeta", "temp", etc.
vor.plot('zeta', legend=True)




#### =========== OPTION 2  (medium) ===============
#### plot variable on elements (after straight avging over nodes)
#els = gpd.read_file('spatial_files/deb_elems.shp')
#isQuad=els['shape']=='quad'
#isTri=els['shape']=='tri'
## find node indices for each element
#nodeIdx = np.array([list(map(int, els['nodes'].values[i].split(','))) for i in range(els.shape[0])])
#nodeIdx[nodeIdx!=-99] = nodeIdx[nodeIdx!=-99] - 1 # switch to 0 based indexing
## associate var with els; use nodeIdx to average over tris (3 nodes) and quads (4 nodes)
#els['varname'] = -99.0
#els['varname'][isTri] = np.mean(plot_var[nodeIdx[isTri,:3]], axis=1)
#els['varname'][isQuad] = np.mean(plot_var[nodeIdx[isQuad,:4]], axis=1)
#els.plot('varname', legend=True)
