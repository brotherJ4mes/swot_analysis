#!/usr/bin/python
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import pandas as pd
plt.ion()


# ====================================================================
# ================    USER CONTROLS    ===============================
# ====================================================================

date_str = '2026-07-01 12:00' # USE 'YYYY-MM-DD HH:MM' format
lake = 'eri' # moniker for lake: sup, mih, eri, ont

# ====================================================================
# ====================================================================
# ====================================================================

igld_shift = dict(sup=183, mih=176, eri=174, ont=75)

# read in netCDF
nc=netCDF4.Dataset('data/'+lake+'_'+date_str[:4]+'.nc')
eta = nc['zeta'][:] + igld_shift[lake] #  [ num_nodes x time_steps]

dts = pd.to_datetime(nc['Itime'][:]*1.0 + nc['Itime2'][:]/1e3/3600/24, unit='D') # read in dates
dts = dts.round('h') # round to hour to fix fractional second
t_sel = dts.get_loc(date_str) # find desired datetime index

date = dts[t_sel]
eta = eta[t_sel,:]

### plot variable on node
nds = gpd.read_file('data/spatial/'+lake+'_nodes.kml')
nds['eta'] = eta
nds.plot('eta', legend=True)

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
