import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# 1. Open the netCDF file (xarray uses netCDF4 under the hood)
ds = xr.open_dataset("SWOT_L2_LR_SSH_Basic_032_132_20250502T171047_20250502T180215_PIC2_01.nc")

# 2. Convert to a pandas DataFrame
#    SWOT files typically have dims: (swath, time) and coords: latitude, longitude
df = ds.to_dataframe()

# 3. Reset the MultiIndex so lat/lon become regular columns
df = df.reset_index()

# 4. Build Point geometries from the lat/lon columns
geom = [Point(lon, lat) for lon, lat in zip(df["longitude"], df["latitude"])]

# 5. Create the GeoDataFrame (set CRS to WGS84)
gdf = gpd.GeoDataFrame(df, geometry=geom, crs="EPSG:4326")
