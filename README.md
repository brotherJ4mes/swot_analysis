
#### contents:

`earth_data_download.py` 
gets data from full pass)

`find_sub_indices.R`   
- determines the relevant indices in num_lines to clip netCDF to just the Great Lakes domain
- can then be executed like` /usr/bin/ncks -d num_lines,2156,2630 input.nc output.nc`

`geos_swot.R`
- reads in SWOT_L2_LR_SSH data and derives  d_eta/dx and d_eta/dy and then geostrophic flow
- [ ] low pass filter SSH before derivation
- [ ] removal of seiches and other modes?

`geos_fvcom.R`
- work in progress... will do the same as `geo_swot.R`



#### timing/data:

CalVal Orbit Start: 2023-01-15 
Science Orbit Start: 2023-07-11
(see table 4.1 in manual for more details)

#### terminology and file_naming

```
SSHA
=
(SSH
−mean_sea_surface_height
−solid_earth_tide_height
−ocean_tide_height
−load_tide_height
−internal_tide_height
−pole_tide_height
−dynamic_atmosphere_correction).
```




####  geminis suggestions:
If you do compute geostrophic currents $u_g = -\frac{g}{f}\frac{\partial \eta}{\partial y},\; v_g = \frac{g}{f}\frac{\partial \eta}{\partial x}$ from L2 `ssha_karin` (or `ssh_karin`):

1. **You must denoise/smooth before differentiating:** At 2 km posting, raw L2 KaRIn noise ($\sim 1\text{--}2\text{ cm}$ RMS) produces unrealistic velocity noise on the order of $0.5\text{--}1+\text{ m/s}$ under direct finite differences. At minimum, apply a 2D spatial low-pass filter (e.g., Gaussian or Boxcar/Lanczos with a length scale of $\sim 10\text{--}20\text{ km}$) before taking gradients.
2. **Don't differentiate across the nadir gap or bad flags:** The ~20 km gap in the center of the KaRIn swath (and masked pixels from `ssha_karin_qual != 0`) will create severe artifacts if an interpolation or smoothing kernel bridges across missing swaths without masking.
3. **Add back Mean Dynamic Topography (MDT) for *total* currents:** `ssha_karin` gives geostrophic velocity *anomalies*. For total geostrophic currents, use `ssh_karin` minus the geoid (or `ssha_karin` + `mean_dynamic_topography` provided in the Expert/L3 files).
> 


