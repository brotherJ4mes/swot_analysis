library(stars)
library(dplyr)

fin = '../swot_data/SWOT_L2_LR_SSH_Basic_563_022_20230626T020202_20230626T025218_PGD0_01.nc'
var_name ='ssh_karin' # must this be defined?


s = read_stars(fin, sub='ssh_karin', curvilinear=c('longitude','latitude'))

# all great lakes
#min_lon = -93
#max_lon = -75
#min_lat = 41
#max_lat = 49

# calval pass for ontario
min_lon = -93
max_lon = -75
min_lat = 41
max_lat = 49


crds = st_coordinates(s) # this is a 2XM vector of all (unraveled) points within the object s

# find the linear indices for points that fall within your bounding box
linear_inds = which(between(crds[,1]-360, min_lon, max_lon)  & between(crds[,2], min_lat, max_lat))

# convert to array indices
ij = arrayInd(linear_inds, .dim=dim(s))

# find extents
my_bbox = apply(ij, 2, range)

#print('full extents of desired region:')
#print(my_bbox)

print('now call ncks with `-d num_lines` using the following indices')
print(ncol(s) - rev(my_bbox[,2]) - 1)




