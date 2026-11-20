import earthaccess as ea


# shortname was determined by searching for a product on earthdata.nasa.gov 
# find granule > "view details" > "view more collection details" > find next to product
# can also find the string by navigating to the DOI for the product.  (e.g. doi.org/10.5067/SWOT-SSH-D) 
# note that the string is NOT the same as the string in the DOI (you have to go there and find it!)


ea.login()


res = ea.search_data(
    #short_name='SWOT_L2_LR_SSH_D',  
    #short_name='SWOT_L2_LR_SSH_D', 
    doi='10.5067/SWOT-SSH-D',
    #temporal=('2023-03-29','2023-07-10'),  # calval orbit period 
    #temporal=('2023-07-10','2027-01-01'),  # science orbit 
    temporal=('2023-07-10','2027-01-01'),  # science orbit 
    bounding_box=(-93, 41, -75, 49),
    granule_name='*Basic*' # this filtering might only work with one of short_name or doi (i forget which)
)





print(f'Found {len(res)} files.  Downloading now...')
files = ea.download(res, "../download", show_progress=True)



