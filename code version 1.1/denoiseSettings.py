#input file to denoise
in_file_name="ex12.tif"

output_folder=".\\output"

#output (filtered) tiffs name (None or file path)
filtered_file_name="data_filtered.tif"

#output masked tiffs  (None or file path)
masked_data_name="data_masked.tif"

#cleared background tiffs (None or file path)
cleared_background_name="cleared_background.tif"

#output masked binary mask  (None or file path)
binary_mask_name="binary_mask.tif"

# threshold value to discard signal below it (range 0-1)
relative_threshold=0.47

#xy_scale (nanometers)
xy_scale=5.0

#z_scale (nanometers)
z_scale=30.0

# gaussian filter size- applied to output image
gaussian_filter_image=5.0

# gaussian filter size- used in compenent extraction
gaussian_filter_for_extraction=10.0

#number of objects to extract (1-3)
include_obj=1

#extra paddig [pixels] to add
padding=15