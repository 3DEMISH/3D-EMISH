#input file name
input_file='ex12.cmap'

#output directory
output_dir=".\\output-analysis"

#xy scale [nanometers]
xy_scale=5

#z scale [nanometers]
z_scale=30

#granularity cott-off[nanometers]
minGranularity=135

#RGB colors for domains
colors=[[1.,0,1.],[0,1.,0],[0,1.,1.],[0.7,0.2,0.7],[1,1,0],[1,0,1]]

#desired cube size [px]of aligned structure
cube_side=400

#show structure projections in principle axis coordinates (True/False)
show_projections=True

#reinterpolate the stack if voxel shape is not cubic
reinterpolate=False

#segmented domain data file output [numpy.matrix]
domain_data="domain_data.npy"

#maximal number of diffusion steps
max_step_nr=200

#stops diffusion when all voxels filled (1), otherwise(0) continues to max_step_nr
stop_diffusion=False

#do not diffuse into backgroung below threshold (0-1)
background_skipped=0.05

#max color intensity in overlay
max_color=0.6
