Example 3D-EMISH image stacks-file and rotating structure- movie files

The folder contains intermediate image data files for three representative chromatin structures (sID 50 for 1 domain, sID 12 for 2 domain, sID 42 for 3 domain). Intermediate image data files (8 for each structures) show image processing steps, from cropped raw EM image stacks to final processed image stacks, identifying sub-domain groups using different colors. In addition, rotational structure movie file is also provided.

For eight different image stacks, Z stack number starts from top to bottom.
Reference point (X=0, Y=0) is Top and Left corner at each plane.
Intermediate eight image data files during 3D-EMISH image processing follow:

[stack_cropped]: Cropped raw image stacks containing a chromatin folding structure.
X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 7x7x50 nm for replicate 1 and 5x5x30 nm for replicate 2.
sID 12, 42, and 50 are all from replicate 2. 

[stack_filtered]: Noise filtered out and resampled in z-direction to obtain the equal voxel size, 5x5x5 nm in XYZ for both replicate 1 and replicate 2.

[stack_cleared]: Unspecific signals in the background were removed. 
Each voxel (or X and Y pixel and Z stack interval ) in the image is 5x5x5 nm.

[binary_mask]: Binary mask of the structure connected components
Each voxel (or X and Y pixel and Z stack interval ) in the image is 5x5x5 nm.

[stack_masked]: Unspecific signals and the background noise were blacked out 
Each voxel (or X and Y pixel and Z stack interval ) in the image is 5x5x5 nm.

[stack_PA_coords]: X, Y, Z coordinates in stack_masked are from image collection. We calculate principal axes of a structure, and realign the image structure using new X, Y, Z coordinates, where X is the primary principal axis, and Y and Z are the secondary and tertiary principal axis, respectively. 
Each voxel (or X and Y pixel and Z stack interval ) in the image is 5x5x5 nm.

[stack_centers_PA_coords]: The mass centers were marked on an image stacks, where X, Y, Z are the principal axes coordinates.
Each voxel (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

[stack_colorized_PA_coords]: The identified subdomains of a structure were distinguished by different color in the principal axes coordinates.
Each voxel (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

