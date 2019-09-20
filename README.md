# 3D-EMISH
Ultra-Resolution Visualization of Sequence-Specific 3D Chromatin Folding Structures

The entire original 3D-EMISH stack file size is about 12TB, not possible to deposit all data files into a public depository space. We first deposit eight intermediate image data files to show our image processing steps, from cropped raw EM image stacks to final processed image stacks, identifying sub-domain groups using different colors. We deposit the three representative chromatin structures (sID 50 for 1 domain, sID 12 for 2 domain, sID 42 for 3 domain). Second, we will deposit cropped raw 3D-EMISH image files for all 229 structures with our image-processing code/script.  

Intermediate eight image data files during 3D-EMISH image processing follow:

1. [stack_cropped]: cropped raw image stacks containing a single structure
   Z stack numbers are aligned by top to bottom
   Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x30 nm

2. [stack_filtered]: image stacks containing a single structure- pixel noise filtered out and resampled in z-direction to obtain the equal voxel sides
   Z stack numbers are aligned by top to bottom
   Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

3. [stack_cleared]: image stacks with unspecific signal in the background removed
   Z stack numbers are aligned by top to bottom
   Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

4. [binary_mask]: binary mask of the structure connected component 
   Z stack numbers are aligned by top to bottom
   Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

5. [stack_masked]: image stacks with unspecific signal and the background blacked out
   Z stack numbers are aligned by top to bottom
   Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

6. [stack_PA_coords]: image stacks with unspecific signal and the background blacked out in principle axis coordinates
  stack aligned according to the structure principle axis  (x-axis), y and z axis are secondary and tercary principle axis
  Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

7. [stack_centers_PA_coords]: the mass centers marked on an image stacks in principle axis coordinates
  stack aligned according to the structure principle axis  (x-axis), y and z axis are secondary and tercary principle axis
  Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm

8. [stack_colorized_PA_coords]: the colorized subdomains in principle axis coordinates
  stack aligned according to the structure principle axis  (x-axis), y and z axis are secondary and tercary principle axis
  Reference point (X=0, Y=0) is Top and Left corner at each plane
   X, Y, Z voxel size (or X and Y pixel and Z stack interval ) in the image is 5x5x5nm
