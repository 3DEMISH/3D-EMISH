import h5py,numpy,scipy,time,os
import scipy.ndimage
import scipy
from PIL import Image
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize
from scipy.misc import imsave as imsave
from skimage.feature import peak_local_max
from skimage import data, img_as_float
from PIL import ImageDraw
import matplotlib.animation as animation
from operator import itemgetter, attrgetter
import denoiseSettings

in_file_name=denoiseSettings.in_file_name
relative_threshold=denoiseSettings.relative_threshold
gaussian_filter_px=denoiseSettings.gaussian_filter_image/denoiseSettings.xy_scale
gaussian_saving_px=denoiseSettings.gaussian_filter_image/denoiseSettings.xy_scale
scale_ratio=denoiseSettings.z_scale/denoiseSettings.xy_scale
aspect_ratio=1.0
include_obj=denoiseSettings.include_obj
try:
	os.makedirs(denoiseSettings.output_folder)
except:
	pass

def get_limits(plane,a,offset=denoiseSettings.padding):
    v1=numpy.max(plane,axis=a)
    for i,x in enumerate(v1):
        if x!=0:
            l=i
            break
            
    for i in range(len(v1)-1,-1,-1):
        x=v1[i]
        if x!=0:
            h=i
            break
    return([max([0,l-offset]),min([h+offset,len(v1)])])

def gauss2d(data,sigma=1.5):
    out=[]
    for plane in data:
        out.append(scipy.ndimage.filters.gaussian_filter(plane,sigma=sigma))
    return(numpy.array(out))
    
def normalize(data):
    data_max=numpy.max(data)
    data_min=numpy.min(data)
    delta=data_max-data_min
    data=numpy.add(-data_min,data)
    data=numpy.multiply(1.0/delta,data)
    return(data)

def save_3dtiff(image,output,dtype="L"):
        plane_nr=image.shape[0]
        plane_y=image.shape[1]
        plane_x=image.shape[2]
        for i in range(plane_nr):
            plane=numpy.array(image[i,:,:],dtype=numpy.uint8)
	    im = Image.frombytes(dtype,(plane_x,plane_y),plane.tostring())
            im.save(output+"_"+str(i)+".tif")
	    
def openTiff3D(stack):
    im =Image.open(stack)
    image_sequence=[]
    print "OPENING STACK..."
    frame_index = 0
    try:
      while 1:
        im.seek(frame_index)
        frame_index = frame_index + 1
        if (im!=None):
	    frame=im
            image_sequence.append(scipy.misc.fromimage(frame,flatten=True))
    except EOFError:
       print "%d FRAMES OPENED..." %frame_index
       im.seek(0)
       pass
    return(numpy.array(image_sequence,dtype=numpy.uint16))

def openCmap(filename):
    f = h5py.File(filename, 'r')
    a_group_key = f.keys()[0]
    data = list(f[a_group_key])
    data=numpy.array(f["Chimera"]["image1"]["data_zyx"].value)
    f=None
    a_group_key=None
    return(data)
    
def openImage(filename):
    ext=os.path.splitext(filename)[1]
    if ext==".cmap":
        return(openCmap(filename))
    if ext==".tif":
        return(openTiff3D(filename))
	
def find_largest_component(labels,num_features):
    out_size=0
    out_label=None
    results=[]
    for i in range(num_features):
        objectpixels=(labels==i)
        object_size=numpy.sum(objectpixels)
        if (object_size>out_size)&(i>0):
            out_size=object_size
            out_label=i
	if i>0:
		results.append([object_size,i])
    results=sorted(results, key=itemgetter(0),reverse=True)
    print "largest object size /vox/=",out_size,out_label
    return(results)

#opening image
print "Input file name=",in_file_name
out_file_name=os.path.splitext(in_file_name)[0]+".cmap"
print "Output file name=",out_file_name
data_orig=openImage(in_file_name)

#filtering and interpolating image
print "FILTERING AND INTERPOLATING..."
data_orig=scipy.ndimage.interpolation.zoom(data_orig,[scale_ratio,1.0,1.0])
data_orig=numpy.multiply(-1.0,data_orig)
data=gauss2d(data_orig,gaussian_filter_px)
data_orig_inv_raw=normalize(-data_orig)
data_orig=gauss2d(data_orig,sigma=gaussian_saving_px)
data_orig=normalize(data_orig)
data_orig_inv=normalize(-data_orig)
data=normalize(data)
z_proj_max=numpy.max(data,axis=0)
x_proj_max=numpy.max(data,axis=1)
y_proj_max=numpy.max(data,axis=2)

#thresholding 
print "THRESHOLDING IMAGE..."
blobs = data > relative_threshold * numpy.max(data)
if denoiseSettings.filtered_file_name!=None:
    save_3dtiff(numpy.multiply(255.0/numpy.max(data_orig_inv),data_orig_inv),os.path.join(denoiseSettings.output_folder,denoiseSettings.filtered_file_name))
z_proj_thr=numpy.max(blobs,axis=0)
x_proj_thr=numpy.max(blobs,axis=1)
y_proj_thr=numpy.max(blobs,axis=2)
s=numpy.ones([3,3,3])
labeled_array, num_features = scipy.ndimage.measurements.label(blobs,structure=s)
print "Number of thresholded objects=",num_features

#extracting objects
print "EXTRACTING CONNECTED COMPONENTS..."
largest_label=find_largest_component(labeled_array,num_features)
z_proj_labels=numpy.max(labeled_array,axis=0)
x_proj_labels=numpy.max(labeled_array,axis=1)
y_proj_labels=numpy.max(labeled_array,axis=2)
if include_obj==1:
	largest_object_mask=((labeled_array==largest_label[0][1]))
if include_obj==2:
	largest_object_mask0=((labeled_array==largest_label[0][1]))
	largest_object_mask1=((labeled_array==largest_label[1][1]))
	largest_object_mask=largest_object_mask0|largest_object_mask1	
if include_obj==3:
	largest_object_mask0=((labeled_array==largest_label[0][1]))
	largest_object_mask1=((labeled_array==largest_label[1][1]))
	largest_object_mask2=((labeled_array==largest_label[2][1]))
	largest_object_mask=largest_object_mask0|largest_object_mask1|largest_object_mask2
largest_object_mask=scipy.ndimage.morphology.binary_fill_holes(largest_object_mask)
z_proj_largest=numpy.max(largest_object_mask,axis=0)
x_proj_largest=numpy.max(largest_object_mask,axis=1)
y_proj_largest=numpy.max(largest_object_mask,axis=2)

#getting and applying structure mask
print "APPLYING MASK"
masked_array = numpy.ma.masked_array(data_orig,(numpy.subtract(1.0,largest_object_mask)))
min_masked=masked_array.min()
data_masked=numpy.multiply(largest_object_mask,numpy.add(-min_masked,data_orig),dtype=numpy.float)
data_masked_inv=numpy.multiply(largest_object_mask,data_orig_inv)
data_masked_inv=numpy.multiply(255.0/numpy.max(data_masked_inv),data_masked_inv)
data_masked_background=numpy.multiply(255.0,numpy.logical_not(largest_object_mask))
data_masked_inv=numpy.add(data_masked_inv,data_masked_background)
data_masked=normalize(data_masked)
z_proj_masked=numpy.max(data_masked,axis=0)
x_l,x_h=get_limits(z_proj_masked,1)
y_l,y_h=get_limits(z_proj_masked,0)
x_proj_masked=numpy.max(data_masked,axis=1)
z_l,z_h=get_limits(x_proj_masked,1)
y_proj_masked=numpy.max(data_masked,axis=2)

#showing image projection
fig2, axs2 = plt.subplots(2,3)
fig3, axs3 = plt.subplots(2,3)
fig4, axs4 = plt.subplots(2,3)
ax1=axs2[0][0]
ax2=axs2[0][1]
ax3=axs2[0][2]
ax1a=axs2[1][0]
ax2a=axs2[1][1]
ax3a=axs2[1][2]
ax1.imshow(z_proj_max,cmap='gray')
ax2.imshow(x_proj_max,cmap='gray',aspect=aspect_ratio)
ax3.imshow(y_proj_max,cmap='gray',aspect=aspect_ratio)
ax1a.imshow(z_proj_thr,cmap='gray')
ax2a.imshow(x_proj_thr,cmap='gray',aspect=aspect_ratio)
ax3a.imshow(y_proj_thr,cmap='gray',aspect=aspect_ratio)
ax10=axs3[0][0]
ax11=axs3[0][1]
ax12=axs3[0][2]
ax20=axs3[1][0]
ax21=axs3[1][1]
ax22=axs3[1][2]
ax10.imshow(z_proj_labels,cmap='Spectral')
ax11.imshow(x_proj_labels,aspect=aspect_ratio,cmap='Spectral')
ax12.imshow(y_proj_labels,aspect=aspect_ratio,cmap='Spectral')
ax20.imshow(z_proj_largest,cmap='gray')
ax21.imshow(x_proj_largest,cmap='gray',aspect=aspect_ratio)
ax22.imshow(y_proj_largest,cmap='gray',aspect=aspect_ratio)
ax30=axs4[0][0]
ax31=axs4[0][1]
ax32=axs4[0][2]
ax40=axs4[1][0]
ax41=axs4[1][1]
ax42=axs4[1][2]
ax30.imshow(z_proj_masked[x_l:x_h,y_l:y_h],cmap='gray')
ax31.imshow(x_proj_masked[z_l:z_h,y_l:y_h],cmap='gray',aspect=aspect_ratio)
ax32.imshow(y_proj_masked[z_l:z_h,x_l:x_h],cmap='gray',aspect=aspect_ratio)
ax40.imshow(z_proj_max,cmap='gray')
ax41.imshow(x_proj_max,cmap='gray',aspect=aspect_ratio)
ax42.imshow(y_proj_max,cmap='gray',aspect=aspect_ratio)

#saving binary mask and masked data
if denoiseSettings.masked_data_name!=None:
	save_3dtiff(data_masked_inv,os.path.join(denoiseSettings.output_folder,denoiseSettings.masked_data_name))
if denoiseSettings.binary_mask_name!=None:
	save_3dtiff(numpy.multiply(255.0/numpy.max(largest_object_mask),largest_object_mask),os.path.join(denoiseSettings.output_folder,denoiseSettings.binary_mask_name))

#getting inspecific grains
print "FILTERING BACKGROUND"
inspecific_grains=numpy.bitwise_xor(blobs,largest_object_mask)
background_mask=numpy.logical_not(blobs)
inspecific_grains_vol=numpy.sum(inspecific_grains)
background_mask_vol=numpy.sum(background_mask)
print "inspecific_grains_vol=",inspecific_grains_vol
print "background_mask_vol=",background_mask_vol
background_cum_signal=numpy.sum(numpy.multiply(background_mask,data_orig_inv_raw))
background_avg=(background_cum_signal/background_mask_vol)
background_dev=numpy.multiply(background_mask,numpy.subtract(data_orig_inv_raw,background_avg))
background_dev2=numpy.multiply(background_dev,background_dev)
background_std=(numpy.sqrt(numpy.sum(background_dev2)/background_mask_vol))
print "background_avg=",background_avg
print "background_std=",background_std

#clearing background
inspecific_grains_inv=numpy.logical_not(inspecific_grains)
inspecific_cleared=numpy.multiply(numpy.random.normal(background_avg,0.1,inspecific_grains.shape),inspecific_grains)
data_orig_inv_cleared=numpy.add(numpy.multiply(data_orig_inv_raw,inspecific_grains_inv),inspecific_cleared)
data_orig_inv_cleared=gauss2d(data_orig_inv_cleared,sigma=gaussian_saving_px)
if denoiseSettings.cleared_background_name:
	save_3dtiff(numpy.multiply(255.0/numpy.max(data_orig_inv_cleared),data_orig_inv_cleared),os.path.join(denoiseSettings.output_folder,denoiseSettings.cleared_background_name))

#saving chimera file
f_out = h5py.File(out_file_name, "w")
grp = f_out.create_group("Chimera/image1")
data_masked=numpy.array(data_masked[z_l:z_h,x_l:x_h,y_l:y_h])
dset = grp.create_dataset("data_zyx",data=data_masked)
f_out.close()
print "OUTPUT SAVED..."

plt.show()