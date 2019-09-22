import h5py,numpy,scipy,time,os
import scipy.ndimage
import scipy
from PIL import Image
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize
#import extensions2 as tiffextensions
from scipy.misc import imsave as imsave
from skimage.feature import peak_local_max
from skimage import data, img_as_float
from PIL import ImageDraw
import matplotlib.animation as animation
from operator import itemgetter, attrgetter
#import extensions2 as tiffextensions

in_file_name="12_nuc12syg1.tif"
in_file_name="50_nuc7syg1.tif"
relative_threshold=0.47
gauss_do_rozmycia=2.0
gauss_do_zapisu=1.0
scale_ratio=6.0
include_obj=1

def save_3dtiff(image,output,dtype="L"):
        plane_nr=image.shape[0]
        plane_y=image.shape[1]
        plane_x=image.shape[2]
        for i in range(plane_nr):
            plane=numpy.array(image[i,:,:],dtype=numpy.uint8)
	    #print numpy.max(plane)
        #tiffextensions.saveTiffMultipageFromSeq(list(self.image_sequence), output, rescaleSeqTo8bit=False,rgbOrder="RGB")
            #imsave(output+"_"+str(i)+".tif",plane)
            im = Image.frombytes(dtype,(plane_x,plane_y),plane.tostring())
            im.save(output+"_"+str(i)+".tif")

out_file_name=os.path.splitext(in_file_name)[0]+".cmap"
print out_file_name

def openTiff3D(stack):
    im =Image.open(stack)
    image_sequence=[]

    print "OPENING FRAMES..."
    frame_index = 0
    try:
      while 1:
        im.seek(frame_index)
        
        frame_index = frame_index + 1
        if (im!=None):
            
	    frame=im
            image_sequence.append(scipy.misc.fromimage(frame,flatten=True))
    except EOFError:
       print "eof after %d frames" %frame_index
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



class classify:
  def __init__(self,data,scale,filename=""):
    ratio=1.0*scale[2]/scale[0]
    self.scale=scale
    #data=openImage("test4.tif")
    #ratio=4.0
    self.filename=filename
    data=scipy.ndimage.interpolation.zoom(data,[ratio,1.0,1.0])
    
    #CENTER OF MASS
    
    print "Data size=",data.shape
    cube_side=int(1.1*numpy.sqrt(numpy.sum(numpy.multiply(data.shape,data.shape))))
    print "cube_side=",cube_side
    pad_size=map(lambda x: [-int(x/2),-int(x/2)],numpy.add(-cube_side,data.shape))
    print "pad_size=",pad_size
    data=numpy.pad(data,pad_size,'constant')
    z_cm,y_cm,x_cm=scipy.ndimage.measurements.center_of_mass(data)
    self.data=data
    print "Resized ata size=",data.shape
    print "Center of mass[z,y,x]=",z_cm,y_cm,x_cm
    self.z_cm=z_cm
    self.y_cm=y_cm
    self.x_cm=x_cm
    self.calculate()

  def tensor_core(self,r,r2,k,l):
            out=0
            #if k==l:
            #	out+=numpy.sum(numpy.multiply(r,r))
            #out-=r[k]*r[l]
            return(out)
            
  def inertia_tensor(self,data,origin):
            tensor=numpy.zeros([3,3])
            t0=time.time()
            coords_cm=[]
            for i in range(3):
                coords_cm.append(numpy.add(-origin[i],range(data.shape[i])))
            
            z_cm,y_cm,x_cm=numpy.meshgrid(coords_cm[0],coords_cm[1],coords_cm[2],indexing="ij")
            mx=numpy.multiply(data,x_cm)
            my=numpy.multiply(data,y_cm)
            mz=numpy.multiply(data,z_cm)
            mxx=numpy.sum(numpy.multiply(mx,x_cm))
            myy=numpy.sum(numpy.multiply(my,y_cm))
            mzz=numpy.sum(numpy.multiply(mz,z_cm))
            mxy=numpy.sum(numpy.multiply(mx,y_cm))
            mxz=numpy.sum(numpy.multiply(mx,z_cm))
            myz=numpy.sum(numpy.multiply(my,z_cm))
            #zz
            tensor[0,0]=myy+mxx
            #yy
            tensor[1,1]=mzz+mxx
            #xx
            tensor[2,2]=mzz+myy
            #xy
            tensor[2,1]=-mxy
            tensor[1,2]=-mxy
            #xz
            tensor[2,0]=-mxz
            tensor[0,2]=-mxz
            #yz
            tensor[1,0]=-myz
            tensor[0,1]=-myz
            
            
            
            return(tensor)
    
    
  def two_norm(self,x, *args):
        m1, m2, s1, s2, k1, k2, = args
        ret = k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1)
        ret += k2*scipy.stats.norm.pdf(x, loc=m2 ,scale=s2)
        return ret
  def save_tiff(self,output):
        plane_nr=self.image_sequence.shape[0]
        plane_y=self.image_sequence.shape[1]
        plane_x=self.image_sequence.shape[2]
        for i in range(plane_nr):
            plane=self.image_sequence[i,:,:]
        #tiffextensions.saveTiffMultipageFromSeq(list(self.image_sequence), output, rescaleSeqTo8bit=False,rgbOrder="RGB")
            #imsave(output+"_"+str(i)+".tif",plane)
            im = Image.frombytes('I;16',(plane_x,plane_y),plane.tostring())
            im.save(output+"_"+str(i)+".tif")
        
  def calculate(self):
        textsize=3
        filename=self.filename
        fig, axs = plt.subplots(2,2)
        fig2, axs2 = plt.subplots(1,3)
        fig3, axs3 = plt.subplots(1,3)
        ax1=axs2[0]
        ax2=axs2[1]
        ax3=axs2[2]
        ax10=axs3[0]
        ax11=axs3[1]
        ax12=axs3[2]
        ax1a=axs[0,0]
        ax2a=axs[0,1]
        ax3a=axs[1,0]
        ax20=axs[1,1]
        z_proj_max=numpy.max(self.data,axis=0)
        x_proj_max=numpy.max(self.data,axis=1)
        y_proj_max=numpy.max(self.data,axis=2).transpose()
        ax1.imshow(-z_proj_max,cmap='gray')
        ax2.imshow(-x_proj_max,cmap='gray')
        ax3.imshow(-y_proj_max,cmap='gray')
        #plt.show()
        #ORTHOGONALIZATION
        z_cm,y_cm,x_cm=[self.z_cm,self.y_cm,self.x_cm]
        tensor=self.inertia_tensor(self.data,[z_cm,y_cm,x_cm]) 
        print "Inertia tensor=",tensor
        e_values,e_vectors=numpy.linalg.eig(tensor)
        order=numpy.argsort(-e_values)
        print "order=",order
        print "order=",list(order)
        print "Eigenvalues=",e_values[order]
        print "Eigenvectors=",e_vectors[:,order].transpose()
        M=e_vectors.transpose()
        #M=[[1,0,0],[0,0.707,-0.707],[0,0.707,0.707]]
       
        order=list(order)
        M=numpy.array([M[order[0]],M[order[1]],M[order[2]]])
        M=M.transpose()
        shift=numpy.subtract([z_cm,y_cm,x_cm],numpy.dot(M,[z_cm,y_cm,x_cm]))
        print "M=",M
        print "origin=",[z_cm,y_cm,x_cm]
        print "shift=",shift


        #z_p,y_p,x_p=numpy.meshgrid(range(data.shape[0]),range(data.shape[1]),range(data.shape[2]),indexing="ij")


        data_new=scipy.ndimage.affine_transform(self.data,M,offset=shift)
        data_new_max=numpy.max(data_new)
        data_new_min=numpy.min(data_new)
        print "data_new_min,data_new_max=",data_new_min,data_new_max
        data_saved=numpy.multiply(65535.0/(data_new_max-data_new_min),numpy.add(-data_new_min,data_new))
        self.data_normalized=numpy.multiply(1.0/(data_new_max-data_new_min),numpy.add(-data_new_min,data_new))
        data_saved=numpy.array(data_saved,dtype=numpy.uint16)
        self.image_sequence=data_saved
        tensor=self.inertia_tensor(data_new,[z_cm,y_cm,x_cm]) 
        print "Tensor in new coordinates=",tensor
        print "data_new.shape=",data_new.shape
        #print data_new  
        self.output_shape=data_new.shape
        z_proj_max_a=numpy.max(data_new,axis=0)
        y_proj_max_a=numpy.max(data_new,axis=1)
        x_proj_max_a=numpy.max(data_new,axis=2).transpose()
        ax1a.imshow(-z_proj_max_a,cmap='gray')
        ax2a.imshow(-y_proj_max_a,cmap='gray')
        ax3a.imshow(-x_proj_max_a,cmap='gray')

        zx_hist=numpy.mean(z_proj_max_a,axis=0,dtype=numpy.float32)
        x=range(len(zx_hist))
        y=zx_hist
        y_max=numpy.max(zx_hist)
        print "y_max=",y_max
        params = [0.35*len(zx_hist), 0.65*len(zx_hist), 0.7*y_max, 0.7*y_max, 0.2*len(zx_hist), 0.2*len(zx_hist)]
        #params = [100, 250, 0.7*y_max, 0.7*y_max, 50, 50]
        fitted_params,_ = scipy.optimize.curve_fit(self.two_norm,x, y, p0=params)
        xx = np.linspace(np.min(x), np.max(x), 1000)
        print "fitted_params=",fitted_params
        print "********************"
        fitted_curve=self.two_norm(xx, *fitted_params)
        fitted_curve_diff=numpy.diff(fitted_curve)
        nr_of_extrema=0
        for i in range(len(fitted_curve_diff)-1):
            if ((fitted_curve_diff[i]>0)&(fitted_curve_diff[i+1]<=0))|((fitted_curve_diff[i]<0)&(fitted_curve_diff[i+1]>=0)):
                nr_of_extrema=nr_of_extrema+1
        print "nr_of_extrema=",nr_of_extrema
        domain_ratio=fitted_params[4]/fitted_params[5]
        if domain_ratio>1.0:
            domain_ratio=fitted_params[5]/fitted_params[4]
        print "filename=",filename
        fig.text(0.01,0.95,filename,size=textsize)
        fig.text(0.01,0.90,"nr_of_extrema="+str(nr_of_extrema),size=textsize)
        print "domain_ratio=",domain_ratio
        fig.text(0.01,0.85,"domain_ratio= %.4f" %domain_ratio,size=textsize)
        dist_centers=numpy.abs(fitted_params[0]-fitted_params[1])*self.scale[0]
        print "distance between centers [nm] = %.2f" %dist_centers
        fig.text(0.01,0.8,"distance between centers [nm] = %.2f" %dist_centers,size=textsize)
        print "domain widths [nm]",fitted_params[2]*self.scale[0],fitted_params[3]*self.scale[0]
        fig.text(0.01,0.75,"domain widths [nm]=%.2f , %.2f " %(fitted_params[2]*self.scale[0],fitted_params[3]*self.scale[0]),size=textsize)
        print "domain centers [nm]",fitted_params[0]*self.scale[0],fitted_params[1]*self.scale[0]
        fig.text(0.01,0.7,"domain centers [nm]=%.2f , %.2f " %(fitted_params[0]*self.scale[0],fitted_params[1]*self.scale[0]),size=textsize)
        #fig2 = plt.figure()
        #ax20 = fig2.add_subplot(111)
        ax20.plot(zx_hist,ls="-",marker=None)
        ax20.plot(xx, self.two_norm(xx, *fitted_params),ls="-",color="green")
        
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax3.set_xticks([])
        ax3.set_yticks([])
        
        
        
        
        
        ax1a.set_xticks([])
        ax1a.set_yticks([])
        ax1a.set_xlabel("x",x=0)
        ax1a.set_ylabel("y",y=0)
        ax2a.set_xticks([])
        ax2a.set_yticks([])
        ax2a.set_xlabel("x",x=0)
        ax2a.set_ylabel("z",y=0)
        ax3a.set_xticks([])
        ax3a.set_yticks([])
        ax3a.set_xlabel("z",x=0)
        ax3a.set_ylabel("y",y=0)
        ax1a.axhline(y=0.5*data_new.shape[2],c="red",linewidth=0.5,ls="dashed")
        ax1a.axvline(x=0.5*data_new.shape[0],c="blue",linewidth=0.5,ls="dotted")
        ax2a.axhline(y=0.5*data_new.shape[1],c="red",linewidth=0.5,ls="dashed")
        ax3a.axvline(x=0.5*data_new.shape[0],c="blue",linewidth=0.5,ls="dotted")
        self.fig1=fig
        self.fig2=fig2
        if nr_of_extrema==1:
            self.category=1
            
        ax10.plot(numpy.mean(z_proj_max_a,axis=0),ls="-",marker=None)
        ax11.plot(numpy.mean(y_proj_max_a,axis=0),ls="-",marker=None)
        ax12.plot(numpy.mean(x_proj_max_a,axis=0),ls="-",marker=None)

  def masscenters(self,data):
        par_size=27
        msize=int(par_size*5.0/self.scale[0])
        fig, axs = plt.subplots(1,3)
        image_max = scipy.ndimage.maximum_filter(data, size=msize, mode='constant')
        image_max_g=scipy.ndimage.filters.gaussian_filter(image_max,sigma=msize/2)
        coordinates = peak_local_max(image_max_g,threshold_rel=0.65,min_distance=msize)
        print coordinates
        self.masscenters=coordinates
        self.masscenters_nr=coordinates.shape[0]
        data_shape=data.shape
        im=[]
        im.append(numpy.max(data,axis=0))
        im.append(numpy.max(data,axis=1))
        im.append(numpy.max(data,axis=2).transpose())
        axs[0].imshow(-im[0],cmap='gray')
        axs[1].imshow(-im[1],cmap='gray')
        axs[2].imshow(-im[2],cmap='gray')
        #axs[1].imshow(numpy.max(image_max,axis=0),cmap='gray')
        #axs[2].imshow(numpy.max(image_max_g,axis=0),cmap='gray')
        axs[0].plot(coordinates[:, 2], coordinates[:, 1], 'r.',ms=10)
        axs[1].plot(coordinates[:, 2], coordinates[:, 0], 'r.',ms=10)
        axs[2].plot(coordinates[:, 0], coordinates[:, 1], 'r.',ms=10)
        plt.show()




        #axs[0].imshow(data,cmap='gray')
        #axs[2].imshow(data,cmap='gray')
        #axs[1].imshow(image_max,cmap='gray')

        #image_max_g=scipy.ndimage.filters.gaussian_filter(image_max,sigma=20)
        #coordinates = peak_local_max(image_max_g,threshold_rel=0.5)
        for i in range(3):
            shape=im[i].shape
            #axs[i].set_xlim(0,shape[1])
            #axs[i].set_ylim(0,shape[0])
            axs[i].set_xticks([])
            axs[i].set_yticks([])
        self.fig3=fig
  def draw_centers(self,output):
      img=numpy.zeros(self.output_shape)
      plane_nr=img.shape[0]
      frames=[]
      for p in range(plane_nr):
          plane=img[p,:,:]
          plane_y=img.shape[1]
          plane_x=img.shape[2]
          frame = Image.frombytes('P',(plane_x,plane_y),plane.tostring())
          frames.append(frame)
      for coordinate in self.masscenters:
          print coordinate
          self.draw_ball(frames,plane_nr,coordinate[2],coordinate[1],coordinate[0],5)
      for p in range(plane_nr):
        frames[p].save(output+"_"+str(p)+".tif")
	#print output+"_"+str(p)+".tif"
  def draw_ball(self,frames,plane_nr,x,y,z,r_max):
	for p in range(plane_nr):
		r=numpy.sqrt(r_max**2-(p-z)**2)
		draw = ImageDraw.Draw(frames[p])
		if r>0:
			print x-r,y-r,x+r,y+r
			draw.ellipse([x-r, y-r, x+r, y+r], fill=155)
        
class diffuse:
    def __init__(self,rho,coords,color=[0,0,0]):
        
        self.rho=rho
        self.array=numpy.zeros(rho.shape,dtype=numpy.double)
        self.color=color
        self.draw_gauss(coords,sigma=3.5)
    def diffuse_onetimestamp(self):
        gradient=numpy.gradient(self.array)
	for i in range(3):
		derivative=numpy.gradient(numpy.multiply(self.rho,gradient[i]),axis=i)
                numpy.add(self.array,derivative,self.array)
    def get_color(self):
        self.narray=numpy.multiply(max_color/numpy.max(self.array),self.array)
        return(numpy.stack((numpy.multiply(self.color[0],self.narray),numpy.multiply(self.color[1],self.narray),numpy.multiply(self.color[2],self.narray)),axis=-1))
    def draw_gauss(self,coords,sigma):
        sigma2=sigma**2
	m=range(-3*int(sigma),3*int(sigma)+1)
        for i in m:
            for j in m:
                for k in m:
                    self.array[coords[0]+i,coords[1]+j,coords[2]+k]=numpy.exp(-(i**2+j**2+k**2)/sigma2)
def updatefig(i,*args):
    global t
    print t
    domains=args
    #print domains
    t=t+1
    for domain in domains:
        domain.diffuse_onetimestamp()
    result=numpy.max(numpy.add(background,numpy.sum(map(lambda x:x.get_color(),domains),axis=0)),axis=0 )
    img.set_array(result)
    
    #tiffextensions.saveTiffMultipageFromSeq(result, "output"+str(t)+".tif", rescaleSeqTo8bit=False,rgbOrder="RGB")
    return img,
    



def calculate_wholedir(input_dir,output_dir):
    import os
    output_nr=0
    for root, dirs, files in os.walk(input_dir):
        for filename in files:
            
            if filename.endswith(".cmap"):
                 scale=0
                 filename_raw=os.path.splitext(filename)[0]
                 filename=(os.path.join(root, filename))
                 if filename.count("5x5x30"):
                         scale=[5,5,30]
                 if filename.count("7x7x50"):
                         scale=[7,7,50]
                 if scale!=0:
                    print "calculating for ",filename 
                    data=openImage(filename) 
                    output=classify(data,scale,filename)
                    output.masscenters(output.image_sequence)
                    output_nr=output_nr+1
                    output.fig1.savefig(os.path.join(output_dir,"fig"+str(output_nr)+".tif"),dpi=300)
                    output.fig3.savefig(os.path.join(output_dir,"centers"+str(output_nr)+".tif"),dpi=300)
                    dirname=str(output_nr)+"_tif"
                    try:
                        os.mkdir(dirname)
                    except:
                        pass
                    centers_dirname=str(output_nr)+"_centers"
                    try:
                        os.mkdir(centers_dirname)
                    except:
                        pass
                    output.save_tiff(os.path.join(dirname,str(output_nr)))
                    output.draw_centers(os.path.join(centers_dirname,str(output_nr)))
                 else:
                    print "skipping ",filename

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
    
ratio=6.0
data_orig=openImage(in_file_name)
data_orig=scipy.ndimage.interpolation.zoom(data_orig,[ratio,1.0,1.0])
#data_orig=openCmap("x.cmap")

data_orig=numpy.multiply(-1.0,data_orig)
data=gauss2d(data_orig,gauss_do_rozmycia)
#data=data_orig
data_orig_inv_raw=normalize(-data_orig)
data_orig=gauss2d(data_orig,sigma=gauss_do_zapisu)
data_orig=normalize(data_orig)
data_orig_inv=normalize(-data_orig)
data=normalize(data)
fig2, axs2 = plt.subplots(2,3)
fig3, axs3 = plt.subplots(2,3)
fig4, axs4 = plt.subplots(2,3)


ax1=axs2[0][0]
ax2=axs2[0][1]
ax3=axs2[0][2]
ax1a=axs2[1][0]
ax2a=axs2[1][1]
ax3a=axs2[1][2]
z_proj_max=numpy.max(data,axis=0)
x_proj_max=numpy.max(data,axis=1)
y_proj_max=numpy.max(data,axis=2)

print z_proj_max
ax1.imshow(z_proj_max,cmap='gray')
ax2.imshow(x_proj_max,cmap='gray',aspect=scale_ratio)
ax3.imshow(y_proj_max,cmap='gray',aspect=scale_ratio)

blobs = data > relative_threshold * numpy.max(data)
save_3dtiff(numpy.multiply(255.0/numpy.max(data_orig_inv),data_orig_inv),".\\output\\data_filtered.tif")

z_proj_thr=numpy.max(blobs,axis=0)
x_proj_thr=numpy.max(blobs,axis=1)
y_proj_thr=numpy.max(blobs,axis=2)
ax1a.imshow(z_proj_thr,cmap='gray')
ax2a.imshow(x_proj_thr,cmap='gray',aspect=scale_ratio)
ax3a.imshow(y_proj_thr,cmap='gray',aspect=scale_ratio)
s=numpy.ones([3,3,3])
labeled_array, num_features = scipy.ndimage.measurements.label(blobs,structure=s)

print num_features

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
    print "results=",results
    return(results)

largest_label=find_largest_component(labeled_array,num_features)
ax10=axs3[0][0]
ax11=axs3[0][1]
ax12=axs3[0][2]
ax20=axs3[1][0]
ax21=axs3[1][1]
ax22=axs3[1][2]
z_proj_labels=numpy.max(labeled_array,axis=0)
x_proj_labels=numpy.max(labeled_array,axis=1)
y_proj_labels=numpy.max(labeled_array,axis=2)
ax10.imshow(z_proj_labels,cmap='spectral')
ax11.imshow(x_proj_labels,aspect=scale_ratio,cmap='spectral')
ax12.imshow(y_proj_labels,aspect=scale_ratio,cmap='spectral')
#(plt.show()

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
ax20.imshow(z_proj_largest,cmap='gray')
ax21.imshow(x_proj_largest,cmap='gray',aspect=scale_ratio)
ax22.imshow(y_proj_largest,cmap='gray',aspect=scale_ratio)

masked_array = numpy.ma.masked_array(data_orig,(numpy.subtract(1.0,largest_object_mask)))
min_masked=masked_array.min()
data_masked=numpy.multiply(largest_object_mask,numpy.add(-min_masked,data_orig),dtype=numpy.float)
data_masked_inv=numpy.multiply(largest_object_mask,data_orig_inv)
data_masked_inv=numpy.multiply(255.0/numpy.max(data_masked_inv),data_masked_inv)
data_masked_background=numpy.multiply(255.0,numpy.logical_not(largest_object_mask))
data_masked_inv=numpy.add(data_masked_inv,data_masked_background)
data_masked=normalize(data_masked)
ax30=axs4[0][0]
ax31=axs4[0][1]
ax32=axs4[0][2]
ax40=axs4[1][0]
ax41=axs4[1][1]
ax42=axs4[1][2]

def get_limits(plane,a,offset=15):
    v1=numpy.max(plane,axis=a)
    for i,x in enumerate(v1):
        if x!=0:
            print i,x
            l=i
            break
            
    for i in range(len(v1)-1,-1,-1):
        x=v1[i]
        if x!=0:
            print i,x
            h=i
            break
    return([max([0,l-offset]),min([h+offset,len(v1)])])


z_proj_masked=numpy.max(data_masked,axis=0)
x_l,x_h=get_limits(z_proj_masked,1)
y_l,y_h=get_limits(z_proj_masked,0)

x_proj_masked=numpy.max(data_masked,axis=1)
z_l,z_h=get_limits(x_proj_masked,1)
y_proj_masked=numpy.max(data_masked,axis=2)
ax30.imshow(z_proj_masked[x_l:x_h,y_l:y_h],cmap='gray')
ax31.imshow(x_proj_masked[z_l:z_h,y_l:y_h],cmap='gray',aspect=scale_ratio)
ax32.imshow(y_proj_masked[z_l:z_h,x_l:x_h],cmap='gray',aspect=scale_ratio)
ax40.imshow(z_proj_max,cmap='gray')
ax41.imshow(x_proj_max,cmap='gray',aspect=scale_ratio)
ax42.imshow(y_proj_max,cmap='gray',aspect=scale_ratio)




save_3dtiff(data_masked_inv,".\\output\\data_masked.tif")

save_3dtiff(numpy.multiply(255.0/numpy.max(largest_object_mask),largest_object_mask),".\\output\\binary_mask.tif")
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
inspecific_grains_inv=numpy.logical_not(inspecific_grains)
inspecific_cleared=numpy.multiply(numpy.random.normal(background_avg,0.1,inspecific_grains.shape),inspecific_grains)
data_orig_inv_cleared=numpy.add(numpy.multiply(data_orig_inv_raw,inspecific_grains_inv),inspecific_cleared)
#data_orig_inv
#(background_avg,inspecific_grains)
data_orig_inv_cleared=gauss2d(data_orig_inv_cleared,sigma=gauss_do_zapisu)
save_3dtiff(numpy.multiply(255.0/numpy.max(data_orig_inv_cleared),data_orig_inv_cleared),".\\output\\data_cleared.tif")

f_out = h5py.File(out_file_name, "w")
grp = f_out.create_group("Chimera/image1")
data_masked=numpy.array(data_masked[z_l:z_h,x_l:x_h,y_l:y_h])
print data_masked.shape

dset = grp.create_dataset("data_zyx",data=data_masked)
f_out.close()
plt.show()