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
import get_gradient
import matplotlib.animation as animation
from Misc import *
from SegmentDomain import *
import shutil
import AnalysisParameters

colors=AnalysisParameters.colors
sorted_colors=[]

class classify:

    def __init__(self,data,scale,filename=""):
        self.ratio=1.0*scale[2]/scale[0]
        self.scale=scale
        self.filename=filename
        print "Input file name=",filename
        if AnalysisParameters.reinterpolate:
            self.data=scipy.ndimage.interpolation.zoom(data,[self.ratio,1.0,1.0])
	else:
	    self.data=data
        self.center_structure()

    def _two_norm(self,x, *args):
        m1, m2, s1, s2, k1, k2, = args
        ret = k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1)
        ret += k2*scipy.stats.norm.pdf(x, loc=m2 ,scale=s2)
        return ret

    def _save_tiff(self,output):
        plane_nr=self.image_sequence.shape[0]
        plane_y=self.image_sequence.shape[1]
        plane_x=self.image_sequence.shape[2]
        try:
            os.mkdir(os.path.join(self.output_dir,output))
        except:
            pass
        for i in range(plane_nr):
            plane=self.image_sequence[i,:,:]
            im = Image.frombytes('I;16',(plane_x,plane_y),plane.tostring())
            im.save(os.path.join(self.output_dir,output,output+"_"+str(i)+".tif"))

    def _plot_pa_projections(self):
	filename=self.filename
        fig, axs = plt.subplots(2,2)
	fig.suptitle('projections in PA coords /avg. density/', fontsize=16)
        fig_z = plt.figure(frameon=False)
	fig.suptitle('Orthogonal projections in PA coords /avg. density/', fontsize=16)
        ax_z = plt.Axes(fig_z, [0., 0., 1., 1.])
        ax_z.set_axis_off()
        fig_z.add_axes(ax_z)
        ax1a=axs[0,0]
        ax2a=axs[0,1]
        ax3a=axs[1,0]
        ax20=axs[1,1]
        ax1a.imshow(-self.z_proj_mean_pa,cmap='gray')
        ax_z.imshow(-self.z_proj_mean_pa,cmap='gray')
        ax2a.imshow(-self.y_proj_mean_pa,cmap='gray')
        ax3a.imshow(-self.x_proj_mean_pa,cmap='gray')
        zx_hist=numpy.mean(self.z_proj_max_pa,axis=0,dtype=numpy.float32)
        ax20.plot(zx_hist,ls="-",marker=None)
        fig.text(0.01,0.95,filename,size=textsize)
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
        ax1a.axhline(y=0.5*self.data_PAcoords.shape[2],c="red",linewidth=0.5,ls="dashed")
        ax1a.axvline(x=0.5*self.data_PAcoords.shape[0],c="blue",linewidth=0.5,ls="dotted")
        ax2a.axhline(y=0.5*self.data_PAcoords.shape[1],c="red",linewidth=0.5,ls="dashed")
        ax3a.axvline(x=0.5*self.data_PAcoords.shape[0],c="blue",linewidth=0.5,ls="dotted")
        self.fig_pa_projections=fig
        self.fig_pa_projections_z=fig_z
        self.PA_axes=[ax1a,ax2a,ax3a]

    def _join_spots_with_domains(self,domains):
        spots_with_domains=[]
        for line in self.spots_coordinates:
            domain_intensities=[]
            for domain in domains:
                intensity=domain.diffusion_data[line[0],line[1],line[2]]
                domain_intensities.append(intensity)
            domain_nr=numpy.argmax(domain_intensities)
            spots_with_domains.append(list(line)+[domain_nr])
        self.spots_with_domains=spots_with_domains


    def _write_spots_coordinates(self):
        try:
            data=self.spots_with_domains
        except:
            self._join_spots_with_domains(self.domains)
            data=self.spots_with_domains
        out_file=open(os.path.join(self.output_dir,'spots_data.txt'), 'w')
        for line in data:
            out_file.write(str(line)+"\n")
        out_file.close()

    def _remove_double_centers(self,centers):
        #removes multiple mass centers at the same loaction
        included_centers=[]
        for center in centers:
            include_flag=1
            for included_center in included_centers:
                dist=self.norm2(center,included_center)
                if dist<5:
                    include_flag=0
            if include_flag:
                included_centers.append(center)
        return(numpy.array(included_centers))

    def _put_masscenters(self,axs):
        coordinates=self.masscenters
        axs[0].plot(coordinates[:, 2], coordinates[:, 1], 'r.',ms=10)
        axs[1].plot(coordinates[:, 2], coordinates[:, 0], 'r.',ms=10)
        axs[2].plot(coordinates[:, 0], coordinates[:, 1], 'r.',ms=10)

    def _draw_centers(self,output):
        #draws mass centers
        print "DRAWING MASS CENTERS"
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
            self._draw_ball(frames,plane_nr,coordinate[2],coordinate[1],coordinate[0],5)
        for p in range(plane_nr):
            frames[p].save(output+"_"+str(p)+".tif")

    def _draw_ball(self,frames,plane_nr,x,y,z,r_max):
        for p in range(plane_nr):
            r=numpy.sqrt(r_max**2-(p-z)**2)
            draw = ImageDraw.Draw(frames[p])
            if r>0:
                print x-r,y-r,x+r,y+r
                draw.ellipse([x-r, y-r, x+r, y+r], fill=155)

    def center_structure(self):
        #Calculating center of mass and centering
        ratio=self.ratio
        data=self.data
        print "STRUCTURE CENTERING..."
        print "Data size=",data.shape
        cube_side=AnalysisParameters.cube_side
        print "Desired cube_side=",cube_side
        pad_size=[]
        for s in data.shape:
            left_pad=int(cube_side-s)/2
            right_pad=cube_side-s-left_pad
            pad_size.append([left_pad,right_pad])
        print "Extra pad size=",pad_size
        data=numpy.pad(data,pad_size,'constant')
        z_cm,y_cm,x_cm=scipy.ndimage.measurements.center_of_mass(data)
        self.data=data
        print "Resized stack size=",data.shape
        print "Center of mass coordinates[z,y,x]=",z_cm,y_cm,x_cm
        self.z_cm=z_cm
        self.y_cm=y_cm
        self.x_cm=x_cm
        self.align()

    def inertia_tensor(self,data,origin):
        print "CALCULATING INERTIA TENSOR..."
        tensor=numpy.zeros([3,3])
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

    def align(self):
        #Structure aligning
        print "ALIGNING STRUCTURE..."
        textsize=3
        filename=self.filename
        #Orthogonalization
        print "ORTHOGONALIZATION..."
        z_cm,y_cm,x_cm=[self.z_cm,self.y_cm,self.x_cm]
        tensor=self.inertia_tensor(self.data,[z_cm,y_cm,x_cm])
        e_values,e_vectors=numpy.linalg.eig(tensor)
        order=numpy.argsort(-e_values)
        M=e_vectors.transpose()
        order=list(order)
        M=numpy.array([M[order[0]],M[order[1]],M[order[2]]])
        M=M.transpose()
        shift=numpy.subtract([z_cm,y_cm,x_cm],numpy.dot(M,[z_cm,y_cm,x_cm]))
        print "Inertia tensor=",tensor
        print "Axis order=",list(order)
        print "Eigenvalues=",e_values[order]
        print "Eigenvectors=",e_vectors[:,order].transpose()
        print "M=",M
        print "origin=",[z_cm,y_cm,x_cm]
        print "shift=",shift
        data_PAcoords=scipy.ndimage.affine_transform(self.data,M,offset=shift)
        data_PAcoords=numpy.multiply(data_PAcoords,numpy.array(data_PAcoords>=0,dtype=numpy.float))
        data_PAcoords_max=numpy.max(data_PAcoords)
        data_PAcoords_min=numpy.min(data_PAcoords)
        print "data_PAcoords_min,data_PAcoords_max=",data_PAcoords_min,data_PAcoords_max
        print "lowest brightness=",data_PAcoords[0][0][0]
        data_saved=numpy.multiply(65535.0/(data_PAcoords_max-data_PAcoords_min),numpy.add(-data_PAcoords_min,data_PAcoords))
        self.data_normalized=numpy.multiply(1.0/(data_PAcoords_max-data_PAcoords_min),numpy.add(-data_PAcoords_min,data_PAcoords))
        data_saved=numpy.array(data_saved,dtype=numpy.uint16)
        self.image_sequence=data_saved
        print "SWITCHING TO NEW COORDINATES..."
        tensor=self.inertia_tensor(data_PAcoords,[z_cm,y_cm,x_cm])
        print "Tensor in new coordinates=",tensor
        print "data_PAcoords cube=",data_PAcoords.shape
        self.output_shape=data_PAcoords.shape
        self.z_proj_max_pa=numpy.max(data_PAcoords,axis=0)
        self.y_proj_max_pa=numpy.max(data_PAcoords,axis=1)
        self.x_proj_max_pa=numpy.max(data_PAcoords,axis=2).transpose()
        self.z_proj_mean_pa=numpy.mean(data_PAcoords,axis=0)
        self.y_proj_mean_pa=numpy.mean(data_PAcoords,axis=1)
        self.x_proj_mean_pa=numpy.mean(data_PAcoords,axis=2).transpose()
        fig2, axs2 = plt.subplots(1,3)
	fig2.suptitle('orthogonal projections in original coords /min. density/', fontsize=16)
        fig3, axs3 = plt.subplots(1,3)
	fig3.suptitle('linear projections in PA coords /avg. density/', fontsize=16)
        ax1=axs2[0]
        ax2=axs2[1]
        ax3=axs2[2]
        ax10=axs3[0]
        ax11=axs3[1]
        ax12=axs3[2]
        z_proj_max=numpy.max(self.data,axis=0)
        x_proj_max=numpy.max(self.data,axis=1)
        y_proj_max=numpy.max(self.data,axis=2)
        ax1.imshow(-z_proj_max,cmap='gray')
        ax2.imshow(-x_proj_max,cmap='gray')
        ax3.imshow(-y_proj_max,cmap='gray')
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax10.plot(numpy.mean(self.z_proj_mean_pa,axis=0),ls="-",marker=None)
        ax11.plot(numpy.mean(self.x_proj_mean_pa,axis=1),ls="-",marker=None)
        ax12.plot(numpy.mean(self.x_proj_mean_pa,axis=0),ls="-",marker=None)
        self.fig2=fig2
        self.data_PAcoords=data_PAcoords
        self._plot_pa_projections()

    @staticmethod
    def save_3dtiff(image,output):
        plane_nr=image.shape[0]
        plane_y=image.shape[1]
        plane_x=image.shape[2]
        for i in range(plane_nr):
            plane=image[i,:,:]
            im = Image.frombytes('RGB',(plane_x,plane_y),plane.tostring())
            im.save(output+"_"+str(i)+".tif")

    @staticmethod
    def norm2(x,y):
        delta=numpy.subtract(x,y)
        delta2=numpy.multiply(delta,delta)
        norm2=numpy.sum(delta2)
        return(norm2)

    def get_spots_coords(self,data):
        par_size=1.0
        msize=int(par_size)
        threshold=numpy.max(data)*0.2
        self.immax=numpy.max(data)
        image_max = scipy.ndimage.maximum_filter(data, size=msize, mode='constant')
        image_max_g=scipy.ndimage.filters.gaussian_filter(image_max,sigma=msize/2)
        coordinates = peak_local_max(data,threshold_rel=0.65,min_distance=msize)
        self.spots_coordinates=coordinates

    def plot_spots(self,data=None):
        global sorted_colors
        if data==None:
            data=self.image_sequence
        try:
            coordinates=self.spots_with_domains
        except(AttributeError):
            self._join_spots_with_domains(self.domains)
            coordinates=self.spots_with_domains
        fig_proj, axs_proj = plt.subplots(1,1)
        axs_proj.imshow(-self.z_proj_mean_pa,cmap='gray')

        for plane_nr in range(data.shape[0]):
            plane=data[plane_nr,:,:]
            fig, axs = plt.subplots(1,1)
            axs.imshow(plane,vmin=0, vmax=self.immax,cmap='gray')
            for line in coordinates:
                if line[0]==plane_nr:
                    axs.plot([line[2]], [line[1]], marker='.',ms=1.0,)
                    axs_proj.plot([line[2]], [line[1]],  marker='.',ms=1.0)
            plt.close(fig)
        fig_proj.savefig(os.path.join(self.output_dir,"spots_projection"),dpi=300)
        plt.close(fig_proj)

    def masscenters(self,data):
        #mass centers visualization
        data_filtered=scipy.ndimage.filters.gaussian_filter(data,0.0)
        msize=int(par_size*5.0/self.scale[0])
        fig, axs = plt.subplots(1,3)
        image_max = scipy.ndimage.maximum_filter(data_filtered,size=msize, mode='constant')
        image_max_g=scipy.ndimage.filters.gaussian_filter(image_max,sigma=msize/2)
        coordinates = peak_local_max(image_max_g,threshold_rel=0.65,min_distance=msize)
        coordinates=self._remove_double_centers(coordinates)
        print "mass centers coordinates=",coordinates
        self.masscenters=coordinates
        self.masscenters_nr=coordinates.shape[0]
        data_shape=data.shape
        im=[]
        im.append(numpy.max(data,axis=0))
        im.append(numpy.max(data,axis=1))
        im.append(numpy.max(data,axis=2).transpose())
        axs[0].imshow(numpy.multiply(-1.0,im[0]),cmap='gray')
        axs[0].autoscale(False)
        axs[1].imshow(numpy.multiply(-1.0,im[1]),cmap='gray')
        axs[1].autoscale(False)
        axs[2].imshow(numpy.multiply(-1.0,im[2]),cmap='gray')
        axs[2].autoscale(False)
        self._put_masscenters(axs)
        self._put_masscenters(self.PA_axes)
        for i in range(3):
            shape=im[i].shape
            axs[i].set_xticks([])
            axs[i].set_yticks([])
        self.fig3=fig




def perform_diffusion(i,*args):
    global nr_of_steps_nochange
    global t
    global sorted_colors
    global background
    domains=args[0]
    output_dir=args[1]
    t=t+1
    print "step:",t

    for domain in domains:
        domain.get_territory()
        change=0
        if (type(domain.territory_old)!=type(None))&(type(domain.territory)!=type(None)):
            
            d_change=numpy.sum(numpy.abs(numpy.logical_xor(domain.territory_old,domain.territory)))
            if d_change>0:
                change=1
        domain.territory_old=domain.territory
    for i,domain1 in enumerate(domains):
        for j,domain2 in enumerate(domains):
            if i!=j:
                domain1.allowed=domain1.allowed&domain2.territory
    for domain in domains:
        domain.diffuse_onetimestamp()


    if change==0:
        nr_of_steps_nochange=nr_of_steps_nochange+1
    else:
        nr_of_steps_nochange=0

    if (t==AnalysisParameters.max_step_nr)|(AnalysisParameters.stop_diffusion&(nr_of_steps_nochange>5)):

        image_domain_index=[]
        for domain in domains:
            image_domain_index.append(domain.diffusion_data)
        image_domain_index=numpy.argmax(image_domain_index,axis=0)

        for i,domain in enumerate(domains):
            domain.colorrange=(image_domain_index==i)

        domain_masses=[]
        for i,domain in enumerate(domains):
            domain.getsize()
            domain_masses.append(domain.size)

        sorted_domains=numpy.argsort(domain_masses)[::-1]

        for i,d  in enumerate(sorted_domains):
            domain=domains[d]
            domain.color=colors[i]

        data3d=numpy.multiply(background,numpy.sum(map(lambda x:x.get_color_range(),domains),axis=0))
        fig = plt.figure(frameon=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        data3d_proj=numpy.max(data3d,axis=0)
        print "numpy.max(data3d_proj)=",numpy.max(data3d_proj)
        data3d_proj=numpy.multiply(1.0/numpy.max(data3d_proj),data3d_proj)
        ax.imshow(data3d_proj)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.set_frame_on(False)
        fig.savefig(os.path.join(output_dir,"result.tif"), bbox_inches='tight', pad_inches=0)
        data3d=numpy.multiply(1.0/numpy.max(data3d),data3d)
        classify.save_3dtiff(numpy.array(numpy.multiply(255,data3d),dtype=numpy.uint8),os.path.join(output_dir,"output"))
        return(0)
    return(1)


def analyze_structure(filename,scale,output_dir):
    global t
    global nr_of_steps_nochange
    global background
    nr_of_steps_nochange=0
    t=0
    data=openImage(filename)
    data=scipy.ndimage.filters.gaussian_filter(data,sigma=0)
    output=classify(data,scale,filename)
    output.output_dir=output_dir
    try:
        os.mkdir(output.output_dir)
        time.sleep(0.5)
        shutil.copyfile(filename, os.path.join(output.output_dir,os.path.basename(filename)))
    except:
        pass
    output.masscenters(output.image_sequence)
    output.get_spots_coords(output.image_sequence)
    output.fig_pa_projections.savefig(os.path.join(output.output_dir,"orth_proj_avg.tif"),dpi=300)
    output.fig_pa_projections_z.savefig(os.path.join(output.output_dir,"z_proj_avg.tif"),dpi=300)
    output.fig3.savefig(os.path.join(output.output_dir,"mass_centers.tif"),dpi=300)
    output._save_tiff("pa_coords")
    output_dir=output.output_dir
    f_prefix="diffused"
    if AnalysisParameters.show_projections:
        plt.show()
    fig7, ax7 = plt.subplots()
    background=numpy.multiply((1.0-AnalysisParameters.max_color),numpy.stack((output.data_normalized,output.data_normalized,output.data_normalized),axis=-1))
    rho_diffuse=numpy.multiply(0.7,numpy.multiply(output.data_normalized,output.data_normalized>AnalysisParameters.background_skipped))
    domains=[]
    log=LogFile(output_dir)
    log.write("scale=",scale)
    log.write("filename",filename)
    log.write("domains_nr",output.masscenters_nr)
    for domain_nr in range(output.masscenters_nr):
        diffusion=diffuse(output.data_normalized,rho_diffuse,output.masscenters[domain_nr],colors[domain_nr])
        domains.append(diffusion)
    for i in range(AnalysisParameters.max_step_nr+1):
        if perform_diffusion(0,domains,output.output_dir)==0:
            break
    for i,domain in enumerate(domains):
        log.write("domain_%s_size"%(str(i)),domain.getsize())
        numpy.save(os.path.join(output.output_dir,"domain"+str(i)+AnalysisParameters.domain_data), domain.domain_density)
    output.domains=domains
    output.plot_spots()
    output._write_spots_coordinates()
    del log
    print "RESULTS SAVED..."


analyze_structure(AnalysisParameters.input_file,[AnalysisParameters.xy_scale,AnalysisParameters.xy_scale,AnalysisParameters.z_scale],AnalysisParameters.output_dir)
