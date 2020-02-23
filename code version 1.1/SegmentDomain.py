import numpy
import os
from Misc import *
import get_gradient
import AnalysisParameters

class LogFile:
    def __init__(self,dirname):
        f=open(os.path.join(dirname,"logfile.txt"),"wb")
        self.file=f
    def write(self,par_name,par_value):
        self.file.write(par_name+"\t"+str(par_value)+"\n")
    def __del__(self):
        self.file.close()

class diffuse:
    def __init__(self,input_data,rho,coords,color=[0,0,0]):
        self.input_data=input_data
        self.rho=rho
        self.diffusion_data=numpy.zeros(rho.shape,dtype=numpy.double)
        self.derivative=numpy.zeros(rho.shape,dtype=numpy.double)
        self.allowed=numpy.ones(rho.shape,dtype=numpy.bool)
        self.color=color
        self.draw_gauss(coords,sigma=6.0)
        self.territory=None
        self.territory_old=None
    def diffuse_onetimestamp(self):
        get_gradient.gradient(self.diffusion_data,self.rho,self.derivative)
        #d=numpy.gradient(self.diffusion_data)
        #e=numpy.sum([numpy.gradient(d[0],axis=0),numpy.gradient(d[1],axis=1),numpy.gradient(d[2],axis=2)],axis=0)
        numpy.add(self.diffusion_data,self.derivative,self.diffusion_data)
        #numpy.add(self.diffusion_data,numpy.multiply(self.allowed,self.derivative),self.diffusion_data)

    def getsize(self):
        domain_support=numpy.multiply(self.diffusion_data>0,self.colorrange)
        size=numpy.sum(domain_support)
        self.size=size
        print "d_size=",size,numpy.sum(self.diffusion_data>0),numpy.sum(self.diffusion_data)
        self.mass=numpy.sum(numpy.multiply(self.diffusion_data>0,numpy.multiply(self.rho,self.colorrange)))
        return(size)

    def get_color(self):
        self.narray=numpy.multiply(AnalysisParameters.max_color/numpy.max(self.diffusion_data),self.diffusion_data)
        return(numpy.stack((numpy.multiply(self.color[0],self.narray),numpy.multiply(self.color[1],self.narray),numpy.multiply(self.color[2],self.narray)),axis=-1))

    def get_color_thresholded(self):
        #self.narray1=numpy.multiply(self.allowed,self.diffusion_data>0)
        self.diffusion_range=numpy.multiply(self.diffusion_data>0,self.allowed)
        self.domain_density=numpy.multiply(self.rho,self.diffusion_range)
        self.narray=numpy.multiply(self.domain_density,AnalysisParameters.max_color/numpy.max(self.domain_density))
        output=(numpy.stack((numpy.multiply(self.color[0],self.narray),numpy.multiply(self.color[1],self.narray),numpy.multiply(self.color[2],self.narray)),axis=-1))
        print output.shape
        return(output)

    def get_color_range(self):
        self.domain_density=numpy.multiply(self.diffusion_data>0,numpy.multiply(self.rho,self.colorrange))
        self.narray=numpy.multiply(self.domain_density,AnalysisParameters.max_color/numpy.max(self.domain_density))
        output=(numpy.stack((numpy.multiply(self.color[0],self.narray),numpy.multiply(self.color[1],self.narray),numpy.multiply(self.color[2],self.narray)),axis=-1))
        print output.shape
        return(output)

    def draw_gauss(self,coords,sigma):
        print "Marking mass centers..."
        sigma2=sigma**2
        m=range(-3*int(sigma),3*int(sigma)+1)
        for i in m:
            for j in m:
                for k in m:
                    self.diffusion_data[coords[0]+i,coords[1]+j,coords[2]+k]=numpy.exp(-(i**2+j**2+k**2)/sigma2 )
        self.diffusion_data=numpy.multiply(self.diffusion_data,(self.input_data>AnalysisParameters.background_skipped))
    def get_territory(self):
        self.territory=(self.diffusion_data==0)
