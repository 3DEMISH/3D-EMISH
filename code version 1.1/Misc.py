from PIL import Image
from pylab import *
import h5py,numpy,scipy,time,os
import AnalysisParameters

par_size=AnalysisParameters.minGranularity/AnalysisParameters.xy_scale
textsize=6

def normalize(M):
	M=numpy.multiply(M,numpy.array(M>=0,dtype=numpy.float))
	Mmin=M.min()
	Mmax=M.max()
	d=Mmax-Mmin
	print "Mmin=",Mmin
	print "Mmax=",Mmax
	res=numpy.multiply((1.0/d),numpy.add(-Mmin,M))
	return (res)

def saveTiffMultipageFromSeq(arrseq, fn, rescaleSeqTo8bit=False, rgbOrder="rgba", **params):
    """
    arrseq can be an iterator that yield 2D(grey) or 3D(color) image

    extension to PIL save TIFF

    if rescaleSeqTo8bit: scale each section (separately!) to 0..255
        (ratios between colors are unchanged)
    **params is directly forwarded to PIL save function
    """
#     if arr.ndim == 4:
#         if arr.shape[1] not in (1,2,3,4):
#             raise ValueError, "can save 4d arrays (color) only with second dim of len 1..4 (RG[B[A]])"
#     elif arr.ndim != 3:
#         raise ValueError, "can only save 3d (grey) or 4d (color) arrays"

    fp = open(fn, 'w+b')

    ifd_offsets=[]

#     if rescaleTo8bit:
#         mi,ma = float(arr.min()), float(arr.max())
#         ra = ma-mi

    params["_debug_multipage"] = True
    for z,a in enumerate(arrseq):
        if rescaleSeqTo8bit:
            mi,ma = float(a.min()), float(a.max())
            ra = ma-mi
            a=(a-mi)*255./ra
            ii = array2image(a.astype(N.uint8), rgbOrder=rgbOrder)
        else:
            ii = array2image(a, rgbOrder=rgbOrder)

        fp.seek(0,2) # go to end of file
        if z==0:
            # ref. PIL  TiffImagePlugin
            # PIL always starts the first IFD at offset 8
            ifdOffset = 8
        else:
            ifdOffset = fp.tell()

        ii.save(fp, format="TIFF", **params)
        
        if z>0: # correct "next" entry of previous ifd -- connect !
            ifdo = ifd_offsets[-1]
            fp.seek(ifdo)
            ifdLength = ii._debug_multipage.i16(fp.read(2))
            fp.seek(ifdLength*12,1) # go to "next" field near end of ifd
            fp.write(ii._debug_multipage.o32( ifdOffset ))

        ifd_offsets.append(ifdOffset)
    fp.close()

def save_3dtiff(image,output):
        plane_nr=image.shape[0]
        plane_y=image.shape[1]
        plane_x=image.shape[2]
	#saveTiffMultipageFromSeq(list(self.image_sequence), output, rescaleSeqTo8bit=False,rgbOrder="RGB")
        for i in range(plane_nr):
            plane=image[i,:,:]
        #tiffextensions.saveTiffMultipageFromSeq(list(self.image_sequence), output, rescaleSeqTo8bit=False,rgbOrder="RGB")
            #imsave(output+"_"+str(i)+".tif",plane)
            im = Image.frombytes('RGB',(plane_x,plane_y),plane.tostring())
            im.save(output+"_"+str(i)+".tif")

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
            if (im.mode!="L"):
                frame=im.convert("L")
            else:
                frame=im
            image_sequence.append(scipy.misc.fromimage(frame,flatten=True))
    except EOFError:
       print "eof after %d frames" %frame_index
       im.seek(0)
       pass
    return(numpy.array(image_sequence))

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

