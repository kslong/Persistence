#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Find the peaks in an image of a fits file or datset


Command line usage (if any):

		usage:	peaks.py  filename
			peaks.py  dataset 

	
			 

Description:  

	The purpose of these routines is to locate some
	regions of a persistence image that can be used
	to evaluate how well persistence subtraction 
	has been done.

	As called from the command line the program reads
	a fits file or a dataset name.  It decides based
	on whether the name contains the sequence 'fit'
	whether it is a fits file or a dataset. 
	
	If it is a dataset it locates the persist file for
	the dataset by reading the list file.  In this case it will
	assume the normal direcgtory structure needed
	for pipleine operation of the code
	
	In either case smooths file using boxcar
	function.  It then atttemps to find the peaks in the smoothed
	image.

	The results are written to a file to be used
	in the next stage of the process.  A number of
	figures are also produced



Primary routines:

Notes:
									   
History:

101207 ksl Coding begun

'''

import sys
import os
import string
import numpy
# This is part of the stsci modules; at one point it was part of scipy but not today
from stsci import convolve 
import pylab
from matplotlib import patches
import per_list
import per_fits


def smooth(image,size):
	'''
	Smooth a 2d image array using 
	a boxcar of size given by size
	'''
	# print 'smooth',numpy.max(image),numpy.min(image),numpy.median(image)
	z=convolve.boxcar(image,(size,size),mode='constant',cval=0.0)
	# print 'smooth',numpy.max(z),numpy.min(z),numpy.median(z)
	return z

def find_peaks(image,mask='none',sep=20,maxno=20):
	'''
	Find brightest peaks in an image that are separated by at least
	"sep" pixels.  

	A list with y,x,peak counts is returned

	mask is normally an array that has the same shape as image, but 
	one can find the peaks when there is no mask.  The way to indicate
	there is no mask is to set mask to 'none'

	The number of peaks returned is limited to maxno.  

	110203	ksl Fixed the way mask is handled.
	111011	ksl Modify the logic to try to do a better job eliminating
		    persistence peaks near bright sources.  The basic change
		    was to first find the peaks regardless of the mask (e.g
		    the original science image, and then to eliminate peaks
		    later.  What had been happening is that we were finding
		    subsidiary peaks in near a bright object.
	'''


	ysize,xsize=image.shape

	# create arrays containing indices into the persistence image
	y,x=numpy.mgrid[0:ysize,0:xsize]

	y=numpy.ravel(y)
	x=numpy.ravel(x)
	im=numpy.ravel(image)
	# print 'test',min(im),max(im)
	zz=numpy.argsort(im)

	# At this point zz provides an index into the array im
	# in order of increasing value in im

	sep2=sep*sep   # This is the closest we are going to allow two "sources"
	
	# We are now going to march through the array trying
	# to find places where there is persistence which
	# are not lying under bright stars in the current image
	# and which are separated from one another.


	sources=[]
	i=len(zz)

	# First find all of the candidate peaks regardless of the mask
	while i>0.8*len(zz) and len(sources)<maxno:
		i=i-1
		ii=zz[i]
		new_source='yes'
		
		# Do not allow sources near the edges of the detector
		if 50>x[ii] or x[ii]>xsize-50:
			continue
		if 50>y[ii] or y[ii]>ysize-50:
			continue

		for source in sources:
			# print source
			yy=source[0]-y[ii]
			xx=source[1]-x[ii]
			if xx*xx+yy*yy < sep2:
				# print 'not a source',source[0],source[1]
				new_source='no'
				break
		if new_source=='yes':
			# print 'Adding ',x[ii],y[ii],im[ii]
			sources.append([y[ii],x[ii],im[ii],mask[y[ii],x[ii]]])
			# print 'Trial source at %3d %3d %8.3f %8.3f'% (x[ii],y[ii],image[y[ii],x[ii]],mask[y[ii],x[ii]])

	# Each record of source  contains an x,y position, the value in the persistence image and the value in the mask
	if mask=='none':
		print 'peaks: No mask! Found are %d peaks in smoothed persistence image with flux > %f' % (len(sources),im[i])
		return sources

	# print 'test: At this point we have %d soures' % len(sources)

	# Note that this is slightly dangerous since we are using the model to exclude points as well as the mask
	# sources contains, the x and y position, and the flux in the image and the mask
	# 111129 - ksl - Attempted to improve the selection of good sources to compare persistence to
	xmed=numpy.median(mask)
	# print 'shape',numpy.shape(mask)
	xsources=[]
	isize=5
	for one in sources:
		# delta=one[3]-xmed
		ixmin=one[1]-isize
		ixmax=one[1]+isize
		iymin=one[0]-isize
		iymax=one[0]+isize
		# print ixmin,ixmax,iymin,iymax
		stamp=mask[iymin:iymax,ixmin:ixmax]
		stamp=numpy.ravel(stamp)
		stamp=numpy.sort(stamp)
		ilen=len(stamp)
		itop=int(0.9*ilen)
		ibot=int(0.1*ilen)
		zmed=numpy.median(stamp)
		qqq=numpy.std(stamp[ibot:itop])
		if zmed <3. and (qqq < 10.*one[2]):
			xsources.append(one)
		
	
	print 'peaks: Found %d peaks in smoothed persistence image with flux > %f' % (len(xsources),im[i])
	return xsources


def circle(x,y,plothandle):
	'''
	plot a circle on a figure that has been displayed with imshow
	'''

	circle = patches.Circle((x, y), 10, facecolor='none',edgecolor=(1,0,0), linewidth=3, alpha=0.5)
	plothandle.add_patch(circle)
		
regionheader='''# Region file format: DS9 version 4.1
# Filename: new.fits
global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image
'''


def make_reg(sources,regionfile='sources.reg'):
	'''
	Make a simple region file from a set of
	source positions, labelling them by
	source number.

	'''

	# g=open(regionfile,'w')
	# os.chmod(regionfile,0770)
	# g.write(regionheader)

	g=per_list.open_file(regionfile)

	i=1
	for source in sources:
		g.write('circle(%6.2f,%6.2f ,5) # text={%03d}\n' % (source[1]+1,source[0]+1,i))
		i=i+1
	
	g.close()
	return



def doit(filename='foo.fits',mask_file='none',box=3,maxno=25,history_file='history.txt',local='no'):
	'''
	Main routine which handles finding of the peaks and creating plots which indicate
	how well the persistence has been subtracted.  It attemps to locate region with 
	a significant amount of persistence, but to avoid regions where the flux from
	stars in the image under consideration is large.
	
	The variables are as follows:
		filename	The persist model
		mask_file	The original science image which here is used to
				identify regions which are contaminated by stars
		box		The size of the boxcar function used to smooth the persist
				image.  
		maxno		The maximum number of peaks to find
		history_file	The place to append history information
		local		The standard way in which to indicate where output files
				should be stored.
		
	
	The routine reads file name foo, and smooths it with boxcar function with kernal of size
	box  It also reads the mask_file, the original science image. This is used in find_peaks
	to exclude regions of the image that have very high count rates, since it will be impossible
	to see persistence agains a bright source.  So that bad pixels pixels are also excluded 
	a large number is added to the mask when a pixel has bad data quality.

	The routine also creates a plot for each of the peaks that indicates where
	the peaks are.  The individual subtractions are not evaluated here, see
	subtract_eval instead.

	111011	ksl	Tried to clean this routine up and make it run a little better
	130204	ksl	Modified the data quality criterion to elimate the possibility
			that the TDF had been set which causes all of the pixels to
			be declared bad.
	130913  ksl     Restored the use of DQ flags to since the problem with loss of 
			lock, which in certain cases, particularly darks, caused all
			pixels to have bad DQ.


	'''

	history=open(history_file,'a')
	history.write('# Starting peaks on file %s using  %s as mask \n' % (filename,mask_file))

	# determine from the filename where we want
	# to put the ouputs

	work_dir=per_list.set_path(filename,'yes',local)  # This will be the Persist directory for the dataset
	fig_dir=work_dir+'/Figs/'



	# Get a name for the root of any output files
	try:
		x=filename.split('/')
		name=x[len(x)-1]
		# print 'name',name
		j=string.rindex(name,'.')
		outroot=name[0:j]
		# print 'OK got',outroot
	except ValueError:
		outroot=filename
		print 'File had no extensions, using entire name'

	# Read the persistence file
	x=per_fits.get_image(filename,1)
	if len(x)==0:
		print 'Error: peaks:  Nothing to do since no data returned for %s' % filename
		history.write('Peaks: Nothing to do since no data returned for %s\n' % filename)
		return 'NOK'

	# Now read the flt file, which will be used to create a mask file in an attempt 
	# to remove from considerations regions of the image where there is persistence which
	# would be overwhelmed by light in the current image.
	# Read both the image and the data quality
	# extension.  The data quality extension is used to help assure that we are tracking
	# persistence from light in earlier exp sures.  For these pixels, we set the value
	# in the persitence image to zero, before the persistence image is smoothed.
	# Then smooth the persitence image  so that peak finding is easier

	if mask_file!='none':
		mask=per_fits.get_image(mask_file,1,rescale='e/s')
		xmask=smooth(mask,box)
		dq=per_fits.get_image(mask_file,3,rescale='no')
		z=numpy.select([dq>0],[0],default=x)
		z=smooth(z,box)
		history.write('Peaks: Using %s to mask exposure\n' % mask_file)
	else:
		history.write('Peaks: No mask from earlier undithered exposures\n')
		z=smooth(x,box)
		mask='none'

	# After having set everything up call the routine that locates the peaks to use to see
	# how well the subtraction works

	sources=find_peaks(z,mask,maxno=maxno)

	# Write out the results in a text file containing
	# the x and y positions

	outfile=work_dir+outroot+'.peaks.dat'

	# 	print 'Error: Could not open %s' % outfile
	
	g=per_list.open_file(outfile)

	g.write('# peaks in %s\n' % filename)
	for source in sources:
		xflux=xmask[source[0],source[1]]
		g.write('%5d %5d  %6.3f %6.3f\n' % (source[1],source[0],source[2],xflux))
	g.close()

	# Write out a set of ds9 region files of the sources
	make_reg(sources,work_dir+outroot+'.peaks.reg')


	# Plot the results. This creates two images, one with the regions selectd for detailed 
	# anaylsis superposed.

	xmed=numpy.median(x)
	pylab.figure(1,(8,12))
	pylab.clf()
	pylab.subplot(211)
	pylab.imshow(x,origin='lower',cmap=pylab.cm.gray,vmin=xmed-0.05,vmax=xmed+0.05)
	plothandle=pylab.subplot(212)
	ysize,xsize=x.shape
	pylab.imshow(z,origin='lower',extent=[0,xsize,0,ysize],cmap=pylab.cm.gray,vmin=xmed-0.05,vmax=xmed+0.05)

	# For reasons I don't quit understand, one needs to add 1 to both axes
	# To get the circles to line up, as if the elements were being counted
	# from one as in IRAF.
	for source in sources:
		circle(source[1]+1,source[0]+1,plothandle)

	# Note that this really needs to be in the Figs subdirectory

	figure_name=fig_dir+outroot+'.peaks.png'
	if os.path.isfile(figure_name):
		os.remove(figure_name)
	pylab.savefig(figure_name)
	os.chmod(figure_name,0770)

	# Generation of this summary plot is now complete.


	# 110325 - This next line was generating errors with the backend I was using.  The desired behavior
	# is that there be no window created when the program is run from the command line.  
	# The error seems to occur # when you try to close the figure without first drawing it. For now 
	# I have deleted the next line.  The behaviou I am currently seeing is taht the window
	# appears as if one had issued a draw command, i.e. the window stays up, but the program
	# continues ulike show()
	# pylab.close('all')

	history.write('# Finished peaks on file %s\n' % filename)

	return  'OK'

def do_dataset(dataset,fileroot='observations',local='no'):
	'''
	This is a routine to set up the peak finding routines
	based on a dataset name and a list file.  

	The variable local indicates where the output files should be
	stored in the same sense as all of the persist programs
	
	This is the way the routine is called for pipline processing

	101214	ksl	Added so that we could run peaks from the top level directory
	110511	ksl	Modified so if there is extenal persistence associated with the dataset
			this is uesed for finding the peaks
	111011	ksl	There was an error in the modiffication above until at least now
	'''
	# Get infomation about this dataset, principally the path so that
	# we can identify where material is.

	record=per_list.read_ordered_list_one(fileroot,dataset)
	try:
		mask_file=record[0]
	except IndexError:
		print 'NOK: dataset %s does not exist in %s.ls' % (dataset,fileroot)
		return 'NOK: dataset %s does not exist in %s.ls' % (dataset,fileroot)



	work_dir=per_list.set_path(record[0],'no',local)  # This will be the Persist directory for the dataset
	fig_dir=work_dir+'/Figs'               # This will be the directory where figures are stored
	history_file=work_dir+dataset+'.txt'

	extper_file=work_dir+dataset+'_extper.fits'
	per_file=work_dir+dataset+'_persist.fits'

	# Fixed 111011
	if os.path.exists(extper_file):
		print 'Using external persistence file for dataset %s' % dataset
		string=doit(per_file,mask_file,history_file=history_file,local=local)
	elif os.path.exists(per_file):
		print 'Using internal persistence file for dataset %s' % dataset
		string=doit(per_file,mask_file,history_file=history_file,local=local)
	else:
		print 'Dataset %s does not have either an external or internal persistence file' % dataset


	print 'Finished peaks for dataset %s' % dataset
	return 'OK: peaks.do_dataset: Finished for dataset %s' %  dataset




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)==3:
		doit(sys.argv[1],sys.argv[2],local='yes')
	elif len(sys.argv)>1:
		name=sys.argv[1]
		if name.count('fit'):  # Assume its a file
			doit(sys.argv[1],'none',local='no')
		else:
			do_dataset(name,local='no')

	else:
		print 'usage: peaks.py  datasetname'

