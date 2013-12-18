#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This provides evaluation images for persistence
subtraction.  

The basic assumption is that this routine is being run from
the directory where the observations.ls file exists, and
that the appropriate files have been generated in earlier
stages of the persistence pipeline.  If the files do
not exist the routine will exit with an error


Command line usage (if any):

		usage: subtract_eval.py  dataset'

Description:  

Primary routines:

Notes:
	The evaluation images are generated here, but
	they are not put into an html file
									   
History:

101214 ksl Coding begun

'''

import sys
import os
import numpy
import glob
import pylab
import per_list
import per_fits

def read_peaks(file_xy):
	'''
	read the file containing the positions of
	the peaks in the persistence image (or
	any simple ascii file with an xy position
	in each row

	110722 - Added code to read the brightness of the predicted persistence

	'''
	try:
		f=open(file_xy,'r')
		lines=f.readlines()
		f.close()
	except IOError:
		print 'Error: read_peaks: file %s not found' % file_xy
		return []

	coords=[]
	for line in lines:
		# print line
		line=line.strip()
		word=line.split()
		if len(word)>1 and word[0][0]!='#':
			x=eval(word[0])
			y=eval(word[1])
			z=eval(word[2])
			coords.append([x,y,z])
	return coords



def do_dataset(dataset='ib6v19bzq',radius=50,local='no'):
	'''
	Examine how well one has succeeded in subtracting persistence
	from a single dataset.  Assuming that all the actual sububraction
	has been done and the peaks file is in place.

	radius here is the half-size of the box that is plotted.

	110107	Changed the name of the output figure files so that it 
		would be easier to keep track of the files that had
		been created.  Also removed some text from figure.
	110203	ksl	Added local switch so testing would be easier
	'''

	# Read information about this dataset from the observations.ls file
	# Note that the file name is hardcoded here
	record=per_list.read_ordered_list_one('observations',dataset)


	# Set up the path, and open the history file

	path=per_list.set_path(record[0],'no',local)
	fig_path=path+'/Figs/'
	history=open(path+dataset+'.txt','a')
	history.write('Start subtract_eval for dataset %s\n' % dataset)

	# Get complete names for all of the files that are to be used.
	file_flt,ext,all=per_fits.parse_fitsname(record[0])
	file_persist=path+dataset+'_persist.fits'
	file_cor=path+dataset+'_flt_cor.fits'
	file_stim=path+dataset+'_stim.fits'
	file_xy=path+dataset+'_persist.peaks.dat'

	# Check that all of these files exist

	ready=True

	if os.path.exists(file_flt) == False:
		print 'Error: subtract_eval.do_dataset: %s does not exist' % file_flt
		ready=False

	if os.path.exists(file_persist) == False:
		print 'Error: subtract_eval.do_dataset: %s does not exist' % file_persist
		ready=False

	if os.path.exists(file_cor) == False:
		print 'Error: subtract_eval.do_dataset: %s does not exist' % file_cor
		ready=False

	if os.path.exists(file_stim) == False:
		print 'Error: subtract_eval.do_dataset: %s does not exist' % file_stim
		ready=False

	if os.path.exists(file_xy) == False:
		print 'Error: subtract_eval.do_dataset: %s does not exist' % file_xy 
		ready=False
	
	if ready==False:
		return 'Error: subtract_eval.do_dataset: Some files are missing'

	# At this point we know all of the necessary files exist

	# Since we are ready we should now delete all the figures from previous
	# runs.  Note that this command is dangerous.

	png_files=glob.glob('%s/%s.peak*png' % (fig_path,dataset))
	# print 'Test',png_files
	for one in png_files:
		os.remove(one)

	# Read the xy positions from file (produced by peaks.py)

	xy=read_peaks(file_xy)

	# Read all of the images
	flt=per_fits.get_image_ext(file_flt,1)
	per=per_fits.get_image_ext(file_persist,1)
	cor=per_fits.get_image_ext(file_cor,1)
	stim=per_fits.get_image_ext(file_stim,1)


	all_orig=[] # This is a place to store histograms of the original data
	all_corr=[] # This is a place to store histograms of the corrected data


	# Set up several arrays that are used in histogram creation (as the x axes)
	# set up the stim array

	# stim_hist=[1.,3.8]
	stim_hist=[1.,4.5] # 30,000 electrons
	# x=4.
	dx=0.2
	x=stim_hist[1]+dx  # We treat everything blow 30,000 as background
	while x<=7:
		stim_hist.append(x)
		x=x+dx

	stim_hist=numpy.array(stim_hist)

	# print 'stim_hist',stim_hist
	stim_hist=10**stim_hist
	# print 'stim_hist',stim_hist

	all_sorig=[] # This is a place to store histograms of the orignal data as a function of stimulus               
	all_scorr=[] # This is a place to store histograms of the corrected data as a function of stimulus               

	# Set up the per_hist array
	per_hist=[]
	qper=0
	dper=0.02
	while qper<=0.3:
		per_hist.append(qper)
		qper=qper+dper



	source_no=0
	for one in xy:  # Main loop for each point
		source_no=source_no+1
		x=one[0]
		y=one[1]
		# print 'test x y',x,y

		# Make the stamps that are needed for each file
		xmin=one[0]-radius
		xmax=one[0]+radius
		ymin=one[1]-radius
		ymax=one[1]+radius

		ysize,xsize=numpy.shape(flt)
		if ymin<1 or xmin<1 or xmax >xsize or ymax>ysize:
			continue

		xflt=flt[ymin-1:ymax,xmin-1:xmax]
		xper=per[ymin-1:ymax,xmin-1:xmax]
	        xcor=cor[ymin-1:ymax,xmin-1:xmax]
	        xstim=stim[ymin-1:ymax,xmin-1:xmax]

		# OK at this point we have all the stamps; now flatten them

		xxflt=numpy.ravel(xflt)
		xxcor=numpy.ravel(xcor)
		xxper=numpy.ravel(xper)
		xxstim=numpy.ravel(xstim)

		med_flt=numpy.median(xxflt)
		max_per=numpy.max(xxper)
		zmin=med_flt-0.05
		zmax=med_flt+0.1

		fig_root=path+'Figs/%s.peak.%03d_%03d.' % (dataset,x,y)

		# Create figure containing the original image, the persistence image and the subtracted 
		# image surrouding a narrow region.  This is the 4 panel figure that appears in the 
		# summary html file for each ragion

		pylab.figure(11,[8,8])
		pylab.clf()
		pylab.subplot(221)
		pylab.imshow(xflt,origin='lower',cmap=pylab.cm.gray,vmin=zmin,vmax=zmax)
		pylab.title('Original')
		pylab.subplot(222)
		# pylab.imshow(xper,origin='lower',cmap=pylab.cm.gray,vmin=0.0,vmax=0.1)
		pylab.imshow(xper,origin='lower',cmap=pylab.cm.gray,vmin=-0.05,vmax=0.1)
		pylab.title('Model')
		pylab.subplot(223)
		pylab.imshow(xcor,origin='lower',cmap=pylab.cm.gray,vmin=zmin,vmax=zmax)
		pylab.title('Corrected')


		# Plot the figure that shows the observed rate as a function of estimate persistence

		pylab.subplot(224)
		pylab.plot(xxper,xxflt,'.',color='green')
		pylab.plot(xxper,xxcor,'.',color='yellow')

		

		# This constructs a histograms of median value of the original and subtracted
		# pixels as a function of the estimated persistence

		orig=[]
		corr=[]

		ii=0
		while ii<len(per_hist)-1:
			value=get_stats(xxper,xxflt,per_hist[ii],per_hist[ii+1])
			orig.append(value)
			value=get_stats(xxper,xxcor,per_hist[ii],per_hist[ii+1])
			corr.append(value)
			ii=ii+1

		# Append the results for this particular point to one for all of the points
		# This is used in the summary slide for the entire field

		# print 'test orig',orig

		all_orig.append(orig)
		all_corr.append(corr)

		# Note that per_hist has one more element than the other arrays so must allow for this
		pylab.plot(per_hist[0:len(per_hist)-1],orig,ls='steps-post',lw=4,color='red')
		pylab.plot(per_hist[0:len(per_hist)-1],corr,ls='steps-post',lw=4,color='blue')

		

		pylab.axis([0,max_per+0.01,med_flt-0.2,med_flt+0.3])
		pylab.xlabel('Est. Persistence (e/s)')
		pylab.ylabel('Flux (e/s)')

		# Finished with this histogram; now write out the figure

		figure_name='%s%d.png' % (fig_root,1)
		if os.path.isfile(figure_name):
			os.remove(figure_name)
		pylab.savefig(figure_name)
		os.chmod(figure_name,0770)


		# print 'test Finished the 4 panel figure'

		# Plot the original and subtracted pixels as a function of distance from a center positions
		# Create an array that contains the distance from the center for each pixel

		z=numpy.arange(-radius,radius+1,1)
		xx,yy=numpy.meshgrid(z,z)
		zz=xx*xx+yy*yy
		zzz=numpy.sqrt(zz)
		zzzz=numpy.ravel(zzz)

		pylab.figure(13,[6,6])
		pylab.clf()
		pylab.plot(zzzz,xxflt,'o')
		pylab.plot(zzzz,xxcor,'o')


		# This constructs a histograms of median value of the original and subtracted
		# pixels as a function of distance from the source
		# Note thatthe size here that is plotted is not determined by radius, but is
		# hardcoded to be 20 pixels. This is typical smaller than radius because 
		# we want to see how well as single star is subtracted.

		meds=[]
		med_corr=[]
		rr=[]
		r=0
		dr=3
		rmax=20
		while r<rmax:
			value=get_stats(zzzz,xxflt,r,r+dr)
			meds.append(value)
			value=get_stats(zzzz,xxcor,r,r+dr)
			med_corr.append(value)
			rr.append(r+0.5*dr)
			r=r+dr

		pylab.plot(rr,meds,ls='steps-mid',lw=3)
		pylab.plot(rr,med_corr,ls='steps-mid',lw=3)
		

		pylab.axis([0,rmax,med_flt-0.2,med_flt+0.3])
		pylab.xlabel('Radius (pixels)')
		pylab.ylabel('Flux (e/s)')



		figure_name='%s%d.png' % (fig_root,3)
		if os.path.isfile(figure_name):
			os.remove(figure_name)
		pylab.savefig(figure_name)
		# print 'test savefig ',figure_name
		os.chmod(figure_name,0770)

		# print 'OK',figure_name


		# 110622 - Elimaated to fix a problem on linux
		# pylab.close('all')

		# next section gathers information about eveything as a function of the stimulus

		i=0
		sorig=[]
		scorr=[]
		while i<len(stim_hist)-1:
			value=get_stats(xxstim,xxflt,stim_hist[i],stim_hist[i+1])
			sorig.append(value)
			value=get_stats(xxstim,xxcor,stim_hist[i],stim_hist[i+1])
			scorr.append(value)
			i=i+1
		all_sorig.append(sorig)
		all_scorr.append(scorr)

		# This ends the main loop for each data point.



	# Now make the first summary figure which is a plot of flux as a function of the model
	# stimulus

	fig_root=path+'Figs/%s.sum1' % (dataset)
	pylab.figure(14,[6,6])
	pylab.clf()

	i=0
	xmax=numpy.max(per_hist)
	ymax=(-1000)
	ymin=1000
	per_hist=numpy.array(per_hist)
	per_hist=per_hist+0.5*dper

	# print 'test len(all_corr)',len(all_corr)
	while i<len(all_corr):
		corr=numpy.array(all_corr[i])
		orig=numpy.array(all_orig[i])

		# print 'test orig',orig

		corr=corr-corr[0]
		orig=orig-orig[0]
		k=0
		while k<len(orig):
			if orig[k]<-900:
				break
			k=k+1
		k=k-1
		if k>0:

			pylab.plot(per_hist[0:k],orig[0:k],'ro-',lw=2)
			pylab.plot(per_hist[0:k],corr[0:k],'bo-',lw=2)
			zmin=numpy.min(corr[0:k])
			if zmin<ymin:
				ymin=zmin
			zmax=numpy.max(orig[0:k])
			if zmax>ymax:
				ymax=zmax
		else:
			print 'Error: subtract_eval.do_dataset: there is a problem, because k=0'
		i=i+1
	pylab.axis([0,xmax+0.05,ymin-0.05,ymax+0.05])
	pylab.xlabel('Est. Persistence (e/s)')
	pylab.ylabel('Flux (e/s)')

	figure_name=fig_root+'.png'
	if os.path.isfile(figure_name):
		os.remove(figure_name)
	pylab.savefig(figure_name)
	os.chmod(figure_name,0770)
	# print 'OK',figure_name,pylab.axis()


	# 110622 - Eliminated to fix a problem on linux
	# pylab.close(14)


	# Now make the second summary figure showing the flux as a function of the stimulus.  
	# Note that because these occurred at different times you will get a range here at
	# any stimulus


	# Construct the x axis for this figure
	i=0
	xstim_hist=[]
	while i<len(stim_hist)-1:
		xstim_hist.append(0.5*(stim_hist[i]+stim_hist[i+1]))
		i=i+1

	fig_root=path+'Figs/%s.sum2' % (dataset)
	pylab.figure(15,[6,6])
	pylab.clf()


	# Now go through each row in the arrays and plot the results
	i=0
	while i<len(all_corr):
		corr=numpy.array(all_scorr[i])
		orig=numpy.array(all_sorig[i])
		corr=corr-corr[0]
		orig=orig-orig[0]
		# print 'orig',orig
		# print 'corr',corr
		k=0
		while k<len(orig):
			if orig[k]<-900:
				break
			k=k+1
		k=k-1
		if k>0:
			pylab.semilogx(xstim_hist[0:k],orig[0:k],'ro-',lw=2)
			pylab.semilogx(xstim_hist[0:k],corr[0:k],'bo-',lw=2)
			zmin=numpy.min(corr[0:k])
			if zmin<ymin:
				ymin=zmin
			zmax=numpy.max(orig[0:k])
			if zmax>ymax:
				ymax=zmax
		i=i+1

	pylab.xlabel('Stimulus (e)')
	pylab.ylabel('Flux (e/s)')
	pylab.axis([3e4,1e7,-0.1,0.3])

	figure_name=fig_root+'.png'
	if os.path.isfile(figure_name):
		os.remove(figure_name)
	pylab.savefig(figure_name)
	os.chmod(figure_name,0770)
	# print 'OK',figure_name,pylab.axis()


	# 110622 - Eliminated to fix a problem on linux
	# pylab.close(15)

	


	


	history.write('End subtract_eval for dataset %s\n' % dataset)
	history.close()



def get_stats(x,y,xmin,xmax):
	'''
	Get statistics in the array y corresponding to elements
	in x between xmin and xmax

	Nntes:


	This is intended as a general tool which takes two arrays
	of the same shape x and y.  
	
	It then selects elements in x between xmin and xmax, and 
	then computes the median, and possibly other statistics of 
	the same array elements in y, 

	If there are no elements within the desired interval it
	returns -999

	110327	ksl	Coded to improve the plots produced by do_dataset

	'''

	# print 'test get_stats',x.shape,xmin,xmax
	# print 'test get_stats x',x
	# print 'test get_stats y',y
	z=numpy.ma.masked_outside(x,xmin,xmax)
	zmask=numpy.ma.getmask(z)
	yy=numpy.ma.array(y,mask=zmask)
	yy=yy.compressed()
	if len(yy)>0:
		return numpy.median(yy)
	else: 
		return -999



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# do_dataset(int(sys.argv[1]))
		do_dataset(sys.argv[1])
	else:
		print 'usage: .subtract_eval.py  dataset'

