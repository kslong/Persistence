#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This program is intended to help with the evaluation of persistance in the images
from WFC3 IR


Description:  

	At present the program simply looks for all of the flt.fits imaage in
	any subdirctory and loads up ds9 with these imaes in time order.

Primary routines:

Notes:
									   
History:

091026 ksl	Coding begun
091215 ksl	Updated to prevent reading UVIS data and put in an option to only
		look at full arrays

'''

import sys
import os
import numpy
import pyraf
import time
import string
import math
import scipy

# This section was imported from ghost.py
import pyfits
from pyraf import iraf
import pylab

#ksl_fits

import ksl_fits
from per_list import *
from ds9 import *
from per_fits import *


# Utilities

def get_scale(filename):
	'''
	Determine whether the image has been exposure corrected and if so
	return the exposure time.  If not return 1.
	'''
	xname=parse_fitsname(filename)

	xxxx=iraf.hselect(xname[2],'exptime,unitcorr','yes',Stdout=1)
	xxxx=xxxx[0].split('\t')
	# print xxxx
	if xxxx[1]=='COMPLETE':
		return eval(xxxx[0])
	else:
		return 1

def get_exposure(filename):
	'''
	Return the expousre time and information on whether the image is
	already exposure time corrected.

	The routine always returns the exposure time.  It also returns 1 if the
	image has been exposure corrected, 0 if not
	'''
	xname=parse_fitsname(filename)

	xxxx=iraf.hselect(xname[2],'exptime,unitcorr','yes',Stdout=1)
	xxxx=xxxx[0].split('\t')

	# print xxxx
	if xxxx[1]=='COMPLETE':
		return eval(xxxx[0]),1
	else:
		return eval(xxxx[0]),0

def get_times(filename):
	'''
	Return the beginning and end times for an image
	'''
	xname=parse_fitsname(filename)

	xxxx=iraf.hselect(xname[2],'expstart,expend','yes',Stdout=1)
	xxxx=xxxx[0].split('\t')

	# print xxxx
	return xxxx



def get_statistics(filename,cuts=[0,50000,100000,300000,1000000]):
	'''
	Get information about the total electron distribution in an image
	'''
	expcorr=get_scale(filename)
	xmin=0
	xmax=1e7/expcorr
	x=[]
	y=[]

	xname=parse_fitsname(filename)
	zz=iraf.imhistogram(xname[2],z1=xmin,z2=xmax,nbins=1000,listout='yes',Stdout=1)
	for word in zz:
		q=word.split()
		x.append(eval(q[0])*expcorr)
		y.append(eval(q[1]))
	
	i=0
	cuts=numpy.array(cuts)
	sums=numpy.zeros(len(cuts))


	while i<len(x):
		j=0
		while j<len(cuts):
			if x[i]>cuts[j]:
				sums[j]=sums[j]+y[i]
			j=j+1
		i=i+1
	print sums
	return sums

def eval_positions(records):
	'''
	Evaluate how positions are changing in a set of observations
	'''

	print 'eval_positions: check0'
	positions=[]
	for record in records:
		x=pyraf.iraf.hselect('%s' % record[0],'$I,expend,ra_targ,dec_targ,postarg1,postarg2,ra_aper,dec_aper,,filter,exptime,targname','yes',Stdout=1)
		x=x[0].replace('\t',' ')
		x=x.split()
		# print x
		positions.append(x)
	
	ra_start=0.0
	dec_start=0.0
	t_start=0
	for position in positions:
		if t_start==0:
			ra_start=eval(position[6])
			dec_start=eval(position[7])
			t_start=eval(position[1])
		ra=eval(position[6])
		dec=eval(position[7])
		t=eval(position[1])

		delta_ra=ra-ra_start
		delta_dec=dec-dec_start
		delta_t=(t-t_start)*24.

		if math.fabs(delta_ra)>0.1 or math.fabs(delta_dec)>0.1:
				print '# Major Move of %10.6f %10.6f' % (delta_ra,delta_dec)
				ra_start=ra
				dec_start=dec
				delta_ra=0
				delta_dec=0

		print '%70s %8.2f %10.6f %10.6f' % (position[0],delta_t,delta_ra,delta_dec)

def get_positions(dataset='all',interval=[-6,6],fileroot='observations'):
	'''
	Get position information from a set of records, and calculate the offsets

	Note that this would normally be for an interval of time
	'''
	if dataset=='all':
		dataset='first'
		interval=[0,1e10]
	

	print 'Check1'
	records=read_ordered_list2(fileroot='observations',dataset='first',interval=interval,outroot='none')
	print 'Check2'

	eval_positions(records)

def get_prog_positions(prog_id=11649,fileroot='observations'):
	'''
	Get position information for a program
	'''
	
	records=prog_search(prog_id)

	if len(records)==0:
		print 'No records found for prog_id %d' % prog_id
		return

	eval_positions(records)


def make_mask(filename,xmin,electrons='yes'):
	'''
	Read and make a mask image for portions
	of the image that could have afterglow

	100302	- Modified to allow one to construct masks without rescaling in those cases where
		the image had already been converted to elec/s
	'''

	print 'make_mask :',filename

	scale=1
	if electrons=='yes':
		scale=get_scale(filename)
		print 'Scale factor ',scale
	else:
		print 'Not rescaling to electons'

	xname=parse_fitsname(filename)

	f=pyfits.open(xname[0])
	data=f[xname[1]].data
	mask=numpy.select([data>xmin/scale],[1],default=0)
	f.close()
	return mask

def get_image(filename):
	'''
	read a science image using pyfits.  Note that the name has to
	include the extension
	'''

	xname=parse_fitsname(filename)

	f=pyfits.open(xname[0])
	data=f[xname[1]].data

	f.close
	return data

def get_dq(filename):
	'''
	read the dataquality.  Note that here we ignore the extension
	in the filename if it exists and assume dq is the 3rd extension

	'''
	xname=parse_fitsname(filename)

	print xname

	f=pyfits.open(xname[0])

	data=f[3].data

	f.close
	return data

# End of utilities
# Beginning of a section intended to measure the decay of the afterglow
# in broad bins


def find_dithers(records):
	'''
	Locate the place where the pointing changed
	'''
	dither=[]
	ndither=0
	ra_old=records[0][12]
	de_old=records[0][13]
	print 'Start RA and DEC',ra_old,de_old
	for record in records:
		ra=record[12]
		de=record[13]
		if de_old!=de or ra_old != ra:
			ndither=ndither+1
		dither.append(ndither)
	print 'There were %d movements in this sequence ' % ndither
	print dither
	return dither




def time_history(dataset='last',delta=3., fileroot='observations',xmin=100000,xmax=200000):
	'''
	Try to find the time history of the afterflow

	Functionall THIS routine works.  Practically one needs to figure out where the back
	ground regions are in each image.  It is possible one could do this by carrying out
	DAOPHOT on each of the images and blotting the source regions.  The size of the blots
	would likely need to depend on the strength of the point source.  
	'''

	# Select the records that run up to a particular dataset
	records=read_ordered_list(fileroot,dataset,delta)


	# Find the places there were dithers
	dithers=find_dithers(records)

	i=1
	while i < len(records):
		# Get the start time and the image in which persistence will be measured
		xtime=eval(records[i][6])
		science=get_image(records[i][0])

		# Make the image x and get the associated times needed for estimating the persistence
		x,t=make_badpix2(records[0:i],outroot='none')

		# print 'chomp',len(numpy.unique(x)),len(numpy.unique(t))


		t=86400.*(xtime-t) # Convert to to delta time to this observation in hours
		dtimes=numpy.unique(t) # Get the uniue values of time in the time array.
		print 'Delta times',dtimes

	


		# Create a mask for x that selects only pixles between xmin and xmax
		mask1=numpy.select([x >= xmin],[0],default=1)
		mask2=numpy.select([x <= xmax],[0],default=1)
		mask=mask1+mask2

		# print 'Mask',numpy.unique(mask)

		for dtime in dtimes:
			tmask=numpy.select([t==dtime],[0],default=1) # Select the pixels with dtime
			xmask=mask+tmask

			# print 'xMask',numpy.unique(xmask)

			# Make the masked science array and find its median
			# print len(science),len(xmask)

			xscience=numpy.ma.array(science,mask=xmask)
			xxscience=xscience.compressed()
			xxxscience=numpy.ravel(xxscience)
			# print len(xscience),len(xxscience),len(xxxscience)

			med_persist=numpy.median(xxxscience)
			med_science=numpy.median(science)

			# print len(xxxscience),med_persist,med_science

			dpersist=med_persist-med_science

			print 'Time %6.2f  Persist %7.3f Pixels %d Background %7.3f' % (dtime,dpersist,len(xxxscience),med_science)
		i=i+1

			








# End of a section to measure the persistence as a function of time

def plot_afterglow(first,second,px_max=200000,py_max=1,new='yes',semi='no'):
	'''
	plot the second image as a function of the first where 
		px_max is the maximum value on the x-axis to plot
		py_max is the maximum value on the y axis to plot
	
	The plot has the x-axis in electrons and y axis in e/s

	100303 Reflected a small change to plot the afterglow in electrons/s
	'''

	# Get rid of the zeros
	mask=numpy.select([first>0],[0],default=1)
	x=numpy.ma.array(first,mask=mask)
	xx=x.compressed()
	y=numpy.ma.array(second,mask=mask)
	yy=y.compressed()

	print 'Lengths',len(first),len(x),len(xx)
	print 'Lengths',len(second),len(y),len(yy)

	# OK now can deal with all of this

	pylab.figure(1,(6,6))
	if new=='yes':
		pylab.clf()

	if semi != 'yes':
		pylab.plot(xx,yy,'.')
		pylab.axis([0,px_max,0,py_max])
		pylab.xlabel('First image (electrons)')
		pylab.ylabel('Afterglow image (elec/s)')
	else:
		pylab.semilogx(xx,yy,'.')
		pylab.axis([0,px_max,0,py_max])
		pylab.xlabel('First image (electrons)')
		pylab.ylabel('Afterglow image (elec/s)')
	pylab.draw()
	return

def plot_ratio(first,second,px_max=200000):
	'''
	plot the ratio of values in two images, where

	px_max is the maximum of the x-axis.  

	100309 Note that this was originally called by do_ghost, but has been
	commeneted out.
	'''


	# Get rid of the zeros
	mask=numpy.select([first>0],[0],default=1)
	x=numpy.ma.array(first,mask=mask)
	xx=x.compressed()
	y=numpy.ma.array(second,mask=mask)
	yy=y.compressed()

	print 'Lengths',len(first),len(x),len(xx)
	print 'Lengths',len(second),len(y),len(yy)

	# OK now can deal with all of this more quickly

	ratio=yy/xx
	xratio=numpy.select([ratio>0],[ratio],default=1e-6)

	# ratio=second/first
	# xratio=numpy.select([ratio>0],[ratio],default=1e-6)

	pylab.figure(2,(6,6))
	pylab.clf()
	# pylab.plot(first,ratio,'.')
	# pylab.axis([0,px_max,0,0.01])
	pylab.semilogy(xx,xratio,'.')
	pylab.axis([0,px_max,1e-5,1])
	pylab.xlabel('First image (electrons)')
	pylab.ylabel('Ratio')
	pylab.draw()
	return



def do_ghost(bright='./11927/Visit10/ibce10adq_flt.fits',afterglow='./11927/Visit10/ibce10aeq_flt.fits',xmin=10000.,xxmax=200000,yymax=0.5):
	'''
	The main routine to plot the afterglow from a bright image in a single afterglow image.  
	
	Here
		xmin defines the lower limit to the mask in electrons
		xxmax defines the maximum value to plot on the xaxis
		yymax defines the maximum value to plot on the yaxis
	'''

	# Get a mask image to define the regions of interest

	mask=make_mask(bright,xmin)

	# Get the bright image again

	first_image=get_image(bright)
	first_scale=get_scale(bright)

	second_image=get_image(afterglow)
	second_scale=get_scale(afterglow)

	first_image=first_image*mask*first_scale
	second_image=second_image*mask
	exp,choice=get_exposure(afterglow)
	if choice==0:
		second=second/exp


	first_image=numpy.ravel(first_image)
	second_image=numpy.ravel(second_image)

	# Now plot it all

	# Note that the next line is hardwired.  One should really make it
	# possible to choose whether one wants linear or semilog in
	# the call to this routine

	# plot_afterglow(first_image,second_image,xxmax,yymax)
	plot_afterglow(first_image,second_image,xxmax,yymax,semi='yes')

	# plot_ratio(first_image,second_image)

def do_ghosts(bright='first',dt=1.,xmin=10000,xxmax=200000,yymax=1):
	'''
	Plot the persistance in multiple files.   At present, written essentially assuming
	that the afterglow files are darks.

	Here:
		dt	The time in hours to go forward from the bright image
		xmin	The minimum brightness in the bright image to consider, in electrons
		xxmax   The maximum brightness in the bright image to plot
		yymax   The upper limit in the plot in elec/s

	100303 - Modified so what is plot is First image electrons, vs electrons per second in
		the various output images

	'''
	records=read_ordered_list2(fileroot='observations',dataset=bright,interval=[0,dt])

	if len(records) == 1:
		print 'do_ghosts: there were no observations within %.1f of %s' % (dt,bright)
		return


	bright=records[0][0] # This allows one to have entered a data set name

	mask=make_mask(bright,xmin)
	bright_image=get_image(bright)
	bright_scale=get_scale(bright)

	# We want the bright image in electrons since the relation to full well matters
	bright_image=mask*bright_scale*bright_image
	bright_image=numpy.ravel(bright_image)

	pylab.figure(1)
	pylab.clf()

	i=1
	while i < len(records):
		record=records[i]
		filename=record[0]
		print 'file --- ', filename
		image=get_image(filename)
		# We want the output images in electrons per second
		image=mask*image
		exp,choice=get_exposure(filename)
		if choice==0:
			image=image/exp

		# image=mask*image*get_scale(filename)
		image=numpy.ravel(image)
		print len(bright_image),len(image)
		plot_afterglow(bright_image,image,xxmax,yymax,new='no')
		i=i+1
	pylab.figure(1)
	# pylab.savefig('ghosts.png')
	return

def test_dark(dark='ibcu0yvqq',dt=12,mask_ratio=2.):
	'''
	Attempt to use a dark to define the ghost pixels and then to find which image caused the problem

	Here dt defines the time to look back in time, and
	the value mask ratio is used to flag afterglow pixels.  Pixels values above mask_ratio*median value
	in the image are said to be persistence pixels

	100302
	'''
	records=read_ordered_list2(fileroot='observations',dataset=dark,interval=[-dt,0])
	dark=records[len(records)-1][0]  # this allows one to have entered the dataset name; otherwise would need to enter exact filename

	exposure=get_exposure(dark)
	dark_image=get_image(dark)
	dq=get_dq(dark)

	print 'Exposure info ',exposure
	if exposure[1]==0 and exposure[0]>0:
		dark_image=dark_image/exposure[0]  # Convert to elec/s

	med=numpy.median(dark_image)
	xmin=med*mask_ratio

	print 'Median',med,xmin

	mask=numpy.select([dark_image>xmin],[1],default=0)  # Create a mask of only those pixels greater than xmin
	xmask=numpy.select([dq==0],[1],default=0)  # Create a mask of the good pixels
	mask=mask*xmask

	xtimes=numpy.ones_like(dark_image) # An array filled with zeros that has the same shape as the dark image

	i=0
	while i<len(records):
		record=records[i]
		filename=record[0]
		print 'file --- ', filename
		zz=get_times(filename)
		print 'times ',zz
		endtime=eval(zz[1])
		image=get_image(filename)
		image=image*get_scale(filename) # This converts the orginal image to electrons
		if i==0:
			ximage=numpy.copy(image)
			xtimes=xtimes*endtime  # Initiall all of the times correspond the end time of the first observation
			
		else:
			xtimes=numpy.select([image>=ximage,image<ximage],[endtime,xtimes],0)
			ximage=numpy.maximum(ximage,image)
			
		i=i+1
	

	# Now just worry about the bright regions
	xtimes=endtime-xtimes

	# maked masked images of everything
	xtimes=xtimes*mask
	ximage=ximage*mask
	dark_image=dark_image*mask

	ksl_fits.put_fits('times.fits',xtimes)
	ksl_fits.put_fits('bright.fits',ximage)

	test_plot(dark_image,ximage,xtimes)

	# Find the unique values of the times
	ztimes=numpy.unique(xtimes)
	print 'times',ztimes

	i=1
	# For a masked array, the passed values are for mask values of 0
	while i<len(ztimes):
		zmask=numpy.select([xtimes==ztimes[i]],[1],default=0)
		fmask=numpy.select([ximage>1e5],[1],default=0)
		zmask=mask*zmask*fmask # Mask which is a bright pixel in the right place
		zzmask=zmask-1 # Mask in which the good pixels are 0

		zbright=numpy.ma.array(ximage,mask=zzmask)
		zzbright=zbright.compressed()

		zdark=numpy.ma.array(dark_image,mask=zzmask)
		zzdark=zdark.compressed()

		print 'What  %6.3f %6d  %8.1f %8.3f' % (ztimes[i],len(zzbright),numpy.median(zzbright),numpy.median(zzdark))
		i=i+1




	# print records

def test_plot(dark,bright,time):
	'''

	A routine called by test_dark

	'''

	xdark=numpy.ravel(dark)
	xbright=numpy.ravel(bright)
	xtime=numpy.ravel(time)  

	# Get rid of the zeros
	mask=numpy.select([xbright>1e6],[0],default=1)
	xxdark=numpy.ma.array(xdark,mask=mask)
	xxxdark=xxdark.compressed()

	xxbright=numpy.ma.array(xbright,mask=mask)
	xxxbright=xxbright.compressed()

	xxtime=numpy.ma.array(xtime,mask=mask)
	xxxtime=xxtime.compressed()


	print 'Lengths',len(xdark),len(xxdark),len(xxxdark)
	print 'Lengths',len(xbright),len(xxbright),len(xxxbright)
	print 'Median', numpy.median(xbright),numpy.median(xxxbright)
	print 'Lengths',len(xtime),len(xxtime),len(xxxtime)

	# OK now can deal with all of this

	pylab.figure(1,(6,6))
	pylab.clf()
	pylab.plot(xxxtime,xxxdark,'.')
	pylab.axis([0,0.5,0,1])
	pylab.xlabel('Time')
	pylab.ylabel('Afterglow')
	pylab.draw()
	return

	



def subtract_persist(science_image,bright_image,norm=0.001,e_fermi=75000,kt=20000):
	'''
	Subtract persistence from a science image using a bright_image as a template
	and assuming a fermi-like subtraction function

	The results are displayed in ds9, and plots are generated.

	100307  - A first attempt to subtract persistence using from a single bright image.
	One can play with the values.  At persent this is a stand alone routine and so
	one must read in the images from the command line.

	Note that there is amore developed attempt called subtract_persist2

	'''

	bright=get_image(bright_image)
	bright_scale=get_scale(bright_image)
	bright=bright*bright_scale    # put bright image into electrons

	persist=calc_fermi(bright,norm,e_fermi,kt)

	science=get_image(science_image)
	subtracted=science-persist

	ksl_fits.put_fits('xbright.fits',bright)
	ksl_fits.put_fits('xscience.fits',science)
	ksl_fits.put_fits('xpersist.fits',persist)
	ksl_fits.put_fits('xsubtracted.fits',subtracted)


	start_ds9()

	print 'Display %s' % 'xbright.fits'
	command='xpaset -p ds9 file new %s' % 'xbright.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')

	print 'Display %s' % 'xpersist.fits'
	command='xpaset -p ds9 file new %s' % 'xpersist.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')


	print 'Display %s' % 'xscience.fits'
	command='xpaset -p ds9 file new %s' % 'xscience.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')


	print 'Display %s' % 'xsubtracted.fits'
	command='xpaset -p ds9 file new %s' % 'xsubtracted.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')

	do_ghost(bright_image,science_image,xxmax=1e6)

	xx=numpy.linspace(0,2e6,1000)
	plot_fermi(xx,norm,e_fermi,kt)
	pylab.savefig('persist.png')







def calc_fermi(x,norm=1,e_fermi=75000,kt=20000):
	'''
	calculate the predicted afterglow assuming a Fermi-Dirac like
	distribution and series of x values

	The Fermi distribution is defined as 

	1/(exp((E_f-E)/kT)+1) where E_f is the Fermi energy.  For E>E_f the exponential term has little  
	effect.  kT determines how sharpe the cutoff is.  For E=E_f we are at the distribution value of
	the function is 0.5 regardless of kT

	100305 - Note that this operates on any dimension array
	
	'''

	y=e_fermi-x
	y=y/kt
	y=numpy.exp(y)
	y=1./(y+1.)
	y=norm*y
	return y

def calc_afermi(x):
	'''
	Calculate the distribuition Adam feels describes the distribution.  It is not a fermi distribution
	but has similar characteristics.
	'''

	y=numpy.power(x/80000.,3)
	y=numpy.exp(-y)
	y=4.2/(y+1)
	y=y*numpy.power(x/100000.,0.1)
	y=y-1.9
	y=numpy.clip(y,0,1e5)
	return y

def calc_atime(t):
	'''
	Calculate the temporal distribution Adam feels describes the fall-off in counts
	'''
	y=numpy.exp(-t/275.) - 9e-6*t+ 0.174
	return y


def plot_fermi(x,norm,e_fermi,kt):
	'''
	Plot a fermi distribution given an array or list of x values.  Usually used to overplot a
	distribution elsewehre.
	'''

	y=calc_fermi(x,norm,e_fermi,kt)
	pylab.figure(1)
	z=pylab.axis()
	pylab.plot(x,y,'-',lw=3)
	pylab.axis(z)

# This is especially for DARKS or a specific program that you want to check for observations
# that might have produced problems.  

def prog_search(prog_id=11929):
	'''
	Find all the records in a given program id

	Note that observations.ls must exist.  This is clearly only currenly
	set up for the _flt exposures.
	'''
	records=read_ordered_list0(fileroot='observations')

	prog_records=[]
	#now locate all the records for this program

	for record in records:
		test=eval(record[2])
		if test==prog_id:
			prog_records.append(record)
	# print prog_records
	return prog_records

def prog_retrieve(prog_id=11929):
	'''
	Write a file which will copy all of the files 
	in the observations.ls file from a certain prog_id

	Notes -- the observtions.ls file is hardwired

	100324 - Written so I could look at individual programs
	for persistence.  
	'''

	
	records=prog_search(prog_id)
	if len(records)==0:
		print 'No files from prog_id %d found' % prog_id
		return

	g=open('Get_%d' % prog_id,'w')
	g.write('mkdir x_%d\n' % prog_id)

	
	for record in records:
		xname=record[0]
		xxname=parse_fitsname(xname)
		g.write('cp %s x_%d\n' % (xxname[0],prog_id))
	g.close()


def check_earlier(prog_id=11929,delta_time=1,save_info='yes'):
	'''
	Look for exposures that are earlier than set of records and summarize the results to a file

	where delta_time is in hours

	The information is saved to prog_id.check_earlier.txt, unless save_info is not 'yes'. In this
	case the information is saved to a temporary file.

	Note that this rooutine assumes observations.ls exists

	'''
	prog_records=prog_search(prog_id)

	if save_info=='yes':
		g=open('prog%d.check_earlier.txt' % prog_id,'w')
	else:
		g=open('tmp.check_earlier.txt','w')

	g.write('# Checking for observations prior to %.1f hour prior to a each dataset in a program %s\n' % (delta_time,prog_id))
	g.write('# Note that earlier exposures from the same programs are not recorded\n')

	for prog_record in prog_records:
		isave=0
		print '\nChecking ',prog_record
		records=read_ordered_list(dataset=prog_record[1],delta_time=delta_time)
		if len(records)==0:
			print 'There was an error.  At least one record should be returned'
		elif len(records)==1:
			print 'There are no earlier exposures from other programs'
		else:
			i=0
			while i<len(records)-1:
				record=records[i]
				if record[2]!=prog_record[2]:
					tlast=eval(prog_record[6])
					tnow=eval(record[6])
					dt=(tlast-tnow)*24.
					print 'Found    ',dt,record
					if isave==0: # We want to write this particular record to the file
						g.write('! Dataset   %10s in %s line %s on %s %s with exposure %7.1f of %15s has precursors\n\n' % (prog_record[1],prog_record[2],prog_record[3],prog_record[7],prog_record[8],eval(prog_record[11]),prog_record[12]))
						isave=1
					g.write('  dt %5.1f  %10s in %s line %s on %s %s with exposure %7.1f of %15s\n' % (dt,record[1],record[2],record[3],record[7],record[8],eval(record[11]),record[12]))
					cuts=[0,50000,100000,300000,1000000]
					zzz=get_statistics(record[0],cuts)
					newstring='      Npix'
					k=0
					while k<len(cuts):
						newstring=newstring+' >%-7d = %-7d' % (cuts[k],zzz[k])
						k=k+1
					print newstring
					g.write('%s\n\n' % newstring)


				i=i+1
	g.close()



# This is where the old per_eval begin
def make_jpegs(records):
	'''
	Make the jpegs using ds9

	100128 - Note that this does not really make the jpegs.  Instead
	it loads the images.  It is not entirely obvious what jpeg one wants.
	'''

	print 'starting to load ds9'

	# First check if ds9 is already extant'
	iok=os.system('xpaget ds9 version &> tmp.txt')
	if iok:
		os.system('ds9 &')
		time.sleep(5)  # give system a chance to instandiate ds9

	else:
		os.system('xpaset -p ds9 frame delete all')
	
	os.system('xpaset -p ds9 scale histequ')
	os.system('xpaset -p ds9 width 600')
	os.system('xpaset -p ds9 height 600')

	# open a text file to record results 
	g=open('per_eval.txt','w')
	# lastime=float(records[len(records)-1][8])
	lastime=float(records[len(records)-1][6])
	for record in records:
		filename=record[0]
		print 'Display %s' % filename
		command='xpaset -p ds9 file new %s' % filename
		# os.system('xpaset -p ds9 scale histequ')
		os.system(command)
		os.system('xpaset -p ds9 zoom to fit')
		os.system('xpaset -p ds9 scale limits -5 200')
		# dt=float(record[8])-lastime
		dt=float(record[6])-lastime
		text='%s  %s %10s  dt: %.2f' % (record[2],record[3], record[12],dt)
		color='green'
		if dt==0.0:
			color='yellow'
		write_text_regionfile(text,650,50,color)
		os.system('xpaset -p ds9 regions load ds9.reg')
		g.write('%8.2f %10s %10s %10s %10s %30s\n' % (dt,record[2],record[3],record[12],record[1],record[0]))
	g.close()

def make_badpix(records,full=100000,outroot='persist'):
	'''

	The routine creates a fits file that is zero for all
	of the pixels that never exceeded full well as defined
	by 'full" and is the maximum value of the various
	images otherwise.  It therefore is a very simple
	bad pixel mask.

	Note that this makes a map from all of the records including
	the last one.  The assumption made here is that you will
	only present the records you want to construct the badpix
	image from, so for example one often sets records to

	records[0:len(records]-1] 
	
	which would not pass the last record to this routine

	Note that this uses iraf and not pyfits

	100128 - Modified so it basically records all the pixel values for
	situations were the pixel value (allowing for absorption) has gone
	above full.
	'''

	# First delete any old temporary files
	os.system('rm tmp*.fits')

	# Now check that records is not empty
	if len(records)==0:
		print 'Error: make_badpix had nothing to process, since len(records) was zero'
		return

	print 'Starting to make bad pixel file'
	g=open(outroot+'.txt','w')


	# Now process the first file
	filename=records[0][0]
	exptime=float(records[0][11])
	print 'Exposure ', exptime
	# pyraf.iraf.imcopy(filename,"tmp1.fits")
	xxxx=pyraf.iraf.hselect(filename,'unitcorr','yes',Stdout=1)
	# print xxxx[0]
	if xxxx[0]=='COMPLETE':
		expression="(a>%f) ? a*%f:0 "  % (full/exptime,exptime)
	else:
		expression="(a>%f) ? a:0 "  % (full)

	# print 'Ok',exptime, expression
	zzz=pyraf.iraf.imstatistics(filename,fields = "npix,mean,stddev,min,max,image",lower=full/exptime,Stdout=1)
	string='Exposure %.1f  statistics %s' % (exptime,zzz[1])
	print '%s' % string
	g.write('%s\n' % string)

	foo=pyraf.iraf.imexp(expr = expression,a=filename,output="tmp1.fits",Stdout=1)
	# Now process each of the remaininng files in order
	i=1
	while i < len(records):
		filename=records[i][0]
		xxxx=pyraf.iraf.hselect(filename,'unitcorr','yes',Stdout=1)
		exptime=float(records[i][11])
		# print xxxx[0]
		if xxxx[0]=='COMPLETE':
			expression="(a>%f) ? a*%f:0 "  % (full/exptime,exptime)
		else:
			expression="(a>%f) ? a:0 "  % (full)

		# print 'Ok',exptime, expression
		zzz=pyraf.iraf.imstatistics(filename,fields = "npix,mean,stddev,min,max,image",lower=full/exptime,Stdout=1)
		string='Exposure %.1f  statistics %s' % (exptime,zzz[1])
		print '%s' % string
		g.write('%s\n' % string)

		foo=pyraf.iraf.imexp(expr = expression,a=filename,output="tmp2.fits",Stdout=1)
		foo=pyraf.iraf.imexp(expr = "max(a,b)",a="tmp1.fits",b="tmp2.fits",output="tmp3.fits",Stdout=1)
		os.system('rm tmp2.fits')
		os.system('mv tmp3.fits tmp1.fits')
		i=i+1
	os.system('mv tmp1.fits %s.fits' % outroot)
	g.close()

def make_badpix2(records,outroot='persist'):
	'''
	A new attempt to make a persistence map with times

	The routine returns arrays containing the highest pixel values, x, and the time
	of same.  The times should be in MJD

	100309 - This versions uses pyfits and numpy instead of iraf
	100427 - Added switch to prevent output fits files being written if outroot=='none'
	'''

	if len(records)==0:
		print 'Error: make_badpix2: No records to construct persistence from'
		return [[]],[[]]

	i=0
	for record in records:
		if i==0:
			x=get_image(record[0])
			scale=get_scale(record[0]) 
			if scale != 1:
				x=x*scale
			time=eval(record[6])
			t=numpy.ones_like(x)
			t=t*time
			# print 'Start', scale, numpy.median(x),numpy.median(t)

			i=i+1
		else:
			xx=get_image(record[0])
			scale=get_scale(record[0]) 
			if scale != 1:
				xx=xx*scale
			time=eval(record[6])
			# print time
			t=numpy.select([xx>x],[time],default=t)
			x=numpy.select([xx>x],[xx],default=x)
			# print 'Now', scale, numpy.median(x),numpy.median(t)
	
	if outroot != 'none':
		ksl_fits.put_fits(outroot+'_cts.fits',x)
		ksl_fits.put_fits(outroot+'_time.fits',t)
	return x,t

def get_bounded_median(x,y,xmin=100000,xmax=300000):
	'''
	Get the median of y for pixels which have values between xmin and
	xmax in x
	'''
	mask=numpy.select([xmin<=x],[0],default=1) # pix first time
	tmask=numpy.select([x <= xmax],[0],default=1) # pix first time
	zmask=tmask+mask

	yy=numpy.ma.array(y,mask=zmask)
	yyy=yy.compressed()
	print 'The bounded median is %7.3f with %d pixels' % (numpy.median(yyy),len(yyy))
	return numpy.median(yyy)



def subtract_persist2(dataset='last',delta=3., fileroot='observations',fermi_n=0.1,fermi_e=100000,fermi_kt=20000,t_const=1.,subtract_back='yes',xplot_max=2e5):
	'''
	Bravely attempt to create an image to subtract

	100516 - Added choice regarding whether one subtracts a background from the counts which are plotted to handle the case where we have
		used a tungsten lamp which illuminates most of th field
	'''

	records=read_ordered_list(fileroot,dataset,delta)

	x,t=make_badpix2(records[0:len(records)-1],'persist2')


	record=records[len(records)-1]
	xtime=eval(record[6])
	print xtime,t
	t=24.*(xtime-t)  # Now t is a delta time in hours

	# So at this point time is the time in hours to the dataset we want to correct.
	dtimes=numpy.unique(t)
	print 'delta_times',dtimes           
	tt=numpy.exp(-t/t_const)
	# This is the image which shows the time deecay function
	ksl_fits.put_fits('persist_time2.fits',tt)

	# This is the image which multiples the brightest pixel value by the decay time 
	xx=x*tt
	ksl_fits.put_fits('persist_cts2.fits',xx)

	xx=calc_fermi(xx,fermi_n,fermi_e,fermi_kt)

	# This is the final persistence image
	ksl_fits.put_fits('persist_cts3.fits',xx)

	science=get_image(record[0])
	# Get some statistics of the science image
	med_science=numpy.median(science)

	delta=science-xx

	# Write out the subtracted image
	ksl_fits.put_fits('subtracted.fits',delta)
	if med_science<0.1:
		xmin=-0.1
		xmax=0.1
	else:
		xmin= med_science - 0.1
		xmax= med_science + 0.1
	print 'Scales are ',xmin,xmax, ' for ',med_science

	start_ds9(scale='linear')


	os.system('xpaset -p ds9 frame delete all')

	command='xpaset -p ds9 file new %s' % 'persist_cts.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits -3e5 3e5')

	print 'Display %s' % record[0] 
	command='xpaset -p ds9 file new %s' % record[0]
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits %.1f %.1f' % (xmin,xmax))

	command='xpaset -p ds9 file new %s' % 'persist_cts3.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits  %.1f %.1f' % (xmin-med_science,xmax-med_science))

	command='xpaset -p ds9 file new %s' % 'subtracted.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits  %.1f %.1f' % (xmin,xmax))



	os.system('xpaset -p ds9 match frames image')

	# Now create some plots
	pylab.figure(1,(6,6))
	pylab.clf()

	# Create a mask of that can select the pixels that were bright at some previous time
	mask=numpy.select([x>50000],[0],default=1) 

	xmed=0.0

	for dtime in dtimes:
		tmask=numpy.select([t==dtime],[0],default=1) # pix first time
		# xmask is the set of pixels that were bright at a specific time in the past
		xmask=mask+tmask

		# Create masked arrays that only contain the pixels that were bright at 
		# at certain time in the past.  xxx contains the counts at an ealire time
		# yyy contions the current count rate
		xxx=numpy.ma.array(x,mask=xmask)
		yyy=numpy.ma.array(science,mask=xmask)

		xxxx=xxx.compressed()
		yyyy=yyy.compressed()

		xxxx=numpy.ravel(xxxx)
	 	yyyy=numpy.ravel(yyyy)

		# At this point xxxx and yyyy are 1-d arrays 

		print 'Median %5.2f for %d pixels at dtime %.2f' % (numpy.median(yyyy), len(yyyy), dtime)
		if xmed==0.0:
			xmed = numpy.median(yyyy)

		# print 'whatever',len(xxxx),len(yyyy)

		# Subtract the background level
		# yyyy=yyyy-med_science

		pylab.plot(xxxx,yyyy,'.')
		fit_fermi(xxxx,yyyy)
		# pylab.axis([0,200000,0,2.*xmed])
		pylab.axis([0,200000,-0.1,2.*(xmed-med_science)])
		pylab.draw()
		persist_median=get_bounded_median(xxxx,yyyy,xmin=100000,xmax=300000)
		print 't %6.2f   Persist %7.3f ' % (dtime,persist_median)

	fx=numpy.linspace(0,xplot_max,1000)
	fy=calc_fermi(fx,fermi_n,fermi_e,fermi_kt)
	# pylab.plot(fx,fy,'-',lw=3)
	for dtime in dtimes:
		# print 'ok',dtime/t_const
	 	fyy=fy*math.exp(-dtime/t_const)
		if subtract_back=='yes':
			fyy=fyy+med_science
			print 'Background has been subtracted from the science image'
	 	pylab.plot(fx,fyy,'-',lw=3)

	pylab.axis([0,xplot_max,-0.1,2*(xmed-med_science)])
	pylab.axis([0,xplot_max,-0.1,1])

	pylab.draw()
	pylab.savefig('subtract2.png')






	return


def subtract_persist_ries(dataset='last',delta=3., fileroot='observations'):
	'''
	Bravely attempt to create an image to subtract using Adam's formulation

	100629 - I'm not clear that this ever worked
	'''

	records=read_ordered_list(fileroot,dataset,delta)

	x,t=make_badpix2(records[0:len(records)-1],'persist')

	record=records[len(records)-1]
	xtime=eval(record[6])
	t=86400*(xtime-t)  # Now t is a delta time in seconds 

	# So at this point time is the time in hours to the dataset we want to correct.
	dtimes=numpy.unique(t)
	print 'delta_times',dtimes           
	tt=calc_atime(t)
	# This is the image which shows the time deecay function
	ksl_fits.put_fits('persist_time2.fits',tt)

	# This is the image which multiples the brightest pixel value by the decay time 
	xx=x*tt
	ksl_fits.put_fits('persist_cts2.fits',xx)

	xx=calc_afermi(xx)

	# This is the final persistence image
	ksl_fits.put_fits('persist_cts3.fits',xx)

	science=get_image(record[0])
	# Get some statistics of the science image
	med_science=numpy.median(science)

	delta=science-xx

	# Write out the subtracted image
	ksl_fits.put_fits('subtracted.fits',delta)
	if med_science<0.1:
		xmin=-0.1
		xmax=0.1
	else:
		xmin= med_science - 0.1
		xmax= med_science + 0.1
	print 'Scales are ',xmin,xmax, ' for ',med_science

	start_ds9(scale='linear')


	os.system('xpaset -p ds9 frame delete all')

	command='xpaset -p ds9 file new %s' % 'persist_cts.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits -3e5 3e5')

	print 'Display %s' % record[0] 
	command='xpaset -p ds9 file new %s' % record[0]
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits %.1f %.1f' % (xmin,xmax))

	command='xpaset -p ds9 file new %s' % 'persist_cts3.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits  %.1f %.1f' % (xmin-med_science,xmax-med_science))

	command='xpaset -p ds9 file new %s' % 'subtracted.fits'
	os.system(command)
	os.system('xpaset -p ds9 zoom to fit')
	os.system('xpaset -p ds9 scale limits  %.1f %.1f' % (xmin,xmax))



	os.system('xpaset -p ds9 match frames image')

	# Now create some plots
	pylab.figure(1,(6,6))
	pylab.clf()

	mask=numpy.select([x>50000],[0],default=1)  # Only worry about bright pixels

	xmed=0.0

	for dtime in dtimes:
		tmask=numpy.select([t==dtime],[0],default=1) # pix first time
		xmask=mask+tmask

		# print 'mask',numpy.unique(mask)
		# print 'tmask',numpy.unique(tmask)
		# print 'xmask',numpy.unique(xmask)
	
		xxx=numpy.ma.array(x,mask=xmask)
		yyy=numpy.ma.array(science,mask=xmask)

		xxxx=xxx.compressed()
		yyyy=yyy.compressed()

		xxxx=numpy.ravel(xxxx)
	 	yyyy=numpy.ravel(yyyy)
		# yyyy=boxcar(yyyy,21)  # This is failed attempt to get more meaningful plots

		print 'Median %5.2f for %d pixels at dtime %.2f' % (numpy.median(yyyy), len(yyyy), dtime)
		if xmed==0.0:
			xmed = numpy.median(yyyy)

		# print 'whatever',len(xxxx),len(yyyy)

		pylab.plot(xxxx,yyyy,'.')
		fit_fermi(xxxx,yyyy)
		pylab.axis([0,200000,0,2.*xmed])
		pylab.draw()

	fx=numpy.linspace(0,200000,1000)
	fy=calc_afermi(fx)
	fy=fy
	# pylab.plot(fx,fy,'-',lw=3)
	for dtime in dtimes:
		# print 'ok',dtime/t_const
	 	fyy=fy*calc_atime(dtime)       
		fyy=fyy+med_science
	 	pylab.plot(fx,fyy,'-',lw=3)

	pylab.axis([0,200000,0,2*xmed])
	pylab.axis([0,200000,0,2])

	pylab.draw()





	return

def boxcar(x,xsize=5):
	'''
	Box car smooth and array

	100629 - Not used by other routines.  Delete ?
	'''

	box=numpy.ones(xsize)
	box=box/xsize

	xx=numpy.convolve(x,box,'same')
	return xx
		
def fit_fermi(x,y):
	'''
	This is a program which is intended to actually fit to a fermi distribution

	100310 - This does not seem to work on real data, most likely because I have too
	many outlyers.
	'''

	import scipy.optimize

	# print 'fit_fermi: ',len(x),len(y)
	# z=numpy.ones_like(y)
	# z=z*0.1


	if len(x)<3:
		print 'fit_fermi: len of input arrays too short: ', len(x),len(y)
		return []

	fitfunc = lambda p, x: calc_fermi(x,p[0],p[1],p[2])
	errfunc = lambda p, x, y: fitfunc(p,x) - y
	# errfunc = lambda p, x, y: numpy.select([abs(fitfunc(p,x) - y)<z],[abs(fitfunc(p,x) - y)],[z])
	p0=[0.1,100000,20000]
	p1,success = scipy.optimize.leastsq(errfunc,p0[:],args=(x,y),ftol=1e-10)

	print p1
	print success

	return p1

def test_fit_fermi():
	'''
	This is simply a little test routine designed to test the fit_fermi
	program above.  It should probably be deleted.

	100309
	'''

	x=numpy.linspace(0,200000,1000)
	y=calc_fermi(x,0.15,35000,10000)
	y=y+0.01*(0.5-scipy.rand(1000))

	fit_fermi(x,y)

def doit(dataset='last',delta_time=24,fileroot='observations',init=1,apertures='full',xbad=100000):
	'''
	main routine for generating a first look at a set of files exhibiting persistence
	'''

	if check_for_ds9():
		 print 'Exiting routine becuse some ds9-related software was not found'
		 return

	if init:
		print 'Generating a time ordered list of .flt files. This may take a while!'
		make_ordered_list(fileroot,apertures)
	else:
		print 'Using previously generated list of .flt files'

	time_sorted=read_ordered_list(fileroot='observations',dataset=dataset,delta_time=delta_time)
	# print 'final',time_sorted
	print 'The number of datasets %d (including the last): ' % len(time_sorted)
	if len(time_sorted)<2:
		print 'Error: doit: The length of the time sorted list is smaller than expected'

	make_jpegs(time_sorted)
	eval_positions(time_sorted)
	# Now construct the bad pixel map from the prevous observations
	make_badpix(time_sorted[0:len(time_sorted)-1],full=xbad,outroot=dataset+'_xhi')



helpstring='''usage: per_eval.py [-h] [-p] [-x time_in_hours] dataset
where:
	-h      produces this message
	-p      causes a previously generated observations.ls file to be used
	-a 	collects data for any aperture including subarrays; otherwise only full reads
	-x time sets the time (in hours) before the named dataset for observations 
		to be displayed
	-e electrons  - Sets the number of electrons to be considered for generating the bad pixel map
and
	dataset is the dataset name for which one wants to investigate persistance

The routine displays .flt files of observations that precede and observation that may
have shown persistance.  

The program requires ds9 and the xpa interface to ds9 to be in your path.
'''

if __name__ == "__main__":
	'''
	If we are running the program from the command line we need to parase sys.argv
	'''
	import sys

	# Put the defaults here

	xtime=12
	dname='last'
	generate_obs=1
	xaps='full'
	xxbad=100000


	if len(sys.argv)==1:
		print 'Processing all .flt files in any subdirectories'
		doit(dataset=dname,delta_time=xtime)
		exit(0)

	i=1
	while i<len(sys.argv):
		if sys.argv[i]=='-h':
                        print '%s\n' % helpstring
			exit(0)
		elif sys.argv[i]=='-x':
			try:
				i=i+1
				xtime=eval(sys.argv[i])
			except NameError:
				print "Problem obtaining time from %s" % sys.argv[i]
				exit(1)
		elif sys.argv[i]=='-e':
			try:
				i=i+1
				xxbad=eval(sys.argv[i])
			except NameError:
				print "Problem obtaining electrons  from %s" % sys.argv[i]
				exit(1)
		elif sys.argv[i]=='-a':
			xaps='all'
		elif sys.argv[i]=='-p':
			print 'Using old observations.ls file'
			generate_obs=0
		# Now get the dataset name
		elif i==(len(sys.argv)-1): # Then get the dataset name
				dname=sys.argv[i]
		else:
			print 'Unknown switch (%s) in the command line; see py_eval -h for allowd options' % sys.argv[i]
			exit(1)
		i=i+1

	if dname=='last':
		print 'Using the last observation as the file with potential ghosts. '
		print 'This may not be what you want, unless you are testing'

	print sys.argv
	print dname
	doit(dataset=dname,delta_time=xtime,init=generate_obs,apertures=xaps,xbad=xxbad)
	exit(0)
