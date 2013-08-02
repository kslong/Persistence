#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine is a routine to subtract the persistence 
from one or more images.  

This is version 2 of this routine and includes an extra term
in the model of persistence at any one time.

The routine requetss a file created by per_list.py that contains
a time-order list of flt files

The outputs include one or more  persistence images, subtracted flt files
and txt files that describe what has happened

The outputs are in a subdirectory of the directory containing the flt files


Description:  

This is a routine which implements a subtraction algorithm for
persistence and creates a modified flt file which contains the corrected
image.  Various other fits files are also created which can be used to
analyze how well one is subtracting the data from the image.

In general, per_list.py should have been run previously 

Basic usage is as follows

subtact_persist.py dateset   - process a single dataset
subtract_persist2.py -all      - process all of the datasets in the .ls file
subtract_persist2.py -all [-after time1]  [-before time2] - process the datasets
	that are in the .ls when the observations occured after time1 and/or before time2
	time1 and time2 are either MJD or ISO formated times, e.g '2009-11-12 15:13:24'
subtract_persist2.py -many file.ls - proces a list of dataset names (one dataset per line)

subtract_persist2.py -obslist obsers2.ls  - specifices something other than the normal
	observation.ls file to use.
subtract_persist2.py -ds9  - causes ds9 to start up and various files to be displayed in
	ds9 during processing
subtract_persisty.py -local - causes the output files to be written in a subdirectory
	of the current working directory instead of with the original flt files as
	is normal

Other switches allow you to control the persistence function that is subtracted, e. g.

-gamma	-- The power law time decay
-n	-- The normalization at 1000 s
-e	-- The fermi energy, nominally the place where the fermi distribution reaches half
	    the maximum
-kT	-- The width of the fermi function
-alpha	-- The power law exponent which accounts for the continued increase of flux at 
	   higher flux levels
-xynorm foo.fits -- replaces the default file which corrects for spatial dependence of persistence 
		with a file foo.fits.  The correction file shold be a file which is on average 1, 
		with values larger or smaller than this indicating the sensitivity to persistence.
		this file is multipled with the overall model
-xynorm none -- Do not make a spatially dependent correction.

	

Primary routines:

Notes:
									   
History:

100603	ksl	Coding begun
101014	ksl	Began adding capabilities to use this in a production-like
		environment
121220	ksl	Darks are not normally processed to e/s but are left in counts/s.
		Modified so that all of the files have the same units as the file
		to be analyzed, even though the base calculation of the model in
		in e/s

'''

import os
import sys
import numpy
import math
import string

import date
import per_list
from per_fits import *
from ds9 import *
from date import *



def calc_fermi(x,norm=1.0,e_fermi=80000,kt=20000,alpha=0.2):
	'''
	calculate the shape of the afterglow assuming a Fermi-Dirac like
	distribution and series of x values.  This routine just calculates 
	the shape.

	The Fermi distribution is defined as 

	1/(exp((E_f-E)/kT)+1) where E_f is the Fermi energy.  For E>E_f the exponential term has little  
	effect.  kT determines how sharpe the cutoff is.  For E=E_f we are at the distribution value of
	the function is 0.5 regardless of kT

	For persistence we see a slow rise beyond this and so we have included an extra term

	(x/e_fermi)**alpha

	to take care of this

	100305 - Note that this operates on any dimension array
	100604 - This should be set up so the normalization is at 1000 s
	1103   - Added a power law term to take care of the continued slow rise of persistence
		at higher stimulus
	110327	Added clipping to eliminated negative values which can occur in dark subtracted
		files
	
	'''

	x=numpy.clip(x,1e-10,1e30)  # Eliminate negative numbers
	y=e_fermi-x
	y=y/kt
	y=numpy.exp(y)
	y=1./(y+1.)
	y=norm*y
	y=y*(x/e_fermi)**alpha
	return y

def calc_decay(dt,norm=0.3,gamma=0.8):
	'''
	Calculate the factor describing the decay of persistent
	'''

	tzero=1000.

	factor=norm*math.pow(dt/tzero,-gamma)
	# print 'calc_decay: dt, factor ',dt,factor

	return factor

def fixup(x,dq):
	'''
	Use the dq image, to find and fix values in the bright image that are likely saturated but 
	have low apparent signal.
	'''

	
	# The next line returns a tuple of arrays
	indices=numpy.where(dq==256)
	xindices=numpy.where(x==0)
	# print 'Bad values in dq ', len(indices),len(indices[0])
	# print 'Bad values in  x ', len(xindices),len(xindices[0])

	i=0
	while i<len(indices[0]):
		ix=indices[0][i]
		iy=indices[1][i]
		# print 'Bad ',ix,iy,x[ix,iy]
		i=i+1
		ixmin=max(0,ix-1)
		ixmax=min(len(x[0])-1,ix+1)
		iymin=max(0,iy-1)
		iymax=min(len(x)-1,iy+1)
		# print 'OK - This is the region',ixmin,ixmax,iymin,iymax
		try:
			# z=max(x[ixmin:ixmax+1][iymin:iymax+1])
			# print 'Values',x[ixmin:ixmax,iymin:iymax]
			z=numpy.max(x[ixmin:ixmax,iymin:iymax])
			# print'Gotcha z',z
			x[ix,iy]=z
		except ValueError:
			print 'No Joy with ',ixmin,ixmax+1,iymin,iymax+1
			print 'Lengths ',len(x),len(x[0])
			print x[ixmin][iymin],x[ixmin][iymax],x[ixmax][iymin],x[ixmax][iymax]


	xindices=numpy.where(x==0)
	# print 'Bad values after ', len(xindices),len(xindices[0])

	return x
	


def calc_persist(x,dq,dt=1000,norm=0.3,alpha=0.2,gamma=0.8,e_fermi=80000,kT=20000):
	'''
	Calculate the persistence of a single image at dt in seconds.  The persistence
	image is returned.

	Note that an image is read in and converted to electons if necessary, but an
	array is returned.  Also, there is no attempt to treat pixels with dq flags
	differently, which mens one is dependent on whatever flux there is there.

	Note also that this routine does not allow for spatial variation in persistence.
	That correction is applied at the in do_dataset()

	101207	ksl	Modified so the arrays are passed to calc_persist rather than
			files which must be read  in order to reduce the 
			number of times an image is read and to allow for subarrays
	110926	ksl	Fixed to avoid the problem due to the very occassional absence
			of a dq array


	'''

	# Next steps are to handle the problem associated with the fact that at least in pre-mid 2010 versions
	# of calwf3 there are certain values in the flt file that were set to zero because it was hard to
	# estimate the count rate from them due to the fact that the souce was very bright

	# 110926 - ksl - Avoid fixup if there is no dq array
	if len(dq)>0:
		dq=numpy.bitwise_and(dq,256)
		x=fixup(x,dq)

	# Calculate the persistence image

	# print 'calc_persist: median ',numpy.median(x)

	x=calc_fermi(x,1.,e_fermi,kT,alpha)

	# print 'calc_persist: median ',numpy.median(x)

	decay=calc_decay(dt,norm,gamma)

	# Return the persistence image

	x=x*decay

	# print 'calc_persist: median ',numpy.median(x)

	return x



def how_much(image):
	'''
	Calculate some basic statistics about how much persistence an image might couase

	Also give a qualitative measure of the persistence,
	'''

	q=numpy.ravel(image)
	qq=numpy.sort(q)
	n1=0
	n2=0
	n3=0
	i=len(qq)-1
	while i>0:
		value=qq[i]
		if value>0.1:
			n1=n1+1
		if value>0.03:
			n2=n2+1
		if value>0.01:
			n3=n3+1
		else:
			break
		i=i-1
	i90=int(0.9*len(qq))
	i99=int(0.99*len(qq))
	i50=int(0.50*len(qq))

	values=[len(qq),n1,n2,n3,qq[len(qq)-1],qq[i50],qq[i90],qq[i99]]


	return values


def get_stats(image,saturation=70000):
	'''
	Calculate basic image statistics for an image that is 
	in electrons and return the midpoint and then number of
	pixels greater than saturateion
	'''

	q=numpy.ravel(image)
	qq=numpy.sort(q)
	midpoint=int(0.5*len(qq))

	i=len(qq)-1
	while i>0 and qq[i]>saturation:
		i=i-1
	
	npts=len(qq)-i

	values=[qq[midpoint],npts]
	return values










def do_dataset(dataset='ia21h2e9q',norm=0.3,alpha=0.2,gamma=0.8,e_fermi=80000,kT=20000,fileroot='observations',ds9='yes',local='no',xynorm='persist_corr.fits'):
	'''
	Create a persistence image for this dataset.  This version works by creating using the 
	persistance model from each of the previous images.  It assumes that the only one which
	matters is the image which created the most persistence according toe the modeld

	All returns from this program should be:
		OK:  something
		NOK: something

	History

	100905	ksl	Added disgnostic files to record the time of the stimulus and the
			value of the stimulus
	101015	ksl	Moved the location of where the output files are stored and 
			added a history file, so we would know what had happened
	110122	ksl	Added a switch to turn off displaying with ds9
	110324	ksl	Changed call to accommodate new persistence model
	110602	ksl	Incorporated the possibility of providing a correction file xynorm to account
			for spatial variations across the detector persistence model
	'''

	cur_time=date.get_gmt()

	print '# Processing dataset %s Norm %4.2f alpha %4.2f gamma %4.2f e_fermi %6.0f kT %6.0f ' % (dataset,norm,alpha,gamma,e_fermi,kT)

	# Read the observations file to get the data set of interest

	delta=8  # Hardwired value for consideration of persistence.  Datasets which occurred more than delta hours earlier not considered.

	records=per_list.read_ordered_list2(fileroot,dataset,interval=[-delta,0],outroot='none')

	# Check the length of records


	if len(records)==0:
		string = 'NOK: subtract_persist.do_dataset :There are no records associated with dataset %s.  Check name in %s.ls' % (dataset,fileroot)
		sum_string='NOK - No record associated with this dataset'
		per_list.update_summary(dataset,'ERROR',sum_string,append='no')
		return string

	# So now we have the list that we need.

	

	science_record=records[len(records)-1]  # The science record is the last record
	sci_progid=science_record[2]
	words=science_record[3].split('.')
	sci_visit=words[0]
	sci_fil=science_record[10]
	sci_exp=eval(science_record[11])
	sci_obj=science_record[14]
	
	# Create the Persist directory if it does not already exist
	path=per_list.set_path(science_record[0],'yes',local)
	if path.count('NOK'):  # Then we were not able to create a plausible directory to put the data ink
		return path

	# Open a history file.  Note that one needs the path before one can do this

	history=per_list.open_file(path+science_record[1]+'.txt')

	# history=open(path+science_record[1]+'.txt','w')
	# os.chmod(path+science_record[1]+'.txt',0770)

	history.write('START:  Persistence processing of file %s\n\n' % science_record[1])
	history.write('! Processed: %s\n' % date.get_gmt())
	history.write('! ProgramID: %s\n' % sci_progid)
	history.write('! Visit:     %s\n' % sci_visit)
	history.write('! FltFile:   %s\n' % science_record[1])
	history.write('! Filter:    %s\n' % sci_fil)
	history.write('! Exposure:  %6.1f\n' % sci_exp)
	history.write('! Object:    %s\n' % sci_obj)

	history.write('\n! Persistence model: norm %6.2f alpha %6.2f e_fermi %6.0f kT %6.0f\n' % (norm,alpha,e_fermi,kT)) 


	# Check whether there is anything to do

	if len(records)==1:
		string='subtract_persist: No persistence for this dataset.  No earlier observations within %.1f hours\n' % (delta)
		history.write('%s\n' % string)
		history.write('! Persistence:  None\n')
		string='OK: subtract_persist: -- None'
		print string
		history.close()
		xstring='  %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f' % (0,0,0,0,0,0)
		per_list.update_summary(dataset,'Persist',xstring,append='no')
		return string



	[t1,t2]=get_times(science_record[0])
	tbright=0.5*(t1+t2)

	# establish which record is the last record that comes from a different visit than the current one

	i=0
	while i<len(records)-1:
		cur_progid=records[i][2]
		words=records[i][3].split('.')
		cur_visit=words[0]
		# print 'OK',cur_progid,cur_visit
		if cur_progid==sci_progid and cur_visit==sci_visit:
			# then we are now into earlier exposures of the same visit
			break
		i=i+1

	# if i = 0, all the exposures being evaluated for persistence create self_persistence
	# if i= len(records), all exposures being evaluated for persistence are from other visits

	last_external=i-1  # This is the last record that created external persistence
	ext_values=[]  # A place to store information about the persistence due to other observers
	ext_persist=[] # This is a place holder for storing the extenal persistence

	if xynorm!='':
		xcorr=get_image(xynorm,1)
		if len(xcorr)==0:
			xynorm=''  # This is an error because we were unable to find the file




	# Now we start actually creating the persistence model
	i=0
	while i<len(records)-1:
		record=records[i]
		# print 'xxx ',i,len(records),record[0]
		[t1,t2]=get_times(record[0])
		dt=(tbright-t2)*86400

		cur_progid=record[2]
		words=record[3].split('.')
		cur_visit=words[0]
		cur_sci_fil=record[10]
		cur_sci_exp=eval(record[11])
		cur_sci_obj=record[14]

		x=get_image(record[0],1,'e',fileref=science_record[0])  # Convert this to electrons
		if len(x)==0:
			xstring='NOK: Problem with science extension of %s' % record[0]
			history.write('%s\n' % xstring)
			print xstring
			return xstring
		dq=get_image(record[0],3,fileref=science_record[0])     # Get the dq 
		if len(dq)==0:
			xstring = 'NOK: Problem with dq extension of %s' % record[0]
			history.write('%s\n' % xstring)
			print xstring
			# 110926 - ksl - modified to allow this to process the image even if there was no dq array
			# return xstring

		if i==0:
			persist=calc_persist(x,dq,dt,norm,alpha,gamma,e_fermi,kT)
			values=how_much(persist)
			stimulus=x   # This is an array which contains the maximum counts in a pixel
			xtimes=numpy.ones_like(persist)
			xtimes=xtimes*dt # This is an array containing the delta time at which the stimulus occured
		else:
			xpersist=calc_persist(x,dq,dt,norm,alpha,gamma,e_fermi,kT)
			values=how_much(xpersist)
			stimulus=numpy.select([xpersist>persist],[x],default=stimulus)
			xtimes=numpy.select([xpersist>persist],[dt],default=xtimes)
			persist=numpy.select([xpersist>persist],[xpersist],default=persist)
		

		# Get some elementary statistics on the stimulus
		xvalues=get_stats(x,70000)

		history.write('\nsubtract_persist: Stimulus by %30s from program %s Visit %s\n' % (record[0],cur_progid,cur_visit))
		history.write('\nsubtract_persist: The filter was %s and exposure was %6.1f for target %s\n' % (cur_sci_fil,cur_sci_exp,cur_sci_obj))
		history.write('subtract_persist: The exposure was %8.0f s earlier than the current exposure\n' % dt)
		history.write('subtract_persist: The median value in the stimulus image was %6.1f and the number of saturated pixels was %d\n' % (xvalues[0],xvalues[1]))
		history.write('subtract_persist:   The maximum value for persistence is %f\n'   % values[4])
		history.write('subtract_persist:    The median value for persistence is %f\n'   % values[5])
		history.write('subtract_persist: 90 percent of persist values less than %f\n'   % values[6]) 
		history.write('subtract_persist: 99 percent of persist values less than %f\n'   % values[7]) 
		history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.10 e/s\n' % (values[1],values[1]*100./values[0]))
		history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.03 e/s\n' % (values[2],values[2]*100./values[0]))
		history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.01 e/s\n' % (values[3],values[3]*100./values[0]))
		string = 'subtract_persist: Finished (%2d of %2d) %s' % (i+1,len(records)-1,record[0])
		print string
		history.write('%s\n' % string)

		# Summarize the Stiumulus printing out the filename, prog_id, visit_name, target, dt and the number of pixels above saturation of 70000
		stimulus_summary='Stimulus: %40s %10s %10s %20s %8.0f %3d\n' % (record[0],cur_progid,cur_visit,cur_sci_obj,dt,xvalues[1])
		history.write('! %s\n' % stimulus_summary)

		if i==last_external:
			ext_values=how_much(persist)
			ext_persist=numpy.copy(persist)

		i=i+1
	
	# First report on the external persistence for this file;

	# Now apply the fix to account for spatial variations in persistence
	if xynorm !='' and numpy.shape(xcorr)==numpy.shape(persist):
		persist=persist*xcorr


	if len(ext_values)>0:

		f1=ext_values[1]*100./ext_values[0]
		f2=ext_values[2]*100./ext_values[0]
		f3=ext_values[3]*100./ext_values[0]
		emeasure='%6.2f %6.2f %6.2f' % (f1,f2,f3)

		history.write('\nsubtract_persist: Estimate of persistence from earlier visits\n')
		history.write('subtract_persist:   The maximum value for persistence is %f\n'   % ext_values[4])
		history.write('subtract_persist:    The median value for persistence is %f\n'   % ext_values[5])
		history.write('subtract_persist: 90 percent of persist values less than %f\n'   % ext_values[6]) 
		history.write('subtract_persist: 99 percent of persist values less than %f\n'   % ext_values[7]) 
		history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.10 e/s\n' % (ext_values[1],f1))
		history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.03 e/s\n' % (ext_values[2],f2))
		history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.01 e/s\n' % (ext_values[3],f3))
	else:
		emeasure='%6.2f %6.2f %6.2f' % (0 ,0,0)
		history.write('\nsubtract_persist: This exposure has no persistence from earlier visits.  All persistence is self-induced\n')



	# Now evaluate the total persistence

	values=how_much(persist)


	f1=values[1]*100./values[0]
	f2=values[2]*100./values[0]
	f3=values[3]*100./values[0]
	measure='%6.2f %6.2f %6.2f' % (f1,f2,f3)


	history.write('\nsubtract_persist: Estimate of total persistence for this file\n')
	history.write('subtract_persist:   The maximum value for persistence is %f\n'   % values[4])
	history.write('subtract_persist:    The median value for persistence is %f\n'   % values[5])
	history.write('subtract_persist: 90 percent of persist values less than %f\n'   % values[6]) 
	history.write('subtract_persist: 99 percent of persist values less than %f\n'   % values[7]) 
	history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.10 e/s\n' % (values[1],f1))
	history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.03 e/s\n' % (values[2],f2))
	history.write('subtract_persist: %7d pixels (or %6.3f percent) greater than 0.01 e/s\n' % (values[3],f3))




	history.write('! PersistenceSum: External   %s\n'% emeasure)
	history.write('! PersistenceSum: Total      %s\n'% measure)


	# Now write out the persistence image

	# First find out the units of the science image

	units=get_keyword(science_record[0],1,'bunit')
	print 'test',units
	if units[0]=='COUNTS/S':
		print 'reducing'
		persist/=2.4

	# subtract and write out the corrected image
	science=get_image(science_record[0],1,'no')
	original=numpy.copy(science)
	science=science-persist


	xname=parse_fitsname(science_record[0])

	persist_file=path+dataset+'_persist.fits'
	ext_persist_file=path+dataset+'_extper.fits'

	# Note: Do not use some with an extension like flt.fits because it will interfere with making list
	corrected_file=path+dataset+'_flt_cor.fits'
	stimulus_file=path+dataset+'_stim.fits'
	time_file=path+dataset+'_dt.fits'

	rewrite_fits(xname[0],persist_file,1,persist)
	rewrite_fits(xname[0],corrected_file,1,science)
	rewrite_fits(xname[0],stimulus_file,1,stimulus)
	rewrite_fits(xname[0],time_file,1,xtimes)
	if len(ext_persist)>0:
		rewrite_fits(xname[0],ext_persist_file,1,ext_persist)
		

	# Get statistics on the images 
        xmed=numpy.median(original)
        zmed=numpy.median(persist)
        zmax=numpy.max(persist)


	pylab.figure(1,[12,12])
	pylab.title(dataset)
	pylab.subplot(221)
	pylab.imshow(original,origin='lower',cmap=pylab.cm.gray,vmin=xmed-0.1,vmax=xmed+0.1)
	pylab.title('Original')

	pylab.subplot(222)
	pylab.imshow(science,origin='lower',cmap=pylab.cm.gray,vmin=xmed-0.1,vmax=xmed+0.1)
	pylab.title('Subtracted')

	pylab.subplot(223)
	pylab.imshow(persist,origin='lower',cmap=pylab.cm.gray,vmin=zmed-0.1,vmax=zmed+0.1)
	pylab.title('Total Persistence')

	if len(ext_persist)>0:
		pylab.subplot(224)
		pylab.imshow(ext_persist,origin='lower',cmap=pylab.cm.gray,vmin=zmed-0.1,vmax=zmed+0.1)
		pylab.title('External Persistence')
	else:
		pylab.subplot(224)
		pylab.imshow(stimulus,origin='lower',cmap=pylab.cm.gray,vmin=0.0,vmax=200000)
		pylab.title('Stimulus')

	fig1=path+'Figs/'+dataset+'_subtract.png'


	if os.path.isfile(fig1):
		os.remove(fig1)
	pylab.savefig(fig1)
	os.chmod(fig1,0770)

	# Eliminated to prevent an error on linux having to do with tkinter
	# pylab.close('all')

	if ds9=='yes':
		LoadFrame(science_record[0],1,0,2,'histequ')
		LoadFrame(persist_file,2,0,2,'histequ')
		LoadFrame(corrected_file,3,0,2,'histequ')
		if len(ext_persist)>0:
			LoadFrame(ext_persist_file,4,0,2,'histequ')
		else:
			LoadFrame(stimulus_file,4,0,1e5,'histequ')

	history.write('# Finished Persistence processing of file %s\n' % science_record[1])
	history.close()
	string='%20s %20s' % (emeasure,measure)
	per_list.update_summary(dataset,'Persist',string,fileroot,append='no')
	return string


def steer(argv):
	'''
	This is a steering routine for subtract persist so that options can be exercised from the 
	command line

	100907	ksl	Added to begin to automate the subtraction process
	'''
	i=1
	dataset_list='none'

	norm=0.3
	alpha=0.2
	gamma=0.8
	e_fermi=80000
	kT=20000
	fileroot='observations'
	xynorm='persist_corr.fits'
	words=[]
	mjd_after=0.0    # A amall number for mjd
	mjd_before=1.e6  # A large number for mjd

	ds9='no'
	local='no'

	while i<len(argv):
		if argv[i]=='-h':
			print __doc__
			return    
		elif argv[i]=='-n':
			i=i+1
			norm=eval(argv[i])
		elif argv[i]=='-e':
			i=i+1
			e_fermi=eval(argv[i])
		elif argv[i]=='-kT':
			i=i+1
			kT=eval(argv[i])
		elif argv[i]=='-alpha':
			i=i+1
			alpha=eval(argv[i])
		elif argv[i]=='-gamma':
			i=i+1
			gamma=eval(argv[i])
		elif argv[i]=='-obslist':
			i=i+1
			fileroot=eval(argv[i])
		elif argv[i]=='-many':
			i=i+1
			dataset_list=argv[i]
			print 'OK you want to evaluate a number of datasets in file %s', dataset_list
		elif argv[i]=='-all':
			dataset_list='!All'
			print 'OK you want to evaluate all the records in the obslist'
		elif argv[i]=='-after':
			i=i+1
			z=argv[i]
			try:
				mjd_after=float(z)
			except ValueError:
				mjd_after=iso2mjd(z)
		elif argv[i]=='-before':
			i=i+1
			z=argv[i]
			try:
				mjd_before=float(z)
			except ValueError:
				mjd_before=iso2mjd(z)
		elif argv[i]=='-ds9':
			ds9='yes'
		elif argv[i]=='-local':
			local='yes'
		elif argv[i]=='-xy':
			i=i+1
			xynorm=argv[i]
			if xynorm=='none':
				xynorm=''
		elif argv[i][0]=='-':
			print 'Error: Unknown switch ---  %s' % argv[i]
			return
		else:
			words.append(argv[i])
		i=i+1
	
	if dataset_list=='none': #  Then we are processing a single file
		dataset=words[0]
		do_dataset(dataset,norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,xynorm)
	# Note it is not clear to me what is going on here.  What distiguishes the next
	# to options, and why isn't per_list being used  111019 -- ksl  !!!
	elif dataset_list=='!All': # Then we are working from the obslist
		f=open(fileroot+'.ls','r')
		lines=f.readlines()
		f.close()

		for line in lines:
			x=line.strip()
			word=x.split()
			if len(word)>1 and x[0]!='#':
				mjd=float(word[6])
				print mjd,mjd_before,mjd_after
				if mjd_after  <= mjd and mjd <= mjd_before:
					do_dataset(word[1],norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,xynorm)
	else:
		f=open(dataset_list,'r')
		lines=f.readlines()
		f.close()

		for line in lines:
			x=line.strip()
			if len(x)>0 and x[0]!='#':
				mjd=float(x[6])
				if mjd_after <= mjd and mjd <= mjd_before:
					do_dataset(x,norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,xynorm)
	
	return

	



	 

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		steer(sys.argv)   
	else:
		print 'subtract_persist2.py  -h to get brief  help text'
	print 'OK done'
	sys.exit(0)
