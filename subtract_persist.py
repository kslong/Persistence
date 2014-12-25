#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine is a routine to subtract the persistence 
from one or more images.  

This version supports two types of models:
	those based on a modified ferrmi formula and
	those based on power laws which depend on expousre time and fluence

The routine requets a file created by per_list.py that contains
a time-order list of flt files

The outputs include one or more persistence images, subtracted flt files
and txt files that describe what has happened

The outputs are in a subdirectory of the directory containing the flt files
or alternatively if the local swithc is set in a subdirectory of the
current working direcory


Description:  

This is a routine which implements a subtraction algorithm for
persistence and creates a modified flt file which contains the corrected
image.  Various other fits files are also created which can be used to
analyze how well one is subtracting the data from the image.

In general, per_list.py should have been run previously 

Basic usage is as follows

subtact_persist.py dateset   - process a single dataset

subtract_persist.py -all      - process all of the datasets in the .ls file

subtract_persist.py -all [-after time1]  [-before time2] - process the datasets
	that are in the .ls when the observations occured after time1 and/or before time2
	time1 and time2 are either MJD or ISO formated times, e.g '2009-11-12 15:13:24'

subtract_persist.py -many file.ls - process a list of dataset names (one dataset per line)

subtract_persist.py -obslist obsers2.ls  - specifices something other than the normal
	observation.ls file to use.

subtract_persist.py -ds9  - causes ds9 to start up and various files to be displayed in
	ds9 during processing

subtract_persisty.py -local - causes the output files to be written in a subdirectory
	of the current working directory instead of the directory containing the
	the original flt files.

Other switches allow you to control the persistence function that is subtracted, e. g.

-model  -- 0 for the original fermi-function based formalism
	   1 for the newer purely observational formula  (This is the defualt)
	   2 for a fermi function with modified for different stimulus exposure times
-gamma	-- The power law time decay
-n	-- The normalization at 1000 s
-e	-- The fermi energy, nominally the place where the fermi distribution reaches half
	    the maximum
-kT	-- The width of the fermi function
-alpha	-- The power law exponent which accounts for the continued increase of flux at 
	   higher flux levels
-pf        foo.pf  -- replaces the default parameter file persist.pf with a new one

Calibration files:

	subtract_persist uses certain calibration files that can be located in the directory 
	from which subtract_persist is being run, a subdirectory that has to be names PerCal
	or in a directory defined by an environment variable PERCAL.  

	At present the calibration files that are needed are 

	persist_corr.fits 
	a_gamma.fits
	fermi.fits

Outputs:

	See do_dataset
	

Primary routines:

	The main routine is do_dataset which handles an individual dataset.  This is
	the routine called by run_persist

	If the routine is called from the command line, then the routine steer handles
	the processing of command line switches and queueing up the datasets.

Notes:
									   
History:

100603	ksl	Coding begun
101014	ksl	Began adding capabilities to use this in a production-like
		environment
121220	ksl	Darks are not normally processed to e/s but are left in counts/s.
		Modified so that all of the files have the same units as the file
		to be analyzed, even though the base calculation of the model in
		in e/s
140609	ksl	Began to incoporate a new model into this routine.  For now at least, we will retain
		the old model as well. This is a complication, but it may be some time before this
		works well.
140805	ksl	At this point in time, subtract supports 3 different models, the original one, and
		one based on a gamma style fits to datasets with different exposure times, and one
		based of fermi fits to different times.  As of this date, the calibration files that
		are needed to run subtract can be found either in the local directory, a subdirectory
		called PerCal or in a directory defined by the environment variable PERCAL
140924	ksl	Incorporated my test version  of subtract, called xsubtract into the standard persistence
		package.
'''

import os
import sys
import numpy
import math
import string

import date
import per_list
from per_fits import *
import astropy.io.fits
from ds9 import *
from date import *

from scipy.interpolate import interp1d


# This is the section where the new a gamma model is calculated 
# These are the global variables which hold the models
model_exp=[]
model_stim=[]
model_a=[]
model_g=[]
model_file=''


def read_file(filename,char=''):
	'''
	Read a file and split it into words, eliminating comments
	
	char is an optional parameter used as the delimiter for
	splitting lines into words.  Otherwise white space is
	assumed.

	History:
	
	110729	ksl	Added optional delimiters
	140617  ksl	Included here so that the A gamma models 
			could be read in.  
	'''

	try:
		f=open(filename,'r')
		xlines=f.readlines()
		f.close()
	except IOError :
		print "The file %s does not exist" % filename
		return []   
	
	lines=[]
	
	i=0
	while i<len(xlines):
		z=xlines[i].strip()
		if char=='':
			z=z.split()
		else:
			z=z.split(char)
		if len(z)>0:
			if z[0][0]!='#':
				lines=lines+[z]
		i=i+1
	return lines

def read_models(filename='per_fermi/fermi.fits'):
	'''
	Read in the fits models for A gamma style models


	Notes:

	Fits reading returns numpy.arrays, which unfortunalely
	do not perform identically to lists in the situation
	where one wants to append rows to one another.  The
	equvalent on append for a list, is vstack, but it works
	differently in the sense that one has to say

	y=numpy.vstack((y,z))

	to produce a new y

	A way to avoid this is to make the outside maxtrix a list and then
	convert it to an array at the end.

	140803	ksl	Coded
	140805	ksl	Replaced older read_models routine
	'''

	global model_exp
	global model_stim
	global model_a
	global model_g
	global model_file
	if filename==model_file and len(model_stim)>0:
		return 'OK'

	xfilename=locate_file(filename)


	try:
		x=astropy.io.fits.open(xfilename)
	except IOError:
		print 'read_models: file %s does not appear to exist' % filename
		return 'NOK'

	i=1
	model_exp=[]
	model_a=[]
	model_stim=[]
	model_g=[]
	while i<len(x):
		model_exp.append(x[i].header['exp'])
		tabdata=x[i].data

		one_stim=tabdata['stim']
		model_stim.append(one_stim)

		one_a=tabdata['a']
		model_a.append(one_a)

		one_gamma=tabdata['gamma']
		model_g.append(one_gamma)
		# print 'diag', one_stim.shape
		i=i+1
	
	model_stim=numpy.array(model_stim[0]) # Use only the first row for the stimulus
	model_a=numpy.array(model_a)
	model_g=numpy.array(model_g)

	


	# print 'read_models:',model_stim.shape,model_a.shape,model_g.shape

	model_file=filename

	return



	






def read_models_orig(filename='per_model/models.ls'):
	'''
	Read in the ascii model files for the A gamma models

	Notes:
		Ultimately a better format for the model files hould be writtne at which point
		this routine will need to be rewritten

	140617	ksl	Coded without ksl.io
	140804 	ksl	This is replaced by the new read_models, and SHOULD BE REMOVED ONCE
			I AM SURE THERE IS NO RESIDUAL PROBLEM WITH READING THE FITS FILES.
	'''


	global model_exp
	global model_stim
	global model_a
	global model_g
	global model_file
	if filename==model_file and len(model_stim)>0:
		return 'OK'

	if filename.count('.fits')>0:
		read_models2(filename)
		return 


	model_exp=[]
	model_stim=[]
	model_a=[]
	model_g=[]

	lines=read_file(filename)
	if len(lines)==0:
		return 'NOK'

	i=0
	for line in lines:
		model_exp.append(eval(line[1]))
		xlines=read_file(line[0])
		a=[]
		g=[]
		for xline in xlines:
			if i==0:
				model_stim.append(eval(xline[0]))
			a.append(eval(xline[1]))
			g.append(eval(xline[2]))
		model_a.append(a)
		model_g.append(g)
		i=i+1
	model_exp=numpy.array(model_exp)
	model_stim=numpy.array(model_stim)
	model_a=numpy.array(model_a)
	model_g=numpy.array(model_g)
	# print model_exp.shape,model_stim.shape,model_a.shape,model_g.shape


	model_file=filename
	# print 'read_models: Finished reading ',filename

	print 'test',model_exp

	print model_stim.shape,model_a.shape,model_g.shape

def read_parameter(parameter_file,parameter):
	'''
	Read a single parameter, such as a calibration file name
	from a parameter file

	Notes:

	This simply finds the parameter file, reads it, and
	returns the parameter as a string

	History:

	141225	ksl	Coded to try to make it easier to switch calibration
			files (which are currently hardcoded
	'''

	parameter_file=locate_file(parameter_file)

	if parameter_file=='':
		print 'Error: read_parameter: Could not read parameter becasue could not locate %s' % parameter_file
		return ''

	f=open(parameter_file)
	lines=f.readlines()
	f.close()

	value=''
	for line in lines:
		line=line.split()
		# print parameter, line
		if len(line)>1 and line[0]==parameter:
			value=line[1]
			break
	if value=='':	
		print 'Error: Could not locate %s in %s'  % (parameter,parameter_file)

	if value=='None' or value == 'none':
		value=''
	
	return value

	
def locate_file(filename):
	'''
	This routine locates a file in one of several places in priority order
	and returns the path to the file

	In the directory where a program is being run
	In a directory PerCal underneath the directory where the program is being run
	In a directory defined by an external variable PERCAL

	If the file is not found then the routine returns an empty string

	Notes:

	This was coded so that by defining an environment variable then subtract_persist
	should work everywhere.  



	

	History:

	140805	ksl	Initially coded
	'''

	if os.path.isfile(filename):
		return filename

	xfile='PerCal/'+filename
	if os.path.isfile(xfile):
		return xfile

	try:
		xpath=os.environ['PERCAL']
	except KeyError:
		print 'Error: subtract_persist.locate_file: Environment variable PERCAL not defined and %s not found either in local directory or PerCal subdiretory' % filename
		return ''

	xfile='%s/%s' % (xpath,filename)
	if os.path.isfile(xfile):
		return xfile
	
	print 'Error: subtract_persist.locate_file: %s not found in the local directory, the PerCal subdiretory, or in the direcotry %s defined by PERCAL' % (filename,xpath)
	return ''

def get_persistence(exp=300.,dt=1000.,models='per_model/models.ls'):
	'''
	Calculate the persistence curve using the tabulated A_gamma model

	where:

	exp	the expousre of the stimulus image
	dt	the time since the stimulus image was taken
	models	the file containing times and links to the persitance
		curves

	We need to interpolate both for the actual stimulus and
	for the fact that we don't have a grid of exposures.

	This just returns the persistence curve (on a specific grid) for a given exposure
	time exp and delta time dt.

	Note:

	Although this was initially developed for a phenomenalogical model where an 
	amplitude and a power law index had be calculated but fitting observatonal data, 
	any model with includes a power law decay can be case in the from of the so-called
	A-gama model, including a fermi-type model.

	History:

	140630	ksl	Added a variable models to allow one to read files in any directory
			of interest

	'''

	# if len(model_stim)==0:
	# 	print 'Reading models',models
	# 	read_models(models)
	read_models(models)

	
	# print 'get_persistence:',len(model_exp),len(model_stim),len(model_a),len(model_g)
	
	# Now we need to interpolate so we have a single model given
	# an exposure time
	i=0
	while i<len(model_exp) and model_exp[i]<exp:
		i=i+1

	# print 'get_persitence:',i,len(model_exp),model_exp[i],exp

	if i==0:
		persist=model_a[0]*(dt/1000.)**-model_g[0]
	elif i==len(model_exp):
		i=i-1
		persist=model_a[i]*(dt/1000.)**-model_g[i]
	else:
		frac=(exp-model_exp[i-1])/(model_exp[i]-model_exp[i-1])
		persist1=model_a[i-1]*(dt/1000.)**-model_g[i-1]
		persist2=model_a[i]*(dt/1000.)**-model_g[i]
		persist=(1.-frac)*persist1+frac*persist2
		# Note as written, frac is the frac of model for the longer exposure to include
		# print 'get_persistence: time',exp,model_exp[i-1],model_exp[i]
		# print 'get_persistence: frac',frac,'low',i-1,'hi',i

	return persist

def make_persistence_image(x,exptime=500,dt=300,models='per_model/models.ls'):
	'''
	Make the persistence image for the A gamma model, given 
	
	an image array x,
	an exposure time for the stimulus image exptime.
	a time since the end of the (last) stimulus image
	a file that contains pointers to the individual persistence amplitude curve

	Note that this routine calls get_persistenee, which returns a persistence curve
	for a particular exptime and dt. The persistence curve is on a fixed grid. Here we
	use scipy.inter1d to produce the persisence image.
	'''

	# print 'make_persistence',exptime,dt,models

	persist_curve=get_persistence(exptime,dt,models)

	# print 'OK got to make_persistence_image:',len(model_stim),len(persist_curve)

	# print 'make_persistence: persist_curve:',persist_curve[0:400:10]

	# Now we need to interpolate this curve
	f=interp1d(model_stim,persist_curve,fill_value=0,bounds_error=False)

	persist=f(x)

	return persist


# These are the routines used to calculate the original fermi formula based model

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

	Notes:
	With older versions of the WRC3 pipeline, saturated pixels in the cores of bright stars were often 
	set to zero.  With the cores set to zero, the estimate of persistence for these were zero.  This
	is an attempt to set the value of a pixel to something plusible.  Newer versios place a large value
	in the pixel, which avoids this preblem.  

	Though better than no estimate at all, what's here still leaves a lot to be desired.
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
	


def calc_persist(x,dq,dt=1000,norm=0.3,alpha=0.159,gamma=1.022,e_fermi=83400,kT=18500):
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
	140812  ksl	Changed parameters to reflect the best of a collection of 
			visits from Cycle 21 program (in prep for cal workshop)

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










def do_dataset(dataset='ia21h2e9q',model_type=0,norm=0.3,alpha=0.2,gamma=0.8,e_fermi=80000,kT=20000,fileroot='observations',ds9='yes',local='no',parameter_file='persist.pf'):
	'''
	Create a persistence image for this dataset.  This version works by creating using the 
	persistance model from each of the previous images.  It assumes that the only one which
	matters is the image which created the most persistence according toe the modeld

	model_type=0 is our our orignal model which has a modified fermi distributioon, controlled by the values of norm, etc.
	model_type=1 is our purely describtive model defined in terms of amplitudes and power laws
	model_type=2 for the fermi model that interpolates bewtween cureves for different expsoure times

	All returns from this program should be:
		OK:  something
		NOK: something

	Outputs:
		If there were IR observations that preceded the observation being analyzed
		then this routines creates various fits files:
			rootname_persist.fits  - The model for the total of internal and external persistnce
			rootname_extper.fits   - The model for persistence from earlier visits, aka
			                         external persistence
			rootname_flt_cor.fits  - The corrected flt file
			rootname_stim.fits     - The stimulus that caused the persistence
			rootname_dt.fits       - The time at which the stimulus occurred
		Plots of the various images are also created, as well as a log file for
		the data set
	Notes:
		This is really the main routine of the program.  When run_persist.py calls this modeule. 
		this is the routine that is called.

		In practice, the only difference between model types 1 and 2 is the calibration file that is read in.  The same interpolation
		routine is used for both.  At present the calibration file names for this are hardwired.

	History

	100905	ksl	Added disgnostic files to record the time of the stimulus and the
			value of the stimulus
	101015	ksl	Moved the location of where the output files are stored and 
			added a history file, so we would know what had happened
	110122	ksl	Added a switch to turn off displaying with ds9
	110324	ksl	Changed call to accommodate new persistence model
	110602	ksl	Incorporated the possibility of providing a correction file xynorm to account
			for spatial variations across the detector persistence model
	140606	ksl	Added correction which puts the correct units in the stimulus file.
	140611	ksl	Began to add in new persistnce model, initially just by short-circuiting everything
	140803	ksl	Switched to fits version of data files
	'''

	cur_time=date.get_gmt()

	xfile=locate_file(parameter_file)
	if xfile=='':
		print '# Error: Could not locate parameter file %s ' % parameter_file
	else:
		parameter_file=xfile


	if model_type==0:
		print '# Processing dataset %s with fermi model: Norm %4.2f alpha %4.2f gamma %4.2f e_fermi %6.0f kT %6.0f ' % (dataset,norm,alpha,gamma,e_fermi,kT)
	elif model_type==1:
		print '# Processing dataset %s with A gamma model' % dataset
	elif model_type==2:
		print '# Processing dataset %s with a time-variable fermi model' % dataset
	else:
		print '# Error: run_persist: Unknown model type %d' % model_type
		return 'NOK'

	# Read the observations file to get the data set of interest

	# delta=8  # Hardwired value for consideration of persistence.  Datasets which occurred more than delta hours earlier not considered.
	# Increased to 16 hours 140617
	delta=16  # Hardwired value for consideration of persistence.  Datasets which occurred more than delta hours earlier not considered.

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

	if model_type==0:
		history.write('\n! Persistence dataset %s with fermi model:  norm %6.2f alpha %6.2f e_fermi %6.0f kT %6.0f\n' % (dataset,norm,alpha,e_fermi,kT)) 
	elif model_type==1:
		history.write('\n! Processing dataset %s with A gamma model' % dataset)
	elif model_type==2:
		history.write('\n! Processing dataset %s with time-variable fermi  model' % dataset)


	# Check whether there is anything to do

	if len(records)==1:
		string='subtract_persist: No persistence for this dataset.  No earlier observations within %4.1f hours\n' % (delta)
		history.write('%s\n' % string)
		history.write('! Persistence:  None\n')
		string='OK: subtract_persist: -- None'
		print string
		history.close()
		xstring='  %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f' % (0,0,0,0,0,0)
		per_list.update_summary(dataset,'Persist',xstring,append='no')
		return string



	# Persistence is calculated for the middle of the interval in which the expousre
	# was taking place

	[t1,t2]=get_times(science_record[0])
	tscience=0.5*(t1+t2)  

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

	xynorm=read_parameter(parameter_file,'xynorm')
	if xynorm!='':
		xynorm=locate_file(xynorm)
		xcorr=get_image(xynorm,1)
		if len(xcorr)==0:
			xynorm=''  # This is an error because we were unable to find the file
	else:
		string='Processing without spatially dependent correction'
		print string
		history.write('%s\n' % string )




	# This is the beginning of the loop for calculating the persistence model
	i=0
	while i<len(records)-1:
		record=records[i]
		print 'subtract: %30s %6.1f model_type %d' % (record[0],eval(record[11]),model_type)
		# dt is measured from the end of the stimulus image to the middle of the
		# science image
		[t1,t2]=get_times(record[0])
		dt=(tscience-t2)*86400

		cur_progid=record[2]
		words=record[3].split('.')
		cur_visit=words[0]
		cur_sci_fil=record[10]
		cur_sci_exp=eval(record[11])
		cur_sci_obj=record[14]
		scan=record[4]

		xfile=record[0]
		# Use the ima file, if file calusing persistence is a scan object
		if scan=='scan':
			xfile=xfile.replace('flt','ima')
			print 'Using ima file for ',record[0],xfile,scan


		x=get_image(xfile,1,'e',fileref=science_record[0])  # Convert this to electrons
		if len(x)==0:
			xstring='NOK: Problem with science extension of %s' % record[0]
			history.write('%s\n' % xstring)
			print xstring
			return xstring

		dq=get_image(xfile,3,fileref=science_record[0])     # Get the dq 
		if len(dq)==0:
			xstring = 'NOK: Problem with dq extension of %s' % record[0]
			history.write('%s\n' % xstring)
			print xstring
			# 110926 - ksl - modified to allow this to process the image even if there was no dq array
			# return xstring

		if model_type==0:
			model_persistence=calc_persist(x,dq,dt,norm,alpha,gamma,e_fermi,kT)
		elif model_type==1:
			# print 'Model type is 1'
			xfile=read_parameter(parameter_file,'a_gamma')
			model_persistence=make_persistence_image(x,cur_sci_exp,dt,xfile)
		elif model_type==2:
			# print 'Model type is 2'
			xfile=read_parameter(parameter_file,'fermi')
			model_persistence=make_persistence_image(x,cur_sci_exp,dt,xfile)
		else:
			print 'Error: subtract_persist: Unknown model type %d' % model_type
			return 'NOK'

		values=how_much(model_persistence)

		if i==0:
			persist=model_persistence
			stimulus=x   # This is an array which contains the maximum counts in a pixel
			xtimes=numpy.ones_like(persist)
			xtimes=xtimes*dt # This is an array containing the delta time at which the stimulus occured
		else:
			xpersist=model_persistence
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
		#130909 Removed print statement as unnneessary.  The string still goes to the history file
		# print string
		history.write('%s\n' % string)

		# Summarize the Stiumulus printing out the filename, prog_id, visit_name, target, dt and the number of pixels above saturation of 70000
		scan='No'
		if record[4]=='scan':
			scan='Yes'

		stimulus_summary='Stimulus: %40s %10s %10s %20s %8.0f %3d %3s %6s\n' % (record[0],cur_progid,cur_visit,cur_sci_obj,dt,xvalues[1],record[9],scan)
		history.write('! %s\n' % stimulus_summary)

		if i==last_external:
			ext_values=how_much(persist)
			ext_persist=numpy.copy(persist)

		i=i+1
	
	# This is the end of the loop where the persistence model is calculated

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


	# Now write out all of the new fits files

	# First find out the units of the science image

	units=get_keyword(science_record[0],1,'bunit')
	# print 'test',units
	if units[0]=='COUNTS/S':
		print 'Reducing model to match units for dataset %s to match %s ' % (dataset,units)
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
	# 140606 - Fix added to put stimulus file in the correct units.
	put_keyword(stimulus_file,1,'bunit','ELECTRONS')
	rewrite_fits(xname[0],time_file,1,xtimes)
	if len(ext_persist)>0:
		rewrite_fits(xname[0],ext_persist_file,1,ext_persist)

	# This completes the section which writes out all of the fits files
		
	# Get statistics on the images and make the 4 panel plot 

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

	# Finished plots

	# Finished everything so wrap it up
	history.write('# Finished Persistence processing of file %s\n' % science_record[1])
	history.close()

	# Upadete the summary file
	string='%20s %20s' % (emeasure,measure)
	per_list.update_summary(dataset,'Persist',string,fileroot,append='no')

	return string


def steer(argv):
	'''
	This is a steering routine for subtract persist so that options can be exercised from the 
	command line.  
	
	Notes:
	The switches are discussed in the documentation section at the tope

	100907	ksl	Added to begin to automate the subtraction process
	140617	ksl	Added new switches to account for the purely observational
			alpha gamma models.
	'''
	i=1
	dataset_list='none'

	model_type=1
	norm=0.3
	alpha=0.2
	gamma=0.8
	e_fermi=80000
	kT=20000
	fileroot='observations'
	parameter_file='persist.pf'
	words=[]
	mjd_after=0.0    # A amall number for mjd
	mjd_before=1.e6  # A large number for mjd

	ds9='no'
	local='no'

	while i<len(argv):
		if argv[i]=='-h':
			print __doc__
			return    
		elif argv[i]=='-model':
			i=i+1
			model_type=eval(argv[i])
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
		elif argv[i]=='-pf':
			i=i+1
			parameter_file=argv[i]
			if parameter_file=='none':
				parameter_file=''
		elif argv[i][0]=='-':
			print 'Error: subtract_persist.steer: Unknown switch ---  %s' % argv[i]
			return
		else:
			words.append(argv[i])
		i=i+1
	
	if dataset_list=='none': #  Then we are processing a single file
		dataset=words[0]
		do_dataset(dataset,model_type,norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,parameter_file)
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
					do_dataset(word[1],model_type,norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,parameter_file)
	else:
		f=open(dataset_list,'r')
		lines=f.readlines()
		f.close()

		for line in lines:
			x=line.strip()
			if len(x)>0 and x[0]!='#':
				mjd=float(x[6])
				if mjd_after <= mjd and mjd <= mjd_before:
					do_dataset(x,model_type,norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,parameter_file)
	
	return

	



	 

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		steer(sys.argv)   
	else:
		print 'subtract_persist.py  -h to get brief  help text'
	print 'OK done'
	sys.exit(0)
