#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

These are utilities supporting reading and writing fits files and
headers.  It is part of the suite of python rotines being used
to analyze the amount of persistence in HST WFC3/IR images

Except for testing purposes this routine is not intended to be
called from the command line


Description:  

Primary routines:

Notes:
									   
History:

100513 ksl Coding begun
110811 ksl Switched to standardized way to set file permissions
121220 ksl Removed pyraf dependencies after adding a routine to
	   to get keywords form files
130909 ksl Standardized printed error messages

'''

import sys
import os

# import pyraf
import pyfits
import pylab
import numpy
import string
import permissions



def parse_fitsname(name,ext=1,force_ext='no'):
	'''
	Parse a fits name, which may contain an extension, and return
	the name of the file, the extension only, and the name with the extension
	as three elements of a tuple 

	e.,g

	('iby03apjq_flt.fits', 2, 'iby03apjq_flt.fits[2]')

	
	If the name does not contain an extension OR
	if force_extension is anythin gbut no 
	set the extension to ext


	Notes

	This is intended to allow one to allow one to move easily between
	iraf and pyfits

	100103	ksl	Coded
	100513	ksl	Mofified so that a new extension could be forced, even if the filename already
			contains an extension.  This should simplify reading of extensions
	'''

	if len(name)==0:
		print 'Error: parse_fitsname: name had length 0, nothing to parse'
		return []

	s=string.replace(name,'[',' ')
	s=string.replace(s   ,']','')
	s=s.strip()
	s=s.split()

	if len(s)==1:
		xname=s[0]
		xext=ext
	elif force_ext!='no':
		xname=s[0]
		xext=ext
	else:         
		xname=s[0]
		try:
			xext=eval(s[1])
		except NameError:
			print 'Error: parse_fitsname: Could not parse extension (%s), using %d' % s[1],ext
			xext=ext
	
	return xname,xext,'%s[%d]' % (xname,xext)


def get_ext_type(filename,exten=1):
	'''
	Open and image and find out what type it is

	101203	Modified so filename is parsed in standard fashion
	'''

	type='UNKNOWN'

	name=parse_fitsname(filename,exten)

	try:
		z=pyfits.open(name[0])
	except IOError:
		print 'Error: get_ext_type: file %s not found' % filename
		return 'UNKNOWN'

	try:
		one_ext=z[name[1]]
		type=one_ext.header['EXTNAME']
	except IndexError:
		print 'Error: get_ext_type: %d exceeds the number of extensions in file %s' % (exten,filename)
	except KeyError:
		print 'Error: get_ext_type: EXTNAME is not found in header extension %d in file %s' % (exten,filename)

	z.close()
	return type


def get_image_pixel_info(filename,exten=1):
	'''
	Get the information needed to map one image onto another in
	detector pixel space

	The routine looks for LTV1 and LTV2 which contain the offset of
	the first pixel in the array to the first physical pixel of the
	detector.  The sense of the mapping is that the first (1,1) in the
	detector is located at pixel (1,1) of the image + (LTV1, LTV2).

	For a subarray in the middle of the detector, both LTV1 and LTV2
	will be negative

	The routine returns
	
	[rows,cols,offset_y,offset_x]

	from the image

	History

	101203	Modified so filename is parsed in standard fashion

	'''

	name=parse_fitsname(filename,exten)

	try:
		z=pyfits.open(name[0])
	except IOError:
		print 'Error: get_pixel_info: file %s not found' % filename
		return 'UNKNOWN'

	try:
		one_ext=z[name[1]]
		offset_x=one_ext.header['LTV1']
		offset_y=one_ext.header['LTV2']
		data=one_ext.data
		# print 'offsets',offset_x,offset_y
		shape=numpy.shape(data)
		rows=shape[0]
		cols=shape[1]

	except IndexError:
		print 'Error: get_image_pixel_info: %d exceeds the number of extensions in file %s' % (exten,filename)
	except KeyError:
		print 'Error: get_image_pixel_info: EXTNAME is not found in header extension %d in file %s' % (exten,filename)

	z.close()
	return [rows,cols,offset_y,offset_x]


def get_keyword(filename,exten,keywords='bunit',default='Unknown'):
	'''
	Get a one or more keywords from a fits file and extension

	Filename and extension are required.  (If the filename 
	contains and extension it will be ignored)
	
	keywords is a string which contains the keywords, separated 
	by commas or spaces.  

	The routine returns the results as a list.

	If a keyword is not found then the default value is inserted
	into that position in the list

	If the file is not found then [default] is returned

	Notes:

	Unlike hsel which returns everything as a string, py_fits
	returns whatever the value was intended to be, and as
	a result one can retrieve a list with both numbers an
	strings

	Unlike hsel, which for multi extension files looks both
	in the named extension and in extension 0 for a keyword
	one needs to explicitly carry out this procedure in pyfits

	121220 ksl Added
	'''

	# Next line means that you must always provide the exentsion
	name=parse_fitsname(filename,exten,'yes')

	answer=default

	try:
		z=pyfits.open(name[0])
	except IOError:
		print 'Error: get_keyword: file %s not found' % filename
		return default

	try:
		zero_ext=z[0]
		one_ext=z[name[1]]
	except IndexError:
		print 'Error: get_keyword: %d exceeds the number of extensions in file %s' % (exten,filename)
		return default


	# Create a list from the inmput string
	keywords=keywords.replace(',',' ')
	words=keywords.split()

	# Populate the output list, looking both in the desired extension and in extension 0
	answer=[]
	for word in words:
		try:
			answer.append(one_ext.header[word])
		except KeyError:
			try:
				answer.append(zero_ext.header[word])
			except KeyError:
				print 'Error: get_keyword: %s is not found in header extension %d in file %s' % (word,exten,filename)
				answer.append(default)

	z.close()


	return answer




def get_image(filename='./11740/ib1cc5gwq_flt.fits',exten=1,rescale='no',fill=0,fileref='none',section=[]):
	'''

	Read a specific extension of a fits file and rescale to electrons or counts depending
	on whether rescale is 'no','e','e/s'  

	If rescale contains the word "flush" then the routine assumes the mimimum exposure time is
	2.9 seconds.  This is specifically included for persistence reasons since the IR obsrvtions
	are always preceded by a flush, and the camera has not shutter.  Normally one would use
	something like e_flush for this option

	Optionally map into another image at the same position as physical pixels
	in the original image file.   This is to allow for subarrays in IR observations.

	A numpy array is returned that is the size of fileref and has values from the fits file designated by
	filename.  The portion of the returned array that does not correspond to the subarray will have a fill
	value designated by fill.

	A section of the original array if section is != [].  Note that this does not work when you want
	to fill out the array

	An empty numpy array is returned if the routine fails.

	Notes

	The purpose of this routine is to allow one to track persistence from sub_arrays to full frame images,
	and vice versus.  

	History

	101123	ksl	Began work
	101202	ksl	Verified that the routine appears to work correctly even in situations where the subimage
			is offset.  
	101203	ksl	Modified call so that is very similar to get_image_ext.  It's possible 
			get_image_ext should be eliminated
	120613	ksl	Add parameters to allow one to get a subsection, specified as a list
	'''

	# Do the simple case first where we are not concerned that we need to map the pixels of one image
	# into the pixels of another

	if fileref=='none':
		data=get_image_ext(filename,exten,rescale)
		if len(data)==0:
			print 'Error: get_image: returning empty array for %s extension %d' % (filename,exten)
		elif len(section)==4:
			xmin=section[0]
			xmax=section[1]+1
			ymin=section[2]
			ymax=section[3]+1
			data=data[xmin:xmax,ymin:ymax]
		return data

	# Do the harder case, where we might have subarrays that affect where physical pixels appear in
	# different images

	source_sizes=get_image_pixel_info(filename,exten)
	# print 'Source sizes',source_sizes

	ref_sizes=get_image_pixel_info(fileref,exten)
	# print 'ref sizes',ref_sizes

	source_data=get_image_ext(filename,exten,rescale)
	if len(source_data)==0:
		print 'Error: get_image: No source data for file %s' % filename
		return []

	ref_data=get_image_ext(fileref,exten)
	if len(ref_data)==0:
		print 'Error: get_image: no refdata for file %s' % filename
		return []

	# print 'shapes',numpy.shape(source_data),numpy.shape(ref_data)

	iymin=ref_sizes[2]-source_sizes[2]
	iymax=iymin+source_sizes[0]

	ixmin=ref_sizes[3]-source_sizes[3]
	ixmax=ixmin+source_sizes[1]

	# print ixmin,ixmax,iymin,iymax

	# 111109 Changed to from zeros_like to zeros to avoid numpy problem
	z=numpy.zeros(ref_data.shape,dtype=ref_data.dtype)

	z=z+fill


	# In situations that we care about the source image should be contained in the reference image or vice versus
	# At this point I have not tried to deal with the other logical possibilities

	# print 'test',len(z),iymin,iymax,ixmin,ixmax
	if iymin >= 0 and ixmin>=0:
		z[iymin:iymax,ixmin:ixmax]=source_data
	else:
		iymin=-iymin
		iymax=iymin+ref_sizes[0]
		ixmin=-ixmin
		ixmax=ixmin+ref_sizes[1]

		z=source_data[iymin:iymax,ixmin:ixmax]

		# print iymin,iymax,ixmin,ixmax

	# Next line for testing purposes only
	# rewrite_fits(oldname=filename,newname='test.fits',ext=1,data=z,clobber='yes')

	return z


def get_image_ext(filename,exten=1,rescale='no'):
	'''
	Read a specific extension of a fits file and rescale to electrons or counts depending
	on whether rescale is 'no','e','e/s'  

	The options for rescale are:

		no	No rescaling is done
		e/s	Units are returned as e/s, even if the original image is in e
		e	Units are returned as e, even if originally in some other unit
		ef	Units are returned as e, with the exposure time altered to reflect
			the likelihood that the maximum exposure time was in the flux, which
			is assumed to be 2.9 sec.  This is specifically included for persistence 
			reasons since the IR obsrvtions are always preceded by a flush, 
			and the camera has not shutter.

	Note that this routine should be deprecated in calls from other routines.  
	get_image is more general and may ultimately replace this.

	101213	ksl	Added option to include the flush time in the exposure.  
	121220	ksl	Removed pyraf.iraf dependency
	'''

	# Note 'yes' means that even if the filename includes and extension the name that will 
	# be returned will contain the extension.
	xname=parse_fitsname(filename,exten,'yes')

	# print 'get_image_ext',xname,exten,rescale


	# print 'xname',xname
	try:
		f=pyfits.open(xname[0])
		data=f[xname[1]].data
		f.close()
	except IOError:
		print 'Error: per_fits.get_image_ext: %s does not appear to exist' % filename
		return []

	if data==None:
		print 'Error: per_fits.get_image.ext: %s exists, but returns None for ext.%d' % (filename,exten)
		return []

	if rescale=='no':
		return data
	
	# Code below is old 121220  ksl - remove when satisfied
	# Determine whether the image has been exposure corrected and the exposure time  
	# xxxx=pyraf.iraf.hselect(xname[2],'exptime,unitcorr,bunit','yes',Stdout=1)
	# try:
	# 	xxxx=xxxx[0].split('\t')
	# except IndexError:
	# 	print 'Error: per_fits.get_image_ext: hsel did not return anything to split'
	# 	print 'xxxx was',xxxx
	# 	return []

	# print 'try',xname
	xxxx=get_keyword(xname[0],exten=xname[1],keywords='exptime,unitcorr,bunit',default='Unknown')
	# print 'new',xxxx


	if xxxx[2].lower() == 'electrons/s':
		units='e/s'
	elif xxxx[2].lower()=='counts':
		units='counts'
	elif xxxx[2].lower()=='counts/s':
		units='counts/s'
	else:
		print 'This file yielded ',xxxx,' which suggests all options are not accounted for, assuming electrons'
		units='e'

	
	# Code below is old 121220  ksl - remove when satisfied
	# texp=eval(xxxx[0])
	texp=xxxx[0]

	if units=='counts' and rescale != 'none':
		print 'Converting file %s in counts to electrons' % filename
		data=data*2.4  # Convert to (by the average gain correction
		units='e'
	if units=='counts/s' and rescale != 'none':
		print 'Converting file %s in counts/s to electrons/s' % filename
		data=data*2.4  # Convert to (by the average gain correction
		units='e/s'

	if rescale=='e' and units=='e/s':
		data=data*texp   # Convert to electrons
	elif rescale=='e/s' and units=='e':
		data=data/texp
	elif rescale=='ef' and units=='e/s':
		if texp<2.9 and rescale=='ef':
			print 'Assuming the minimum exposure is 2.9 sec'
			texp=2.9
		data=data*texp   # Convert to electrons
	elif rescale=='ef' and units=='e':
		if texp<2.9 and rescale=='ef':
			data=data*(2.9/texp)
	elif rescale!='none' and rescale != units:
		print 'Not quite sure how this image was to be rescaled ',rescale,xxxx
		print 'Assuming no rescaling was desired'

	return data

dummy=numpy.array([0])

def rewrite_fits(oldname='old.fits',newname='new.fits',ext=1,data=dummy,clobber='yes'):
	'''
	Rewrite a multi extension fits file using the header from
	the old file and updating one of the extension with the values
	contained in data

	If the oldname and newname are identical the file will simply
	be updated, something which is dangerous so be careful.

	100603
	'''

	if oldname==newname:
		x=pyfits.open(oldname,'update')
		x[ext].data=data
		x.close()
		return

	if os.path.exists(newname) and clobber=='yes':
		os.remove(newname)
	x=pyfits.open(oldname)
	x[ext].data=data
	x.writeto(newname)
	x.close()
	permissions.set(newname)
	return



def get_times(filename):
        '''
	Return the beginning and end times for an image as [t1,t2]

	If the file is missing, an empty list is returned

	100513	ksl	This version evaluates the times to floats
	111107	ksl	Fixed so that it would not fail if the file 
			was missing
	121220	ksl	Removed pyraf.iraf dependency
	'''

	xname=parse_fitsname(filename)

	xxxx=get_keyword(xname[0],exten=xname[1],keywords='expstart,expend',default='Unknown')
	if xxxx=='Unknown':
		return []
	return xxxx

	# Old 121220
	# xxxx=pyraf.iraf.hselect(xname[2],'expstart,expend','yes',Stdout=1)
	# print 'get_times', filename,xxxx
	# try:
	# 	xxxx=xxxx[0].split('\t')
	# 	return [eval(xxxx[0]),eval(xxxx[1])]
	# except IndexError:
	# 	print 'get_times: hsel did not work for file %s ' % xname[2]
	# 	return []



# This begins a section on following the history of individual pixels in the raw images.

def get_ext_info(filename,type='SCI'):
	'''
	get the extensions in a file which are of a certain type.  For SCI extensions
	also return the samptime, if that is possible

	A typical return looks like  this for type SCI

	[[1, 352.93951399999997], [6, 327.93899499999998], ... 

	for othe types the SAMPTIME will normally be missing, that is the records
	will simply be
	[[2],[7] ....


	If a file is not found, and empty list will be returned.

	Notes:

	samptime is only present for SCI extensions in WFC3 data


	History

	100308	ksl	Written to begin to deal with raw data files.  Not obvious how
			generally useful it is.
	100525	ksl	 This version opens and then closes the file
	100831	ksl	Modified to work for files other than SCI extensions
			but it is not clear how useful it is in this case.  It 
			might be better to have a routine which gets all of the
			extensions of a certain type and another routine which
			gets a keyword from each of the extensions
	'''
	try:
		z=pyfits.open(filename)
	except IOError:
		print 'Error: get_ext_info: file %s not found' % filename
		return []

	i=1
	ext=[]
	while i<len(z):
		one_ext=z[i]

		line=[]
		try:
			if one_ext.header['EXTNAME']==type:
				line.append(i)
				if type=='SCI':
					try:
						samptime=one_ext.header['SAMPTIME']
						line.append(samptime)
					except KeyError:
						print 'Error: get_ext_info: SAMPTIME keyword not found for ext %d' % i
				ext.append(line)
		except KeyError:
			print 'EXTNAME keyword not found for ext %d' % i



		i=i+1
	z.close()
	return ext

def doit(filename):
	'''
	This is intended to be used as test routine

	121220	ksl	Removed most of routine as part of eliminating pyraf.iraf
			dependencies.  Just a stub at this point
	'''

	xx=get_ext_info(filename)


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		doit((sys.argv[1]))
	else:
		print 'usage: per_fits.py  filename'

