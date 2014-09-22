#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine is intended to contain all of the additional routines necessary
to determine whether an observation is a scanned observation, and if so to
create a processed image that can be used as the input to the persistence
estimation routines


Command line usage (if any):

	usage: scan.py filename

Description:  

Primary routines:

Notes:
									   
History:

140225 ksl Coding begun

'''

import sys
import per_fits
from astropy.io import fits
import wfc3tools


def is_scan(filename='./Visit43/ic9t43j1q_flt.fits[1]'):
	'''
	Determine whether or not a file is associated with an observation that contains 
	a spatial scan

	Return:
		yes is a spatial scan
		no  if the spt file exists, but it is not a spatial san
		noinfo  if the spt file does not exist
	
	Notes:

	The routine looks for the spt file corresponding to the dataset and
	parses the header to find out if it is a scanned observation.
	
	130225  Coded and Debugged
	'''


	# First strip of the extension if any
	xname=per_fits.parse_fitsname(filename,0,'yes')
	xfile=xname[0]
	xfile=xfile.replace('flt','spt')
	xfile=xfile.replace('raw','spt')
	xfile=xfile.replace('ima','spt')
	# Now we should have the name of the file

	# print xfile

	try:
		x=fits.open(xfile)
	except IOError:
		return 'noinfo'

	# print x[0].header

	if x[0].header['SCAN_TYP'] == 'N':
		return 'no'

	return 'yes' 

def cal_scan(filename='./Visit43/ic9t43j1q_flt.fits[1]'):
	'''
	Recalibrate a wfc3 file changing switches to assure
	that only those parts needed for scanned observations are
	carried out
	'''

	# First strip of the extension if any
	xname=per_fits.parse_fitsname(filename,0,'yes')
	xfile=xname[0]
	xfile=xfile.replace('flt','ima')
	xfile=xfile.replace('raw','ima')
	# Now we should have the name of the file

	try:
		x=fits.open(xfile)
	except IOError:
		print 'Error: Could not open %s' % filename
		return 

	extime=x[1].header['EXPTIME']
	data=x[1].data
	data=data*exptime


def compare(filename='./Visit43/ic9t43j1q_flt.fits[1]'):
	'''
	Compare information in the flt files and the ima file

	notes:  This is just test routine
	'''
	xname=per_fits.parse_fitsname(filename,0,'yes')
	xfile=xname[0]
	ima_file=xfile
	ima_file=ima_file.replace('flt','ima')
	flt_file=xfile
	flt_file=flt_file.replace('ima','flt')

	# So now we have the name of the ima and the flt file no matter what

	ima_data=per_fits.get_image(ima_file,1,'e',fileref=flt_file)
	flt_data=per_fits.get_image(flt_file,1,'e')

	print ima_data.shape,flt_data.shape

	per_fits.rewrite_fits(flt_file,'x_flt.fits',1,flt_data)
	per_fits.rewrite_fits(flt_file,'x_ima.fits',1,ima_data)










# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		print 'usage: scan.py filename'
