#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is a routine to look at raw (and possibly ima) files on a
extension by extension basis.


Description:  

Primary routines:

Notes:
	At this point this is just a routine to carry out an imstat on
	a single file
									   
History:

100610 ksl Coding begun

'''


import sys
from per_fits import *
import pyraf

def doit(filename='ibel03e2q_raw.fits'):

	extensions=get_ext_info(filename)
	for ext in extensions:
		xfilename='%s[%d]' % (filename,ext[0])
		zzz=pyraf.iraf.imstatistics(xfilename,fields = "npix,mean,stddev,min,max,image", lower='INDEF', upper='INDEF',nclip=2, lsigma=3.0, usigma=3.0, binwidth=0.1,format='no',Stdout=1)
		print 'ext %2d time %5.1f' % (ext[0],ext[1]),zzz

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		print 'usage: per_raw.py  whatever.fits'

