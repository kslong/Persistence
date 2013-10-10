#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine is intended to make a correction file for a uniform persistence model.
It assumes that the only thing that is different between pixels is the amplitude
of persistence.  The idea is that the correction file is 1 + a small correction file


Command line usage (if any):

	usage: make_corr.py filename [outroot]

	where filename is the list of flt files that should be used to make the correction
	and outroot is an optional rootname for output fits files

Description:  

	Prior to running this program one should run_run persist on a set of darks
	which have already been processed, normally without having included a persistence
	correction file

Primary routines:

	doit is the main routine 

Notes:
									   
History:

110601 ksl Coding begun
120606	ksl	Revised to account for current way in which the persist directories are
		set up, and to allow it to accept a list of flt files.  Also made it 
		possible to determine the root part of the output.

'''

import sys
import os
import per_fits
import numpy
from scipy import signal   # Used to smooth the image

def read_file(filename):
	'''
	Read a file and split it into words, eliminating comments
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
		z=xlines[i].split()
		if len(z)>0:
			if z[0][0]!='#':
				lines=lines+[z]
		i=i+1
	return lines


def doit(filename='files.ls',outroot='persist'):
	'''
	Create an outuput file which contains the correction flat for
	persistence.

	The routine accepts a list of flt files.  It looks for persistence
	models (which should have been constructed without a spatially
	dependend correction, and constructs the correction flat from
	these

	outroot is the root name for all of the files produced

	120606	ksl	Rewrote sections to take account of new standard
			directory structure for the persist programs
	'''
	lines=read_file(filename)
	# print 'what',len(lines)
	models=[]
	residuals=[]
	for one in lines:
		# print one[0]
		z=per_fits.parse_fitsname(one[0])
		words=z[0].split('/') 
		basename=words[len(words)-1]
		# Get path to persist directory
		i=0
		name=''
		while i<len(words)-1:
			name=name+words[i]+'/'
			i=i+1
		name=name+'Persist/'
		xmod=basename.replace('_flt.fits','_persist.fits')
		xres=basename.replace('_flt.fits','_flt_cor.fits')
		model=name+xmod
		residual=name+xres
		# print z,name,basename
		# print 'OK',model,residual
		models.append(model)
		residuals.append(residual)
	# First check that the files we need exist

	icheck=0
	for one in models:
		if os.path.exists(one)==False:
			print 'Error: %s does not appear to exist' % one
			icheck=icheck+1
	for one in residuals:
		if os.path.exists(one)==False:
			print 'Error: %s does not appear to exist' % one
			icheck=icheck+1
	if icheck==0:
		print 'All of the necessary files appear to exist'
	else:
		print '%d files are missing.  Quitting' % icheck
		return
	


	# Now construct the correction flat
	i=0
	while i<len(models):
		model=models[i]
		residual=residuals[i]
		print 'Working on:',model,residual
	
		xmod=per_fits.get_image(model,1)
		xresid=per_fits.get_image(residual,1)

		med=numpy.median(xresid)
		mmed=numpy.median(xmod)
		print 'Model',mmed,numpy.std(xmod),numpy.median(xresid),numpy.std(xresid)
		xresid=xresid-med
		print 'Model',mmed,numpy.std(xmod),numpy.median(xresid),numpy.std(xresid)
		r=xresid/xmod
		print 'Result', numpy.median(r),numpy.std(r)
		if i==0:
			template=model
			rr=mmed*r
			norm=mmed
		else:
			rr=rr+mmed*r
			norm=norm+mmed
		i=i+1
	
	corr=1.+ rr/ norm
	print 'Final',numpy.median(corr),numpy.std(corr)

	# Now write the output files
	per_fits.rewrite_fits(template,outroot+'_corr.fits',1,corr)
	xcorr=signal.medfilt2d(corr,5)
	per_fits.rewrite_fits(template,outroot+'_smooth.fits',1,xcorr)
	diff=corr-xcorr
	per_fits.rewrite_fits(template,outroot+'_diff.fits',1,diff)



# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)==2:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	elif len(sys.argv)==3:
		doit(sys.argv[1],sys.argv[2])
	else:
		print __doc__

