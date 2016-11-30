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
from astropy.io import fits as pyfits
from pyraf import iraf
import pylab


# DS9 utilities

def check_for_ds9():
	'''
	Check for the system software needed to use ds9  
	'''

	x=os.system('which ds9')
	if x:
		print 'ds9 is not in your path'
	y=os.system('which xpaset')
	if y:
		print 'xpa software is not in your path.  Modify your path or install from links on XPA in the Reference manual for ds9'
	
	if x or y:
		return 1
	return 0

	
def start_ds9(height=600,width=600,scale='histequ',reset='no'):
	'''
	Start ds9 if it does not exist wiht height, width and scale

	If ds9 already exists then nothing will happen, unless reset is set
	to yes, in which case ds9 will be reinitialized


	'''


	# First check if ds9 is already extant'
	iok=os.system('xpaget ds9 version &> tmp.txt')
	if iok:
		print 'starting  ds9'
		os.system('ds9 &')
		time.sleep(5)  # give system a chance to instandiate ds9
		os.system('xpaset -p ds9 scale %s'% scale)
		os.system('xpaset -p ds9 width %d' % width)
		os.system('xpaset -p ds9 height %d'% height)
	elif reset=='yes':
		os.system('xpaset -p ds9 scale %s'% scale)
		os.system('xpaset -p ds9 frame delete all')
		os.system('xpaset -p ds9 width %d' % width)
		os.system('xpaset -p ds9 height %d'% height)
	
	return
	


def write_text_regionfile(text,x,y,color='green',filename='ds9.reg'):
	'''
	Create a region file to put text specific position in physical coordinates

	'''

	header='''# Region file format: DS9 version 4.0
# Filename: ./dark/Visit47/ibcu47juq_flt.fits[SCI]
global color=%s font="helvetica 14 bold" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source
physical''' % (color)
	
	reg=open(filename,'w')
	reg.write('%s\n' % header)
	reg.write('# text(%s,%s) text={%s}\n' % (x,y,text))
	reg.close()



def LoadFrame(filename,frame=1,xmin=0,xmax=1000,scale='linear'):
	'''
	Load a frame into ds9 with fixed limits on a linear
	or histogram equilized scale.  If scale is something else
	then the scaling will be the inherited scaling in ds9
	'''

	start_ds9(800,800,'linear')


	string='xpaset -p ds9 frame %d' % frame   # Select the frame
	os.system(string)

	string='xpaset -p ds9 frame reset' # Reset the frame to its defaults
	os.system(string)

	string='xpaset -p ds9 file %s' % filename     # Load the frame 
	os.system(string)

	# Next lines allow one to get close to the scaling of interest
	if scale.count('lin'):
		string='xpaset -p ds9 scale linear'
		os.system(string)
	elif scale.count('hist'):
		# print 'Histogram equalized'
		string='xpaset -p ds9 scale histequ'
		os.system(string)
	

	string='xpaset -p ds9 cmap invert'
	os.system(string)

	string='xpaset -p ds9 scale limits %f  %f' % (xmin,xmax)
	os.system(string)

	os.system('xpaset -p ds9 zoom to fit')

	return






