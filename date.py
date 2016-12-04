#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

	Routines to convert times from ISO times to MJD, and back      

	Also includes a routine to get the current GMT time in a reasonable
	format


Description:  

	These routines make use of the time and calendar modules of standard python


Notes:
	The conversion reoutines appear to be accurate to about a second for dates
	near the current date.  
	
	The iso2mjd routine does not have unlimited range, but they
	do cover the period from 1970 to at least 2020.

	The mjd2iso routine gives the correct zero time, i.e. 1858-11-17 produces
	and mjd date of zero
									   
History:

101015	ksl	Coding begun
101220	ksl	Added get_gmt to these routines, not so much because it is associated
		with the other routines specifically but simply because it seemed a 
		good idea to put all of the routines that use the time and calendar
		modules in one place

'''

import time
import calendar


# This utility was added here simply because it seemed useful to put all the time related routines in one place

def get_gmt(format="%Y-%m-%d %H:%M:%S"):
	'''
	Return the current gmt time

	101222	ksl	Changed format so returns ISO time   
	'''
	x=time.strftime(format, time.gmtime())
	return x

mjd1970=40587.  # MJD for 1970.000

def parse_iso(iso='2009-10-21 06:15:09'):
	'''
	Parse an iso-formatted  date and return the time in seconds since Epoch0

	Note the formatting requirments for this are strict.  One can enter just
	a date or a date and at time, but the components of the date must be of
	the form
		2009-10-21
	and the components of the time must be
		06:15:09
	
	One can enter just the date or the date and the time
	'''


	# print 'xxxzzz',iso
	# The formatting for this is quite strict
	try:
		z=time.strptime(iso,'%Y-%m-%d %H:%M:%S')
		# print 'parsed day and hours'
	except ValueError:
		z=time.strptime(iso,'%Y-%m-%d')
		# print 'parsed day only'

	# print 'parsed before timegm',z

	# This is correct to GMT time in seconds since Epoch 0
	# Note that time.mktime is the  local time.  This requires a time sstructue

	zz=calendar.timegm(z)
	# print 'parse_iso output of timegm ',zz

	return zz



def iso2mjd(iso='2009-10-23 04:00:32'):
	'''
	Convert a string containing the ISO time into MJD
	'''

	# Verify that 1970 0 0 is the Epogh
	time_zero=time.gmtime(0.0)  # Get epoch 0

	if time_zero[0] != 1970:
		print('Error iso2mjd - Epoch is not 1970')
	# MJD of 1970.0


	# convert iso time to gmt

	delta=parse_iso(iso)  # This produces a float that gives the time in seconds since EPOCH 0

	delta=delta/86400.  # Convert to days

	mjd=mjd1970+delta  # This is the mjd time in seconds


	# print 'mjd in days',mjd

	return mjd  # In days

def mjd2iso(t=5.512526926e4):
	'''
	Convert from MJD (days) to date in iso format

	ISO format is  2009-10-23 04:00:32
	'''

	# print 'input  in days %8.2f' % t

	x=t-mjd1970

	# print 'offset in days %8.2f' % x

	x=x*86400.

	# print 'offset in s', x

	x=time.gmtime(x)

	# print 'ouput of gmtime',x

	z=time.strftime('%Y-%m-%d %H:%M:%S',x)

	return  z

def test():
	'''
	This is a test of some of the conversion  routines.

	101018	ksl - On my Mac, the forward conversion had an error of 1 sec on the return
	'''

	iso='2009-10-21 06:27:44'

	print('Forward')

	mjd=iso2mjd(iso)

	# print 'mjd after converting',mjd

	xiso=mjd2iso(mjd)

	print('%s --> %6.2f --> %s' % (iso,mjd,xiso))

	print('Reverse')
	mjd=5.512526e4

	iso=mjd2iso(mjd)

	xmjd=iso2mjd(iso)

	print('%6.2f --> %s --> %6.2f' % (mjd,iso,xmjd))

