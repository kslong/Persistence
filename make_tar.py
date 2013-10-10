#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Read the summary file with various switches and construct
tar files of the results

Command line usage :

	make_tar.py 

	-h		print this help infomation
	-prog_id xxxx	construct html file for prog id xxxx
	-all            construct for all datasets that have been processed
	-today          construct for datasets processed today
	-start          construct for datasets taken after a givn date in mjd or iso format
	-stop           construct for datasets taken before date in mjd or iso format

if subtract_sum.py is called with no options.  The program will create tar files for
the day's work

Description:  

	The routine reads the sum file to find the records that have been
	completed with the swithces file and makes a very simple html file
	to allow one to navigate the individual html files easily

Primary routines:

Notes:
									   
History:

111117	ksl	Adapted from subtract_sum
		on average

'''

import sys
import date
import os
import numpy
import subprocess
import time
import date

def read_sum_file(fileroot='observations',status='Complete',prog_id=0,mjd_start=0,mjd_stop=0,proc_start=0,proc_stop=0):
	'''
	Read the summary file and return selected sets of records in the observations.sum with the selection determined by

	prog_id		program_id (a single number)
	mjd_start       datasets taken after mjd_start
	mjd_stop	datasets taken before mjd_stop
	proc_start	datasets processed after proc_start
	proc_stop       datasets processed before proc_stop

	In general, if a value is 0, it is not used in selecting the data set

	111011	ksl	Modified the way the routine works so that to select something on the basis of status,
			one just needs to include the string that is in status.  This is to accommodate versioning.
			NOTE -- It is not obvious that one should check for status at all, since one usually wants
			to see what happened to a serios of files. What is more relevebant is the program ID or
			whether one attempted to process it in a certain period.


	'''

	filename=fileroot+'.sum'

	try:
		f=open(filename,'r')
		xlines=f.readlines()
		f.close()
	except IOError :
		print "The file %s does not exist" % filename
		return []   
	
	lines=[]

	# print 'OK0',len(xlines)
	
	i=0
	while i<len(xlines):
		z=xlines[i].strip()
		z=z.split()
		if len(z)>0:
			if z[0][0]!='#':
				# print z[5].count(status),z
				if z[5].count(status) or status=='All':
					lines.append(z)
		i=i+1
	
	# print 'OK1',len(lines)
	
	# At this point we have split the records into words and we have chosen according to status

	if prog_id != 0:
		xlines=lines
		lines=[]
		for z in xlines:
			if int(z[1])==int(prog_id):
				lines.append(z)
	# At this point our list has been reduced to those of the status we want and the prog_id we want

	# print 'OK2',len(lines)

	if mjd_stop > 0 and mjd_start <  mjd_stop:
		xlines=lines
		lines=[]
		for z in xlines:
			xtime=eval(z[2])
			if mjd_start <= xtime and xtime <= mjd_stop:
				lines.append(z)

	# Convert proc start and stop time to MJD
	# print 'OK3',len(lines)

	if proc_start!=0:
		try:
			proc_start=float(proc_start)
		except ValueError:
			proc_start=date.parse_iso(proc_start)

	if proc_stop!=0:
		try:
			proc_stop =float(proc_stop )
		except ValueError:
			proc_stop =date.parse_iso(proc_stop )

	if proc_stop != 0 or proc_start <  proc_stop:
		xlines=lines
		lines=[]
		for z in xlines:
			# print z
			# Construct the process time
			proc_time='%s %s' % (z[3],z[4])
			# print 'Starting',proc_time
			proc_time=date.parse_iso(proc_time)  # Convert to Mjd
			# print proc_time, 'xxx',proc_start,proc_stop
			if proc_start <= proc_time and proc_time<=proc_stop:
				lines.append(z)
	# print 'OK4',len(lines)

	return lines

def doit(fileroot='observations',status='Complete',prog_id=0,mjd_start=0,mjd_stop=0,proc_start=0,proc_stop=0,censor='yes'):
	'''
	This is the main routine for the subtract_asum program

	It calls the routine to read the main parts of the summary file and then reformats the results
	so a table can be made.  It then calls the routine that makes the html file, and another routine
	that attempts to locate examples of the worst persistence, and another to make a figure.

	If censor=='yes', then certain calibration programs and all Tungsten lamp exposures are
	excluded from analysis

	
	'''
	lines=read_sum_file(fileroot,status,prog_id,mjd_start,mjd_stop,proc_start,proc_stop)
	xstart=time.time()
	print 'xstart',time
        cur_time=date.get_gmt()


	# Get programs, visits, and locations

	location=[]
	for one in lines:
		file=one[12]
		words=file.split('/')
		try:
			x=file[0:file.rindex('/')]
			location.append(x)
		except ValueError:
			print 'Error: no / in %s\n' % file
	
	x=set(location)
	x=list(x)

	ntot=len(x)
	print 'There are %d directories to process' % ntot

	n=1
	for one in x:
		print 'Starting %s' % one
		g=open('TarEm','w')
		word=one.split('/')
		tar_dir='long_%s/%s_%s' % (word[1],word[2],word[3])
		src_dir=one
		string=        'rm -r %s\n' % tar_dir
		string=string +'mkdir %s\n' % tar_dir
		string=string +'ln %s/*persist.fits %s\n' % (src_dir,tar_dir)
		string=string +'ln %s/*cor.fits %s\n' % (src_dir,tar_dir)
		string=string +'tar czf long_%s/%s.%s.tar.gz %s\n' % (word[1],word[2],word[3],tar_dir)
		string=string +'rm -r %s\n' % tar_dir
		g.write('%s\n' % string)
		g.close()
		proc=subprocess.Popen('source TarEm',shell=True,stdout=subprocess.PIPE)
		x=proc.communicate()[0]
		dt=time.time()-xstart
		print dt
		print '# Completed dataset %d of %d. Elapsed time is %0.1f s (Ave %0.1f)' % (n,ntot,dt,dt/n)
		n=n+1
	print x


		

		





def steer(argv):
	'''
	Steering routine for this routine which is intended to make it easier
	to inspect the results of persistence subtraction.
	'''

	fileroot='observations'
	status='Complete'
	prog_id=0
	mjd_start=0
	mjd_stop=0
	proc_start=0
	proc_stop=0

	i=1
	while i<len(argv):
		if argv[i]=='-h':
			print __doc__ 
			return
		elif argv[i]=='-prog_id':
			i=i+1
			prog_id=int(argv[i])
		elif argv[i]=='-all':
			status='All'
		elif argv[i]=='-status':
			i=i+1
			status=argv[i]
		elif argv[i]=='-start':
			i=i+1
			z=argv[i]
			try:
				mjd_start=float(z)
			except	ValueError:
				mjd_start=date.iso2mjd(z)
				print 'Start',z,mjd_start
		elif argv[i]=='-stop':
			i=i+1
			z=argv[i]
			try:
				mjd_stop=float(z)
			except	ValueError:
				mjd_stop=date.iso2mjd(z)
				print 'Stop',z,mjd_stop
		elif argv[i]=='-today':
			# today=date.get_gmt('%Y-%m-%d')
			today=date.get_gmt()
			today=date.parse_iso(today)
			proc_start=today-86400
			proc_stop=today
		else:
			print 'Unknown switch %s' % argv[i]
			return

		i=i+1
	
	# parse_iso returns a unix time in seconds since Epoch0
	
	if len(argv)==1:
		today=date.get_gmt()
		# print 'now',today
		today=date.parse_iso(today)
		proc_start=today-3600
		proc_stop=today
		print 'Creating summary file for last records processes in last hour'
		print proc_start, proc_stop



	doit(fileroot,status,prog_id,mjd_start,mjd_stop,proc_start,proc_stop)
	return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		steer(sys.argv)   
	else:
		print 'subtract_sum.py -h to get help'
		steer(sys.argv)   
