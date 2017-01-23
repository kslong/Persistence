#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Read the summary file with various switches and construct
an html version of the summary file so that the results
can be easily inspected.


Command line usage :

	subtract_sum.py 

	-h		print this help infomation
	-prog_id xxxx	construct html file for prog id xxxx
	-all            construct for all datasets that have been processed
	-today          construct for datasets processed today
	-start          construct for datasets taken after a givn date in mjd or iso format
	-stop           construct for datasets taken before date in mjd or iso format

if subtract_sum.py is called with no options.  The program will crate a summary file for
the last hour's work

Description:  

	The routine reads the sum file to find the records that have been
	completed with the swithces file and makes a very simple html file
	to allow one to navigate the individual html files easily

Primary routines:

Notes:
									   
History:

110104	ksl	Coding begun
110201	ksl	help changed
110428	ksl	Rewrote using stand alone regions, not markup.py to allow for makeing a table
111115	ksl	Began adding routines to include a figure of how bad the persistence is
		on average

'''


# Added so this can run in the background without a $DISPLAY 
# environment variable.
import matplotlib
matplotlib.use('Agg')

import sys
import date
import os
import per_list
import numpy
import xhtml
import pylab

def link2file(link,word=''):
	'''
	Make the html to include a link
	'''
	if word == '':
		word=link
	
	string='<a href="file:%s" > %s </a>'  % (link,word) 
	return string


def make_html(lines,filename='observations.html'):
	'''

	This routine makes the summary html file.  
	
	110428 	ksl	This routine was writtend because markup.py did not seem to be to handle
			a simple table.

	'''

	header='''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
	<html lang="en">
	<head>

	<meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
	  <title>Summary Evaluation Page of Persistence</title>
	  </head>
	  <body>
	  <p>This page contains links to the individual html files for each
	  dataset</p>
	  <hr size="3" width="100%">
	  <table border="1" cellpadding="2" cellspacing="2" width="100%">
	  '''

	string=header

	for line in lines:
		row='<tr>\n'
		for word in line:
			row=row+'<td> %s </td>\n' % word
		row=row+'<tr>\n'
		string=string+row

	trailer='''
	    </tbody>
	  </table>
	  <hr size="3" width="100%">
	  </body>
	  </html>
	  '''

	string=string+trailer


	g=per_list.open_file(filename)


	g.write('%s\n' % string)
	return

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
		print("The file %s does not exist" % filename)
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

	# 130103 - Added to fix bug where start time is provided but no end time
	if mjd_start>0 and mjd_stop==0:
		mjd_stop=mjd_start+1e6  # A large number

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

def worst(lines,records,nmax=500):
	'''
	Find the datasets with the most persitence and write information regarding these
	to the file subtract_eval_worst.txt.  

	The selection is made on the fraction of pixels with persistence EXT2.

	111110	ksl	Modified to eliminate certain programs that are ringers, like ksl's persistence testing programs
	111116  ksl	Moved censoring to the main routine, and changed inputs to reflect this
	'''


	xlines=lines

	xx=[]
	for one in xlines:
		xx.append(eval(one[7]))
	xx=numpy.array(xx)

	order=numpy.argsort(xx)
	
	n=0
	i=len(order)-1
	g=open('subtract_eval_worst.txt','w')
	table=[]
	while i>=0:
		# line=per_list.read_ordered_list_one(dataset=xlines[order[i]][0])
		record=records[i]
		# print record
		string1= '%50s %10s %8s %8s %10s %5s %10s %20s %20s ' % (record[0],record[1],record[2],record[3],record[9],record[10],record[11],record[14],record[16])
		one=xlines[order[i]]
		# print records[order[i]]
		string2= '%8s %8s %8s %8s %8s %8s %30s' % (one[6],one[7],one[8],one[9],one[10],one[11],one[12])
		string=string1+string2
		# print string
		g.write('%s\n' % string)
		n=n+1
		if n>nmax:
			break
		tline=[]
		tline.append(html.link(record[1],one[12]))
		tline.append(record[2])
		tline.append(record[3])
		tline.append(record[9])
		tline.append(record[10])
		tline.append(record[11])
		tline.append(record[14])
		tline.append(record[16])
		tline.append(one[6])
		tline.append(one[7])
		tline.append(one[8])
		tline.append(one[9])
		tline.append(one[10])
		tline.append(one[11])
		table.append(tline)


		i=i-1
	g.close()

	string=html.begin('Worst Affected by Persistence')
	string=string+html.table(table)
	string=string+html.end()

	g=open('Worst.html','w')
	g.write(string)
	g.close()


	return

def make_fig(lines):
	'''
	Create a figure with showing how common various levels of persistence are
	'''
	values=[]
	for one in lines:
		values.append([one[6:12]])
	values=numpy.array(values)
	values=numpy.asfarray(values)  # convert everything to floats
	values=numpy.transpose(values)

	label=['Ext >  0.1 e/s','Ext >0.03 e/s','Ext >0.01 e/s','Int >  0.1 e/s','Int >0.03 e/s','Int >0.01 e/s']
	
	# print values[0]
	# print len(values)
	i=0

	pylab.figure(1,(6,6))
	pylab.clf()
	while i<len(values):
		y=numpy.ravel(values[i])
		# print y
		z=numpy.sort(y)
		z=numpy.flipud(z)
		x=numpy.linspace(0.,100.,num=len(z))
		# print x
		# print z
		pylab.plot(x,z,label=label[i])
		pylab.xlabel('% of images')
		pylab.ylabel('% of pixels')
		pylab.legend(loc='upper right')
		pylab.axis([0,100,0,10])
		i=i+1

	pylab.draw()
	pylab.savefig('persist_stats.png')


	# Just plot external persistence
	pylab.figure(2,(6,6))
	pylab.clf()
	i=0
	while i<3:
		y=numpy.ravel(values[i])
		# print y
		z=numpy.sort(y)
		z=numpy.flipud(z)
		x=numpy.linspace(0.,100.,num=len(z))
		# print x
		# print z
		pylab.plot(x,z,label=label[i])
		pylab.xlabel('% of images')
		pylab.ylabel('% of pixels')
		pylab.legend(loc='upper right')
		pylab.axis([0,25,0,10])
		i=i+1
	pylab.draw()
	pylab.savefig('persist_ext_stats.png')


	return


def doit(fileroot='observations',status='Complete',prog_id=0,mjd_start=0,mjd_stop=0,proc_start=0,proc_stop=0,censor='yes'):
	'''
	This is the main routine for the subtract_asum program

	It calls the routine to read the main parts of the summary file and then reformats the results
	so a table can be made.  It then calls the routine that makes the html file, and another routine
	that attempts to locate examples of the worst persistence, and another to make a figure.

	If censor=='yes', then certain calibration programs and all Tungsten lamp exposures are
	excluded from analysis

	
	'''
	# print 'start,stop',mjd_start,mjd_stop
	# lines=read_sum_file(fileroot,status,prog_id)
	lines=read_sum_file(fileroot,status,prog_id,mjd_start,mjd_stop,proc_start,proc_stop)
	xlines=[]
	for line in lines:
		if line[5].count('Compl'):
			xlines.append(line)
	lines=xlines

	print('Records to process (before censoring):',len(lines))

	if len(lines)==0:
		print('There is no lines to process, so returning')
		return


	# Now get the corresponding lines in the observations.ls file
	records=[]
	for line in lines:
		one=per_list.read_ordered_list_one(dataset=line[0])
		records.append(one)
	
	# So now we have two parallel lists one containing the summary file and one containing the .ls file

	xlines=[]
	title=['Rootname','Prog Id','Date_Proc','Time_Proc','Status','Ext1','Ext2','Ext3','Int1','Int2','Int3']
	xlines.append(title)
	for line in lines:
		print(line)
		xline=[link2file(line[12],line[0])]
		xline.append(line[1])
		xline.append(line[3])
		xline.append(line[4])
		xline.append(line[5])
		xline=xline+line[6:12]
		xlines.append(xline)




	# print 'How many',len(xlines)
	make_html(xlines)


	# Now censor the list to produce the worst

	ignore=[11423,11915,11930,12351,12338,11931,11696,12283]
	if censor=='yes':
		xlines=[]
		xrecords=[]
		i=0
		while i<len(lines):
			ok='yes'
			prog_id=lines[i][1]
			for value in ignore:
				# print value, prog_id
				if value==prog_id:
					ok='no'
					break
			if records[i][14]=='TUNGSTEN':
				ok='no'
			if ok=='yes':
				xlines.append(lines[i])
				xrecords.append(records[i])
			i=i+1
		lines=xlines
		records=xrecords
	

	



		
	print('Records to process (after censoring) :',len(lines))
	# print lines[0]
	# print lines[len(lines)-1]

	worst(lines,records)
	make_fig(lines)
	return





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
			print(__doc__) 
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
				print('Start',z,mjd_start)
		elif argv[i]=='-stop':
			i=i+1
			z=argv[i]
			try:
				mjd_stop=float(z)
			except	ValueError:
				mjd_stop=date.iso2mjd(z)
				print('Stop',z,mjd_stop)
		elif argv[i]=='-today':
			# today=date.get_gmt('%Y-%m-%d')
			today=date.get_gmt()
			today=date.parse_iso(today)
			proc_start=today-86400
			proc_stop=today
		else:
			print('Unknown switch %s' % argv[i])
			return

		i=i+1
	
	# parse_iso returns a unix time in seconds since Epoch0
	
	if len(argv)==1:
		today=date.get_gmt()
		# print 'now',today
		today=date.parse_iso(today)
		proc_start=today-3600
		proc_stop=today
		print('Creating summary file for last records processes in last hour')
		print(proc_start, proc_stop)



	doit(fileroot,status,prog_id,mjd_start,mjd_stop,proc_start,proc_stop)
	return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		steer(sys.argv)   
	else:
		# sys.argv.append('-h')
		steer(sys.argv)   
