#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

The purpose of this routine is to identify observations that 
would have been bad actors had an observation occurred after
them.  

It uses the observations.ls file as the file list to
explore


Basic ussage is as follows:

	bad_actor.py dataset           = Carry out a bad actor evaluation of this data set
	bad_actor.py -prog_id  program - Carry out a bad actor evaluation for this program
	bad_actor.py -all              - Carry out a bad actor evaluation for everything.  This will take 
					a while
	bad_actor.py -all [-start time1]  [-stop time2] - process the datasets
	        			that are in the .ls when the observations occured after time1 and/or before time2
		        		time1 and time2 are either MJD or ISO formated times, e.g '2009-11-12 15:13:24'


Description:  

	Three files are produced:
		bad_actor.txt     	Contains the results of analyze, a sorted list of the 
			programs from most egregious bad actor to least
		bad_actor_all.txt	Contains results for each file that was analyzed
		Problems.txt		Identifies files which may have been moved since
			the observations.ls file was created by per_list.py

Primary routines:

Notes:
									   
History:

111019 ksl Coding begun

'''

import sys
import numpy
import pylab
import date
import per_list
import per_fits

def do_dataset(fileroot='observations',dataset='ia21h2e9q'):
	'''
	Determine for a single dataset how much persistence this object is likely 
	to cause

	The routine returns a record for each dataset contining
		The total number of pixels
		The medium value
		Fraction of pixels with more than 2*saturation
		Fraction of pixels with more than saturation
		Fraction of pixels with more than half saturation

	
	111025	ksl	Covert so results are provided as persentages
	111031	ksl	Modified so that if the image is not read then the routine
			returns an list of length 0
	'''
	record=per_list.read_ordered_list_one(fileroot,dataset)

	if len(record)==0:
		print 'Could not find dataset %s' % dataset
		return

	# print 'xxxx',record
	x=per_fits.get_image(record[0],1,'e')  # Convert this to electrons
	if len(x)==0:
		print 'Skipping dataset %s because it was not read' % record[0]
		return []

	q=numpy.ravel(x)
	mpts=len(q)
	qq=numpy.sort(q)     # So now we have all sorts of things we can do with the mean image
	midpoint=int(0.5*len(qq))
	values=[mpts,qq[midpoint]]

	sat=70000

	i=len(qq)-1

	levels=[2*sat,sat,0.5*sat]  # These need to be in reverse order
	j=0
	while j<len(levels):
		saturation=levels[j]
		while i>0 and qq[i]>saturation:
			i=i-1
		npts=len(qq)-i
		values.append(float(npts)/mpts)
		j=j+1


	print dataset,values
	return values

def get_progs(records):
	'''
	Get a sorted list of the program ids from a
	standard set of records
	'''

	progs=[]
	for record in records:
		progs.append(record[2])
	if len(progs)<1:
		print 'Houston, there were no records for get_progs'
		return []
	progs.sort() # Sort the program ids in place
	progs=frozenset(progs) # Get the unique elments of the set
	progs=list(progs) # Turn it back into a list
	return progs


def analyze(records,results):
	'''
	Produce plots etc of the results

	Plot the fraction with less than a certain number of pexels
	'''
	results=numpy.array(results)

	progs=get_progs(records)
	print progs

	# The quickest way to proceed is going to be to process 
	# The reords only once, but this means we need to set up 
	# initalized lists for each row
	aaa=[]
	info=[]
	for prog in progs:
		aaa.append([0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) # nimages, average for each level, worst fo each level
		info.append(['UnkownPi','UnknownWorstTarget'])
	i=0
	while i<len(records):
		j=0
		while records[i][2] != progs[j]:
			j=j+1
		if aaa[j][0]==0:
			info[j][0]=records[i][16]

		aaa[j][0]=aaa[j][0]+1
		aaa[j][1]=aaa[j][1]+results[i][2]
		aaa[j][2]=aaa[j][2]+results[i][3]
		aaa[j][3]=aaa[j][3]+results[i][4]
		if aaa[j][4]<results[i][2]:
			aaa[j][4]=results[i][2]
		if aaa[j][5]<results[i][3]:
			aaa[j][5]=results[i][3]
			info[j][1]=records[i][14]
		if aaa[j][6]<results[i][4]:
			aaa[j][6]=results[i][4]
		i=i+1
	
	xsort=[]
	for one in aaa:
		one[1]=one[1]/one[0]
		one[2]=one[2]/one[0]
		one[3]=one[3]/one[0]
		xsort.append(one[2])  # This is the average fraction of pixels that exceed saturation
	for one in aaa:
		print one

	# Now find the worse bad actors programs. Sort on the average number of pixels that exceed
	# saturation

	xsort=numpy.array(xsort)
	xindex=numpy.argsort(xsort)
	xindex=numpy.flipud(xindex)

	g=open('bad_actor.txt','w')
	g.write('# Program     PI        WorstTarget    Nimages   2xSat     Sat   1/2XSat   Worst2x   WorstSat  Worst1/2Sat\n')
	for i in xindex:
		one=aaa[i]
		string='%10s %20s %20s %7d ' % (progs[i],info[i][0],info[i][1],one[0])
		j=1
		while j<len(one):
			string=string+' %6.3f ' % (one[j]*100.)
			j=j+1
		print string
		g.write('%s\n' % string)
		# print progs[i],aaa[i],info[i]
	g.close()




	# Now make some plots
	x=numpy.transpose(results)  # Transpose columns to rows for easier access

	sat=x[2]*100 
	npts=float(len(sat))
	number=numpy.linspace(0.,100.,npts)
	sat=numpy.sort(sat)
	sat=numpy.flipud(sat) # Reverse so it is the percentage with more than this saturation

	pylab.figure(1,[5,5])
	pylab.clf()
	print len(number),len(sat)
	pylab.plot(number,sat,'-',linewidth=3,label='%Pix > 2x    Saturation')
	sat=x[3]*100 
	sat=numpy.sort(sat)
	sat=numpy.flipud(sat) # Reverse so it is the percentage with more than this saturation
	pylab.plot(number,sat,'-',linewidth=3,label='%Pix > 1x    Saturation')
	sat=x[4]*100 
	sat=numpy.sort(sat)
	sat=numpy.flipud(sat) # Reverse so it is the percentage with more than this saturation
	pylab.plot(number,sat,'-',linewidth=3,label='%Pix > 0.5x Saturation')
	pylab.legend(loc='upper right')
	pylab.ylabel('% Pixels')
	pylab.xlabel('% images')
	pylab.axis([0,20,0.,3.0])
	pylab.savefig('bad_actor.png')
	pylab.draw()
	return
	



def steer(argv):
	'''
	This is a steering routine for bad_actor so that options can be exercised from the 
	command line

	111019	ksl	Adapted from the same routine in subtract_persist3
	'''
	i=1
	dataset_list='none'
	fileroot='observations'
	words=[]
	mjd_after=0.0    # A amall number for mjd
	mjd_before=1.e6  # A large number for mjd
	prog_id=0


	while i<len(argv):
		if argv[i]=='-h':
			print __doc__
			return    
		elif argv[i]=='-obslist':
			i=i+1
			fileroot=eval(argv[i])
		elif argv[i]=='-many':
			i=i+1
			dataset_list=argv[i]
			print 'OK you want to evaluate a number of datasets in file %s', dataset_list
		elif argv[i]=='-all':
			dataset_list='All'
			print 'OK you want to evaluate many, possibly all, of the records in the obslist'
		elif argv[i]=='-start':
			dataset_list='All'
			i=i+1
			z=argv[i]
			try:
				mjd_after=float(z)
				print 'OK you want records after %s' % z
			except ValueError:
				mjd_after=date.iso2mjd(z)
		elif argv[i]=='-stop':
			dataset_list='All'
			i=i+1
			z=argv[i]
			try:
				mjd_before=float(z)
				print 'OK you want records before %s' % z
			except ValueError:
				mjd_before=date.iso2mjd(z)
		elif argv[i]=='-prog_id':
			dataset_list='All'
			i=i+1
			prog_id=int(argv[i])
		elif argv[i][0]=='-':
			print 'Error: Unknown switch ---  %s' % argv[i]
			return
		else:
			words.append(argv[i])
		i=i+1
	
	print dataset_list
	if dataset_list=='none': #  Then we are processing a single dataset
		dataset=words[0]
		do_dataset(fileroot,dataset)
	elif dataset_list=='All': # Then we are handling multiple datasets
		g=open('Problems.txt','w')
		f=open('bad_actor_all.txt','w')
		records=per_list.read_ordered_list_progid(fileroot,prog_id,mjd_after,mjd_before)
		print 'There are %d images to process' % len(records)
		results=[]
		records_out=[]
		for word in records:
			xxxx=do_dataset(fileroot,word[1])
			if len(xxxx)>0:
				records_out.append(word)
				results.append(xxxx)
				string='%10s %10s %10s %10s %10s %20s %20s' % (word[1],word[2],word[3],word[10],word[11],word[14],word[16])
				scale=100.
				string2=' %8.1f %8.3f %8.3f %8.3f ' % (xxxx[1],xxxx[2]*scale,xxxx[3]*scale,xxxx[4]*scale)

				f.write('%s\n' % (string+string2))
			else:
				print('Error: Ignoring %s' % word[1])
				g.write('%30s %30s %30s\n' % (word[0],word[1],word[2]))
		# print 'OK',len(records),len(records_out),len(results)
		g.close()
		f.close()
		analyze(records_out,results)

	else:
		print "Don't know what to do"
	
	return

	



	 

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		steer(sys.argv)   
	else:
		print 'bad_actor.py  -h to get brief  help text'
	print 'OK done'
	sys.exit(0)
