#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine is the steering routine for persistence subtraction,
and evaluation.  It is effectively the user interfaace for 
evaluating persistence for one or more datastes

The routine requiress a file created by per_list.py that contains
a time-order list of flt files that are to be considered.


Description:  

This is a routine to run the other persistence subtraction,
evaluation and logging routines.

It calls other routines, such as subtract_persist to carry
out the actual work.  See the documentaion on these routines
to learn more

Basic usage is as follows

run_persist.py dateset   
	process a single dataset

-all     
	process all of the datasets in the .ls file

-all [-start time1]  [-stop time2] - process the datasets
	that are in the .ls file when the observations occured after time1 and/or before time2
	time1 and time2 are either MJD or ISO formated times, e.g '2009-11-12 15:13:24'

-many file.ls 
	process a list of dataset names (one dataset per line).  The assumption here
	is that you know exactly which datasets you want and so the time limits
	should not work.  If the file has multiple columns then the assumption is
	the dataset name is the first column

[-start time1]  [-stop time2] -prog_id  xxxxx 
	Evaluate the persistence for all of the visits of
	a program id xxx  with optional start and stop times

-obslist obsers2  
	specifices something other than the normal observation.ls file to use.
	as the file containing the list of possible flt files

-dryrun
	Indicates in a file what datasets would be processed

-clean
	Indicates that instead of processing datasets that the Persist directories should
	be removed from the datasets specified instead of processed.

-local
	Indicates that instead of putting all of the output files in directories beneath
	the visit directories, all the output files will be put in a subdirectory
	Persist of the current working direcgtory

Other switches allow you to control the persistence function that is subtracted, e. g.

Note that all of the numbers should be positive numbers!

-model  -- 0 for the orginal fermi-function based formalism
	   1 for the newer purely observational a gammma model
	   2 for the variable fermi function with differrent prescriptions depending on the
	     length of the exposre
-n	-- The normalization at 1000 s
-e	-- The fermi energy, nominally the place where the fermi distribution reaches half
	    the maximum
-kT	-- The width of the fermi function
-alpha	-- The power law to model very oversatruated events.
-gamma	-- The power law decay with time
-pf foo.pf  -- Replace the default parameter file contianing links to all of the files read
		by do persist to foo.pf
	
Outputs:

	a history.log file, which is a high level record of the processing 
		and errors, and 
	a command.log which is a record of the command sent to the
		steering routine for reference
	


Primary routines:

	steer  - This routine parses the input command line, and determines what datasets are to be
		processed
	doit  -  This routine controls the processing of a single dataset. 

Notes:
									   
History:

100603	ksl	Coding begun
101014	ksl	Began adding capabilities to use this in a production-like
		environment
110323	ksl	Modified to incorporate our new fitting function for persistence
		which requires and extra variable associated with the power law
		rise at very high stimulus levels.
110721	ksl	Added a crude form of versioning
121220	ksl	Updated documentation and cleaned up cals
140923	ksl	Changed version number to 3.0 to reflect the fact that this will
		be a major change in the software to allow for different types of
		models
141225	ksl	Minor changes in parameters the program accepts to reflect the
		use of a parameter file to contain names of files needed by the
		persistence model

'''

import os
import sys
import numpy
import math
import string
import time

import date

# Next two lines set the graphics backend to one that is noninteractive
import matplotlib
matplotlib.use('AGG')

import per_list
import subtract_persist
import peaks
import subtract_eval
import subtract_html


# Change this number when ever significant revisions of the S/W are made

import config
VERSION=config.version
	

def log(string,filename='history.log',option='a'):
	'''
	Log anything to a file and also print it to
	the screen

	Notes:
	121220 It might be better to assume that the string
	had no trailing carraibe return and add it
	'''
	history=open(filename,option)
	os.chmod(filename,0770)
	history.write(string)
	history.close()
	# Remove carriage returns from the string
	string=string.replace('\n',' ')
	print string
	return


def do_dataset(dataset='ia21h2e9q',model_type=1,norm=0.3,alpha=0.2,gamma=0.8,e_fermi=80000,kT=20000,fileroot='observations',ds9='yes',local='no',pffile='persist.pf'):
	'''

	Run the persistence 'pipeline' for a single dataset.

	where	dataset is the dataset name
		norm,alpha,gamma, e_fermi, kT define the model
		fileroot is the rootname of the file produced by per_lsist
		ds9 indicates whehter on should try to show resulst in ds9
		local=='yes'  imples that the output files are created directly below
			the place the program is being run, insteand of in the
			Visit directories
		pffile is the parameter file which contains the names of files needed
			to create models of persistence
		


	History

	101215	ksl	Begain coding to provide a mechanism to control the running of
			persistence suhtraction software for larger numbers of datasets
	110103	ksl	Added local switch for testing
	140929	ksl	Added new variable model_type to do_dataset, which is simply passed to subtract_persist.do_dataset.
	'''



	cur_time=date.get_gmt()

	print '# Processing dataset %s at %s' % (dataset,cur_time)


	log('# Starting dataset %s at %s\n' % (dataset,cur_time))

	record=per_list.read_ordered_list_one(fileroot,dataset)
	if len(record)>0:
		log('run_persist: flt file is %s\n' % record[0])
	else:
		log('run_persist: dataset %s not found in %s.ls\n' % (dataset,fileroot))
	        log('NOK: Processing aborted for dataset %s\n' %dataset)
		log('# Finished dataset %s at %s\n' % (dataset,cur_time))
		return 'Error: Dataset not found at %s' % cur_time

	# Carry out peristence subtraction for this dataset
	string=subtract_persist.do_dataset(dataset,model_type,norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,pffile)

	log('%s\n' % string)
	if string[0:3]=='NOK':
		sum_string='%s %s' % (cur_time,string)
		log('NOK: Processing aborted for dataset %s\n' %dataset)
		cur_time=date.get_gmt()
		log('# Finished dataset %s at %s\n' % (dataset,cur_time))
		return

	print string

	# Carry out peaks identification for this dataset
	string=peaks.do_dataset(dataset,fileroot,local=local)
	log('%s\n' % string)
	if  string[0:2]!='OK':
	                log('NOK: Processing aborted for dataset %s\n' %dataset)
			cur_time=date.get_gmt()
			log('# Finished dataset %s at %s\n' % (dataset,cur_time))
	                return
	
	# Evaluate the results at the positions identified in peaks

	string=subtract_eval.do_dataset(dataset,local=local)
	log('%s\n' % string)

	# Make an html file for the dataset

	string=subtract_html.do_dataset(dataset,local=local)
	log('%s\n' % string)
	if string[0:2]!='OK':
		log('NOK: Processing aborted for dataset %s\n' %dataset)
		cur_time=date.get_gmt()
		log('# Finished dataset %s at %s\n' % (dataset,cur_time))
		return


	words=string.split()

	# Now update the summary file
	per_list.update_summary(dataset,'Complete_%s'% VERSION,words[2])


	cur_time=date.get_gmt()
	log('# Finished dataset %s at %s\n' % (dataset,cur_time))
	return

def steer(argv):
	'''
	This is a steering routine for subtract persist so that options can be exercised from the 
	command line.  See the top level documentaion for details

	100907	ksl	Added to begin to automate the subtraction process
	101215	ksl	Moved to a separate routine so that one would be able to split various portions
			of persistence subtraction and evaluation of the results into multiple files
	111114	ksl	Added a command.log to keep track of all of the run_persist commands
	140924	ksl	Updated to allow for varioua model types
	'''

	log('# Start run_persist  %s\n' % date.get_gmt())
	xstart=time.clock()
	cur_time=date.get_gmt()

	i=1
	dataset_list='none'

	model_type=1
	norm=0.3
	alpha=0.174
	gamma=1.0
	e_fermi=90000
	kT=20000
	fileroot='observations'
	words=[]
	mjd_start=0.0    # A amall number for mjd
	mjd_stop=1.e6  # A large number for mjd
	dryrun='no'
	clean='no'
	ds9='no'
	local='no'
	pffile='persist.pf'

	switch='single'

	while i<len(argv):
		if argv[i]=='-h':
			print __doc__
			return    
		elif argv[i]=='-model':
			i=i+1
			model_type=int(argv[i])
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
			switch='many'
			print 'OK you want to evaluate a number of datasets in file %s', dataset_list
		elif argv[i]=='-all':
			switch='all'
			print 'OK you want to evaluate all the records in the obslist'
		elif argv[i] =='-prog_id':
			i=i+1
			prog_id=int(argv[i])
			switch='prog_id'
			print 'OK, you want to evaluate all the records for program %d' % prog_id
		elif argv[i]=='-start':
			i=i+1
			z=argv[i]
			try:
				mjd_start=float(z)
			except ValueError:
				mjd_start=date.iso2mjd(z)
			if switch !='prog_id':
				switch='all'
		elif argv[i]=='-stop':
			i=i+1
			z=argv[i]
			try:
				mjd_stop=float(z)
			except ValueError:
				mjd_stop=date.iso2mjd(z)
			if switch !='prog_id':
				switch='all'
		elif argv[i]=='-dryrun':
			dryrun='yes'
			print 'OK, This will be a dry run!'
		elif argv[i]=='-clean':
			dryrun='yes'
			clean='yes'
			print 'OK. This run will clean out various Persist directories, and revise the .sum file'
		elif argv[i]=='-ds9':
			ds9='yes'
		elif argv[i]=='-local':
			local='yes'
		elif argv[i]=='-pf':
			i=i+1
			pffile=argv[i]
		elif argv[i][0]=='-':
			print 'Error: Unknown switch ---  %s' % argv[i]
			return
		else:
			words.append(argv[i])
		i=i+1

	# At this point all of the options have been processed and we can
	# begin the processing of individual datasets

	# Check that the listfile actually exists, and if not exit with a stern warning

	listfile=fileroot+'.ls'
	if os.path.exists(listfile)==False:
		print 'Error: run_persist.steer - No %s file in run directory.  EXITING!!!' % listfile
		return
	
	# At this point, we are able to determine exactly what datasets to process
	log('# Starting at %s\n' % (cur_time),'command.log')
	log('# Starting at %s\n' % (cur_time))
	string=''
	for one in argv:
		string=string+'%s ' % one
	log('Command:  %s\n' % string,'command.log')
	log('Command:  %s\n' % string)


	datasets=[]
	
	if switch=='single': #  Then we are processing a single file
		datasets.append(words[0])
	elif switch=='all': # Then we are working from the obslist
		records=per_list.read_ordered_list_mjd(fileroot,mjd_start,mjd_stop)
		for record in records:
			datasets.append(record[1])
	elif switch=='many':  # Then we are reading a file with rootnames of the files we want to process
		f=open(dataset_list,'r')
		lines=f.readlines()
		f.close()
		for line in lines:
			x=line.strip()
			if len(x)>0 and x[0]!='#':
				xx=x.split()  #  Take the first word to be the dataset name
				datasets.append(xx[0])
	elif switch=='prog_id':
		records=per_list.read_ordered_list_progid(fileroot,prog_id,mjd_start,mjd_stop)
		for record in records:
			datasets.append(record[1])
	else:
		print 'Error: run_persist: Unknown switch %s'% switch

	# Ok, now, unless this is a dryrun we actually process the data
	ntot=len(datasets)
	print 'There are %d datasets to process' % ntot          

	dry=[]
	if dryrun=='yes':
		for one in datasets:
			record=per_list.read_ordered_list_one(fileroot,one)
			dry.append(record)
		per_list.write_ordered_list('dryrun',dry)
		if clean=='yes':
			# xclean=open('CleanFiles','w')
			# os.chmod('CleanFiles',0770)

			xclean=per_list.open_file(Cleanfiles)

			for one in dry:
				xxx=per_list.set_path(one[0])
				xclean.write('rm -r -f %s\n' % xxx)
			xclean.close()
		# return
	else:
		n=1
		for one in datasets:
			do_dataset(one,model_type,norm,alpha,gamma,e_fermi,kT,fileroot,ds9,local,pffile)
			print '# Completed dataset %d of %d. Elapsed time is %0.1f s' % (n,ntot,time.clock()-xstart)
			n=n+1


	
	dtime=time.clock()-xstart
	cur_time=date.get_gmt()
	log('# End  %s  (Elapsed time for %d datasets was %.1f (or %.1f per dataset)\n' % (date.get_gmt(),ntot,dtime,dtime/ntot))
	cur_time=date.get_gmt()
	log('# End  %s  (Elapsed time for %d datasets was %.1f (or %.1f per dataset)\n' % (date.get_gmt(),ntot,dtime,dtime/ntot),'command.log')

	log('# Finished at %s\n' % (cur_time))
	log('# Finished at %s\n' % (cur_time),'command.log')
	return

	



	 

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		steer(sys.argv)   
	else:
		print 'run_persist.py  -h to get brief  help text'
	print 'OK done'
	sys.exit(0)
