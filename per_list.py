#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Usage (as a standalone program: 
	per_list.py 

or 

per_list.py [various options]  fileroot  

where 	fileroot     the root name of the outputfiles
	aperture is 'all' or 'full' with the latter indicating one creates a list containing only the full IR frames, no subarrays

	Options:
	-h		Print this help
	-all		Create a new .ls file and a new summary file
	-new_sum 	Create a new summary file
	-daily		Create a new .ls but just update the old summary file
	-file_type	Instead of the defult flt files, create a file continaing
			a different file type, e.g. -file_type raw  to get the
			raw data files

without any arguments the call is effectively

	per_list.py  -all -file_type flt observations 


	per_list.py
	
Synopsis:  


If this is run as a standalone program, the routine creates two files.  The first
is set of ascii records that describe all the WFC3/IR files of a certain type 
that are in this directory and any subdirectories.  The second is an observations.sum
file that really is only pertinent to the persistence processing software

More generally this module contains routines both to make an ordered list of the records
and then to read them back, as well as a number of utilities to handle matters 
related to where the Persist directories are located.

Description:  

	The routine uses 'find' to locate all of the files of a certain type,
	aand then orders the list in time.  If there are duplicate files then
	it uses the creation date as the one to put in the list.  The dupliccates
	are in the list also, but are commented out. 

Primary routines:

	steer 			Parses the command line and runs the program
	mske_ordered_list	Creates .ls file
	make_sum_file		Creates/updates the summary file

Notes:
	This should probably be rewritten to remove the iraf dependencies
									   
History:

091026 ksl	Coding begun
091215 ksl	Updated to prevent reading UVIS data and put in an option to only
		look at full arrays
1006	ksl	Split out from other routines so this routine could be used as
		a library of routines for dealing with ascii lists
100608	ksl	Added the ability to create lists from the command line by running 
		per_list.py
100906	ksl	Added what amounts to a standalone utility find_latest to locate
		the last version of a file in any subdirectory if there is no
		such file in the current directory.  This is here in this file
		becasue per_list.py contains all the file locating routines
		not because it is really part of the rest of per_list
110103	ksl	Modified various portions of the way with which summary files are
		handled and introduced a steering routine as the input options
		became more complicated.
120330	ksl	Modified the help so that it reflects how the routine actually 
		works
130909 ksl	Standardized error printouts
150317	ksl	Remove iraf/pyraf from per_list
151022	ksl	Restored iraf/pyraf for performance reasons
160104	ksl	Changes to lock and unlock files so that parallel processing
		will not corrupt the obervations.sum file

'''

import sys
import os
import numpy
import time
import string
import math
import scipy
import subprocess
import pylab
import shutil
import per_fits
from astropy.io import fits
import pyraf
from multiprocessing import Pool


# Utilities
def backup(filename,tformat='%y%m%d',force='no'):
	'''
	Backup a file with a standard way of naming the file.  
	
	The original file is copied to the backup unless the
	backup file already exists or force is 'yes'.
	
	The format can be any valid time.strftime format, but
	probably should have no spaces.

	The routine returns the name of the backup file, or 
	'None' if there was no backup made

	Notes:  The idea of this is that one wants to keep
	a history of a file over time, but not at infinite
	time resolution.

	110105	ksl	Coded as part of effort to be able
			to use per_list as aprt of a crontab
	'''

	curtime=time.strftime(tformat,time.localtime())
	if os.path.exists(filename)==False:
		print 'Warning: No file %s to backup' % filename
		return 'None'

	backup_file='%s.%s.old' % (filename,curtime)

	if os.path.exists(backup_file)==False or force=='yes':
		shutil.copy(filename,backup_file)
		return backup_file
	else:
		return 'None'

def open_file(filename,permiss=0770):
	'''
	Open a file for writing and set its permissions

	This was written in an attempt to get all of the
	files to have a common set of permissions
	no mattter who writes them

	'''

	if os.path.isfile(filename):
		os.remove(filename)
	g=open(filename,'w')
	try:
		os.chmod(filename,permiss)
	except OSError:
		print 'Error: open_file: OSerror trying to set permissions'
	return g



def set_path(name,mkdirs='yes',local='no'):
	'''
	Check that the directory Persist where the outputs are to go exists, and if not
	create it, and also create an underlying directory .../Persist/Figs.  

	The input is the complete filename, including the path.  The routine 
	requires that there be a '/' in the complete name.  If not, it will
	create Persist below the current working directory 

	The routine returns the path include the directory name Persist

	Note:

	Normally this routine is given either the complete name to a file in the
	visit directory, e. g. 
		./D/VisitAJ/ib0lajdyq_flt.fits[1] or 
	a name to a file in the  Persist directory, e. g.
		./D/VisitAJ/Persist/filename
	This accounts for the somewhat convoluted attempt to locate where the Persist
	directory should be

	mkdirs must be 'yes' for directores to be made.  If it is for example, 'no', 
	then the path will be return but there will be not attempt to make the
	directory

	if local is anything but 'no', then the path will be set to ./Persist

	

	101014	ksl	Added
	101214 	ksl	Moved into per_list since this is mostly used in conjunction
			with other routines that are there. It is not obvious that
			this is the correct place for this.
	101215	ksl	Added creation of the Figs directory
	110105	ksl	Made creation of the directories an option
	110121	ksl	Made change to control permissions of the directories
	110203	ksl	Made change to allow the directories to be written beneath
			the current working directory, a change that is primarily
			for testing
	110811	ksl	The directory permissions are set so that the group name 
			should be inherited.
	110811 	ksl	Added commands to set the group name, but only if we are
			in the Quicklook2 directory structure.  These commands
			would need to change if the group names change or directory
			names change.  They should have no effect outsdide the standard
			structure
	'''

	# 110120  - I am not sure why but I had to reimport string
	import string

	if len(name)==0:
		print 'Error: set_path: name had length 0, nothing to parse'
		return ''


	# Determine where we want Persist to be located.  
	if local=='no':
		try:
			i=string.rindex(name,'Persist')
			path=name[0:i]
		except ValueError:
			# Find the last / in the name

			try:
				i=string.rindex(name,'/')
				path=name[0:i]
			except ValueError:
				print 'Warning: set_path: Assuming the intent was to work in the current directory'
				path='./'
	else:
		path='./'

	# Check whether the parent directory for Persits exists
	if os.path.exists(path)==False:
		string='set path: The directory %s contained in %s does not exist' % (path,name)
		print 'Error:  %s' % (string)
		return 'NOK %s ' % string


	path=path+'/Persist/'


	if mkdirs!='yes':  
		return path

	# If one has reached this point then you want to make any necessary directories
	# set the chmod of the diretory so those within the group can delete the directory
	# when necessary
	
	# Note that the group name assignments below only are valid with the current
	# Quicklook2 directory structure and group names



	if os.path.exists(path)==False:
		try:
			os.mkdir(path)
			os.chmod(path,02770)
			if path.count('QL_GO'):
				os.chown(path,-1,6047)
			elif path.count('Quicklook'):
				os.chown(path,-1,340)


		except OSError:
			print 'Error: set_path: Could not create %s for %s' % (path,name)
			return 'NOK Could not create %s' % path


	# Add a figs directory if it does not exist as well
	figs=path+'/Figs'
	if os.path.exists(figs)==False:
		os.mkdir(figs)
		os.chmod(figs,02770)
		if figs.count('QL_GO'):
			os.chown(figs,-1,6047)
		elif figs.count('Quicklook'):
			os.chown(figs,-1,340)

	# print 'set_path',figs,os.path.exists(figs)
	
	return path

def parse_dataset_name(name):
	'''
	Check if we have been given the name of the fits file instead of the 
	dataset name, and if so try to determine and return the dataset name

	100103	ksl	Coded because it is natural in some cases to give the filename
	'''
	xname=name
	if string.count(xname,'.') > 0 or string.count(xname,'/') > 0:
		# This must be a filename
		i=string.rindex(xname,'/')
		xname=xname[i+1:i+10]
	
	return xname

def parse_creation_time(xtime='2010-07-10T19:11:52'):
	'''
	Parse the time string for the date at which a fits file
	was created and return a time that can be used to choose
	which of two fits files was created last

	The time returned approximates the number of days since
	1990 but does not account for the fact that there are
	different numbers of days in a month.  

	The routine returns 0.0 and raises an IndexErro exception 
	if the time cannot be parsed

	Note the use of int instead of eval because leading zeros
	caused eval to interprete the numbers as octal

	100817	Added error checking which arose because the values
		being passed to routine were actually not the
		creation date

	'''

	xtime=xtime.replace('-',' ')
	xtime=xtime.replace('T',' ')
	xtime=xtime.replace(':',' ')
	xtime=xtime.split()

	try:
		day_sec=int(xtime[3])*3600.+int(xtime[4])*60+int(xtime[5])
		day_frac=day_sec/86400
		# Next section for day is really not accurate, but we don't care
		year_frac=(int(xtime[1])+int(xtime[2])/31.)/12.
		day=365.*((int(xtime[0])-1990.)+year_frac)
	except IndexError:
		raise IndexError
		return 0.0

	pseudo_time=day+day_frac

	return pseudo_time

def check4duplicates(records):
	'''
	Check the list 4 duplicate records, and choose the one that was
	created last if that is possible

	The routine returns a list that contains 'ok' for files that 
	are the ones to use and 'nok' for those that are duplicates

	100817	Added checks to trap problems parsing times.
	110120	This is a new attempt to find the duplicates and select the last one
	'''


	ok=[]
	for record in records:
		ok.append('ok')
	
	i=0
	while i<len(records):

		# if we have already determined that this record is nok, then don't investigate
		# further.
		if ok[i]!='ok':
			i=i+1
			continue

		one=records[i]
		hold=[i]

		j=i+1
		while j<len(records):
			two=records[j]
			if one[1]==two[1]:
				hold.append(j)
			j=j+1

		if len(hold)>1:  # The length of hold must be greated than 1 to have to worry about duplicates

			times=[]
			for one in hold:
				times.append(records[one][17])

			times=numpy.array(times)
			order=numpy.argsort(times) # Give me the order of the times

			last=order[len(order)-1]

			k=0
			while k<len(hold):
				if k==last:
					ok[hold[k]]='ok'
				else:
					ok[hold[k]]='nok'
				k=k+1
		# Now go on to the next record.  
		i=i+1
	return ok
	




def find_latest(filename='foo'):
	'''

	This is a simple routine to locate a specific version of a file.  If
	the file is in the current directory that is the file name that
	will be returned.  If it is in one of the subdirectories of the
	current directory then the one that was modified most recently
	will be returned. 

	Notes:

	This routine uses the unix utility find.  It is therefore
	likely to be slow in large directory structures

	This uses subprocess which handles stdin and stdout, unlike
	os.system

	History:

	100906 ksl Coding begun.  There is a standalone version of this 
		   called find.py (in my normal py_progs/scripts directory


	'''

	proc=subprocess.Popen('find . -follow -name %s -print ' % filename,shell=True,stdout=subprocess.PIPE)
	# Before
	lines=proc.stdout.readlines()
	if len(lines)==0:
			'Warning: find_best: No versions of %s found' % filename
			return ''
	tbest=0
	for line in lines:
		fname=line.strip()
		if fname.count('/') <= 1:
			fbest=fname
			break
		else:
			time=os.path.getmtime(fname)
			if time>tbest:
				fbest=fname
				tbest=time
	return fbest







def check4scan(filename='./Visit43/ic9t43j1q_flt.fits[1]'):
	'''
	Determine whether or not a raw, ima, or flt file is associated with an observation involving  
	a spatial scan.  

	Return:
		scan	if a spatial scan
		stare   if the spt file exists, but it is not a spatial san
		unknown	if the spt file does not exist
	
	Notes:

	The routine looks for the spt file corresponding to the dataset and
	parses the header to find out if it is a scanned observation.
	
	130225  Coded and Debugged
	130307  Replaced routine using astropy.fits with iraf because astropy.fits
		was very slow
	'''

	xscan='unknown'

	# First strip of the extension if any
	xname=per_fits.parse_fitsname(filename,0,'yes')
	xfile=xname[0]
	xfile=xfile.replace('flt','spt')
	xfile=xfile.replace('raw','spt')
	xfile=xfile.replace('ima','spt')

	# Now we should have the name of the spt file
	if os.path.exists(xfile)==True:
		# Revert to iraf/pyraf for performance reasons
		# xx=per_fits.get_keyword(xfile,0,'SCAN_TYP')
		xx=pyraf.iraf.hselect(xfile+'[0]','SCAN_TYP','yes',Stdout=1)
		xx=xx[0].split('\t')
	else: 
		return 'no_spt'



	if xx[0]=='N':
		xscan='stare'
	else: 
		xscan='scan'


	return xscan  
	

# End utilities

# Routines for making and updating the summary file

def update_summary(dataset,status_word='Unknown',results='Whatever you want',fileroot='observations',append='yes'):
	'''
	Update a record in the summary file, where the inputs have the following meaning
	
	dataset		rootname of the file, usually the flt file, being processed
	status_word	one word, e.g. Complete or Error, to define the current status of the processing
	results		A string which can be appended to the current set of results or which can replace it
	fileroot	The rootname of the summarry file, almost always obsevations
	append		if 'yes', then append the new results to the old, else replace the old results with the
			new


	Notes:

	Aside from the dataset name everything is free format in terms of results


	It is not really obvious that appending to a line is the best approach.  If I werre writing this
	today, I might opt for a fixed format from the beginning, but this is part of a larger question
	of whether one should use astropy. tables, or some kind of database system - ksl - 151006

	Before 1601 this routine modified the observation.sum file directly, but now it just records the
	information needed to update the summary file in a directory tmp_sum.  This change was to prevent
	multiple processing writing to the summary file at the same file.  A new routine fixup_summary_file
	below now inserts the results into the data.  


	101221	ksl	Coded as part of effort to get a way to check what files
			had been processed with the persistence software
	110103	ksl	Added a better way to keep enough summary files that one
			might be able to roll back
	110117	ksl	Dealt with the situation where there were not old_results
	110721	ksl	Improved the description of the routine
	151006	ksl	Fixed Error #553, which arose because the routine assumes one
			often wants to append to a summary line.  Because we use index
			to find where in the line we want to start replacing data, one
			needs to be careful that there is no possibility that one is
			indexing to the wrong position.
	160105	ksl	Modified for multiprocessing

	'''


	if os.path.isdir('tmp_sum')==False:
		os.mkdir('tmp_sum')
	

	summary_file=fileroot+'.sum'
	gmt=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())


	try:
		f=open(summary_file,'r')
		lines=f.readlines()
		f.close()
	except IOError:
		print 'Error: update_summary: File %s does not exist' % summary_file
		return


	# this opens the files for writing and sets the permissions
	xfile='tmp_sum/%s.txt' % dataset
	if os.path.isfile(xfile)==True:
		gg=open(xfile,'r')
		old_results=gg.readline()
		gg.close()
	else:
		i=0
		while i <len(lines):
			line=lines[i].split()
			if line[0]==dataset:
				old_results=lines[i]
				break
			i=i+1
	
	# OK at this point I have the old results

	results=results.strip()

	old_results=old_results.strip()
	line=old_results.split()

	if len(line)>6:
		k=old_results.index(line[5])+len(line[5])
		old_results=old_results[k:len(old_results)] # This should be the bits that were append previously
		old_results=old_results.strip()

	if append=='yes' and len(old_results)>0:
		results='%s %s' % (old_results,results)
	string='%-10s %5s %20s  %20s %-20s %s' % (line[0],line[1],line[2],gmt,status_word,results)

	gg=open_file('tmp_sum/%s.txt' % dataset)
	gg.write('%s\n' % string)
	gg.close()

	return



def  fixup_summary_file(datasets,fileroot='observations'):
	'''
	Update the observations.sum file from the information stored in the 'tmp_sum' directory

	where 
		datasets is a list of the datasets which have been processed
		fileroot is the root name of the .sum file to be updated
	
	Notes:
		
		The directory where the data is stored is hardwired

		The routine simple finds the line in the observations.sum file associated with
		the dataset and replaces that line with the line that is in the .txt file
	
	History:

	160105	ksl	Coded as part of the effort to add multiprocessing

	'''


	# Get the data
	good=[]
	data=[]
	for dataset in datasets:
		try:
			xname='tmp_sum/%s.txt' % dataset
			f=open(xname,'r')
			line=f.readline()
			data.append(line)
			good.append(dataset)
		except IOError:
			print 'Error: fixup_summary_file: %s does not appear to exist' % xname
	
	# OK at this point we have found all of the good files

	summary_file=fileroot+'.sum'

	try:
		f=open(summary_file,'r')
		lines=f.readlines()
		f.close()
	except IOError:
		print 'Error: update_summary: File %s does not exist' % summary_file
		return


	i=0
	while i<len(lines):
		line=lines[i].split()
		j=0
		while j<len(good):
			if line[0]==good[j]:
				lines[i]=data[j]
				print lines[i].strip()
				break
			j=j+1
		i=i+1
	

	g=open_file('tmp.sum')
	for one in lines:
		g.write(one)




	backup(summary_file)
	proc=subprocess.Popen('mv %s %s' % ('tmp.sum',summary_file),shell=True,stdout=subprocess.PIPE)

	return





	




def make_sum_file(fileroot='observations',new='no'):
	'''
	Make the file that will record the results of persistence processing.  If
	it already exists, give the user the option of merging new records into
	the old file or creating a new file

	Note that each new output line should have the following values

	dataset	prog_id MJD  ProcessDate Unprocessed


	101221	ksl	Coded as part of effort to put a recording mechanism
			in place
	110103	ksl	Modified the outputs when new records are inserted
	110103	ksl	Changed so the routine itself reads the observation.ls
			file
	110119	ksl	Rewrote to assure that there is exactly one summary 
			file line for each observation file line
	'''



	# Read the entire observations.ls file
	print '# (Re)Making summary file'
	records=read_ordered_list0(fileroot)
	print '# The number of records in the old per_list file is %d' % (len(records))

	gmt=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

	summary_file=fileroot+'.sum'


	if os.path.exists(summary_file)==False or new=='yes':
		print '# Making a pristine summary file'

		g=open_file(summary_file)

		for record in records:
			string='%-10s %5s %20s  %20s %-20s' % (record[1],record[2],record[6],gmt,'Unprocessed')
			g.write('%s\n' % string)
		g.close()
	else:
		print '# Merging new records into old list'
		try:
			f=open(summary_file,'r')
			summary=f.readlines()
		except  IOError:
			print 'Error: make_sum_file: Could not read summary file %s ' % summary_file

		# Get the dataset names in the summary file
		datasets=[]
		for one in summary:
			word=one.split()
			datasets.append(word[0])
			

		g=open_file('tmp.sum')

		# The next portion is not very smart since it always goes through the summary information
		# from the beginning.  But is should be "sure". 
		for record in records:
			i=0
			while i<len(datasets):
				if datasets[i]==record[1]:
					break
				i=i+1

			if i==len(datasets):
				string='%-10s %5s %20s  %20s %-20s' % (record[1],record[2],record[6],gmt,'Unprocessed')
				g.write('%s\n' % string)
			else:
				g.write(summary[i])

		g.close()

		ndup=check_sum_file('tmp.sum',summary_file)
		if ndup>0:
			print 'Error: make_sum_file: Since there were duplicates in th tmp file, not moving to %s'  % (summary_file)
			return 'NOK  - Since there were duplicates in th tmp file, not moving to %s'  % (summary_file)

		# Now move the files around.  Note that the next 3 lines need to be the same as in update_summary above
		# gmt=time.strftime("%y%m%d.%H%M", time.gmtime())  # Create a string to use to name the updated file. As written a new file a minute
		# proc=subprocess.Popen('mv %s %s.%s.old' % (summary_file,summary_file,gmt),shell=True,stdout=subprocess.PIPE)
		backup(summary_file)
		proc=subprocess.Popen('mv %s %s' % ('tmp.sum',summary_file),shell=True,stdout=subprocess.PIPE)
	return


def check_sum_file(new='tmp.sum',old='none'):
	'''
	Check a summary file for unique dataset names

	Description

	The routine checks the new summary file to assure 
	that there is only one line with a given dataset name.
	It returns the number of duplicate records

	If a file name is given in old, then the routine checks
	whether there are any records that were in the old file, that
	are not in the new file.  This can happen and not be an error,
	but the results are recorded in a file with the extension problem

	Notes:

	There should be one and only one record  for each
	record in observaitions.ls

	110119 - Expanded the questions that were asked 
	150127 ksl Rewrote so that it was much faster by using sets
	'''

	# Read the new summary file
	summary_file=new
	if os.path.exists(summary_file)==False:
		print 'Error: check_sum: There is no summary file names %s.sum to check' % fileroot
		return
	f=open(summary_file,'r')
	lines=f.readlines()
	f.close()

	
	# Open a diagnostic file to record any problems
	g=open_file(summary_file+'.problem')
	string='# Checking new (%s) against (%s) old summary file' % (new,old)
	print string
	g.write('%s\n' % string)

	# Check for duplicate datasets in new summary file

	xstart=time.time()
	names=[]
	for line in lines:
		x=line.strip()
		word=x.split()
		names.append(word[0])

	unique=set(names)
	if len(unique)==len(names):
		dups=0
	else:  # There is a problem and so the duplicate check must be carried out
		print 'Error: There are duplicate datasets in the the newly created summary file %s' % new
		print 'A detailed check is now underway'
		dup_names=[]
		dups=0
		for one in unique:
			jj=names.count(one)
			if jj>1:
				g.write('Duplicate: %5d %s\n' % (jj,one))
				dup_names.append(one)
				dups=dups+1

	string='# Check of %s revealed %d duplicates' % (new,dups)
	print string
	g.write('%s\n' % string)
	delta_t=time.time()-xstart
	print 'check_sum: time to search for duplicate records in the new .sum file: ',delta_t

	# Check for datasets that were in the .sum file before but are missing now

	if old!='none':
		print 'Checking for datasets that were in the last .sum file but are now missing'
		xstart=time.time()
		old_names=[]

		try:
			f=open(old,'r')
			old_lines=f.readlines()
			f.close()

			for one in old_lines:
				 x=one.strip()
				 word=x.split()
				 old_names.append(word[0])
		except IOError:
			old_lines=[]

		nlost=0
		for one in old_names:
			if names.count(one)==0:
				g.write('Missing: %s\n' % (one))
				nlost=nlost+1

		string='# Check of %s revealed %d datasets no longer in the .sum file' % (old,nlost)
		if nlost>0:
			print 'Warning: Lost datasets may not be a problem, but are noted as a warning of something possibly worth invesigation'
		print string
		g.write('%s\n' % string)
		print '# check_sum_file: The check for missing records took %s s' % delta_t







	# xstart=time.time()
	# names=[]
	# dups=0
	# for line in lines:
	# 	x=line.strip()
	# 	word=x.split()
	# 	j=0
	# 	while j<len(names):
	# 		name=names[j]
	# 		if word[0]==name:
	# 			g.write('%5d %s' % (j,line))
	# 			dups=dups+1
	# 			break
	# 		j=j+1
	# 	if j==len(names):
	# 		names.append(word[0])
# 	
	# string='# Check of %s revealed %d duplicates' % (new,dups)
	# print string
	# g.write('%s\n' % string)
	# delta_t=time.time()-xstart
	# print '# check_sum_file: The check for duplicates took %s s' % delta_t


	# if old!='none':
	# 	xstart=time.time()
	#
	#	try:
	#		f=open(old,'r')
	#		old_lines=f.readlines()
	#		f.close()
	#	except IOError:
	#		old_lines=[]
#
	#	nlost=0
	#	for one in old_lines:
	#		x=one.strip()
	#		word=x.split()
	#		j=0
	#		while j<len(names):
	#			name=names[j]
	#			if word[0]==name:
	#				break   # Then we have a match
	#			j=j+1
	#		if j==len(names): # There was no match
	#			print 'Lost record:',x
	#			g.write('%5d %s' % (j,one))
	#			nlost=nlost+1
#
	#	string='# Check of %s revealed %d lost records' % (old,nlost)
	#	print string
	#	g.write('%s\n' % string)
	#	print '# check_sum_file: The check for missing records took %s s' % delta_t

	g.close()


	return dups





def get_info(lines,apertures,filetype):
	'''
	Get all of the keyword information for a set of files and return this

	Notes:

	This section of the old make_ordered list was put into a separate function
	so that it could be run in parallel 

	History

	160118	ksl	Added
	'''

	records=[]
	times=[]

	if len(lines)==0:
		print 'There were no %s files in the directory structure' % filetype
		return []
	else:
		print 'There are %d datasets to process' % len(lines)

	i=0
	for line in lines:

		line=line.strip()
		line=line.split()

		# 100807 Added date of file creation so could handle non-unique data sets
		# This is the old version using pyraf, which has been put back for perfomance reasons
		xfile='%s[1]' % line[0]
		if pyraf.iraf.imaccess(xfile):
			x=pyraf.iraf.hselect(xfile,'$I,rootname,proposid,linenum, instrume,detector,expstart,date-obs,time-obs,aperture,filter,exptime,crval1,crval2,targname,asn_id,pr_inv_L,date','yes',Stdout=1)
			x=x[0].split('\t')
			# Kluge for raw data files which have two ROOTNAME keywords for unknown reasons
			if x[1]==x[2]:
				x.pop(2)

			x[16].replace(' ','-')  # Get rid of spaces in PI names
			# Another kludge for raw files.  The is no 'date' field in the first extension as there is for flt and ima files, but DATE does exist in extension 0
			if filetype=='raw':
				xname=per_fits.parse_fitsname(xfile,0,'yes')
				xx=pyraf.iraf.hselect(xname[2],'$I,DATE','yes',Stdout=1)
				xx=xx[0].split('\t')
				x.append(xx[1])
			x[0]=line[0]

#		# Replaced upcoming lines with iraf/pyraf for performance reasons
#		xfile=line[0]
#		if os.path.isfile(xfile) == True:
#			x=per_fits.get_keyword(xfile,1,'rootname,proposid,linenum, instrume,detector,expstart,date-obs,time-obs,aperture,filter,exptime,crval1,crval2,targname,asn_id,pr_inv_L')
#			x=[xfile]+x

#			# for raw files the date is in extension 0, but for the others it is in 1.
#			if filetype=='raw':
#				date=per_fits.get_keyword(xfile,0,'date')
#			else:
#				date=per_fits.get_keyword(xfile,1,'date')
#			x.append(date[0])



			scan=check4scan(xfile)

			if x[5]=='IR':
				x[4]=scan
				j=string.count(x[9],'SUB')
				if j == 0 or apertures != 'full':
					records.append(x)
					times.append(float(x[6]))
		else: 
			print 'File %s does not really exist' %  xfile
		i=i+1
		if i%100 == 1:
			print 'Inspected %6d of %6d datasets --> %6d IR datasets' % (i,len(lines),len(records))
	print 'Inspected %6d of %6d datasets --> %6d IR datasets' % (i,len(lines),len(records))
	return times,records


def info_helper(args):
	'''
	This repackages the argments so that they can be used in parallel processing

	Notes:

	The issue here is that you only can pass a single argument via the Pool mechanism.
	This simple allows that.
	'''

	return get_info(*args)



def make_ordered_list(fileroot='observations',apertures='full',filetype='flt',new='no',np=1):
	'''
	find all of the observations in all subdiretories and make
	a time ordered list of the observations from files of a
	given filetype, e. g. flt.  If apertures == 'full' only
	put full frames into the list.  Otherwise, include subarrays



	Note - This could be done as a true database, but it's 
	simpler for now just to make it a file

	This routine still uses iraf.hselect

	100308 - Added option to deal with other types of files aside from flt files.  The 
		assumption made is that the filetype is part of the filename
	100427 - Added crval1 and crval2 to output list because want to use this to determine
		whether there was a dither or move between observations
	100817 - Added checks to verify that a file found by find actually exists and has the
		appropriate image extension (namely 1)
	110103	- Added the switch new as a stepping stone to making the routine more general. AT
		Some point we are going to need to merge new files into old ones and that is
		partially the reason for this.  Right now it is simply to avoid being asked
		if one wants to make a new file
	110104	Split off the actual writing of the file into a separate routine called
		write_time_sorted in order to make sublists.  Any change in the information 
		that is gathered by make_ordered list, needs to be relected in that routine
	110121	Added a line to assure that PI names had no spaces.  This was because running
		split on records caused with spaces int he PI name was failing.  The alternative
		which was to change the separated character to a tab seemed more trouble
	120330	Added a kludge to add a date for a raw data file.  The prupose of this date is 
		resolve situations in which multiple versions of the same file are in the
		directory structure.  The choice is only accurate to the day
	130820	Revised the way in which searched for the date in raw data files.
	140306	Add capability to check for scanned observations.
	160118  Added the possibility of running multiple processors
	160119  Switched from time.clock to time.time so that everything would be in wall clock
		time
	'''

	if filetype!='flt':
		fileroot=fileroot+'_'+filetype


	backup(fileroot+'.ls')
	xstart=time.time()
	os.system('find . -follow -name \*%s.fits  -print 1>files.ls 2>/tmp/foo; rm /tmp/foo' % filetype)
	search_time=dtime=time.time()-xstart
	print '# Found all of the files to read in %f s' % dtime
	xstart=time.time()


	f=open('files.ls','r')
	lines=f.readlines()
	f.close()

	if np<=1 or len(lines)<np:
		times,records=get_info(lines,apertures,filetype)
	else:
		inputs=[]
		idelta=len(lines)/np +1
		i=0
		while i < len(lines):
			imin=i
			imax=i+idelta
			if imax>len(lines):
				imax=len(lines)
			xinputs=lines[imin:imax]
			inputs.append([xinputs,apertures,filetype])
			i=i+idelta

		# i=0
		# while i<len(inputs):
		# 	print i,inputs[i]
		# 	i=i+1

		p=Pool(np)  

		times=[]
		records=[]
		for one in p.map(info_helper,inputs):
			times=times+one[0]
			records=records+one[1]
	

	# print 'times',len(times)
	# print 'records',len(records)

	key_time=dtime=time.time()-xstart
	print 'Inspected all of the files in %f s' % dtime
	

	if len(times)==0:
		print 'There were no IR observations to consider'
		return []
	# Now sort this all on the time
	# This returns an index of the order of the lines
	xstart=time.time()
	order=numpy.argsort(times)
	sort_time=dtime=time.time()-xstart
	print '# Time to sort the records %s s' % dtime

	lastime=float(records[len(order)-1][6])
	
	time_sorted=[]
	for index in order:
		time_sorted.append(records[index])

	# Now check for uniqueness files



	# Now write the time sorted file
	xstart=time.time()
	write_ordered_list(fileroot,time_sorted)
	write_time=time.time()-xstart
	print '# Time to write the .ls file  %s s' % write_time

	print '# Completed creating %s.ls' % (fileroot)

	print '# make_ordered_list: times',search_time,key_time,sort_time,write_time

	return time_sorted
	
def write_ordered_list(fileroot='observations',records=[]):
	'''
	Write the ordered_list file from a set of records

	110104	ksl	Split from make_ordered list in order to ease the creation of sublists
	'''

	# Now write it all to an ascii file
	# f=open(fileroot+'.ls','w')
	# os.chmod(fileroot+'.ls',0770)

	f=open_file(fileroot+'.ls')

	# Note that the duplication check here is awkward, but we want to
	# keep track of duplicate records in the output file
	ok=check4duplicates(records)

	# Now write out each record one by one
	i=0
	while i<len(records):
		record=records[i]
		# print 'xxx',record
		if ok[i]=='ok':
			xstring='%-50s ' % record[0]  # File name
		else:
			xstring='# %-48s ' % record[0]  # File name
		xstring=xstring+'%-5s ' % record[1]  # dataset name
		xstring=xstring+'%-5s ' % record[2]  # Propid
		xstring=xstring+'%-5s ' % record[3]  # Lineno
		xstring=xstring+'%-5s ' % record[4]  # Intrum
		xstring=xstring+'%-5s ' % record[5]  # detector
		# Next line needed for pyraf
		xstring=xstring+'%14.7f ' % (eval(record[6]))  # expstart
		# Next line removed to switch from astropy to pyraf
		# xstring=xstring+'%14.7f ' % record[6]  # expstart
		xstring=xstring+'%-5s ' % record[7]  # date-obs
		xstring=xstring+'%-5s ' % record[8]  # time-obs
		xstring=xstring+'%-5s ' % record[9]  # aperture
		xstring=xstring+'%-5s ' % record[10] # filter 
		# Next line needed for pyraf
		xstring=xstring+'%7.1f ' % (eval(record[11])) # exptime
		xstring=xstring+'%10.6f ' % (eval( record[12])) # crval1 
		xstring=xstring+'%10.6f ' % (eval(record[13])) # crval2  
		# Next lines removed in swich back to pyraf
		# xstring=xstring+'%7.1f ' % record[11] # exptime
		# xstring=xstring+'%10.6f ' % record[12] # crval1 
		# xstring=xstring+'%10.6f ' % record[13] # crval2  
		xstring=xstring+'%-5s ' % record[14] # targname
		xstring=xstring+'%-5s ' % record[15] # Assn_id 
		xstring=xstring+'%-5s ' % record[16] # PI      
		xstring=xstring+'%-5s ' % record[17] # File creation date
		f.write('%-5s\n' % xstring)
		i=i+1
	f.close()
	return

	

def read_ordered_list0(fileroot='observations'):
	'''
	This simply reads the list and returns all of the records

	Eventually should replace what is in read_ordered_list and
	read_ordered_list2
	'''

	# Read the entire file
	records=[]
	times=[]
	try:
		f=open(fileroot+'.ls','r')
		lines=f.readlines()
		f.close
	except IOError:
		print 'Error: read_ordered_list0: Could not open %s ' % (fileroot+'.ls')
		return []

	for line in lines:
		line=line.strip()
		words=line.split()
		if words[0][0]!='#':
			records.append(words)
	return(records)


def read_ordered_list2(fileroot='observations',dataset='first',interval=[-1,2],outroot='none'):
	'''
	Given a dataset name and an interval, return the datasets that
	were obtained in the interval around the data set.

	if outroot is anything but none, the retrieved information will also be written to
	a file called outroot+.ls

	111019	ksl	Removed lines which read entire file, and replaced with call to read_ordered_list0
	'''

	# Read the entire file
	records=read_ordered_list0(fileroot)

	# Check the dataset name
	dataset=parse_dataset_name(dataset)

	# print 'Looking for ',dataset
	# locate the record with a given dataset name
	izero=0
	while izero< len(records):
		record=records[izero]
		# print record[1],dataset
		if dataset==record[1]:
			# print 'Found %s at record %d' % (dataset,izero)
			break
		izero=izero+1

	if izero==len(records):
		print 'Error: read_ordered_list2: Could not locate record for dataset %s in %s.ls' % (dataset,fileroot)
		# print 'Using first record'
		# izero=0
		# Changed this return 101215 to trap datasets that are not in observations.ls
		return []

	zero_time=eval(records[izero][6])

	# locate the datasets within the interval
	istart=-1
	istop=-1
	i=0
	times=[]
	while i<len(records):
		time=eval(records[i][6])
		dt=(time-zero_time)*24  # convert MJD to hours
		if dt >= interval[0] and istart == -1:
			istart=i
		if istart!=-1 and dt <= interval[1]:
			times.append(dt)
			istop=i
		i=i+1

	# print 'Check lengths: ',len(times),len(records[istart:istop+1])

	xxx=records[istart:istop+1]
	i=0
	while i<len(times):
		# print '%10.1f %50s %10s %10s %10s' % (times[i],xxx[i][0],xxx[i][7],xxx[i][8],xxx[i][12])
		i=i+1


	# Now write it all to an ascii file
	if outroot!='none' and outroot != fileroot:

		# f=open(outroot+'.ls','w')
		# os.chmod(outroot+'.ls',0770)

		f=open_file(outroot+'.ls')

		for record in xxx:
			xstring='%-50s ' % record[0]
			for word in record[1:len(record)]:
				xstring=xstring+' %s' % word
			f.write('%-5s\n' % xstring)

		f.close()

	return records[istart:istop+1]


def read_ordered_list(fileroot='observations',dataset='last',delta_time=24):
	'''
	read the file written by the make_ordered_list routine, where
	fileroot is the rootname of the fileroot+'.ls' file which contains
	the all of the infomration about each dataset, datset is the
	dataset for which one is searchin, and delt_time is the time in hours
	preceding this dataset for which you would like information.

	if dataset='last', then the assumption is that one is looking for the
	very last dataset in the file.


	If the fileroot+'.ls' file is missing, or the dataset is not in the 
	file, the routine returns an empty list

	110121	ksl	Modified so that returns empty array if a dataset is request
			which does not exist 
	111019	ksl	Replaced the section that reads the entire file
	'''

	# Read the entire file
	records=read_ordered_list0(fileroot)

	if dataset!='last':
		# locate the record with a given dataset name
		ilast=0
		while ilast< len(records):
			record=records[ilast]
			if dataset==record[1]:
				# print 'Found %s at record %d' % (dataset,ilast)
				break
			ilast=ilast+1

		if ilast==len(records):
			print 'Error: read_ordered_list: Did not find dataset %s, returning' % dataset
			return []

		end_time=eval(records[ilast][6])
	else:  # Use the last record as the endpoint
		# end_time=eval(records[len(records)-1][8])
		end_time=eval(records[len(records)-1][6])
		ilast=len(records) - 1

	# print 'End time', end_time
	if delta_time>0:
		delta_time=delta_time/24.
		# Locate the first record we want
		ifirst=0
		# dt=end_time-eval(records[ifirst][8])
		dt=end_time-eval(records[ifirst][6])
		while dt > delta_time and ifirst < ilast:
			# dt=end_time-eval(records[ifirst][8])
			dt=end_time-eval(records[ifirst][6])
			ifirst=ifirst+1
		# print 'Returning %d to %d ' % (ifirst,ilast)x

		return records[ifirst:ilast+1]
	else:
		return records[ilast]



def read_ordered_list_mjd(fileroot,mjd_start=0,mjd_stop=0):
	'''
	Read an observation file and return all records between
	mjd_start and mjd_stop.

	If mjd_start=0, start from the beginning of the list.
	if mjd_stop=0,  continue to the end of the list

	101109	ksl	Coded and debugged
	'''

	records=read_ordered_list0(fileroot)
	records_out=[]
	nrec=len(records)

	print 'The total number of records is ',nrec

	i=0
	if mjd_start>0:
		while i < nrec and mjd_start > eval(records[i][6]):
			i=i+1
	ifirst=i

	if mjd_stop==0:
		ilast=nrec
	else:
		while i<nrec and mjd_stop >= eval(records[i][6]):
			i=i+1
		ilast=i

	print 'Retrieving %d records from %d to %d ' % (ilast-ifirst,ifirst,ilast)
	
	return records[ifirst:ilast]

def read_ordered_list_progid(fileroot='observations',prog_id=11216,mjd_start=0,mjd_stop=0):
	'''
	Read a list and return the records for a given program id, and if mjd_start and 
	mjd_stop are not both 0, between these two values

	If both mjd_start and mjd_stop are 0, then no tme check is done and all of the 
	records for that program will be returned.

	111019	ksl	Modified so that if prog_id is 0 or less then everything is returned
			In this case the routine behaves exactly like red_ordered_list_mjd
	'''

	records=read_ordered_list0(fileroot)
	records_out=[]
	nrec=len(records)


	if prog_id<=0:
		xrec=records
	else:
		xrec=[]
		for record in records:
			if int(record[2])==int(prog_id):
				xrec.append(record)
	
	if len(xrec)==0:
		print 'No records from prog_id %d found in %s.ls' % (prog_id,fileroot)
		return []

	if mjd_start==0 and mjd_stop==0:
		return xrec

	zrec=[]
	for record in xrec:
		xtime=eval(record[6])
		# print xtime
		if mjd_start<=xtime and xtime <= mjd_stop:
			zrec.append(record)
	
	if len(zrec)==0:
		print 'Although %d were found for prog_id %d, none between mjd %e and %e' % (len(xrec),prog_id,mjd_start,mjd_stop)
	
	return zrec




def read_ordered_list_one(fileroot='observations',dataset='last'):
	'''
	return a single record corresponding to the dataset name from 
	the ordered list whose root name is given by fileroot

	Note that if the dataset is not found, or the observaions.ls
	file is missing, the routine returns an empty list.

	101214	ksl	Added because it is needed in situations where
			one wants to work on a single dataset, but we
			need information about where everything is stored
	'''
	record=read_ordered_list(fileroot,dataset,delta_time=0)
	return record


def steer(argv):
	'''
	This is a steering program for per_list

	110103	ksl	Added as the options for per_list became more
			complicated, see top level for details
	120330	ksl	Added a fix so make_sum_file would receive
			a rootname that it could use in instances
			where the ftype was not flt

	'''

	ftype='flt'
	aperture='all'
	new_ls_file='no'
	new_summary_file='no'
	root='observations'
	np=1


	i=1

	while i<len(argv):
		if argv[i]=='-h':
			print __doc__ 
			return
		elif argv[i]=='-all':
			new_ls_file='yes'
			new_summary_file='yes'
		elif argv[i]=='-new_sum':
			new_summary_file='yes'
		elif argv[i]=='-new_ls':
			new_ls_file='yes'
		elif argv[i]=='-daily':  # This is the standard switch crontab generation
			new_ls_file='yes'
			new_summary_file='no'
		elif argv[i]=='-file_type':
			i=i+1
			ftype=argv[i]
		elif argv[i]=='-np':
			i=i+1
			np=int(argv[i]) 
		else:
			if i != len(argv)-1:
				print 'Could not understand argument %d :%s' % (i,argv[i])
				return
			else:
				root=argv[i]
		i=i+1

	# At this point we have fully parsed the observation list

	xstart=time.time()
	make_ordered_list(root,aperture,ftype,new_ls_file,np)
	dtime=time.time()-xstart
	print '# Time to make the .ls file: ',dtime

	# 120330 - kludge of a change to get to work for file types other than flt.  It's not 
	# obvious to me what would make this easier though as make ordered list has to do
	# a file search for files of a specific type.
	if ftype!='flt':
		root=root+'_'+ftype

	dtime=time.time()-xstart
	make_sum_file(root,new_summary_file)
	dtime=time.time()-xstart
	print '# Time to make the .sum file: ',dtime

	return

		


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":

	import sys

	steer(sys.argv)
	# print 'Got back here to main.  Obscure IRAF related error occurs here on some systems'
	
