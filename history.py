#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Standard routine for logging a string to a file

Command line usage (if any):

	None.  Normally this simple module should be
	import as 
		from history import log

	

Description:  
	This is a simple set of logging routines that
	may grow into something more general over time

Primary routines:

	open_log(filename)
	log(string,filename='history.log',option='a'):

Notes:
	130909:  As written one needs to be careful about
		how this is used since you can only have
		one history file open at the same time. 
									   
History:

120410 ksl Coding begun

'''

import sys

current_name='history.log'


def open_log(filename,option='w'):
	'''
	Set the hame of the history file
	'''
	global current_name

	if filename!='':
		current_name=filename
	log('# %s\n' % filename,option=option)
	return

def log(string,filename='',option='a'):
	'''
	Log a string  to a file.  

	Once you have used log for the
	first time.  It will always
	write to the same file unless
	you give it a new name.
	'''
	global current_name

	# Change the name of the current_filename if 
	# the current_name has not been initialized
	# or the filename is different than before

	if filename!='' and filename!=current_name:
		current_name=filename


	history=open(current_name,option)
	history.write(string)
	history.close()

