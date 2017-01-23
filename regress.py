#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Carry out a regression test on the persistence software, assuming
at most one regression test per day will be carried out.


Command line usage (if any):

	usage: regress.py filename

Description:  

Primary routines:

Notes:
									   
History:

140922 ksl Coding begun

'''

import sys
import time
import os
import subprocess

def run_command(one_command='ls -l'):
	'''
	run a command and trap the outputs
	'''

	proc=subprocess.Popen(one_command,shell=True,stdout=subprocess.PIPE)
	x=proc.communicate()[0]
	print(x)

	return 

def doit(dirname='foo',data_dir='D'):
	'''

	Setup and run a regression test

	140924	ksl	Modified so that the files one wants to be
			tested should be in a directory D, unless data_dir
			is set to something else

	'''

	if dirname=='':
		x=time.strftime('%y%m%d')
		dirname=x
	
	print('Running regression test in ',dirname)

	# Create a directory for running the regression test

	try:
		os.stat(dirname)
		print('Directory %s already exists, removing and restarting' % dirname)
		run_command('rm -r %s' % dirname)
	except OSError: 
		print('Created directory %s ' % dirname)

	os.mkdir(dirname)
	
	os.chdir(dirname)

	run_command('rm D')
	run_command('ln -s ../%s D' % data_dir)
	run_command('ln -s ~/WFC3/Calfiles/persist_xcorr.fits  persist_corr.fits')

	run_command('per_list.py')

	run_command('run_persist.py -all -local')

	run_command('subtract_sum.py -all')

	







# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		doit('')
