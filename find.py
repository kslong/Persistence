#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is a simple routine to locate a specific version of a file.  If
the file is in the current directory that is the file name that
will be returned.  If it is in one of the subdirectories of the
current directory then the one that was modified most recently
will be returned. 


Command line usage (if any):

		usage: find.py  filename'

Description:  

Primary routines:

Notes:
									   
History:

100906 ksl Coding begun

'''

import os
import sys
import subprocess



def find_best(filename='foo'):
	'''
	Find the best version of a file, defined to be the one in
	the current directory or the most recent one in any of the 
	subdirectories.
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




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		x=find_best(sys.argv[1])
		print x
	else:
		print 'usage: find.py  filename'

