#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is a routine to handle the control permissions.
for the persistence software


Description:  

	set(name)  

	The routine has no effect if one is not 
	operating on files in the Quicklook and
	QL_GO file hierarchy.  

	If you are in this file structure though
	it attempts set the group name to either
		wfc3   if in Quicklook or
		wfc3ql if in QL_GO

	If the name refers to a directory then the
	permissons are set to 02770 (which is intended
	to make the group name inherited and to allow
	other users to be able to delete files in
	directories below.

	If the name refers to a file then the
	permissions for the file are set to 0770

Primary routines:

Notes:
	gr_ids are hardwired.  They can be obtained using the 
	unix command id.  If these change, the program will 
	break

	The program will also fail if the Quicklook2 data
	structure is changed.

	This routine  needs to work withing the production environment, and
	it also needs to work for the average user so that we don't
	have to rewrite the soltware if it is exported.
									   
History:

110811 ksl Coding begun

'''

import sys
import os

def open_permiss(name):
	'''
	Open a file and set its permissions
	in a robust fashion
	'''

def set(name):
	'''
	Set the permissions for a directory or file

	The routine assumes that QL_GO are to be
	restricted to the wfc3ql group, wherease
	the routines in QuickLook can seen by
	the wfc3 group as a whole

	'''

	# First verify that this is a directory, and if not return

	if os.path.isdir(name)==True:
		isdir=True 
	elif os.path.isfile(name)==True:
		isdir=False
	else:
		print 'Error: set: %s does not exist' % name     
		return False

	# Next get the absolute path name of the file or directory

	name=os.path.abspath(name)


	# Now check whether this is a Quicklook or QL_Go directory
	# If yes then set the group name of the directory to wfc3 
	# or wfc3_ql as appropriate

	if name.count('QL_GO'):
		group='wfc3ql'
		os.chown(name,-1,6047)
	elif name.count('Quicklook'):
		group='wfc3'
		os.chown(name,-1,340)
	else:
		group='unknown'
		return True

	# Now set the permission of the directory or file
	if isdir:
		os.chmod(name,02770)
	else:
		os.chmod(name,0770)

	return True
