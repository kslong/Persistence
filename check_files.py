#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

The purpose of this routine is to check the integrity of
the file structure assumed by the persist routines.  It
checks both the obs*ls file and the obs*sumfile


Command line usage (if any):

	usage: check_files.py [rootname]

Description:  

Primary routines:

Notes:
									   
History:

111107 ksl Coding begun

'''

import sys
import os
import per_list
import per_fits

def read_summary_file(fileroot='observations'):
	'''
	Read the summary file.  

	This routine assumes there are no comment
	lines in the summary file
	'''

	summary_file=fileroot+'.sum'

	f=open(summary_file,'r')
	summary=f.readlines()
	f.close()

	records=[]
	for record in summary:
		record=record.split()
		records.append(record)
	return(records)




def doit(fileroot='observations'):
	'''
	Main routine which locates any missing files
	and writes their names to the screen and 
	to a file called 'Missing.txt'

	Note: It is possible this should be included 
	as part of per_list
	'''
	records=per_list.read_ordered_list0(fileroot)
	sums=read_summary_file(fileroot)

	print len(records),len(sums)

	g=open('Missing.txt','w')

	i=0
	nmissing=0
	while i<len(records):
		record=records[i]
		sum=sums[i]

		# Get the flt file name, eliminating the extension that is part of the name in the .ls file
		flt_name=record[0]
		flt_name=per_fits.parse_fitsname(flt_name)
		flt_name=flt_name[0]

		persist_name='None'
		if sum[5].count('Complet'):
			try:
				persist_name=sum[12]
			except IndexError:
				print 'Incomplete obs.sum record: ',sum
		# print flt_name,persist_name

		flt='yes'
		if os.path.exists(flt_name)==False:
			flt='no'
		persist='yes'
		if persist_name!='None' and os.path.exists(persist_name)==False:
			persist='no'

		if persist=='no' or flt=='no':
			string= '%10s %50s %10s  %60s %10s' % (record[1],flt_name,flt,persist_name,persist)
			print string
			g.write('%s\n' % string)
			nmissing=nmissing+1

		i=i+1
	g.close()
	print 'The number of missing files was %d.  See Missing.txt' % nmissing
	return




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		doit('observations')
