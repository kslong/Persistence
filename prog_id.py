#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This routine creates html files  that allow one to explore the html outputs
associated with persistence by program ID.  

A summary page is produced for all the programs, and this links to an
individual page for each program.  The individual page has links to the
pages for the individual data sets page

It is intended to produce a set of html files that can be viewed by the community


Command line usage (if any):

	usage: prog_id.py [fileroot] 

		where fileroot is the rootname of the file that is needed by 
		run_persist.py, normally observations.ls

	 	If fileroot is omitted 'observations' is assumed

Description:  

	The routine obtains all its information the observations.ls file and 
	the observations.sum file.  The observations.ls and observations.sum
	file are produced by per_list.py and the observations.sum file is 
	updated as run_persist.py is run.

	The routine reads these two files and picks out all of the program 
	Ids.  For each progam ID it creates an html file that contains information
	about persistence for the various datasets in that program.

	It also creates a Summary.html file that is just a listing of all the
	program IDs and PIs, with links to the individual html files for each
	program ID


Primary routines:

Notes:
									   
History:

110722 ksl Coding begun
110811  ksl     Switched to standardized way to set file permissions

'''

import sys
import html
import date
import os
import per_list
import subtract_sum
import permissions

def prog_html(sum,records,filename):
	'''
	Make an html file for a specfic program id, having
	already read the summary file and ls file

	111003	ksl	Fixed small error in the title 
	'''

	# first get the information we would like to put in the table

	i=0

	title=[['Rootname','Visit.Line','Obs. Date', 'Obs. Time','Target','Date_Proc','Status','Ext1','Ext2','Ext3','tot1','Tot2','Tot3']]
	lines=title
	title=[['','','YYYY-MM-DD', 'HH:MM:SS','','YYYY-MM-DD','','%>0.1 e/s','%>0.03 e/s','%>0.01 e/s','%>0.1 e/s','%>0.03 e/s','%>0.01 e/s']]
	lines=lines+title
	while i<len(sum):
		one_sum=sum[i]
		one_record=records[i]
		root=one_record[1]
		visit=one_record[3]
		obsdate=one_record[7]
		obstime=one_record[8]
		targ=one_record[14]
		xline=[root,visit,obsdate,obstime,targ]
		xline=xline+[one_sum[3]]+one_sum[5:12]
		lines.append(xline)
		i=i+1
	
	page=html.begin('Summary Page: Persistence in  program %s with PI %s' % (records[0][2],records[0][16]))

	string='''The table belows gives and indication of whether persistence is likely to be a problem in any of the datasets
	obtained to date in this program.  The columns of the table are as follows:'''

	page=page+html.paragraph(string)

	comment=[]
	comment.append('Rootname: Rootname of the specific exposure')
	comment.append('Visit.Line: Visit and line number in phase II proposal for this exposure')
	comment.append('Obs. Date: Observation date for this exposure')
	comment.append('Obs. Time: Start time for this exposure')
	comment.append('Target: Target of this exposure')
	comment.append('Date Proc: Date on which persistence processing was carried out')
	comment.append('Status: Status message about the persistence processing.  Normally Complete with version number')
	comment.append('Ext1: Percentage of pixels which external persistence is estimated to be greater than 0.1 e/s')
	comment.append('Ext2: Percentage of pixels which external persistence is estimated to be greater than 0.03 e/s')
	comment.append('Ext3: Percentage of pixels which external persistence is estimated to be greater than 0.01 e/s')
	comment.append('Tot1: Percentage of pixels which total persistence is estimated to be greater than 0.1 e/s')
	comment.append('Tot2: Percentage of pixels which total persistence is estimated to be greater than 0.03 e/s')
	comment.append('Tot3: Percentage of pixels which total persistence is estimated to be greater than 0.01 e/s')

	page=page+html.add_list(comment)

	string='''Note: External persistence refers to persistence caused by visits that preceed this dataset. This is usually the type of
	persistence that causes the most trouble in terms of data analysis since it can appear anywhere in the image.  Internal persistence, in our terminolgy is
	persistence that is due to earlier exposures in the same visit.  This is usually, though not always, less of a worry than external
	persistence because unless the patterns and dithers that were used in the visit will be smalll, and the persistence will mainly
	occur near bright objects.  If one has used large dithers or has multiple points within a single visit, one needs to take that into account
	when evaluating problems associated with persistence.'''

	page=page+html.paragraph(string)
	page=page+html.table(lines)
	page=page+html.end()

	g=open(filename,'w')
	g.write(page)
	g.close()





	

def doit(fileroot='observations'):
	'''
	Do something useful
	'''
	records=per_list.read_ordered_list0(fileroot)

	Dir='./Summary'

	progs=[]
	for record in records:
		progs.append(record[2])
	
	progs=set(progs)
	progs=sorted(progs)


	
	if os.path.exists(Dir)==False:
		try:
			os.mkdir(Dir)
			permissions.set(Dir)
		except OSError:
			print('!NOK Could not create %s' % (Dir))

	# Create each of the individual html files

	pi=[]
	for prog in progs:
		print(prog)
		sum=subtract_sum.read_sum_file(fileroot,'All',prog)
		records=per_list.read_ordered_list_progid(fileroot,prog)
		pi.append(records[0][16])

		for one in sum:
			print(one)
			# subtract_sum.make_html(sum,Dir+'/prog_%s.html'% prog)
			prog_html(sum,records,Dir+'/prog_%s.html'% prog)
		

	# Now create the master summary file

	page=html.begin('Persistence Summary by Program ID')
	string='''This page is intended to help those who have IR images evaluate whether they need to be 
	concerned about persistence.  It provides a link to each program with IR images.  Simply follow the
	link to specific infomation about your program'''
	page=page+html.h3(string)
	# Now create a list of lines, we will put in a list

	sum_lines=[]
	i=0
	while i< len(progs):
		prog=progs[i]
		one_pi=pi[i]
		xstring=html.link(prog,'prog_%s.html'% prog)+'-- %s' % one_pi  
		sum_lines.append(xstring)
		i=i+1
	
	page=page+html.add_list(sum_lines)

	page=page+html.end()

	g=open(Dir+'/Summary.html','w')
	g.write(page)


	

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		doit('observations')
