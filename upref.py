#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

The purppse of this routine is to check whether there are calibration files that should be updated in wfc3
data before reprocessing them It presumes you have use StarView to create a csv file that contains the
current best ref files


Command line usage (if any):

		usage: upref.py  bestref.csv fitsfile

Description:  

	The routine reads a comma separated csv file produced by StarView's best ref screen,
	It then reads the header keywords of a (raw) fits file.  It uses the dataset name to 
	match a line the .csv file to the keywords and at present simply prints out a comparison.  
	It would be very easy to add a switch
	to cause an update of the file, either directly using pyfits or using pyraf.


Notes:
	The routine currently forces everthing to lower case.  Ideally StarView would return
	everything in lower case.


									   
History:

101110	ksl Coding begun
101227	ksl Updated to handle both IR and UVIS files

'''

import sys
import pyfits

cross_ref=[
		['Dataset' , 'Dataset Name'     ],
		['BPIXTAB' , 'w3r_best_bpixtab' ],
		['CCDTAB'  , 'w3r_best_ccdtab'  ],
		['ATODTAB' , 'w3r_best_atodtab' ],
		['OSCNTAB' , 'w3r_best_oscntab' ],
		['BIASFILE', 'w3r_best_biasfile'],
		['FLSHFILE', 'w3r_best_flshfile'],
		['CRREJTAB', 'w3r_best_crrejtab'],
		['SHADFILE', 'w3r_best_shadfile'],
		['DARKFILE', 'w3r_best_darkfile'],
		['NLINFILE', 'w3r_best_nlinfile'],
		['PFLTFILE', 'w3r_best_pfltfile'],
		['DFLTFILE', 'w3r_best_dfltfile'],
		['LFLTFILE', 'w3r_best_lfltfile'],
		['GRAPHTAB', 'w3r_best_graphtab'],
		['COMPTAB' , 'w3r_best_comptab' ],
		['IDCTAB'  , 'w3r_best_idctab'  ],
		['DGEOFILE', 'w3r_best_dgeofile'],
		['MDRIZTAB', 'w3r_best_mdriztab'],
		]

# Next line was header line of the csv file that came from Tim's second version of the BestRef screen
# ir_sv='"w3r_best_bpixtab"', '"w3r_best_ccdtab"', '"w3r_best_atodtab"', '"w3r_best_oscntab"', '"w3r_best_biasfile"', '"w3r_best_flshfile"', '"w3r_best_crrejtab"', '"w3r_best_shadfile"', '"w3r_best_darkfile"', '"w3r_best_nlinfile"', '"w3r_best_pfltfile"', '"w3r_best_dfltfile"', '"w3r_best_lfltfile"', '"w3r_best_graphtab"', '"w3r_best_comptab"', '"w3r_best_idctab"', '"w3r_best_dgeofile"', '"w3r_best_mdriztab"\r\n']

def read_csv_file(filename):
	'''
	Read the comma separated file produced by StarView and split lines
	which are not commented out into words.  In the process strip off 
	double quotes, etc.
	'''

	try:
		f=open(filename,'r')
	except IOError :
		print "Error: The file %s does not appear to exist!" % filename
		return []   
	xlines=f.readlines()
	
	lines=[]
	
	i=0
	while i<len(xlines):
		z=xlines[i].split(',')
		j=0
		while j<len(z):
			z[j]=z[j].strip()
			z[j]=z[j].strip('"')
			j=j+1
		if len(z)>0:
			if z[0][0]!='#':
				lines=lines+[z]
		i=i+1
	return lines

def get_bestref_values(dataset,bestref_file):
	'''
	Associate keywords with those in records and return the values of the keywords  for a single
	dataset.

	Here best_ref_file is the name of the .csv file produced by StarView, and dataset is
	the dataset name.
	'''

	# Read the csv file produced by StarView best ref
	records=read_csv_file(bestref_file)

	if len(records)==0:
		return []

	# Index the results of StarView to the kewword values we need

	index=[]
	for one in cross_ref:
		# print one
		header=records[0]
		j=0
		while j<len(header):
			if one[1]==header[j]:
				# print 'adding',j
				index.append(j)
				break
			j=j+1
		if j==len(header):
			print 'Could not associate ',one
			index.append(-1)
	
	# print index

	values=[]


	i=1
	xdataset=dataset.lower()
	while i < len(records):
		record=records[i]
		zdataset=record[index[0]].lower()
		# print xdataset,zdataset
		if xdataset==zdataset:
			j=0
			while j<len(index):
				values.append(record[index[j]])
				j=j+1
			break

		i=i+1

	if i==len(records):
		print 'Could not locate dataset %s in %s' % (dataset,bestref_file)
		return []

	return values

def get_current_values(fits_file):
	'''
	Get all of the data processing keyword values for a dataset
	using pyfits

	Notes:

	The dataset name is not necessarily the rootname of a file.  It is
	insteand usually, but not always the assoc_id.  Only if the assoc_id
	is NONE should you use the root name or alternatively the exposure
	name


	101227	ksl	Fix problem with getting the dataset name from
			the fits files.
	'''

	try:
		z=pyfits.open(fits_file)
	except IOError:
		print 'Error get_current_values: file %s not found' % fits_file
		return []

	ext=z[0]  # Choose extension 0

	# print 'gotcha',ext.header['ASN_ID'],ext.header['ROOTNAME']

	try:
		assn=ext.header['ASN_ID']
		rootname=ext.header['ROOTNAME']
	except KeyError: 
		print 'Error: Could not identify one of either ASN_ID or ROOTNAME in %s' % fits_file
		return []

	if assn=='NONE':
		dataset=rootname
	else:
		dataset=assn

	# print 'Got dataset  xx%sxx' % dataset



	values=[dataset]

	i=1
	while i<len(cross_ref):
		one=cross_ref[i]
		try:
			answer=ext.header[one[0]]
			values.append(answer)
		except KeyError:
			# print 'Keyword %s not found in %s' % (one[0],fits_file)
			values.append('Unknown')
		i=i+1

	return values




def do_one(bestref='new.csv',fits_file='ibcj04jbq_flt.fits'):
	'''
	Read and compare the current and best values for the 
	calibration keywords for a single fits file

	101227	ksl	Added more error checking and the ability to
			create an hedit file
	'''

	# Open both the txt output file and the iraf command file
	g=open('upref.txt','a')
	h=open('upref.cl','a')

	g.write('# Check best reference files for %s\n' % fits_file)
	h.write('# Updates for file %s\n' % fits_file)
	

	current=get_current_values(fits_file)

	if len(current)==0:
		string= 'Error: could not get information for fits file  %s' % fits_file
		print string
		g.write('%s\n' % string)
		h.write('# !!! %s\n' % string)
		return

	dataset=current[0]
	best=get_bestref_values(dataset,bestref)

	if len(best)==0:
		string= 'Error: could not get best_ref information from %s' % bestref 
		print string
		g.write('%s\n' % string)
		h.write('# !!! %s\n' % string)
		return


	if len(current)!=len(best):
		string= 'Error: current (%d) and best (%d)  do not have the same number of elements ' % (len(current),len(best))
		print string
		g.write('%s\n' % string)
		h.write('# !!! %s\n' % string)
		return

	# At this point we should have the correct information to process this file

	i=0
	while i<len(current):
		best[i]=best[i].lower()
		current[i]=current[i].lower()
		if current[i]!='unknown'and best[i]!='n/a':
			status='OK'
			if current[i].count(best[i]):
				status='OK'
				format = '%20s %20s %30s  %30s OK'
			else:
				status='X'

			# Write the results to the upref.txt file
			string = '%20s %20s %30s  %30s %s' % (cross_ref[i][0],cross_ref[i][1],current[i],best[i],status)
			print string
			g.write('%s\n' % string)

			# Add to the upref.cl file if there was a change
			if status=='X':
				# Get the logical directory name where the files should be stored
				try:
					k=current[i].rindex('$')
					prefix=current[i][0:k+1]
				except ValueError:
					prefix=''
				# Now write out the hedit command
				h.write('hedit ("%s[0]", "%s", "%s%s", add=no, addonly=no, delete=no, verify=no, show=yes, update=yes)\n' % (fits_file,cross_ref[i][0],prefix,best[i]))



		i=i+1
	g.close()
	h.close()



# Next lines permit one to run the routine from the command line
helpstring='''
		usage: upref.py bestref.csv fitsfile ...

where the first argument refers to a comma separated file of the best reference files obtained 
from the StarView  (See http://starview.stsci.edu/web/# and the WFC3 Best Reference Files form)
and the remaining arguments are the names of fits files (usually raw fits files) containg data
obtained with wfc3.  

The routine will check the files designated by the calibration file keywords in the fits files 
against the currently recommended calibration files contained in the bestref.csv file.

The routine prints out the results of its comparison to the screen and to a text file upref.txt.
It also creates a file upref.cl that can be used to update the headers from within iraf or
pyraf with the following command: cl <upref.cl

'''
if __name__ == "__main__":
	import sys
	if len(sys.argv)>=3:
		i=2
		while i<len(sys.argv):
			do_one(sys.argv[1],sys.argv[i])
			i=i+1
	else:
		print helpstring

