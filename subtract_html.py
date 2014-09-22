#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is a routine to that creates one or more html files
that summarize the results of the persistence subtraction


Command line usage (if any):

		usage: subtract_html.py  dataset

Description:  

Primary routines:

Notes:
									   
History:

100903	ksl	Coding begun
110505	ksl	Switched from markup to my own html module

'''

import os
import sys
import glob
import per_fits
import numpy
import pylab
import time
import per_list
import html

def read_history(history_file):
	'''
	Read and get the most important lines of the history file.
	'''

	f=open(history_file,'r')
	
	line=f.readline()
	keep=[]
	table_stimulus=[['Filename','Program Id','Visit','Target','dt (sec)','# Saturated Pixels','Aperture','Scan']]
	table_sum=[['Type of Persistence','% Pix > 0.1 e/s','% Pix >0.03 e/s','% Pix >0.01 e/s']]
	while line !='':
		line=line.strip()
		if len(line)>0 and line[0]=='!':
			line=line.strip()
			line=line[2:len(line)]
			words=line.split()
			if words[0]=='Stimulus:':
				table_stimulus.append(words[1:len(words)])
			elif words[0]=='PersistenceSum:':
				table_sum.append(words[1:len(words)])
			else:
				keep.append(line)
		line=f.readline()
	return keep,table_stimulus,table_sum

def make_html(coords=[[100,200],[300,400]],fileroot='ib2v09kzq'):
	'''
	Make an html file that gathers all of the images that have been made
	'''

	# page=markup.page()

	title='Persistence Removal Evaluation for %s' % fileroot
	page=html.begin(title)

	# page.init(title='Persistence Removal Evaluation for %s' % fileroot)
	# page.h1('Persistence Removal Evaluation for %s' % fileroot)
	# page.p('''This page contains images for the evaluation of how well persistence has been removed from an image''')
	page=page+html.paragraph('''This page contains images for the evaluation of how well persistence has been removed from an image''')

	# page.hr(size='3',width='100%')
	page=page+html.hline(size='3',width='100')

	for coord in coords:

		# page.p('Images for positions: %3d  %3d' % (coord[0],coord[1]))
		page=page+html.paragraph('Images for positions: %3d  %3d' % (coord[0],coord[1]))
		fig1='Figs/Fig_%s_%04d_%04d_1.png' % (fileroot,coord[0],coord[1])
		fig2='Figs/Fig_%s_%04d_%04d_2.png' % (fileroot,coord[0],coord[1])
		fig3='Figs/Fig_%s_%04d_%04d_3.png' % (fileroot,coord[0],coord[1])

		# page.img( src=fig2, width=900, height=300, alt="Thumbnails" )

		page=page+html.image(f2,width=900, height=300, alt="Thumbnails" )

		# page.p('Left: Original, Center: Persistence model, Right: Subtracted')
		page=page+html.paragraph('Left: Original, Center: Persistence model, Right: Subtracted')

		# page.img( src=fig1, width=400, height=400, alt="Thumbnails" )
		page=page+html.image(fig1,width=400, height=400, alt="Thumbnails" )
		# page.img( src=fig3, width=400, height=400, alt="Thumbnails" )
		page=page+html.image(fig3,width=400, height=400, alt="Thumbnails" )

		# page.p('Left: Original and subtracted data as a function of the estimated persistence, Right: Original and subtracted data as a function of distance from x,y')
		page=page+html.paragraph('Left: Original and subtracted data as a function of the estimated persistence, Right: Original and subtracted data as a function of distance from x,y')
		# page.hr(size='3',width='100%')

		page=page+hline(size='3',width='100%')

	
	# Write the page to a file
	name=xpath+'Persist_%s.html' % fileroot


	g=per_list.open_file(name)
	# g=open(name,'w')
	# os.chmod(name,0770)

	g.write('%s' % page)
	g.close()
	return fileroot



def do_dataset(dataset='ia21h2eaq',fileroot='observations',local='no'):
	'''
	Make html files for a single dataset

	110203	ksl	Added local swithch which controls where the
			real working directory is to make testing
			easier
	140307	ksl	Added information about scans and subarray observations
	'''

	record=per_list.read_ordered_list_one(fileroot,dataset)
	if len(record)==0:
		return 'NOK: make_html failed becaouse could not find dataset %s' % dataset

	work_dir=per_list.set_path(record[0],'no',local)  # This will be the Persist directory for the dataset
        fig_dir=work_dir+'/Figs/'               # This will be the directory where figures are stored

	html_filename=work_dir+dataset+'_persist.html'

	# page=markup.page()
	title='Persistence Removal Evaluation for dataset %s' % dataset
	page=html.begin(title)

	# page.init(title='Persistence Removal Evaluation for dataset %s' % dataset)
	# page.h1('Persistence Removal Evaluation for %s' % dataset)

	# page.p('''This page contains images for the evaluation of how well persistence has been removed from an image''')
	page=page+html.paragraph('''This page contains images for the evaluation of how well persistence has been removed from an image''')

	# Look for the history file for this dataset

	history_file=dataset+'.txt'
	
	if os.path.exists(work_dir+history_file):
		string='''The history file for the processing of this dataset is '''
		string=string+html.link("here",href=history_file)
		page=page+html.paragraph(string)

		# read history simply returns all of the lines in the history file that begin with !
		# And so any processing of these lines still has to be done
		lines,table1,table2=read_history(work_dir+history_file)
		for line in lines:
			page=page+html.paragraph(line)
		if len(table1)>0:
			page=page+html.h2('Earlier exposures that could affect this image')
			page=page+html.table(table1)
		if len(table2)>0:
			page=page+html.h2('External and total persistence for this image')
			string='''External persistence is persistance from previous visits; internal persistence
			is persistence induced from exposures in this vist.  Total persistence includes both
			internal and external persistence.  . Generally, self-induced or internal persistence is  
			only important if the dithers larger than the psf have been used within the visit'''
			page=page+html.paragraph(string)
			page=page+html.table(table2)

	else:
		page=page+html.paragraph(''' The history file for this dataset appears to be missing.  Check that the file has been processed''')


	page=page+html.hline(size='3',width='100')

	string='''The next 4-panel image shows the original flt image (upper left), the corrected flt image (upper right), 
	the persistence model (lower left) and the stimulus (lower right).  The stimulus is simply the image constructed
	maximum value in electrons of any of the images that went into the stimulus model'''


	# Look for the summary image

	xname=dataset+'_subtract.png'
	if os.path.exists(fig_dir+xname):
		# page.img(src='Figs/'+xname,width=600,height=600,alt="Thumbnails")
		page=page+html.image(image='Figs/'+xname,width=600,height=600,alt="Thumbnails")
	else:
		# page.p('''The summary image is missing''')
		page=page+html.paragraph('''The summary image is missing''')


	# page.hr(size='3',width='100%')
	page=page+html.hline(size='3',width='100')

	# Now include the evaluation images

	string='''As a qualitative indicator of how well the persistence correction has worked, some of the regions with
	the highest predicted persistence have been examined. 
	The next two images give an indication of how well the persistence has been subtracted from the images.
	Both images have the original data in red and the persistence-subtracted data in blue.  The first image is
	a plot of flux vs the persisence model, the second is flux as a function of the stimulus. Ideally the blue 
	curves would all center around 0. The utility of these plots depends on how isolated the presistence peaks
	are from stars in the image. If these plots are empty, no good regions for evaluation persistence were found.'''

	page=page+html.paragraph(string)

	xname=dataset+'.sum1.png'
	if os.path.exists(fig_dir+xname):
		# page.img(src='Figs/'+xname,width=300,height=300,alt="Thumbnails")
		page=page+html.image('Figs/'+xname,width=300,height=300,alt="Thumbnails")
	else:
		# page.p('''The first evaluation image showing the subtraction is missing''')
		page=page+'''The first evaluation image showing the subtraction is missing'''




	xname=dataset+'.sum2.png'
	if os.path.exists(fig_dir+xname):
		# page.img(src='Figs/'+xname,width=300,height=300,alt="Thumbnails")
		page=page+html.image('Figs/'+xname,width=300,height=300,alt="Thumbnails")
	else:
		# page.p('''The second evaluation image showing the subtraction is missing''')
		page=page+html.paragraph('''The second evaluation image showing the subtraction is missing''')



	# page.hr(size='3',width='100%')
	page=page+html.hline(size=3,width=100)

	# Look for the peaks summary

	string='''This figures indicates what regions were selected for evaluation. The two panels are
	identical except the regions selected are indicated in the lower panel. '''

	page=page+html.paragraph(string)

	xname=dataset+'_persist.peaks.png'
	if os.path.exists(fig_dir+xname):
		# page.img(src='Figs/'+xname,width=600,height=1000,alt="Thumbnails")
		page=page+html.image('Figs/'+xname,width=900,height=900,alt="Thumbnails")
	else:
		# page.p('''The summary figure for peak identification is missing''')
		page=page+html.paragraph('''The summary figure for peak identification is missing''')



	# Now find all of the individual peak files:

	searchstring=fig_dir+dataset+'.peak.*.1.png'
	print searchstring

	try:
		peaks_file=work_dir+dataset+'_persist.peaks.dat'
		p=open(peaks_file,'r')
		lines=p.readlines()
		p.close
	except IOError:
		print 'Warning: %s not found' % peaks_file
		lines=[]
	
	

	xlines=[]
	for one in lines:
		one=one.strip()
		if one[0]!='#'and len(one)>0:
			xlines.append(one)

	if len(xlines)>0:
		string='''The results for individual regions are shown below. The four panels are a subsection of the original flt file, the predicted persistence in that region, the persistence subtracted flt file, and a plot of pixel values as a function of predicted persistence in the region. Green points are the original values; yellow point are the corrected values. The red and blue lines show the mean values in the original and corrected and corrected images, respectively.'''
		page=page+html.paragraph(string)
		page=page+html.hline(size='3',width='100')

		for one in xlines:
			word=one.split()
			x=int(word[0])
			y=int(word[1])
			z=eval(word[2])
			zz=eval(word[3])
			# page.p('Persistence at x = %3d, y=%3d' %(x,y))
			page=page+html.paragraph('Persistence at x = %3d, y=%3d is about %6.3f e/s compared to science image flux of %6.3f e/s' %(x,y,z,zz))
			xname='%s.peak.%03d_%03d.1.png' % (dataset,x,y)
			if os.path.exists(fig_dir+xname):
				# page.img(src='Figs/'+xname,width=400,height=400,alt="Thumbnails")
				page=page+html.image('Figs/'+xname,width=400,height=400,alt="Thumbnails")
			else:
				# page.p('Figure %s not present' % (work_dir+xname))
				page=page+html.paragraph('Figure %s not present' % (work_dir+xname))
			# page.hr(size='3',width='100%')
			page=page+html.hline(size='3',width='100')
	else:
		string='''Unfortunately, no good regions for evaluating persistence were found.'''
		page=page+html.paragraph(string)
		page=page+html.hline(size='3',width='100')


	
		

	page=page+html.end()


	# Open the html file with the appropriate permissions, and then write it
	g=per_list.open_file(html_filename)
	g.write('%s' % page)
	g.close()


	return 'OK: subtract_html: %s' % html_filename 




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		do_dataset(sys.argv[1])
	else:
		print 'usage: subtract_html.py  dataset  '
