#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

These are routines which are specifically intended to be used for analyzing
the Persist II calibraton data.  A lot of the routines seem to be only partially
completed.


Description:  

Primary routines:

	doit - Reads flt files from a visit, uses imstat to get statistics on the image, and then performs a power
		law fit to the data

Notes:

	100831 - Its not obvious that anythin here deals with the problem of
	what to do with ima files, despite having one of the programs looking
	as if it would do this.  It may be that the routine for that was not
	finished.
									   
History:

100513	ksl	Coding begun
100831	ksl	Removed routines that are in per_fits; per_list had been removed
		previously

'''

import sys
import pyraf
import pyfits
import pylab
import numpy
from per_list import *
from per_fits import *

from scipy import optimize

def plot_rows(filename,ext=1,nrow=500,dnrow=50,new='yes'):
	'''
	plot the median value of a number of rows in an image as a function of
	column.  

	The image is defined by a filename and extension (ext).  The extension for a simple fits file would be 0
	The center of the region plottted is definded by nrow and the number of rows to consider by dnrow
	If new is 'yes', then a new plot is created.  Otherwise calls will be overwritten.
	

	101012	Coded and debugged
	'''

	x=get_image_ext(filename,ext)

	nmin=int(nrow-0.5*dnrow)
	nmax=int(nrow+0.5*dnrow)
	if nmin<0:
		nmin=0
	if nmax>len(x):
		nmax=len(x)
	
	ncol=len(x[nmin])

	print nmin,nmax

	xx=[]
	yy=[]
	n=0
	while n<ncol:
		zz=x[nmin:nmax,n]
		z=numpy.median(zz)
		xx.append(n)
		yy.append(z)
		n=n+1
	
	# print xx
	# print yy

	pylab.figure(1,(6,6))
	if new=='yes':
		pylab.clf()

	pylab.plot(xx,yy,'-')

	lims=pylab.axis()
	pylab.axis((0,1024,0,lims[3]))
	print lims

	pylab.xlabel('Column')
	pylab.ylabel('Rate (e/s)')
	
	pylab.draw()

def make_fig_rows(row=500):
	'''
	This is essentially a hard-coded figure for the ISR
	associated with the Cycle 17 data

	101012 - ksl
	'''
	plot_rows('ibel04l8q_flt.fits',1,row,100)
	plot_rows('ibel03e4q_flt.fits',1,row,100,new='no')
	plot_rows('ibel11l7q_flt.fits',1,row,100,new='no')
	plot_rows('ibel02g8q_flt.fits',1,row,100,new='no')
	plot_rows('foo.fits',0,row,100,new='no')
	pylab.savefig('fig_rows.png')

def plot_pixel_dist(filename,ext=1,dq=0):
	'''
	plot the pixel distribution in an image, eliminating badpixels if the dq is non-zero,
	and also print out some statatistics associated with the image, namely the mean value
	and the 97.5% limits.

	Notes 

	101011  - Written

	This was created in an attempt to better describe variations in the amount of the persistence
	in an image, although it will do this for any image.   Once written, it was not obvious that
	it was extremely useful, in part because readout noise is an issue when talking about persistence
	of individual pixels.

	The plot needs better limits to make the routine useful
	
	'''

	# Check the extensions'
	type1=get_ext_type(filename,ext)

	x=get_image_ext(filename,ext)
	if dq!=0:
		type1=get_ext_type(filename,dq)
		if type1 == 'DQ':
			dq=get_image_ext(filename,dq)
			# print dq
			xx=numpy.ma.array(x,mask=dq)
			# print 'test',len(xx)
			xxx=xx.compressed()
			# print 'test',len(xxx)
			# print xxx
		else:
			print 'plot_pixel_dist: Extension %d in %s is not DQ, but %s' % (ext,filename,type1)
			xxx=x
	else:
		xxx=x
	
	x=numpy.ravel(xxx)
	xx=numpy.sort(x)

	# Get some statistics associated with the image
	xmin=xx[0]
	xmax=xx[len(xx)-1]

	nlow=int(0.025*len(xx))
	nmid=int(0.5*len(xx))
	nhigh=int(0.975*len(xx))

	print 'There are %d points in the image' % len(xx)
	print 'The min and maximum valuese are:' ,xmin,xmax
	print '0.025, median, and 0.975 vlaues are: ',xx[nlow],xx[nmid],xx[nhigh]

	# Next section is commented out since it is a way to generate histograms, but this is done also within
	# pylab.hist
	# z=numpy.linspace(xmin,xmax,1000)
	# y,z=numpy.histogram(xx,bins=1000)


	# Set a reasonable range to plot

	delta=2*(xx[nhigh]-xx[nlow])
	zrange=(xx[nmid]-delta,xx[nmid]+delta)

	pylab.figure(1,(6,6))
	pylab.clf()

	# pylab.hist(xx,bins=10000,range=(-1,5),histtype='step',log=True,normed=True)
	# pylab.hist(xx,bins=10000,range=(-1,5),histtype='step',cumulative=True,normed=True)

	pylab.hist(xx,bins=10000,range=zrange,histtype='step',log=True,normed=True)
	pylab.hist(xx,bins=10000,range=zrange,histtype='step',cumulative=True,normed=True)

	pylab.draw()


	
	




def get_sci_stats(filename):
	'''
	Get the samptime and mean all of the SCI extensions in a image, subtracting
	the bias values if the samptime of the last SCI extension is 0.

	The routine is intended to work on both raw and ima files

	This routine use iraf.imstat and so does not account for bad pixels, The units
	are those of the files read.  That is there is not an attempt to convert to
	a common set of units.  

	The routine returns a list, in which each element has a time and the mean value.
	Time is the time in the header

	100901	ksl	Coded

	'''

	ext=get_ext_info(filename,'SCI')

	if len(ext)==0:
		print 'There were no SCI extensions in %s' % filename
		return 

	# Read the extensions in reverse order, that is,  so that they will be in time order
	summary=[]
	i=len(ext)-1
	while i>=0:
		one_ext=ext[i]
		complete_name=filename+'[%s]' % one_ext[0]
		zzz=pyraf.iraf.imstatistics(complete_name,fields = "npix,mean,stddev,min,max,image", lower='INDEF', upper='INDEF',nclip=2, lsigma=3.0, usigma=3.0, binwidth=0.1,format='no',Stdout=1)
		zzzz=zzz[0].split()
		# print one_ext[1],zzzz
		summary.append([one_ext[1],float(zzzz[1])])
		i=i-1
	# print summary
	if summary[0][0]!=0:
		print 'For some reason the initial value does not look like a bias, returning all times and means'
		return summary
	bias=summary[0][1]
	xsummary=[]
	i=1
	while i<len(summary):
		xsummary.append([summary[i][0],summary[i][1]-bias])
		i=i+1
	# print xsummary
	return xsummary


def get_rates(filename):
	'''
	Get the time and mean rates in the SCI extensions of an ima or raw image
	The time returned is the middle of the interval.  

	The routine returns 2 one-d arrays, one containing times, and the other containing rates
	
	The extensions are differenced so that the rates are the rates measrured in two subsequent extensions

	This does not eliminate bad pixels, since get_sci_stats uses imstat (at present).


	'''

	z=get_sci_stats(filename)

	print 'stats ', z

	t=[]
	rate=[]

	if filename.count('ima') > 0:
		t.append(z[0][0]/2.)
		rate.append(z[0][1])
		i=1
		while i<len(z):
			xrate=z[i][0]*z[i][1]-z[i-1][0]*z[i-1][1]
			xrate=xrate/(z[i][0]-z[i-1][0])
			tt=0.5*(z[i][0]+z[i-1][0])
			rate.append(xrate)
			t.append(tt)
			i=i+1
	elif  filename.count('raw')>0:
		t.append(z[0][0]/2.)
		rate.append(z[0][1]/z[0][0])
		i=1
		while i<len(z):
			xrate=z[i][1]-z[i-1][1]
			xrate=xrate/(z[i][0]-z[i-1][0])
			tt=0.5*(z[i][0]+z[i-1][0])
			rate.append(xrate)
			t.append(tt)
			i=i+1
	else:
		print 'Do not know how to operate on file %s ',filename
	
	print 'Results'
	print t
	print rate
	return t,rate


def plot_rate(filename,new='yes'):
	'''
	Plot rates using imstats for each SCI extension within a file

	This uses get_rates above
	'''

	x,y=get_rates(filename)



	pylab.figure(1,figsize=(6,6))
	if new=='yes':
		pylab.clf()



	# pylab.loglog(x,y,'o',label='test')
	pylab.plot(x,y,'o',label='test')

	pylab.xlabel('Time (s) ')
	pylab.ylabel('Afterglow (e/s)')

	pylab.savefig('Persistence_sum.png')
	pylab.draw()

	return

def plot_rates(fileroot='observations_raw',bright='ibel04l3q'):
	'''
	Plot the rates as established from the individual SCI extensions
	of a series of files.   It was intended for looking at ima or raw
	files.

	The routine use get_rates, so the rates are the rates beween
	two subsequent science extensions.  Also, since, get_rates uses imstat
	bad pixels are not explicitly excluded.
	'''
	pylab.figure(1,figsize=(6,6))
	pylab.clf()

	records=read_ordered_list2(fileroot,dataset=bright,interval=[0,2],outroot='none')

	bright=records[0][0]
	bright_time=get_times(bright)
	bright_time=bright_time[1]

	i=1
	while i<len(records):
		record=records[i]
		xname=parse_fitsname(record[0])
		print xname
		tt=get_times(xname[2]) # Assumption hre is that 1 is a science record
		offset=(tt[0]-bright_time)*86400
		print 'xxxx',offset,tt,bright_time
		times,rates=get_rates(xname[0])
		xtimes=[]
		for time in times:
			xtimes.append(time+offset)

		pylab.loglog(xtimes,rates,'o',label='test')
		i=i+1


	pylab.xlabel('Time (s) ')
	pylab.ylabel('Afterglow')

	pylab.savefig('Persistence.png')
	pylab.draw()



	

def get_pixel_value(filename,row=500,col=500,type='SCI'):
	'''
	get the value of a pixel for all of the extensions of one type,
	so for example a raw or ima file will likely return

	Really only makes sense as written for SCI data

	100525 - Coded from plot_pixel
	100831	ksl	Moved back into this routine because not generic and
			only used here.
	'''
	sci_ext=get_ext_info(filename,type)

	if len(sci_ext) == 0:
		print 'Error: No extension of type %s for %s' % (type,filename)
		return [],[]

	print 'test',sci_ext
	z=pyfits.open(filename)

	# The data are inverted in time, so the last extension has the
	# first data
	xx=[]
	yy=[]


	j=len(sci_ext)-1
	print 'j',j

	zero=z[sci_ext[j][0]].data
	print 'test', sci_ext[j][0], zero[row,col]

	j=j-1
	while j>=0:
		try:
			scidata=z[sci_ext[j][0]].data
			delta=scidata[row,col] - zero[row,col]
			print j, sci_ext[j][0], sci_ext[j][1],scidata[row,col], zero[row,col],delta 
			xx.append(sci_ext[j][1])
			yy.append(delta)
		except KeyError:
			break
		j=j-1
	
	xx=numpy.array(xx)
	yy=numpy.array(yy)
	return xx,yy

def plot_pixel(filename,row=500,col=500,new='yes',offset=0):
	'''
	Plot a single pixel from all the science extensions in an image
	Row and col have the same meaning as in ds9 or iraf

	Written with WFC3/IR in mind.  The zero read is subtracted.
	Called by pixel_history below

	100308 - Coded
	100525 - Modified to split getting the pixel values from the actual plotting
	'''

	xx,yy=get_pixel_value(filename,row,col)

	if len(xx)==0:
		print 'Error: plot_pixel: no pixel values to plot'
		return xx,yy
	
	xx=xx+offset

	pylab.figure(1,(6,6))
	if new=='yes':
		pylab.clf()
	pylab.plot(xx,yy,'.')
	pylab.xlabel('Time (s)')
	pylab.ylabel('Electrons')
	pylab.draw()
	return xx,yy

def pixel_history_fit(fileroot='observations_raw',row=500,col=500,bright='ibel01p1q'):
	'''
	Retrieve the pixel history for a single pixel in a set of images and 
	fit a pwoer law to the results
	'''
	records=read_ordered_list2(fileroot,dataset=bright,interval=[0,2],outroot='none')
	xx=[]
	yy=[]
	i=0
	for record in records:
		xname=parse_fitsname(record[0])
		if i==0:
			print xname[0]
			word=xname[0]
			word=word.replace('ima','flt')
			print word
			tt=get_times(word+'[1]')
			tzero=tt[1]
			dt=0
			i=i+1
		else: 
			x,y=get_pixel_value(xname[0],row,col,type='SCI')
			time=eval(record[6])
			dt=86400*(time-tzero)
			x=x+dt
			xx.append(list(x))
			yy.append(list(y))

	xx=numpy.array(xx)
	yy=numpy.array(yy)
	xx=numpy.ravel(xx)
	yy=numpy.ravel(yy)
	print xx
	print yy

	pylab.figure(2,(6,6))
	pylab.clf()
	p1=fit_power(xx,yy)
	pylab.axis([0,6000,0,2.5])
	print 'Fit results',p1
	pylab.xlabel('Time (s)')
	pylab.ylabel('Electrons')
	pylab.draw()
	return







def pixel_history(fileroot='observations_raw',row=500,col=500,cont='no'):
	'''
	Plot the history of a pixel for a series of raw data frames

	If cont is anything but no, the history will be plotted as 
	a continuous function of time.  no means to start the time
	over for each image.

	100308 - Coded to begin to look at whether one could use the short-term values
	contained in a raw data file to determine how the persistence decayed with time.

	100526 - Copied here from per_eval, as I begin to modify it to fit ima files

	'''
	records=read_ordered_list2(fileroot,dataset='first',interval=[0,2],outroot='none')
	i=0
	for record in records:
		xname=parse_fitsname(record[0])
		if i==0:
			tzero=eval(record[6])
			print tzero
			plot_pixel(xname[0],row,col,new='yes',offset=0)
			i=i+1
		else: 
			time=eval(record[6])
			dt=86400*(time-tzero)
			print 'times',time, dt
			if cont=='no':
				plot_pixel(xname[0],row,col,new='no',offset=0)
			else:
				plot_pixel(xname[0],row,col,new='no',offset=dt)
	return

# End utilities

# Section which is intended to characterize the persistance from a tungsten image and a 
# follow-on image.

def imshow(xx,yy,xmax=0,ymax=0):
	'''

	Called by examine
	
	'''
	from matplotlib.ticker import FuncFormatter


	if xmax==0:
		xmax=2*numpy.median(xx)
	if ymax==0:
		ymax=2*numpy.median(yy)

	isize=100	

	i=0
	image=numpy.zeros((isize,isize))
	while i<len(xx):
		x=int(xx[i]/xmax*isize)
		y=int(yy[i]/ymax*isize)
		x=min(x,isize-1)
		y=min(y,isize-1)

		x=max(x,0)
		y=max(y,0)

		# print x,y


		image[y][x]=image[y][x]+1
		i=i+1
	
	# inshow displays image in iraf order e. g. row column
	z=numpy.max(image)
	print 'The maximum value in the image is ',z

	# image[25][75]=z

	image=image+1
	image=numpy.log(image)
	pylab.figure(2,figsize=(6,6))
	pylab.clf()
	pylab.imshow(image,origin='lower')
	pylab.xlabel('Stimulus (e)')
	pylab.ylabel('Persistence (e/s)')

	j=0
	xticks=[]
	xxticks=[]

	if xmax>20000:
		idelta=max(int(xmax/200000) % 5 ,1)*200000
	elif xmax>50000:
		idelta=max(int(xmax/50000) % 5 ,1)*50000
	else:
		idelta=max(int(xmax/10000) % 5 ,1)*10000



	while j<xmax:
		x=(j/xmax*isize)
		xticks.append(x)
		xxticks.append(j)
		# j=j+50000
		j=j+idelta



	if ymax>0.5:
		jdelta=0.2
	elif ymax > 0.1:
		jdelta=0.04
	elif ymax > 0.05:
		jdelta-0.02
	else:
		jdelta=0.005

	j=0
	yticks=[]
	yyticks=[]
	while j<ymax:
		y=(j/ymax*isize)
		yticks.append(y)
		yyticks.append(j)
		# j=j+0.5
		j=j+jdelta


	pylab.xticks(xticks,xxticks)
	pylab.yticks(yticks,yyticks)


	pylab.draw()

def construct_persist(x,y,ncomb=1000,filename='none'):
	'''
	Construct a histogram curve of all of the values in x and y

	Also writes the persistence to a file
	'''

	print 'Starting construct_persist'

	# sort the array
	index=numpy.argsort(x)
	xx=[]
	yy=[]

	i=0
	while i+ncomb<len(x):
		j=0
		zx=0
		zy=0
		while j<ncomb:
			ii=index[i]
			zx=zx+x[ii]
			zy=zy+y[ii]
			i=i+1
			j=j+1
		zx=zx/ncomb
		zy=zy/ncomb
		xx.append(zx)
		yy.append(zy)

		print 'OK ',zx,zy

	pylab.figure(3,figsize=(6,6))
	pylab.clf()
	pylab.plot(xx,yy,'.')
	zz=pylab.axis()
	zz=list(zz)
	zz[0]=0
	zz[2]=0
	pylab.axis(zz)
	print 'axis',zz
	pylab.xlabel('Stimulus (e)')
	pylab.ylabel('Persistence (e/s)')

	pylab.draw()

	if filename=='none':
		return

	g=open(filename,'w')
	g.write('# Stimulus(e/s)   Persistence(counts/s)\n')
	i=0
	while i<len(xx):
		g.write('%8.0f %8.4f\n' % (xx[i],yy[i]))
		i=i+1
	g.close()
	return



	

def examine(bright='ibel01p1q',persist='1',fileroot='observations'):
	'''
	Examine the persistence in a single image, and plot the observed
	persistence as a function of the stimulus

	101011 - It's not clear this was left in a completed state.

	The imputs to the routine appear to be the image bright that 
	causes persistetnce and persist, which can either be a number
	in which case the intent was to read an offset, or a filename.
	'''

	zbright=bright
	zpersist=persist
	#Read the list of observations that are associated with this particular visit
	records=read_ordered_list2(fileroot,dataset=bright,interval=[0,3],outroot='none')

	bright=records[0][0]

	bright_image=get_image_ext(bright,1,'e')  #Get the data from the bright image in electrons


	#Read the image containing thing the persistence

	try:
		number=eval(persist)
		persist=records[number][0]
	except NameError:
		records=read_ordered_list2(fileroot,dataset=persist,interval=[0,3],outroot='none')
		persist=records[0][0]
	
	print 'Persistence file is ',persist

	persist_image=get_image_ext(persist,1,'e/s')

	dq=get_image_ext(persist,3)
	dq=get_image_ext(bright,3)

	# Seems like one would want to do somethig with these two dataquaility arrays, such as logical or them.

	mask=numpy.select([bright_image>30000],[0],default=1)

	mask=dq

	xxx=numpy.ma.array(bright_image,mask=mask)
	yyy=numpy.ma.array(persist_image,mask=mask)

	xxxx=xxx.compressed()
	yyyy=yyy.compressed()

	xxxx=numpy.ravel(xxxx)
	yyyy=numpy.ravel(yyyy)

	xmed=numpy.median(xxxx)
	ymed=numpy.median(yyyy)

	xmax=2*xmed
	ymax=3*ymed

	# So xxxx and yyyy are the stimulus and persistence images with the dq associated with the bright pixel removed.
	i=0
	fx=[]
	fy=[]
	while i<len(xxxx):
		fx.append(xxxx[i])
		fy.append(yyyy[i])
		i=i+1000



	print 'There are %d datapoints to plot  %8.1f  %8.1f' % (len(xxxx),xmed,ymed)

	pylab.figure(1,figsize=(6,6))
	pylab.clf()
	pylab.plot(xxxx,yyyy,'.')
	pylab.plot(fx,fy,'.')
	pylab.axis([0,xmax,0,ymax])
	pylab.xlabel('Stimulus (e)')
	pylab.ylabel('Persistence (e/s)')

	# Write a file for Adam

	name='Persistence_dist_%s_%s.txt' % (zbright,zpersist)

	pylab.draw()

	imshow(xxxx,yyyy,xmax,ymax)

	construct_persist(xxxx,yyyy,filename=name)






# End of this section

# def sum_power(names=['ibel03e2q.txt','ibel04l3q.txt','ibel01p1q.txt','ibel02g7q.txt','ibel05s2q.txt','ibel12gbq.txt','ibel11l5q.txt']):
def sum_power(names=['ibel01p1q.txt','ibel02g7q.txt','ibel03e2q.txt','ibel04l3q.txt','ibel05s2q.txt','ibel11l5q.txt','ibel12gbq.txt']):
	'''
	Create summary power law fits and a single plot of the data written to a text file by  doit.  Note that although
	the fit is calculated the model is not contained in the plot.

	This assumes that the data one wants to fit has been obtained by some other means.  It's
	really intended to generate a figure and re-do the fits that were obtained previously.

	Plot here is a linear plot
	'''

	pylab.figure(1,figsize=(6,6))
	pylab.clf()

	for name in names:
		f=open(name,'r')
		lines=f.readlines()
		f.close

		x=[]
		y=[]

		for line in lines:
			line=line.split()
			xx=eval(line[0])
			yy=eval(line[1])
			x.append(xx)
			y.append(yy)



		# pylab.plot(x,y,'bo')
		p1=fit_power(x,y) 
		print 'Fit results',p1

	pylab.xlabel('Time')
	pylab.ylabel('Afterglow (e/s)')

	pylab.savefig('Persistence_sum.png')
	pylab.draw()

	return

def sum_same():
	'''
	make a log log plot of the 3 examples that should be the same exposure
	'''
	sum_log_log(['ibel03e2q.txt','ibel12gbq.txt','ibel11l5q.txt'])
	sum_power(['ibel03e2q.txt','ibel12gbq.txt','ibel11l5q.txt'])
	return


def sum_Visit41():
	'''
	Compare visit 41 to 3 examples in Cycle 17  that should be the same exposure
	'''
	sum_log_log(['ibmf41lrq.txt','ibel01p1q.txt','ibel12gbq.txt','ibel11l5q.txt'])
	sum_power(['ibmf41lrq.txt','ibel01p1q.txt','ibel12gbq.txt','ibel11l5q.txt'])
	return

def sum_log_log(names=['ibel01p1q.txt','ibel02g7q.txt','ibel03e2q.txt','ibel04l3q.txt','ibel05s2q.txt','ibel11l5q.txt','ibel12gbq.txt'],labels=[]):
	'''
	Create summary log_log plot of the persisenced as  produced by doit (or some other means)


	Notes:
		This has been "customized" to produces plots for the ISR, and therefore may need to be rewritten for general use, if that is
		desirable

		This routine does not do any model fitting
	'''

	pylab.figure(2,figsize=(8,8))
	pylab.clf()

	plotno=0
	for name in names:
		f=open(name,'r')
		lines=f.readlines()
		f.close()

		x=[]
		y=[]

		for line in lines:
			line=line.split()
			xx=eval(line[0])
			yy=eval(line[1])
			x.append(xx)
			y.append(yy)



		if len(labels)==0:
			if plotno==6:
				plotno=11
				pylab.loglog(x,y,'o',label='%02d' % (plotno+1))
		else:
			pylab.loglog(x,y,'o',label=labels[plotno])
		plotno=plotno+1

	pylab.xlabel('Time (s)')
	pylab.ylabel('Persistence (e/s)')
	pylab.legend(numpoints=1)

	pylab.savefig('Persistence_sum_loglog.png')
	pylab.draw()

	return

def fit_alpha_only(names=['ibel01p1q.txt','ibel02g7q.txt','ibel03e2q.txt','ibel04l3q.txt','ibel05s2q.txt','ibel11l5q.txt','ibel12gbq.txt'], zmin=0.05):
	'''
	This is a routine intended to find the best alpha for a number of persistence files.  The idea is that the normalization
	might change but we want a single alpha.  The files need to to consist of times and persistence rates., such as those which 
	are produced by doit.  

	zmin is a value below which we do not think the persistence measurement is valid 

	The idea is that we will do individual fits first, and then renormalize each of the individual files based on the fits.  At that point
	one can do a global fit.  It's a kluge but can be done quickly

	101013	Coding begun
	'''

	pp=[]
	all_times=[]
	all_values=[]

	# First get all of the data
	for name in names:
		f=open(name,'r')
		lines=f.readlines()
		f.close()

		x=[]
		y=[]

		for line in lines:
			line=line.split()
			xx=eval(line[0])
			yy=eval(line[1])
			if yy>zmin:
				x.append(xx)
				y.append(yy)
		if len(y)>4:
			all_times.append(x)
			all_values.append(y)
		else:
			print 'Did not include %s' % name
	# So at this point we can obtain the fits and rescale 
	# to 1 at 1000 seconds

	# Note tthat all times and all values do not have the same number of records in each
	# row so it can not be converted to a numpy.array



	print all_times


	i=0
	while i<len(all_times):
		time=numpy.array(all_times[i])
		values=numpy.array(all_values[i])
		p=fit_power(time,values)
		pp.append(p)
		values=values/p[0]
		if i==0:
			new_times=numpy.copy(time)
			new_values=numpy.copy(values)
		else:
			new_values=numpy.append(new_values,values)
			new_times=numpy.append(new_times,time)
		i=i+1

	print new_times
	print new_values
	print 'OK now determine the global fits with all scaled together'
	p=fit_power(new_times,new_values,new='yes')

	print p



	

	return 





def fit_power(x,y,new='no'):
	'''
	Fit a power law to the data, and plot the result.  If new is anything but 'no'
	then the figure will be cleared before the plot is made
	'''

	x=numpy.array(x)
	y=numpy.array(y)


	fitfunc = lambda p, x: p[0]*pow(x/1000.,-p[1])  # powerlaw
	errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
	p0 = [1,1000., 1.] # Initial guess for the parameters
	p0 = [1, 1.] # Initial guess for the parameters
	p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))


	time = numpy.linspace(x.min(), x.max(), 1000)

	if new!='no':
		pylab.clf()
	pylab.plot(x, y, "o", time, fitfunc(p1, time), "-") # Plot of the data and the fit
	pylab.xlabel('Time (s) ')
	pylab.ylabel('Afterglow (e/s)')

	# Calculate the error

	zzz=errfunc(p1,x,y)
	zzz=abs(zzz)
	print 'Fit error: typical  %0.3f max %0.3f' % (numpy.median(zzz),numpy.max(zzz))

	zzz=zzz*zzz
	q=numpy.sqrt(zzz)

	print 'RMS error',numpy.average(q)

	return p1

def doit_ima(bright='ibel01p1q',fileroot='observations'):
	'''

	This version is expected to work on the ima files

	bright is the dataset name of the tungsten exposure
	fileroot is the rootname for the list of observations files

	This routine reads the bright image and the follwing images.  It then runs
	imstats on the images and fits a power law to them.

	It's not obvious that this does anything special for ima files


	'''

	xbright=bright  # Kluge

	#Read the list of observations that are associated with this particular visit
	records=read_ordered_list2(fileroot,dataset=bright,interval=[0,3],outroot='none')

	bright=records[0][0]

	bright_image=get_image_ext(bright,1,'e')  #Get the data from the bright image in electrons
	bright_time=get_times(bright)
	bright_time=bright_time[1]

	x=[]
	y=[]

	i=1
	while i<len(records):
		record=records[i]
		filename=record[0]

		times=get_times(filename)
		dtime=0.5*(times[0]+times[1])-bright_time
		dtime=dtime*86400
		print dtime



		x.append(dtime)
		zzz=pyraf.iraf.imstatistics(filename,fields = "npix,mean,stddev,min,max,image", lower='INDEF', upper='INDEF',nclip=2, lsigma=3.0, usigma=3.0, binwidth=0.1,format='no',Stdout=1)
		zzzz=zzz[0].split()
		# zzzz=zzz[1].split('\t')
		print zzzz
		y.append(eval(zzzz[1]))
		i=i+1
	
	print x
	print y
	g=open(xbright+'.txt','w')
	i=0
	while i<len(x):
		g.write('%6.1f %6.3f\n' % (x[i],y[i]))
		i=i+1
	g.close

	pylab.figure(1,figsize=(6,6))
	pylab.clf()

	pylab.plot(x,y)

	pylab.xlabel('Time (s) ')
	pylab.ylabel('Afterglow (e/s)')

	p1=fit_power(x,y)
	print 'Fit results',p1

	pylab.savefig('Persistence_'+xbright+'.png')
	pylab.draw()

	return


def do_sect(bright='ibel01p1q',ncol=500,nrow=500,nsize=100,fileroot='observations'):
	'''
	Using the bright image as a reference, construct a the history of persistence
	in a subsection of the following images, and fit this to a power law.  This
	version does not use imstats and hence bad pixels can be eliminated.

	This works on flt files.

	101013 - Work begun
	'''


	#Read the list of observations that are associated with this particular visit
	records=read_ordered_list2(fileroot,dataset=bright,interval=[0,3],outroot='none')
	bright=records[0][0]

	# Get the end time for the persistence image
	bright_time=get_times(bright)
	bright_time=bright_time[1]

	# These are arrays to stores times and fluxes
	x=[]
	y=[]

	# Cycle through the files obtaining the necessary statistics
	i=1
	while i<len(records):
		record=records[i]
		filename=record[0]

		# Find the midpt of the time for this image
		times=get_times(filename)
		dtime=0.5*(times[0]+times[1])-bright_time
		dtime=dtime*86400
		x.append(dtime)

		sc=get_image_ext(filename,1)
		dq =get_image_ext(filename,3)

		# Establish the boundaries of the region we want statistics for

		nrows=len(sc)
		ncols=len(sc[0])

		nrow_min=int(nrow-0.5*nsize)
		nrow_max=nrow_min+nsize

		ncol_min=int(ncol-0.5*nsize)
		ncol_max=ncol_min+nsize

		if nrow_min<0:
			nrow_min=0
		if nrow_max>nrows:
			nrow_max=nrows

		if ncol_min<0:
			ncol_min=0
		if ncol_max>ncols:
			ncol_max=ncols

		scx=sc[nrow_min:nrow_max,ncol_min:ncol_max]
		dqx=dq[nrow_min:nrow_max,ncol_min:ncol_max]
		scxx=numpy.ma.array(scx,mask=dqx)
		scxxx=scxx.compressed()
		scxxx=numpy.ravel(scxx)

		value=numpy.median(scxxx)

		y.append(value)
		i=i+1
	
	print x
	print y

	pylab.figure(1,figsize=(6,6))
	pylab.clf()

	pylab.plot(x,y)

	pylab.xlabel('Time (s) ')
	pylab.ylabel('Afterglow (e/s)')

	p1=fit_power(x,y)
	print 'Fit results',p1

	pylab.draw()

	return p1

def do_sections(bright='ibel01p1q',nsize=100,fileroot='observations'):
	'''
	Run do_sect in a grid over the whole array and determine
	how uniform the power law is for each section
	'''

	irow=nsize/2
	norm=[]
	alpha=[]
	while irow < 1000:
		icol=nsize/2
		while icol<1000:
			x=do_sect(bright,icol,irow,nsize,fileroot) 
			norm.append(x[0])
			alpha.append(x[1])
			icol=icol+nsize
		irow=irow+nsize
	
	print 'The Results of Sections'
	print 'Norm',norm
	print 'Alpha',alpha

	norm=numpy.array(norm)
	alpha=numpy.array(alpha)
	print 'Norm - Summary', numpy.mean(norm),numpy.std(norm)
	print 'alpha- Summary', numpy.mean(alpha),numpy.std(alpha)
	

	

def doit(bright='ibel01p1q',fileroot='observations'):
	'''

	This routine reads the bright image and the follwing images.  It then runs
	imstats on the images and fits a power law to them.

	bright is the dataset name of the tungsten exposure. 
	fileroot is the rootname for the list of observations files

	The routine writes the results of imstats to a txt file, and creates
	a plot that compares the fits and the data.  The results of the fits
	only appear on the screen.

	Notes:

	The tungsten image does not figure in the fits except to set time zero, so
	it is not necessary that the tungsten image be correctly fluxed for this
	routine.


	'''

	xbright=bright  # Kluge to avoid the fact that I use bright for something elese later

	#Read the list of observations that are associated with this particular visit
	records=read_ordered_list2(fileroot,dataset=bright,interval=[0,2],outroot='none')
	bright=records[0][0]

	# bright_image=get_image_ext(bright,1,'e')  #Get the data from the bright image in electrons
	bright_time=get_times(bright)
	bright_time=bright_time[1]

	x=[]
	y=[]

	i=1
	while i<len(records):
		record=records[i]
		filename=record[0]

		times=get_times(filename)
		dtime=0.5*(times[0]+times[1])-bright_time
		dtime=dtime*86400
		x.append(dtime)

		zzz=pyraf.iraf.imstatistics(filename,fields = "npix,mean,stddev,min,max,image", lower='INDEF', upper='INDEF',nclip=3, lsigma=3.0, usigma=3.0, binwidth=0.1,format='no',Stdout=1)
		zzzz=zzz[0].split()
		print 'Statistics',zzzz,'at time',dtime
		y.append(eval(zzzz[1]))
		i=i+1
	
	# print x
	# print y
	g=open(xbright+'.txt','w')
	i=0
	while i<len(x):
		g.write('%6.1f %6.3f\n' % (x[i],y[i]))
		i=i+1
	g.close

	pylab.figure(1,figsize=(6,6))
	pylab.clf()

	pylab.plot(x,y)

	pylab.xlabel('Time (s) ')
	pylab.ylabel('Afterglow (e/s)')

	p1=fit_power(x,y)
	print 'Fit results',p1

	pylab.savefig('Persistence_'+xbright+'.png')
	pylab.draw()

	return

def do_one(xcen=600,ycen=800,fileroot='observations'):
	'''
	Analyze a single pixel in detail
	'''

	records=read_ordered_list0(fileroot)

	n=0
	i=0
	istop=[]
	t=[]
	y=[]
	exposure=[]

	istart=[]
	i=0
	while i<len(records):
		record=records[i]
		type=record[14]
		times=get_times(record[0])
		if type=='TUNGSTEN':
			exposure.append(eval(record[11]))
			istart.append(i)
			bright_time=times[1]
			dtime=0
		else:
			dtime=0.5*(times[0]+times[1])-bright_time
			dtime=dtime*86400

		image=get_image_ext(record[0],1)
		t.append(dtime)
		y.append(image[xcen,ycen])
		i=i+1
	
	i=0
	while i < len(istart)-1:
		istop.append(istart[i+1])
		i=i+1

	istop.append(len(records))
	

	print istart
	print istop
	print len(t),len(y)

	pylab.figure(1,figsize=(6,6))
	pylab.clf()
	
	i=0
	while i<len(istart):
		tt=t[istart[i]+1:istop[i]]
		yy=y[istart[i]+1:istop[i]]
		# pylab.plot(tt,yy)
		p1=fit_power(tt,yy)
		print 'Fit results  %6.1f  %6.3f  %6.2f' % (exposure[i],p1[0],p1[1])
		i=i+1

	pylab.xlabel('Time (s) ')
	pylab.ylabel('Afterglow (e/s)')

	

	pylab.draw()
	return

# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		print 'usage: per_cal.py  filename'
