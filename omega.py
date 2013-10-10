#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

This is test routine intended to explore a better fitting function
for persistence using the OmegaCen data


Command line usage (if any):

	usage: omega.py filename

Description:  

Primary routines:

Notes:
									   
History:

110301 ksl Coding begun
120410	ksl	Began to adapt for Cycle 19

'''

import sys
import numpy
import pylab
from scipy import optimize
from history import *



def calc_fermi_pow(x,t,norm=1.0,e_fermi=80000,kt=20000,alpha=0.3,gamma=1):
	'''
	calculate the shape of the afterglow assuming a Fermi-Dirac like
	distribution and series of x values.  This routine just calculates 
	the shape.

	The Fermi distribution is defined as 

	1/(exp((E_f-E)/kT)+1) where E_f is the Fermi energy.  For E>E_f the exponential term has little  
	effect.  kT determines how sharpe the cutoff is.  For E=E_f we are at the distribution value of
	the function is 0.5 regardless of kT

	100305 - Note that this operates on any dimension array
	100604 - This should be set up so the normalization is at 1000 s
	
	'''

	y=e_fermi-x
	y=y/kt
	y=numpy.exp(y)
	y=1./(y+1.)
	y=y*(x/e_fermi)**alpha
	y=norm*y


	# if t>0:
	# 	tfactor=(t/1000.)**(-gamma)
	# else:
	# 	tfactor=1
	tfactor=(t/1000.)**(-gamma)
	y=y*tfactor
	return y



def calc_fermi2(x,norm=1.0,e_fermi=80000,kt=20000,alpha=0.3):
	'''
	calculate the shape of the afterglow assuming a Fermi-Dirac like
	distribution and series of x values.  This routine just calculates 
	the shape.

	The Fermi distribution is defined as 

	1/(exp((E_f-E)/kT)+1) where E_f is the Fermi energy.  For E>E_f the exponential term has little  
	effect.  kT determines how sharpe the cutoff is.  For E=E_f we are at the distribution value of
	the function is 0.5 regardless of kT

	calc_fermi2 allows for a slow rise after saturation, unlike calc_fermi

	100305 - Note that this operates on any dimension array
	100604 - This should be set up so the normalization is at 1000 s
	
	'''

	y=e_fermi-x
	y=y/kt
	y=numpy.exp(y)
	y=1./(y+1.)
	y=y*(x/e_fermi)**alpha
	y=norm*y
	return y


def calc_fermi(x,norm=1.0,e_fermi=80000,kt=20000):
	'''
	calculate the shape of the afterglow assuming a Fermi-Dirac like
	distribution and series of x values.  This routine just calculates 
	the shape.

	The Fermi distribution is defined as 

	1/(exp((E_f-E)/kT)+1) where E_f is the Fermi energy.  For E>E_f the exponential term has little  
	effect.  kT determines how sharpe the cutoff is.  For E=E_f we are at the distribution value of
	the function is 0.5 regardless of kT

	100305 - Note that this operates on any dimension array
	100604 - This should be set up so the normalization is at 1000 s
	
	'''

	y=e_fermi-x
	y=y/kt
	y=numpy.exp(y)
	y=1./(y+1.)
	y=norm*y
	return y



def fit_fermi_pow(x,t,y):
	'''
	Fit the combination of a ferm oermi function and a power loaw 
	time decat to the data


	101303	ksl	Coded as part of effort to cacculate a model
			with a fermi-like function for the persistence
			and a pwoer law decay



	'''

	# To get leasq to produce anything it was necessar to flatten
	x=numpy.ravel(x)
	y=numpy.ravel(y)
	t=numpy.ravel(t)


	fitfunc = lambda p, x, t: calc_fermi_pow(x,t,p[0],p[1],p[2],p[3],p[4])


	errfunc = lambda p, x, t, y: fitfunc(p, x, t) - y # Distance to the target function

	p0 = [1, 80000,20000,0.5,-1] # Initial guess for the parameters


	# Now fit the results

	p1, success = optimize.leastsq(errfunc, p0[:], args=(x, t,y))

	# Calculate the error
	model=fitfunc(p1,x,t)
	print 'Model    : typical  %0.3f max %0.3f' % (numpy.median(model),numpy.max(model))

	zzz=errfunc(p1,x,t,y)
	q=zzz
	zzz=abs(zzz)
	print 'Fit error: typical  %0.3f max %0.3f' % (numpy.median(zzz),numpy.max(zzz))

	zzz=zzz*zzz
	q=numpy.sqrt(zzz)

	print 'RMS error',numpy.average(q)


	return p1



def fit_fermi(x,y):
	'''
	Fit a fermi function to the data

        calc_fermi(x,norm=1.0,e_fermi=80000,kt=20000):
	
	y=p0*x**-p1 to the data.
	
	Return the fit, and also a "continuous" function containing
	the power law fit.  
	
	Note:  This version does not plot the results (as the corresponding
	routine in per_cal.py did.  

	101227	ksl	Adapted from same routine in per_cal

	'''

	x=numpy.array(x)
	y=numpy.array(y)


	# fitfunc = lambda p, x: p[0]*pow(x/1000.,-p[1])  # powerlaw

	# fitfunc = lambda p, x: calc_fermi(x,p[0],p[1],p[2])
	fitfunc = lambda p, x: calc_fermi2(x,p[0],p[1],p[2],p[3])
	errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

	p0 = [1, 80000,20000] # Initial guess for the parameters
	p0 = [1, 80000,20000,0.5] # Initial guess for the parameters
	p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))


	# Calculate the error

	zzz=errfunc(p1,x,y)
	q=zzz
	zzz=abs(zzz)
	print 'Fit error: typical  %0.3f max %0.3f' % (numpy.median(zzz),numpy.max(zzz))

	zzz=zzz*zzz
	q=numpy.sqrt(zzz)

	print 'RMS error',numpy.average(q)

	time = numpy.linspace(x.min(), x.max(), 1000)
	fit=fitfunc(p1, time)

	# pylab.figure(2,(6,6))
	# pylab.clf()
	# pylab.plot(x,q)
	# pylab.draw()


	return p1,[time,fit]


def fit_power(x,y):
	'''
	Fit a power law of the form
	
	y=p0*x**-p1 to the data.
	
	Return the fit, and also a "continuous" function containing
	the power law fit.  
	
	Note:  This version does not plot the results (as the corresponding
	routine in per_cal.py did.  

	101227	ksl	Adapted from same routine in per_cal

	'''

	x=numpy.array(x)
	y=numpy.array(y)


	fitfunc = lambda p, x: p[0]*pow(x/1000.,-p[1])  # powerlaw
	errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
	p0 = [1,1000., 1.] # Initial guess for the parameters
	p0 = [1, 1.] # Initial guess for the parameters
	p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))


	# Calculate the error

	zzz=errfunc(p1,x,y)
	zzz=abs(zzz)
	print 'Fit error: typical  %0.3f max %0.3f' % (numpy.median(zzz),numpy.max(zzz))

	zzz=zzz*zzz
	q=numpy.sqrt(zzz)

	print 'RMS error',numpy.average(q)

	time = numpy.linspace(x.min(), x.max(), 1000)
	fit=fitfunc(p1, time)

	return p1,[time,fit]


def get_data(filename):
	'''
	Read the two column file containing the persistence, elimiating comments
	'''

	try:
		f=open(filename,'r')
		xlines=f.readlines()
		f.close()
	except IOError :
		print "The file %s does not exist" % filename
		return [],[]   
	
	x=[]
	y=[]

	i=0
	while i<len(xlines):
		z=xlines[i].split()
		if len(z)>0:
			if z[0][0]!='#':
				x.append(eval(z[0]))
				y.append(eval(z[1]))

						
		i=i+1
	x=numpy.array(x)
	y=numpy.array(y)
	return x,y


def do_one_fermi_fit(filename='foo'):
	'''
	Fit a single example of persistence in
	an image
	'''
	x,y=get_data(filename)    

	if len(x)==0:
		print 'There is nothing to do'
		return
	i=0
	while i<len(x):
		print x[i],y[i]
		i=i+1

	# pl,fit=fit_power(x,y)
	pl,fit=fit_fermi(x,y)

	pylab.figure(1,figsize=(6,6))
	pylab.clf()

	pylab.plot(x,y,'o')
	pylab.plot(fit[0],fit[1],'-')
	pylab.xlabel('Stimulus')
	pylab.ylabel('Persitence')
	pylab.draw()
	print pl
	return pl

def do_many_fermi(files='files.ls'):
	'''
	Read a list of files, fit each one individually to a fermi function and
	print out the results
	'''

	# Read the names of the files that contain the data
	# The filename is the first word in the file.
	f=open(files,'r')
	lines=f.readlines()
	
	filenames=[]
	for line in lines:
		x=line.strip()
		x=x.split()
		if x[0]!='#':
			filenames.append(x[0])
	
	results=[]
	for one in filenames:
		p=do_one_fermi(one)
		results.append(p)
	
	i=0
	while i < len(filenames):
		print filenames[i],results[i]
		i=i+1


def do_global_fit(fileroot='stat_ghost'):
	'''
	Read a list of files and try to construct a global fit to the data

	This attempts to find a model for both the shape of persitence in individual images
	and the time decay. 

	It was intended for modeling the OmegaCen calibration program in Cycle 18, but it
	should be possible to use it in any situation where there is a stimulus image
	which produces a wide range of exposures in the stimulus image followed by a series
	of darks
	'''

	# Read the names of the files that contain the data
	# The filename is the first word in the file.
	# The delay time associated with the file is the second word
	# in the file

	# This secont also creates a set of data in which earh row is one of the data
	# sets

	f=open(fileroot+'.ls','r')
	lines=f.readlines()
	
	filenames=[]
	i=0
	for line in lines:
		x=line.strip()
		x=x.split()
		if x[0]!='#':
			filenames.append(x[0])
			dt=eval(x[1])
			x,y=get_data(x[0])    

			#120410 - Add individual fits for x and y
			params,dummy=fit_fermi(x,y)
			# print 'xxx',params
			string='# Individual fits: %7.1f' % dt
			for one in params:
				string=string+' %10.3f ' % one
			print string
			log('%s\n' % string)



			t=numpy.ones_like(x)
			t=t*dt
			if i==0:
				xx=numpy.copy(x)
				yy=numpy.copy(y)
				tt=numpy.copy(t)
			else:
				xx=numpy.vstack((xx,x))
				yy=numpy.vstack((yy,y))
				tt=numpy.vstack((tt,t))
				
			i=i+1


	g=open(fileroot+'.txt','w')
	# So at this point we have read in all of the data and are ready to attempt a fit
	params=fit_fermi_pow(xx,tt,yy)

	string=''
	for one in params:
		string=string+' %6.3f ' % one

	print 'The fit params are ', string

	g.write('# Results %s\n' % string)
	log('# Global fits: %20s %s\n' % (fileroot,string))


	

	fitfunc = lambda p, x, t: calc_fermi_pow(x,t,p[0],p[1],p[2],p[3],p[4])

	xxx=numpy.ravel(xx)
	yyy=numpy.ravel(yy)
	ttt=numpy.ravel(tt)

	z=fitfunc(params,xxx,ttt)

	pylab.figure(1,(8,8))
	pylab.clf()
	pylab.plot(xxx,yyy,'.')
	pylab.plot(xxx,z,'.')
	pylab.xlabel('Stimulus (elec)')
	pylab.ylabel('Persistence (elec/s)')
	pylab.draw()
	pylab.savefig(fileroot+'.png')

	print 'Result for %s' % fileroot
	string='# Time  MeanP  MeanDiff  Sigma     max    min AbsMeanDiff   AbsMaxDiff'
	print string
	g.write('%s\n' % string)

	pylab.figure(2,(6,6))
	pylab.clf()
	# Now go through each model and compare it to the data
	i=0
	while i<len(xx):
		# x,y,t are one datasett
		x=xx[i]
		y=yy[i]
		t=tt[i]
		# z is the results of the best fit model for that dataset
		z=fitfunc(params,x,t)

		# delta is the difference between the data and the model
		delta=z-y
		absdelta=abs(delta)
		ratio=delta/z
		amean=numpy.mean(absdelta)
		amax=numpy.max(absdelta)

		ymean=numpy.mean(y)
		mean=numpy.mean(delta)
		std=numpy.std(delta)
		max=numpy.max(delta)
		min=numpy.min(delta)
		string='%7.1f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f' %  (t[0],ymean,mean,std,max,min,amean,amax)
		print string
		g.write('%s\n' % string)
		# pylab.plot(x,delta,'o')
		pylab.semilogx(x,ratio,'-')
		pylab.axis([1e4,1e6,-1,1])
		pylab.draw()
		i=i+1
	
	


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		do_global_fit(sys.argv[1])
	else:
		print 'usage: omega.py filename'
	

