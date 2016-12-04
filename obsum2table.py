#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Convert the current observations.sum file to an astropy
table to work with run_persist


Command line usage (if any):

    usage: obsum2table.py filename

Description:  

Primary routines:

    doit

Notes:
                                       
History:

161203 ksl Coding begun

'''

import sys
from astropy.io import ascii
from astropy.table import Table
import numpy


def read_file(filename,char=''):
    '''
    Read a file and split it into words, eliminating comments
    
    char is an optional parameter used as the delimiter for
    splitting lines into words.  Otherwise white space is
    assumed.

    History:
    
    110729    ksl    Added optional delimiters
    141209    ksl    Reinstalled in my standard startup
            script so there was flexibility to
            read any ascii file
    '''

    try:
        f=open(filename,'r')
        xlines=f.readlines()
        f.close()
    except IOError :
        print ("The file %s does not exist" % filename)
        return []   
    
    lines=[]
    
    i=0
    while i<len(xlines):
        z=xlines[i].strip()
        if char=='':
            z=z.split()
        else:
            z=z.split(char)
        if len(z)>0:
            if z[0][0]!='#':
                lines=lines+[z]
        i=i+1
    return lines




def read_table(filename='foo.txt',format=''):
    '''
    Read a file using astropy.io.ascii and 
    return this 

    Description:

    Notes:

    History:


    '''
    try:
        if format=='':
            data=ascii.read(filename)
        else:
            data=ascii.read(filename,format=format)
    except IOError:
        print ('Error: file %s does not appear to exist' % filename)
        return

    print ('Here are the column names:')
    
    print (data.colnames)

    return data


def doit(filename='observations.sum',outputfile='observations.sum.txt'):
    '''
    Do something magnificent

    Description:

    Notes:

    History:


    '''
    
    


    data=read_file(filename)

    dataset=[]
    prog_id=[]
    date=[]
    xobs=[]
    xtime=[]
    processed=[]
    ext1=[]
    ext2=[]
    ext3=[]
    int1=[]
    int2=[]
    int3=[]
    qhtml=[]
    for one in data:
        dataset.append(one[0])
        prog_id.append(one[1])
        date.append(one[2])
        xobs.append(one[3])
        xtime.append(one[4])
        processed.append(one[5])
        if len(one)==13:
            ext1.append(one[6])
            ext2.append(one[7])
            ext3.append(one[8])
            int1.append(one[9])
            int2.append(one[10])
            int3.append(one[11])
            qhtml.append(one[12])
        else:
            ext1.append(-99)
            ext2.append(-99)
            ext3.append(-99)
            int1.append(-99)
            int2.append(-99)
            int3.append(-99)
            qhtml.append('--')



    x=Table([dataset],names=['Dataset'])
    x['ProgID']=prog_id
    x['ExpStart']=date
    x['Proc-Date']=xobs
    x['Proc-Time']=xtime
    x['ProcStat']=processed
    x['E0.10']=ext1
    x['E0.03']=ext2
    x['E0.01']=ext3
    x['I0.10']=int1
    x['I0.03']=int2
    x['I0.01']=int3
    x['PerHTML']=qhtml


    
    print (x)


    # This format is the easy to read back automatically
    ascii.write(x,outputfile,format='fixed_width_two_line')

    return


    





# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)>1:
        # doit(int(sys.argv[1]))
        doit(sys.argv[1])
    else:
        print ('usage: obsum2table.py filename')
