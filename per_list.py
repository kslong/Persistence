#! /usr/bin/env python 

'''
                    Space Telescope Science Institute

Usage (as a standalone program: 
    per_list.py 

or 

per_list.py [various options]  fileroot  

where     fileroot     the root name of the outputfiles
    aperture is 'all' or 'full' with the latter indicating one creates a list containing only the full IR frames, no subarrays

    Options:
    -h          Print this help
    -np 6       implies to run in parallel mode with 6 processors
    -all        Create a new .ls file and a new summary file
    -new_sum    Create a new summary file
    -daily      Create a new .ls but just update the old summary file
    -file_type  Instead of the defult flt files, create a file continaing
                a different file type, e.g. -file_type raw  to get the
                raw data files

without any arguments the call is effectively

    per_list.py  -all -file_type flt observations 


    per_list.py
    
Synopsis:  


If this is run as a standalone program, the routine creates two files.  The first
is set of ascii records that describe all the WFC3/IR files of a certain type 
that are in this directory and any subdirectories.  The second is an observations.sum
file that really is only pertinent to the persistence processing software

More generally this module contains routines both to make an ordered list of the records
and then to read them back, as well as a number of utilities to handle matters 
related to where the Persist directories are located.

Description:  

    The routine uses 'find' to locate all of the files of a certain type,
    aand then orders the list in time.  If there are duplicate files then
    it uses the creation date as the one to put in the list.  The dupliccates
    are in the list also, but are commented out. 

Primary routines:

    steer             Parses the command line and runs the program
    mske_ordered_list    Creates .ls file
    make_sum_file        Creates/updates the summary file

Notes:

History (Recent):

091026  ksl Coding begun
160104  ksl Changes to lock and unlock files so that parallel processing
            will not corrupt the obervations.sum file
161130  ksl Remove pyraf dependencies because it is going away.
161204  ksl Updated for python3, and to use astropy tables

'''

import sys
import os
import numpy
import time
import math
import scipy
import subprocess
import pylab
import shutil
import per_fits
from astropy.io import fits
from astropy.io import ascii
from multiprocessing import Pool
from astropy.table import Table
from astropy.table import join
from astropy.table import vstack


# Utilities

def read_table(filename='foo.txt',format=''):
    '''
    Read a file using astropy.io.ascii and 
    return this 

    Description:

        This simply reads a table, but it also
        assures that strings have been converted
        to objects.  This can be important if one
        intends to add or manipulate strings, 
        because astropy tends to truncate such
        objects

    Notes:

    History:

    161214  ksl Added in order to handle potential isssues
        with string truncation

    '''
    try:
        if format=='':
            data=ascii.read(filename)
        else:
            data=ascii.read(filename,format=format)
        for col in data.itercols():
            if col.dtype.kind in 'SU':
                data.replace_column(col.name,col.astype('object'))
    except IOError:
        print ('Error: file %s does not appear to exist' % filename)
        raise IOError
        return
    except Exception as e:
        print('Error: Strange exception of file %s' % filename)
        raise IOError
        sys.exit(0)



    return data

def backup(filename,tformat='%y%m%d',force='no'):
    '''
    Backup a file with a standard way of naming the file.  
    
    The original file is copied to the backup unless the
    backup file already exists or force is 'yes'.
    
    The format can be any valid time.strftime format, but
    probably should have no spaces.

    The routine returns the name of the backup file, or 
    'None' if there was no backup made

    Notes:  The idea of this is that one wants to keep
    a history of a file over time, but not at infinite
    time resolution.

    110105    ksl    Coded as part of effort to be able
            to use per_list as aprt of a crontab
    '''

    curtime=time.strftime(tformat,time.localtime())
    if os.path.exists(filename)==False:
        print('Warning: No file %s to backup' % filename)
        return 'None'

    backup_file='%s.%s.old' % (filename,curtime)

    if os.path.exists(backup_file)==False or force=='yes':
        shutil.copy(filename,backup_file)
        return backup_file
    else:
        return 'None'

def open_file(filename,permiss=0o770):
    '''
    Open a file for writing and set its permissions

    This was written in an attempt to get all of the
    files to have a common set of permissions
    no mattter who writes them

    '''

    if os.path.isfile(filename):
        os.remove(filename)
    g=open(filename,'w')
    try:
        os.chmod(filename,permiss)
    except OSError:
        print('Error: open_file: OSerror trying to set permissions')
    return g



def set_path(name,mkdirs='yes',local='no'):
    '''
    Check that the directory Persist where the outputs are to go exists, and if not
    create it, and also create an underlying directory .../Persist/Figs.  

    The input is the complete filename, including the path.  The routine 
    requires that there be a '/' in the complete name.  If not, it will
    create Persist below the current working directory 

    The routine returns the path include the directory name Persist

    Note:

    Normally this routine is given either the complete name to a file in the
    visit directory, e. g. 
        ./D/VisitAJ/ib0lajdyq_flt.fits[1] or 
    a name to a file in the  Persist directory, e. g.
        ./D/VisitAJ/Persist/filename
    This accounts for the somewhat convoluted attempt to locate where the Persist
    directory should be

    mkdirs must be 'yes' for directories to be made.  If it is for example, 'no', 
    then the path will be returned but there will be not attempt to make the
    directory

    If local is 'no' then the path for the Persist directory is determined by the name
    of the file that accompanies the definition.  
    If local is 'yes' then the directory will be names Persist, and it will sit immediatel
        below the run directory
    If local is anything else, then the name will still be below the run directory, but
        but have the name of the local variable.

    

    101014  ksl Added
    101214  ksl Moved into per_list since this is mostly used in conjunction
                with other routines that are there. It is not obvious that
                this is the correct place for this.
    101215  ksl Added creation of the Figs directory
    110105  ksl Made creation of the directories an option
    110121  ksl Made change to control permissions of the directories
    110203  ksl Made change to allow the directories to be written beneath
            the current working directory, a change that is primarily
            for testing
    110811  ksl The directory permissions are set so that the group name 
            should be inherited.
    110811  ksl Added commands to set the group name, but only if we are
                in the Quicklook2 directory structure.  These commands
                would need to change if the group names change or directory
                names change.  They should have no effect outside the standard
                structure
    170314  ksl Modified so that if local is not yes or no, then one will try 
                to create a local directoy of a given name
    '''


    if len(name)==0:
        print('Error: set_path: name had length 0, nothing to parse')
        return ''


    # Determine where we want Persist to be located.  
    if local=='no':
        try:
            i=name.rindex('Persist')
            path=name[0:i]
        except ValueError:
            # Find the last / in the name

            try:
                i=name.rindex('/')
                path=name[0:i]
            except ValueError:
                print('Warning: set_path: Assuming the intent was to work in the current directory')
                path='./'
    else:
        path='./'

    # Check whether the parent directory for Persist exists
    if os.path.exists(path)==False:
        string='set path: The directory %s contained in %s does not exist' % (path,name)
        print('Error:  %s' % (string))
        return 'NOK %s ' % string


    if local=='no' or local=='yes':
        path=path+'/Persist/'
    else:
        path='%s/%s/' % (path,local)


    if mkdirs!='yes':  
        return path

    # If one has reached this point then you want to make any necessary directories
    # set the chmod of the diretory so those within the group can delete the directory
    # when necessary
    
    # Note that the group name assignments below only are valid with the current
    # Quicklook2 directory structure and group names



    if os.path.exists(path)==False:
        try:
            os.mkdir(path)
            os.chmod(path,0o2770)
            if path.count('QL_GO'):
                os.chown(path,-1,6047)
            elif path.count('Quicklook'):
                os.chown(path,-1,340)


        except OSError:
            print('Error: set_path: Could not create %s for %s' % (path,name))
            return 'NOK Could not create %s' % path


    # Add a figs directory if it does not exist as well
    figs=path+'/Figs'
    if os.path.exists(figs)==False:
        os.mkdir(figs)
        os.chmod(figs,0o2770)
        if figs.count('QL_GO'):
            os.chown(figs,-1,6047)
        elif figs.count('Quicklook'):
            os.chown(figs,-1,340)

    
    return path

def parse_dataset_name(name):
    '''
    Check if we have been given the name of the fits file instead of the 
    dataset name, and if so try to determine and return the dataset name

    100103    ksl    Coded because it is natural in some cases to give the filename
    '''
    xname=name
    if xname.count('.') > 0 or xname.count('/') > 0:
        # This must be a filename
        i=xname.rindex('/')
        xname=xname[i+1:i+10]
    
    return xname

def parse_creation_time(xtime='2010-07-10T19:11:52'):
    '''
    Parse the time string for the date at which a fits file
    was created and return a time that can be used to choose
    which of two fits files was created last

    The time returned approximates the number of days since
    1990 but does not account for the fact that there are
    different numbers of days in a month.  

    The routine returns 0.0 and raises an IndexErro exception 
    if the time cannot be parsed

    Note the use of int instead of eval because leading zeros
    caused eval to interprete the numbers as octal

    100817    Added error checking which arose because the values
        being passed to routine were actually not the
        creation date

    '''

    xtime=xtime.replace('-',' ')
    xtime=xtime.replace('T',' ')
    xtime=xtime.replace(':',' ')
    xtime=xtime.split()

    try:
        day_sec=int(xtime[3])*3600.+int(xtime[4])*60+int(xtime[5])
        day_frac=day_sec/86400
        # Next section for day is really not accurate, but we don't care
        year_frac=(int(xtime[1])+int(xtime[2])/31.)/12.
        day=365.*((int(xtime[0])-1990.)+year_frac)
    except IndexError:
        raise IndexError
        return 0.0

    pseudo_time=day+day_frac

    return pseudo_time

def check4duplicates(records):
    '''
    Check the list 4 duplicate records, and choose the one that was
    created last if that is possible

    The routine returns a list that contains 'ok' for files that 
    are the ones to use and 'nok' for those that are duplicates

    100817    Added checks to trap problems parsing times.
    110120    This is a new attempt to find the duplicates and select the last one
    160127  Modified again to speed this up (when there are large numbers of files
        to examine.  If there are no duplicates this is quite quick.  Even
        when there are duplicates this is about 3x faster than the old method.
    161208  Modify so uses tables
    '''


    xstart=time.time()

    # Remove the duplicates file if it exists to avoid confusion
    if os.path.isfile('duplicate_files.txt'):
        os.rename('duplicate_files.txt','duplicate_files.txt.old')

    names=list(records['Dataset'])
    print(names)
    print('What',len(names),len(records))
    unique=set(names)
    if len(unique)==len(records):
        print('check4duplicates: There are no duplicates in the directory structure')
        return records
    else:
        duplicate_names=[]
        hold_all=[]
        records2delete=[]
        print('check4duplicates: Warning: There are %d duplicate datasets in the directory structure' % (len(names)-len(unique)))
        for one in unique:
            icount=names.count(one)
            if icount>1:

                duplicate_names.append(one)  # these are the names of the duplicate datasets

                hold=[]
                j=0
                while j<len(names):
                    if names[j]==one:
                        hold.append(j)
                        if len(hold)==icount:
                            break   # We have them all
                    j=j+1

                # hold contains the record numbers of the records for this duplicate dataset

                hold_all=hold_all+hold # This indexes all of the duplicate dataset

                times=[]
                for one in hold:
                    times.append(records['File-date'][one])
                times=numpy.array(times)
                order=numpy.argsort(times) # Give me the order of the times
                # Order has index of the hold array in time order.  We want to keep the
                # last one, and delete the others
                last=order[len(order)-1]


                k=0
                while k<len(hold):
                    if k!=last:
                        records2delete.append(hold[k])
                    k=k+1
    duplicates=records[hold_all]

    duplicates.write('duplicate_files.txt',format='ascii.fixed_width_two_line')
    records.remove_rows(records2delete)
    print('Check for duplicates in direcory structure took:',time.time()-xstart)

    return records
    




def find_latest(filename='foo'):
    '''

    This is a simple routine to locate a specific version of a file.  If
    the file is in the current directory that is the file name that
    will be returned.  If it is in one of the subdirectories of the
    current directory then the one that was modified most recently
    will be returned. 

    Notes:

    This routine uses the unix utility find.  It is therefore
    likely to be slow in large directory structures

    This uses subprocess which handles stdin and stdout, unlike
    os.system

    History:

    100906 ksl Coding begun.  There is a standalone version of this 
           called find.py (in my normal py_progs/scripts directory


    '''

    proc=subprocess.Popen('find . -follow -name %s -print ' % filename,shell=True,stdout=subprocess.PIPE)
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







def check4scan(filename='./Visit43/ic9t43j1q_flt.fits[1]'):
    '''
    Determine whether or not a raw, ima, or flt file is associated with an observation involving  
    a spatial scan.  

    Return:
        scan    if a spatial scan
        stare   if the spt file exists, but it is not a spatial san
        unknown    if the spt file does not exist
    
    Notes:

    The routine looks for the spt file corresponding to the dataset and
    parses the header to find out if it is a scanned observation.
    
    130225  Coded and Debugged
    130307  Replaced routine using astropy.fits with iraf because astropy.fits
        was very slow
    '''

    xscan='unknown'

    # First strip of the extension if any
    xname=per_fits.parse_fitsname(filename,0,'yes')
    xfile=xname[0]
    xfile=xfile.replace('flt','spt')
    xfile=xfile.replace('raw','spt')
    xfile=xfile.replace('ima','spt')

    # Now we should have the name of the spt file
    if os.path.exists(xfile)==True:
        # Revert to iraf/pyraf for performance reasons
        xx=per_fits.get_keyword(xfile,0,'SCAN_TYP')
        # Elimimate pyraf dependencies because it is going away
        # xx=pyraf.iraf.hselect(xfile+'[0]','SCAN_TYP','yes',Stdout=1)
        # xx=xx[0].split('\t')
    else: 
        return 'no_spt'



    if xx[0]=='N':
        xscan='stare'
    else: 
        xscan='scan'


    return xscan  
    

# End utilities



def update_summary(dataset,status_word,keys=[],values=[],fileroot='observations',reinit='yes'):
    '''
    Update a record in the summary file, where the inputs have the following meaning
    
    dataset         rootname of the file, usually the flt file, being processed
    status_word     one word, e.g. Complete or Error, to define the current status of the processing
    keys            the database values that we want to update, e.g PerHtml
    values          A list of the values we want to put into the output table
    fileroot        The rootname of the summarry file, almost always obsevations
    reinit          if 'yes', then reinitialize to a default state the fields that new values have not been
                    included in the file


    Notes:

    This routine does not actually update the summary file, and so its name is a misnomer. It
    just records the information need for the update into small tables in temp_sum.  This is all
    set up this way to allow parallel processing. A new routine fixup_summary_file
    below now inserts the results into the data. 


    History:

    101221  ksl Coded as part of effort to get a way to check what files
                had been processed with the persistence software
    160105  ksl Modified for multiprocessing
    161212  ksl Recoded to use astropy tables

    '''


    if os.path.isdir('tmp_sum')==False:
        os.mkdir('tmp_sum')

    summary_file=fileroot+'.sum'
    gmt=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    words=gmt.split()


    try:
        lines=read_table(summary_file)
    except IOError:
        print('Error: update_summary: File %s does not exist' % summary_file)
        return


    # Read the old data from tmp_sum if it exists, or failing that take the
    # results from the old obs_sum file
    xfile='tmp_sum/%s.txt' % dataset
    if os.path.isfile(xfile)==True:
        old_results=read_table(xfile)
    else:
        i=0
        while i <len(lines):
            line=lines[i]
            if line['Dataset']==dataset:
                old_results=lines[i]
                break
            i=i+1
    
    # OK at this point I have the old results

    old_results['ProcStat']=status_word
    old_results['Proc-Date']=words[0]
    old_results['Proc-Time']=words[1]

    old_results=Table(old_results)  # Not sure why this is needed

    i=0
    while i<len(keys):
        if keys[i] in set(old_results.colnames):
            old_results[keys[i]]=values[i]
            if type(values[i])==float:
                old_results[keys[i]].format='8.4f'
        else:
            print ('Error: Update_summary: key %s not found' % keys[i])
        i+=1


    # At present this is a row of and old table and so we need to convert to a Table
    # in order to be able to write it out

    old_results=Table(old_results)
    old_results.write('tmp_sum/%s.txt' % dataset,format='ascii.fixed_width_two_line')

    return



def  fixup_summary_file(datasets,fileroot='observations'):
    '''
    Update the observations.sum file from the information stored in the 'tmp_sum' directory

    where 
        datasets is a list of the datasets which have been processed
        fileroot is the root name of the .sum file to be updated
    
    Notes:
        
        The directory where the data is stored is hardwired

        The routine simple finds the line in the observations.sum file associated with
        the dataset and replaces that line with the line that is in the .txt file
    
    History:

    160105  ksl Coded as part of the effort to add multiprocessing
    161204  ksl Update to handle tables.

    '''


    # Get the data
    data=[]
    for dataset in datasets:
        try:
            xname='tmp_sum/%s.txt' % dataset
            xdata=read_table(xname,format='fixed_width_two_line')
            if len(data)==0:
                data=xdata
            else:
                data=vstack([data,xdata])
        except IOError:
            print('Error: fixup_summary_file: %s does not appear to exist' % xname)
    


    gmt=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

    # Now read the old observations.sum file
    summary_file=fileroot+'.sum'
    try:
        xsum=read_table(summary_file,format='fixed_width_two_line')
    except  IOError:
        print('Error: update_summary : Could not read summary file %s ' % summary_file)
        return


    # Assume that there are fewer lines to update than there are in
    # the obsum file  (xsum is the old summary file)

    i=0
    while i<len(data):
        j=0
        while j<len(xsum):
            if data['Dataset'][i]==xsum['Dataset'][j]:
                break
            j+=1
        if j<len(xsum): # Then we have to parse the ith datasrig and insert it
            xsum[j]=data[i]
        i+=1


    backup(summary_file)
    xsum.write(summary_file,format='ascii.fixed_width_two_line')

    return



def make_sum_file(fileroot='observations',new='no'):
    '''
    Make the file that will record the results of persistence processing.  If
    it already exists, give the user the option of merging new records into
    the old file or creating a new file

    Note that each new output line should have the following values

    dataset    prog_id MJD  ProcessDate Unprocessed


    101221    ksl    Coded as part of effort to put a recording mechanism
            in place
    110103    ksl    Modified the outputs when new records are inserted
    110103    ksl    Changed so the routine itself reads the observation.ls
            file
    110119    ksl    Rewrote to assure that there is exactly one summary 
            file line for each observation file line
    '''



    # Read the entire observations.ls file
    print('# (Re)Making summary file')
    x=read_ordered_list0(fileroot)
    print('# The number of records in the old per_list file is %d' % (len(x)))


    gmt=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    words=gmt.split()
    proc_date=words[0]
    proc_time=words[1]

    summary_file=fileroot+'.sum'


    # Create pristine file

    g=x['Dataset','ProgID','ExpStart']
    g['Proc-Date']=proc_date
    g['Proc-Time']=proc_time
    g['ProcStat']='Unprocessed'
    g['E0.10']=-99.
    g['E0.03']=-99.
    g['E0.01']=-99.
    g['I0.10']=-99.
    g['I0.03']=-99.
    g['I0.01']=-99.
    g['PerHTML']='--'

    # The next steps are necessary, whenever the an ascii table is
    # read in, and there is a chance the updates will result in 
    # strings which are longer than the current format.  Withot
    # this the strings will be truncated in the ouput table
    for col in g.itercols():
        if col.dtype.kind in 'SU':
               g.replace_column(col.name, col.astype('object'))




    if os.path.exists(summary_file)==False or new=='yes':
        print('# Making a pristine summary file')
        g.write(summary_file,format='ascii.fixed_width_two_line')

    else:
        print('# Merging new records into old list')
        try:
            xsum=read_table(summary_file,format='fixed_width_two_line')
            if 'ProcStat' in xsum.colnames == False:
                print ('Error: make_summary_file: %s read but does not have correct table format' % summary_file)
                return
        except:
            print('Error: make_sum_file: Could not read summary file %s ' % summary_file)
            return



        j=0
        for one in g:
            i=0

            while i<len(xsum):
                if xsum['Dataset'][i]==one['Dataset']:
                    g[j]=xsum[i]
                    break
                i+=1
            j=j+1

        g.write('tmp.sum.txt',format='ascii.fixed_width_two_line')


        ndup=check_sum_file('tmp.sum.txt',summary_file)

        if ndup>0:
            print('Error: make_sum_file: Since there were duplicates in th tmp file, not moving to %s'  % (summary_file))
            return 'NOK  - Since there were duplicates in th tmp file, not moving to %s'  % (summary_file)

        # Now move the files around.  Note that the next 3 lines need to be the same as in update_summary above
        # gmt=time.strftime("%y%m%d.%H%M", time.gmtime())  # Create a string to use to name the updated file. As written a new file a minute
        # proc=subprocess.Popen('mv %s %s.%s.old' % (summary_file,summary_file,gmt),shell=True,stdout=subprocess.PIPE)
        backup(summary_file)
        proc=subprocess.Popen('mv %s %s' % ('tmp.sum.txt',summary_file),shell=True,stdout=subprocess.PIPE)
    return


def check_sum_file(new='tmp.sum',old='none'):
    '''
    Check a summary file for unique dataset names

    Description

    The routine checks the new summary file to assure 
    that there is only one line with a given dataset name.
    It returns the number of duplicate records

    If a file name is given in old, then the routine checks
    whether there are any records that were in the old file, that
    are not in the new file.  This can happen and not be an error,
    but the results are recorded in a file with the extension problem

    Notes:

    There should be one and only one record  for each
    record in observaitions.ls

    110119 - Expanded the questions that were asked 
    150127 ksl Rewrote so that it was much faster by using sets
    '''

    # Read the new summary file
    summary_file=new
    if os.path.exists(summary_file)==False:
        print('Error: check_sum: There is no summary file names %s.sum to check' % fileroot)
        return
    f=open(summary_file,'r')
    lines=f.readlines()
    f.close()

    
    # Open a diagnostic file to record any problems
    g=open_file(summary_file+'.problem')
    string='# Checking new (%s) against (%s) old summary file' % (new,old)
    print(string)
    g.write('%s\n' % string)

    # Check for duplicate datasets in new summary file

    xstart=time.time()
    names=[]
    for line in lines:
        x=line.strip()
        word=x.split()
        names.append(word[0])

    unique=set(names)
    if len(unique)==len(names):
        dups=0
    else:  # There is a problem and so the duplicate check must be carried out
        print('Error: There are duplicate datasets in the the newly created summary file %s' % new)
        print('A detailed check is now underway')
        dup_names=[]
        dups=0
        for one in unique:
            jj=names.count(one)
            if jj>1:
                g.write('Duplicate: %5d %s\n' % (jj,one))
                dup_names.append(one)
                dups=dups+1

    string='# Check of %s revealed %d duplicates' % (new,dups)
    print(string)
    g.write('%s\n' % string)
    delta_t=time.time()-xstart
    print('check_sum: time to search for duplicate records in the new .sum file: ',delta_t)

    # Check for datasets that were in the .sum file before but are missing now

    if old!='none':
        print('Checking for datasets that were in the last .sum file but are now missing')
        xstart=time.time()
        old_names=[]

        try:
            f=open(old,'r')
            old_lines=f.readlines()
            f.close()

            for one in old_lines:
                 x=one.strip()
                 word=x.split()
                 old_names.append(word[0])
        except IOError:
            old_lines=[]

        nlost=0
        for one in old_names:
            if names.count(one)==0:
                g.write('Missing: %s\n' % (one))
                nlost=nlost+1

        string='# Check of %s revealed %d datasets no longer in the .sum file' % (old,nlost)
        if nlost>0:
            print('Warning: Lost datasets may not be a problem, but are noted as a warning of something possibly worth invesigation')
        print(string)
        g.write('%s\n' % string)
        print('# check_sum_file: The check for missing records took %s s' % delta_t)

    g.close()


    return dups





def get_info(lines,filetype):
    '''
    Get all of the keyword information for a set of files and return this

    Notes:

    This section of the old make_ordered list was put into a separate function
    so that it could be run in parallel 

    History

    160118    ksl    Added
    '''

    times=[]

    if len(lines)==0:
        print('There were no %s files in the directory structure' % filetype)
        return []
    else:
        print('There are %d datasets to process' % len(lines))

    i=0
    filename=[]
    dataset=[]
    propid=[]
    line_no=[]
    instrument=[]
    detector=[]
    expstart=[]
    date_obs=[]
    time_obs=[]
    aperture=[]
    xfilter=[]
    exptime=[]
    ra=[]
    dec=[]
    target=[]
    assn=[]
    investigator=[]
    file_create=[]
    file_mod=[]


    for one in lines:

        xfile=one['File']

        if os.path.isfile(xfile) == True:
            x=per_fits.get_keyword(xfile,1,'rootname,proposid,linenum, instrume,detector,expstart,date-obs,time-obs,aperture,filter,exptime,crval1,crval2,targname,asn_id,pr_inv_L')
            if x[4]=='IR':
            # for raw files the date is in extension 0, but for the others it is in 1.
                if filetype=='raw':
                    date=per_fits.get_keyword(xfile,0,'date')
                else:
                    date=per_fits.get_keyword(xfile,1,'date')
            # Check for scan returns scan,stare, or unkown
                scan=check4scan(xfile)
                x[4]=scan

                filename.append(xfile)
                dataset.append(x[0])
                propid.append(x[1])
                line_no.append(x[2])
                instrument.append(x[3])


                detector.append(x[4])
                expstart.append(x[5])
                date_obs.append(x[6])
                time_obs.append(x[7])


                aperture.append(x[8])
                xfilter.append(x[9])
                exptime.append(x[10])
                ra.append(x[11])
                dec.append(x[12])
                target.append(x[13])
                assn.append(x[14])
                investigator.append(x[15])

                file_create.append(date[0])
                file_mod.append(one['Mod-time'])


        else: 
            print('File %s does not really exist' %  xfile)
        i=i+1
        if i%100 == 1:
            print('Inspected %6d of %6d datasets --> %6d IR datasets' % (i,len(lines),len(dataset)))

    # Now put the results back into a table
    print('Inspected %6d of %6d datasets --> %6d IR datasets' % (i,len(lines),len(dataset)))

    x=Table()
    x['File']=filename
    x['Dataset']=dataset
    x['ProgID']=propid
    x['LineNo']=line_no
    x['Inst']=instrument
    x['Det']=detector
    x['ExpStart']=expstart
    x['Date-Obs']=date_obs
    x['Time-Obs']=time_obs
    x['Aper']=aperture
    x['Filter']=xfilter
    x['Exptime']=exptime
    x['RA']=ra
    x['Dec']=dec
    x['Target']=target
    x['Assoc']=assn
    x['PI']=investigator
    x['File-date']=file_create
    x['Mod-time']=file_mod


    # Fixup the formats so strings are objects
    # This is important to prevent truncation of strings
    for col in x.itercols():
        if col.dtype.kind in 'SU':
            x.replace_column(col.name,col.astype('object'))

    return x


def info_helper(args):
    '''
    This repackages the argments so that they can be used in parallel processing

    Notes:

    The issue here is that you only can pass a single argument via the Pool mechanism.
    This simple allows that.
    '''

    return get_info(*args)



def make_ordered_list(fileroot='observations',filetype='flt',use_old='yes',np=1):
    '''
    find all of the observations in all subdiretories and make
    a time ordered list of the observations from files of a
    given filetype, e. g. flt.  If apertures == 'full' only
    put full frames into the list.  Otherwise, include subarrays



    Note - This could be done as a true database, but it's 
    simpler for now just to make it a file


    100308 - Added option to deal with other types of files aside from flt files.  The 
        assumption made is that the filetype is part of the filename
    100427 - Added crval1 and crval2 to output list because want to use this to determine
        whether there was a dither or move between observations
    100817 - Added checks to verify that a file found by find actually exists and has the
        appropriate image extension (namely 1)
    110103    - Added the switch new as a stepping stone to making the routine more general. AT
        Some point we are going to need to merge new files into old ones and that is
        partially the reason for this.  Right now it is simply to avoid being asked
        if one wants to make a new file
    110104    Split off the actual writing of the file into a separate routine called
        write_time_sorted in order to make sublists.  Any change in the information 
        that is gathered by make_ordered list, needs to be relected in that routine
    110121    Added a line to assure that PI names had no spaces.  This was because running
        split on records caused with spaces int he PI name was failing.  The alternative
        which was to change the separated character to a tab seemed more trouble
    120330    Added a kludge to add a date for a raw data file.  The prupose of this date is 
        resolve situations in which multiple versions of the same file are in the
        directory structure.  The choice is only accurate to the day
    130820    Revised the way in which searched for the date in raw data files.
    140306    Add capability to check for scanned observations.
    160118  Added the possibility of running multiple processors
    160119  Switched from time.clock to time.time so that everything would be in wall clock
        time
    '''

    if filetype!='flt':
        fileroot=fileroot+'_'+filetype


    backup(fileroot+'.ls')
    xstart=time.time()
    os.system('find . -follow -name \*%s.fits  -print 1>files.ls 2>/tmp/foo; rm /tmp/foo' % filetype)
    search_time=dtime=time.time()-xstart
    print('# Found all of the files to read in %f s' % dtime)
    xstart=time.time()

    # lines=ascii.read('files.ls',format='no_header')
    lines=read_table('files.ls',format='no_header')
    lines.rename_column('col1','File')

    # Get the MOdification date for all of the files
    mod_date=[]
    for one in lines:
        x=os.path.getmtime(one['File'])
        mod_date.append(time.strftime('%Y-%b-%d-%H:%M:%S', time.gmtime(x)))
    lines['Mod-time']=mod_date

    old_lines=[] # This just sets up defaults
    new_lines=lines

    if use_old=='yes':
        try:
            # obs_old=ascii.read(fileroot+'.ls')
            obs_old=read_table(fileroot+'.ls')
            # Next couple of lines can be deleted once fully transitioned to new version
            # of per_list
            xset=set(obs_old.colnames)
            if 'ModDate' in xset == False:
                use_old='no'
        except IOError:
            use_old='no'


    if use_old=='yes':
        xjoin=join(lines,obs_old,keys=['File'],join_type='left')
        old=[]
        new=[]
        i=0
        while i<len(xjoin):
            one=xjoin[i]
            if one['Mod-time_1']==one['Mod-time_2']:
                old.append(i)
            else:
                new.append(i)
            i+=1
        if len(old)>0:
            old_lines=xjoin[old]
            old_lines.rename_column('Mod-time_2','Mod-time')
            old_lines.remove_column('Mod-time_1')
        if len(new)>0:
            new_lines=lines[new]
        else:
            new_lines=[]

    #  At this point we know which files need to be read
    #  Note that an explicit assumption here is that the spt files have not changed
    #  Ignore this for now

    print('Of %d files, %d are old, and %d are new' % (len(lines),len(old_lines),len(new_lines)))


    if len(new_lines)==0:
        pass
    elif np<=1 or len(lines)<np:
        records=get_info(new_lines,filetype)
    else:
        inputs=[]
        idelta=len(new_lines)/np +1
        i=0
        while i < len(new_lines):
            imin=i
            imax=i+idelta
            if imax>len(new_lines):
                imax=len(new_lines)
            xinputs=lines[imin:imax]
            inputs.append([xinputs,filetype])
            i=i+idelta


        p=Pool(np)  

        i=0
        # 180811 - one will have length 0 if all
        # of the datasets for a particular processor 
        # have no IR components
        for one in p.map(info_helper,inputs):
            if len(one)>0:
                if i==0:
                    records=one
                else:
               	    records=vstack([records,one])

                i+=1
    
    # At this point, records contains all of the new_records
    if len(old_lines)>0 and len(new_lines)>0:
        records=vstack([old_lines,records])
    elif len(old_lines)>0:
        records=old_lines


    key_time=dtime=time.time()-xstart
    print('Inspected all of the files in %f s' % dtime)
    

    if len(records)==0:
        print('There were no IR observations to consider')
        return []
    # Now sort this all on the time
    # This returns an index of the order of the lines

    records.sort('ExpStart')

    records=check4duplicates(records)

    records.write(fileroot+'.ls',format='ascii.fixed_width_two_line')
    


    # Now write the time sorted file
    write_time=time.time()-xstart
    print('# Time to write the .ls file  %s s' % write_time)

    print('# Completed creating %s.ls' % (fileroot))

    print('# make_ordered_list: times',search_time,key_time,write_time)

    return records

def read_ordered_list0(fileroot='observations'):
    '''
    This simply reads the list and returns all of the records

    Eventually should replace what is in read_ordered_list and
    read_ordered_list2

    This is the table version

    History

    170208  ksl Fixed a problem where the line id was being interpreted
                as a float
    '''

    try:
        x=Table.read(fileroot+'.ls',format='ascii.fixed_width_two_line')
    except:
        print('Error: read_ordered_list0: Could not open %s ' % (fileroot+'ls'))
        return []

    if x['LineNo'].dtype.kind=='f':
        xline=[]
        for one in x['LineNo']:
            xline.append('%05.3f' % one)
        x.replace_column('LineNo',xline)


    return(x)

def read_ordered_list2(fileroot='observations',dataset='ibel01p4q',interval=[-1,2],outroot='none'):
    '''
    Given a dataset name and an interval, return the datasets that
    were obtained in the interval around the data set.

    if outroot is anything but none, the retrieved information will also be written to
    a file called outroot+.ls

    111019    ksl    Removed lines which read entire file, and replaced with call to read_ordered_list0
    '''

    # Read the entire file

    x=read_ordered_list0(fileroot)


    # Check the dataset name
    dataset=parse_dataset_name(dataset)

    # print 'Looking for ',dataset
    # locate the record with a given dataset name
    izero=0
    while izero< len(x):
        record=x[izero]
        if dataset==record['Dataset']:
            zero_time=record['ExpStart']
            break
        izero=izero+1

    if izero==len(x):
        print('Error: read_ordered_list2: Could not locate record for dataset %s in %s.ls' % (dataset,fileroot))
        return []


    # locate the datasets within the interval
    istart=-1
    istop=-1
    i=0
    times=[]
    while i<len(x):
        time=x['ExpStart'][i]
        dt=(time-zero_time)*24  # convert MJD to hours
        if dt >= interval[0] and istart == -1:
            istart=i
        if istart!=-1 and dt <= interval[1]:
            times.append(dt)
            istop=i
        i=i+1

    # print 'Check lengths: ',len(times),len(records[istart:istop+1])

    xxx=x[istart:istop+1]


    if outroot!='none' and outroot != fileroot:
        xxx.write(outroot+'.ls',format='ascii_fixed_width_two_line')


    return xxx 



def read_ordered_list(fileroot='observations',dataset='last',delta_time=24):
    '''
    read the file written by the make_ordered_list routine, where
    fileroot is the rootname of the fileroot+'.ls' file which contains
    the all of the infomration about each dataset, datset is the
    dataset for which one is searchin, and delt_time is the time in hours
    preceding this dataset for which you would like information.

    if dataset='last', then the assumption is that one is looking for the
    very last dataset in the file.


    If the fileroot+'.ls' file is missing, or the dataset is not in the 
    file, the routine returns an empty list

    110121    ksl    Modified so that returns empty array if a dataset is request
            which does not exist 
    111019    ksl    Replaced the section that reads the entire file
    '''

    # Read the entire file
    x=read_ordered_list0(fileroot)


    if dataset!='last':
        # locate the record with a given dataset name
        ilast=0
        while ilast< len(x):
            one=x[ilast]
            if dataset==one['Dataset']:
                # print 'Found %s at record %d' % (dataset,ilast)
                end_time=one['ExpStart']
                break
            ilast=ilast+1
        if ilast==len(x):
            print('Error: read_ordered_list: Did not find dataset %s, returning' % dataset)
            return []
    else:  # Use the last record as the endpoint
        ilast=len(x) - 1
        end_time=x['ExpStart'][ilast]

    # Locate the first reciord we want

    if delta_time>0:
        delta_time=delta_time/24.
        # Locate the first record we want
        ifirst=0
        dt=end_time-x['ExpStart'][ifirst]
        while dt > delta_time and ifirst < ilast:
            dt=end_time-x['ExpStart'][ifirst]
            ifirst+=1
        return x[ifirst:ilast+1]
    else:
        return x[ilast]



def read_ordered_list_mjd(fileroot='observations',mjd_start=0,mjd_stop=0):
    '''
    Read an observation file and return all records between
    mjd_start and mjd_stop.

    If mjd_start=0, start from the beginning of the list.
    if mjd_stop=0,  continue to the end of the list

    101109    ksl    Coded and debugged
    '''

    x=read_ordered_list0(fileroot)

    nrec=len(x)

    print('The total number of records is ',nrec)

    i=0
    if mjd_start>0:
        while i < nrec and mjd_start > x['ExpStart'][i]:
            i=i+1
    ifirst=i

    if mjd_stop==0:
        ilast=nrec
    else:
        while i<nrec and mjd_stop >= x['ExpStart'][i]:
            i=i+1
        ilast=i

    print('Retrieving %d records from %d to %d ' % (ilast-ifirst,ifirst,ilast))
    
    return x[ifirst:ilast]


def read_ordered_list_progid(fileroot='observations',prog_id=11216,mjd_start=0,mjd_stop=0):
    '''
    Read a list and return the records for a given program id, and if mjd_start and 
    mjd_stop are not both 0, between these two values

    If both mjd_start and mjd_stop are 0, then no tme check is done and all of the 
    records for that program will be returned.

    111019    ksl    Modified so that if prog_id is 0 or less then everything is returned
            In this case the routine behaves exactly like red_ordered_list_mjd
    '''


    #xxOLD try:
    #xxOLD     x=Table.read(fileroot+'.ls',format='ascii.fixed_width_two_line')
    #xxOLD except:
    #xxOLD     print('Error: read_ordered_list0: Could not open %s ' % (fileroot+'.ls'))
    #xxOLD     return []


    x=read_ordered_list0(fileroot)

    if prog_id<=0:
        pass
    else:
        select=[]
        i=0
        while i<len(x):
            if x['ProgID'][i]==prog_id:
                select.append(i)
            i+=1
        if len(select)==0:
            return []
        x=x[select]

    if mjd_start==0 and mjd_stop==0:
        return x


    i=0
    nrec=len(x)
    if mjd_start>0:
        while i < nrec and mjd_start > x['ExpStart'][i]:
            i=i+1
    ifirst=i

    if mjd_stop==0:
        ilast=nrec
    else:
        while i<nrec and mjd_stop >= x['ExpStart'][i]:
            i=i+1
        ilast=i

    x=x[ifirst:ilast]
    if len(x)==0:
        print('Although %d were found for prog_id %d, none between mjd %e and %e' % (len(xrec),prog_id,mjd_start,mjd_stop))
        return []
    
    return x




def read_ordered_list_one(fileroot='observations',dataset='last'):
    '''
    return a single record corresponding to the dataset name from 
    the ordered list whose root name is given by fileroot

    Note that if the dataset is not found, or the observaions.ls
    file is missing, the routine returns an empty list.

    101214    ksl    Added because it is needed in situations where
            one wants to work on a single dataset, but we
            need information about where everything is stored
    '''
    record=read_ordered_list(fileroot,dataset,delta_time=0)
    return record


def steer(argv):
    '''
    This is a steering program for per_list

    110103    ksl    Added as the options for per_list became more
            complicated, see top level for details
    120330    ksl    Added a fix so make_sum_file would receive
            a rootname that it could use in instances
            where the ftype was not flt

    '''

    ftype='flt'
    use_old='yes'
    new_summary_file='no'
    root='observations'
    np=1


    i=1

    while i<len(argv):
        if argv[i]=='-h':
            print(__doc__) 
            return
        elif argv[i]=='-new_all':
            use_old='no'
            new_summary_file='yes'
        elif argv[i]=='-new_sum':
            new_summary_file='yes'
        elif argv[i]=='-new_ls':
            new_ls_file='yes'
        elif argv[i]=='-daily':  # This is the standard switch crontab generation
            use_old='yes'
            new_summary_file='no'
        elif argv[i]=='-file_type':
            i=i+1
            ftype=argv[i]
        elif argv[i]=='-np':
            i=i+1
            np=int(argv[i]) 
        else:
            if i != len(argv)-1:
                print('Could not understand argument %d :%s' % (i,argv[i]))
                return
            else:
                root=argv[i]
        i=i+1

    if  use_old=='no':
        print('Rebuilding the .ls file from scratch')
    else:
        print('Rebuilding the .ls file using information the previous .ls file')
    if new_summary_file=='yes':
        print('Rebulding the summary file such that information from earlier persistence processing will be "forgotten"')
    else:
        print('Interpolating old information about persistence processing into the new .sum file')

    # At this point we have fully parsed the observation list

    xstart=time.time()
    make_ordered_list(root,ftype,use_old,np)
    dtime=time.time()-xstart
    print('# Time to make the .ls file: ',dtime)

    # 120330 - kludge of a change to get to work for file types other than flt.  It's not 
    # obvious to me what would make this easier though as make ordered list has to do
    # a file search for files of a specific type.
    if ftype!='flt':
        root=root+'_'+ftype

    dtime=time.time()-xstart
    make_sum_file(root,new_summary_file)
    dtime=time.time()-xstart
    print('# Time to make the .sum file: ',dtime)

    return

        


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":

    import sys

    steer(sys.argv)
    # print 'Got back here to main.  Obscure IRAF related error occurs here on some systems'
    
