#! /usr/bin/env python

'''
This is the python wrapper for the WFC3 IR Persistence software.
'''

from pyql.logging.logging_functions import configure_logging
from pyql.logging.logging_functions import log_info
from pyql.logging.logging_functions import log_fail

from datetime import datetime
from datetime import timedelta
import logging
import os
import subprocess


class Chdir(object):
    '''
    Exception-safe implementation of the chdir() function. 
    '''         
    def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)

    def __del__( self ):
        os.chdir( self.savedPath )


@log_fail
@log_info
def persistence_wrapper_main():
    '''
    The main controller.
    '''
    foo = Chdir('/grp/hst/wfc3a/persistence/workspace')
    logging.info(os.getcwd())

    # Run per_list.py
    logging.info('Running: python /grp/hst/wfc3a/persistence/code/per_list.py')
    subprocess.call(['python','/grp/hst/wfc3a/persistence/code/per_list.py'])
    logging.info('per_list.py complete')

    # Run run_persist.py
    five_days_ago = datetime.strftime((datetime.today() - timedelta(days=5)), '%Y-%m-%d')
    today = datetime.strftime((datetime.today()), '%Y-%m-%d')
    logging_string = 'Running: '
    logging_string += 'python '
    logging_string += '/home/long/WFC3/corvina/py_progs/scripts/run_persist.py '
    logging_string += '-start {} -stop {}'.format(five_days_ago, today)
    logging.info(logging_string)
    subprocess.call(['python', 
                    '/home/long/WFC3/corvina/py_progs/scripts/run_persist.py',
                    '-start', five_days_ago , '-stop', today])
    logging.info('run_persist.py complete') 


if __name__ == '__main__':
    configure_logging('persistence_wrapper')
    persistence_wrapper_main()
