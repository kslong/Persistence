#! /usr/bin/env python

'''
This is the python wrapper for the WFC3 IR Persistence software.
'''

from datetime import datetime
from datetime import timedelta

from pyql.logging.logging_functions import configure_logging
from pyql.logging.logging_functions import log_info
from pyql.logging.logging_functions import log_fail

from subprocess import PIPE
from subprocess import Popen

import logging
import os

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

    code_path = '/grp/hst/wfc3a/software_git/automated_scripts/cal_ir_make_persistence'

    # Run per_list.py
    logging.info('Running: python {}'.format(os.path.join(code_path, 'per_list.py')))
    child_process = Popen(['python', os.path.join(code_path, 'per_list.py')],
                         stdout=PIPE, bufsize=1)
    for line in iter(child_process.stdout.readline, b''):
        print(line,)
        #logging.info(line,)
    child_process.communicate() 
    logging.info('per_list.py complete')

    # Run run_persist.py
    five_days_ago = datetime.strftime((datetime.today() - timedelta(days=5)), '%Y-%m-%d')
    today = datetime.strftime((datetime.today()), '%Y-%m-%d')
    logging_string = 'Running: '
    logging_string += 'python '
    logging_string += '{} '.format(os.path.join(code_path, 'run_persist.py'))
    logging_string += '-start {} -stop {}'.format(five_days_ago, today)
    logging.info(logging_string)
    child_process = Popen(['python', 
                    os.path.join(code_path, 'run_persist.py'),
                    '-start', five_days_ago , '-stop', today], stdout=PIPE, 
                    bufsize=1)
    for line in iter(child_process.stdout.readline, b''):
        logging.info(line,)
    child_process.communicate() 
    logging.info('run_persist.py complete')

    # Run subtract_sum.py
    logging.info('Running: {}'.format(os.path.join(code_path, 'subtract_sum.py')))
    child_process = Popen(['python', os.path.join(code_path, 'subtract_sum.py')])
    for line in iter(child_process.stdout.readline, b''):
        logging.info(line,)
    child_process.communicate()    
    logging.info('subtract_sum.py complete')


if __name__ == '__main__':
    configure_logging('persistence_wrapper')
    persistence_wrapper_main()
