#!/usr/bin/env python
"""
Main program for fitting 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import sys 
from tls import * 

def main():
    args = sys.argv[1:]
    if args[0] == 'angular2D':
        return angular2D(args[1:])
    else:
        raise NameError(args)

def angular2D(args):
    datatype = args[0]
    label = args[1]
    test = option_exists(args, '-t')
    batch = option_exists(args, '-b')
    figname = 'test_fig' 

    if batch:
        jobname = 'fit'
        queue = '1nd'
        cmd = create_batch_cmd(main='./fit.py')
        bashfile = set_file('./', label, figname, '.sh', test=test)
        pwd = os.getcwd()
        pre= 'cd %s' % pwd 
	update_bashfile_cmd(bashfile, cmd, pre=pre, test=test)
        logfile = set_file('./', label, figname, '.log', test=test)
        if jobname and '[' in jobname and ']' in jobname:
            logfile += '.%I'

        test = True 
        bsub_jobs(logfile, jobname, bashfile, test=test, queue=queue)
       
        
if __name__ == '__main__':
    main()

    
