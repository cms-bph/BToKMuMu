"""
Tools for the A_FB analysis scripts. 

"""

__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


import os
import sys
#import atr # tools shoul be independent from attribute 
import subprocess
import json
import logging

#from progressbar import Bar, Percentage, ProgressBar
from filetools import *
import tls

# ----------------------------------------------------------
# Functions
# ----------------------------------------------------------
def main(args): 
    if args[0] == 'convert_bytes':
        function = getattr(tls, args[0])
        bytes_ = args[1] 
        sys.stdout.write('%s bytes = %s\n' %(bytes_, function(bytes_)))
    
    elif args[0] == 'duration':
        function = getattr(tls, args[0])
        seconds = float(args[1])
        sys.stdout.write('%s seconds = %s \n' %( seconds, function(seconds)))
    
    else:
        raise NameError(args)

def bsub_jobs(logfile, jobname, bash_file, test=True, queue='1nd') :
    if os.access(logfile, os.F_OK) and not test and logfile != '/dev/null':
        os.remove(logfile)


    if jobname and '[' in jobname and ']' in jobname:
        logfile += '.%I'

    cmd = 'bsub -q %s -J %s -o %s  %s ' %( 
        queue, jobname, logfile, bash_file) # no need to use "" for subprosess.  

    stdout = proc_cmd(cmd, test)
    if stdout is not None:
        sys.stdout.write(stdout)


def canvas_output(canvas, figname, label):
    canvas.SaveAs(figname)
    return 
    

def check_job_log(stage, datatype, label, job):
    logfile = set_logfile(stage, datatype, label, 'evt', job)
    parsed = LogFile(logfile)

    if stage == 'skm':
        passed = parsed.skimmed
    elif stage == 'sel':
        passed = parsed.selected
    elif stage == 'gen':
        passed = 'N/A'
    elif stage == 'ntp':
        passed = 'N/A'
    else:
        raise NameError(stage)
    
    sys.stdout.write('| %s |  %s | %s | %s |\n' % (
            job, parsed.processed, passed, parsed.duration))

    return job, parsed.processed, passed, parsed.duration


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
        

def convert_bytes(bytes_):
    bytes_ = float(bytes_)
    if bytes_ >= 1099511627776:
        terabytes_ = bytes_ / 1099511627776
        size = '%.2fT' % terabytes_
    elif bytes_ >= 1073741824:
        gigabytes_ = bytes_ / 1073741824
        size = '%.2fG' % gigabytes_
    elif bytes_ >= 1048576:
        megabytes_ = bytes_ / 1048576
        size = '%.2fM' % megabytes_
    elif bytes_ >= 1024:
        kilobytes_ = bytes_ / 1024
        size = '%.2fK' % kilobytes_
    else:
        size = '%.2fb' % bytes_
    return size

def create_batch_cmd(job=None, main='afb'):
    cmd = sys.argv[1:] 
    cmd.remove('-b')
    if job is not None:
        range_idx = cmd.index('submit') + 1
        cmd[range_idx] = str(job)

    cmd = main + ' ' + ' '.join(cmd)
    return cmd

    

def create_bash_file_pyline(label, datatype, pyline,
                            bashname, optcommand='', subdir='', test=False):
    script_dir =  os.path.join(atr.afbpath, 'src', 'sh')

    bash_content =  '''#!/bin/sh
date
hostname

.  ~/.bashrc
setafb 

%s

python -c '%s'

date
''' % (optcommand, pyline)
    
    bash_file = os.path.join(script_dir, bashname)

    verbose = 0
    if test:
        verbose = 1

    f = UserFile()
    f.append(bash_content) 
    f.output(bash_file, verbose=verbose)
    os.chmod(bash_file, 0755)
    return bash_file


def create_jobs_db(stage, datatype, label, test=False, size_max='2G',
                   input_label=None):
    if input_label is None:
        input_label= label
    
    if stage == 'srm':
        paths = atr.srmpath(input_label, datatype)
        lscmd = 'cmsLs' # 'lcg-ls -l'
        # if label in ['B2KstarMuMu/RECO_1M_v2.4']:
        #     lscmd = 'cmsLs'

    elif stage == 'ntp':
        paths = atr.ntppath(input_label, datatype)
        lscmd = 'rfdir'
        if label == '5ifbv2.8.1':
            lscmd = 'ls -l'

    elif stage == 'skm':
        skmpath = atr.skmpath(input_label)
        paths = [os.path.join(skmpath, datatype, input_label)]
        lscmd = 'rfdir'

    else:
        raise NameError(stage)

    dbpath = os.path.join(atr.datpath, stage, 'db', label)
    dbname = 'jobs_%s.db' % datatype
    if test:
        dbname += '.test'
    dbfile = check_and_join(dbpath, dbname)

    db = open(dbfile, 'w')
    path_name_size_list = []

    for p in paths:
        if p.startswith('/castor/'):
            lscmd = 'rfdir'
        cmd = '%s %s' % (lscmd, p)
        output = proc_cmd(cmd)

        if 'lcg-ls' in lscmd:
            f = LcgFile(data=output.split('\n'))
            name_size_list = f.name_size_list

        elif 'rfdir' in lscmd:
            name_size_list = parse_ls(output, ignore='.test')

        elif 'ls -l' in lscmd:
            name_size_list = parse_ls(output, ignore='total')
            print 'done. ', name_size_list

        elif 'cmsLs' in lscmd:
            name_size_list = parse_ls(output, poz=(1,4))

        else:
            raise NameError(lscmd)

        for name, size in name_size_list:
            filename = os.path.join(p, name)
            path_name_size_list.append((filename, size))

    rootfiles = group_files_by_size(path_name_size_list, size_max)

    job = 0
    jobs = {}
    for rootfile in rootfiles:
        job += 1
        jobs[job] = {}
        jobs[job]['rootfile'] = rootfile
        jobs[job]['status'] = 'created'
        
    json.dump(jobs, db)
    db.close()

    check_update_status(dbfile)
    sys.stdout.write('created %s jobs.\n' %len(jobs))
    return dbfile 

    
def duration(seconds):
    seconds = long(round(seconds))
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    years, days = divmod(days, 365.242199)
 
    minutes = long(minutes)
    hours = long(hours)
    days = long(days)
    years = long(years)
 
    duration = []
    if years > 0:
        duration.append('%d year' % years + 's'*(years != 1))
    else:
        if days > 0:
            duration.append('%d day' % days + 's'*(days != 1))
        if hours > 0:
            duration.append('%d hour' % hours + 's'*(hours != 1))
        if minutes > 0:
            duration.append('%d minute' % minutes + 's'*(minutes != 1))
        if seconds > 0:
            duration.append('%d second' % seconds + 's'*(seconds != 1))
    return ' '.join(duration)


def fourvecs_xyzm(px, py, pz, m):
    from ROOT import TLorentzVector
    p4 = TLorentzVector()
    p4.SetXYZM(px, py, pz, m)
    return p4 


def get_name_from_label(label, ver=None):
    if ver is not None:
        label = get_label_by_version(label, ver)
        
    com_name = label.replace('/', '_')
    com_name = com_name.replace('.', '_')
    com_name = com_name.replace('-', '_')
    return com_name


def get_crab_cfg_template(label):
    cfg = UserFile()

    if label in ['5ifbv2.3', '5ifbv2.4']:
        cfg.data = ['[CMSSW]\n', 'total_number_of_lumis  = -1 \n',
                    'number_of_jobs   = \n',
                    'pset             = \n',
                    'datasetpath      = \n',
                    'lumi_mask        = \n',
                    'output_file      = \n',
                    '[USER]\n',
                    'return_data      = 0 \n',
                    'email            = Xin.Shi@cern.ch  \n',
                    'copy_data        = 1 \n',
                    'storage_element  = srm-cms.cern.ch \n',
                    'storage_path     = /srm/managerv2?SFN=/castor/cern.ch \n',
                    'user_remote_dir  = /user/x/xshi/afb/dat/srm/ \n',
                    '[CRAB] \n',
                    'scheduler        = glite \n',
                    'jobtype          = cmssw \n',
                    'use_server       = 1 \n' ]
    else:
         cfg.data = ['[CMSSW]\n', 'total_number_of_lumis  = -1 \n',
                    'lumis_per_job    = \n',
                    'pset             = \n',
                    'datasetpath      = \n',
                    'lumi_mask        = \n',
                    'output_file      = \n',
                    '[USER]\n',
                    'return_data      = 0 \n',
                    'email            = Xin.Shi@cern.ch  \n',
                    'copy_data        = 1 \n',
                    'storage_element  = srm-cms.cern.ch \n',
                    'storage_path     = /srm/managerv2?SFN=/castor/cern.ch \n',
                    'user_remote_dir  = /user/x/xshi/afb/dat/srm/ \n',
                    '[CRAB] \n',
                    'scheduler        = glite \n',
                    'jobtype          = cmssw \n',
                    'use_server       = 1 \n' ]
        
    return cfg 


   
def get_file_size(f, lscmd='rfdir'):
    cmd = '%s %s' % (lscmd, f)
    output = proc_cmd(cmd)
    
    if 'No such file or directory' in output:
        size = 0
    else:
        if lscmd == 'rfdir':
            size = int(output.split()[4])
        elif lscmd == 'cmsLs':
            size = int(output.split()[1])
        else:
            raise NameError(lscmd)
        
    return size

    

def get_files_from_ls(cmd):
    files = []
    if 'cmsLs' in cmd:
        hostname = os.environ['HOSTNAME']
        if 'lxplus' not in hostname:
            sys.stdout.write('Please do this on lxplus. \n')
            sys.exit()
    output = proc_cmd(cmd)
    if 'No such file or directory' in output:
        sys.stdout.write(output)
        sys.stdout.write('Please check cmd: \n')
        sys.stdout.write(cmd)
        sys.stdout.write('\n')
        sys.exit()
        
    lines = output.split('\n')
    for l in lines:
        l = l.strip()
        if l == '':
            continue 
        try: 
            f = l.split()[-1].split('/')[-1]
            if 'cmsLs' in cmd:
                size = l.split()[1]
                if int(size) < 5:
                    sys.stdout.write('File size is too small: \n')
                    sys.stdout.write('%s\n' %f)
        except IndexError:
            sys.stdout.write('Cannot parse file name: \n')
            print "Canafdsfa"+l+"AAA"
            sys.exit()

        files.append(f)

    return files


def get_files_name_size_list_from_path(p, lscmd='ls -l'):
    cmd = '%s %s' % (lscmd, p)
    output = proc_cmd(cmd)

    files_name_size_list = []
    if 'No such file or directory' in output:
        sys.stdout.write(output)
        sys.stdout.write('Please check cmd: \n')
        sys.stdout.write(cmd)
        sys.stdout.write('\n')
        sys.exit()
        
    lines = output.split('\n')
    for l in lines:
        l = l.strip()
        if l == '' or len(l.split()) < 3: 
            continue 

        name = l.split()[-1]
        size = l.split()[-5]

        files_name_size_list.append((name, size))

    return files_name_size_list
    
 
# def get_func_name_from_label(label, ver=None):
#     if ver is not None: 
#         label = get_label_by_version(label, ver)

#     com_name = label.replace('/', '_')
#     com_name = com_name.replace('.', '_')
#     com_name = com_name.replace('-', '_')
#     return com_name

  
def get_jobs(stage, datatype, label, verbose=False):
    dbpath = os.path.join(atr.datpath, stage, 'db', label)
    dbname = 'jobs_%s.db' % datatype
    dbfile = check_and_join(dbpath, dbname)
    
    db = open(dbfile)
    jobs = json.load(db)
    db.close()
    if verbose: 
        sys.stdout.write('Loaded db file: %s \n' % dbfile)
    return jobs


def get_label_by_version(label, ver):
    if '/' in label: 
        labels = label.split('/')
       
    else:
        labels = [label]

    new_labels = []
    for label in labels:
        if ver == 'x':
            part = '.'.join(label.split('.')[:1])
        elif ver == 'x.':
            part = label.split('.')[0]
        elif ver == 'x.x':
            part = '.'.join(label.split('.')[:2])
        else:
            raise NameError(ver)
        new_labels.append(part)

    if len(new_labels) > 1: 
        new_label = '/'.join(new_labels)
    else:
        new_label = new_labels[0]

    return new_label



def get_number_from_label(label):
    if '_1M_' in label:
        number = 1000000
    elif '_10M_' in label:
        number = 10000000
    elif '_100M_' in label:
        number = 100000000
    elif '_5E5_' in label or '-5E5-' in label:
        number = 500000
    elif '_2E7_' in label or '-2E7-' in label:
        number = 20000000
    else:
        raise NameError(label)
    return number


def get_obj_db(dbfile):
    db = open(dbfile)
    obj = json.load(db)
    db.close()
    return obj    


def get_option(args, key):
    idx = args.index(key) + 1
    option = args[idx]
    return option 

    
# def get_options(args, key): 
#     opt = False 
#     if key == 'test':
#         if '-t' in args:
#             opt = True

#     elif key == 'pbar':
#         if '-p' in args:
#             opt = True

#     elif key == 'loglevel':
#         opt = 'INFO'
#         if '-log' in args:
#             value_idx = args.index('-log') + 1
#             opt = args[value_idx]

#     elif key == 'batch':
#         if '-b' in args:
#             opt = True

#     elif key == 'grid':
#         if '-grid' in args:
#             opt = True

#     elif key == 'job':
#         if '-job' in args:
#             value_idx = args.index('-job') + 1
#             opt = args[value_idx]

#     elif key == 'submit' or key == 'status' or \
#         key == 'stager_qry' or key == 'stager_get' and key in args:
#         value_idx = args.index(key) + 1
#         try:
#             val = args[value_idx]
#         except IndexError:
#             sys.stdout.write('Please specify the range or use "all". \n\n')
#             raise

#         opt = get_range_from_str(val)#, start, stop)

#     elif key == 'queue':
#         opt = '1nd'
#         if '-queue' in args:
#             qidx = args.index('-queue') + 1
#             opt = args[qidx]

#     elif key == 'jobname':
#         if '-j' in args:
#             idx = args.index('-j') + 1
#             opt = args[idx]

#     elif key == 'rng':
#         if '-rng' in args:
#             value_idx = args.index('-rng') + 1
#             try:
#                 val = args[value_idx]
#             except IndexError:
#                 sys.stdout.write('Please specify the range or use "all". \n\n')
#                 raise
#             opt = get_range_from_str(val)

#     elif key == 'rmdupl':
#         if '-rmdupl' in args:
#             opt = True

#     elif key == 'seq':
#         if '-seq' in args:
#             value_idx = args.index('-seq') + 1
#             try:
#                 opt = args[value_idx]
#             except IndexError:
#                 sys.stdout.write('Please specify the sequence or use "all". \n\n')
#                 raise

#     elif key == 'stat':
#         if '-stat' in args:
#             opt = True

#     else:
#         raise NameError(key)
    
#     return opt


# def get_options_name(args, key):
#     name = None
#     if key == 'rng':
#         if '-rng' in args:
#             value_idx = args.index('-rng') + 1
#             try:
#                 val = args[value_idx]
#             except IndexError:
#                 sys.stdout.write('Please specify the range or use "all". \n\n')
#                 raise

#             name = val
#             name = name.replace('-', '_')
            
#     else:
#         raise NameError(key)

#     return name 


#def get_progressbar(maxval=10000):
#    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=maxval).start()
#    return pbar


def get_range_from_str(val, start=0, stop=None):
    if val == 'all' :
        return val
        
    def get_range_hypen(val):
        start = int(val.split('-')[0])
        tmp_stop = val.split('-')[1]
        if tmp_stop != '':
            stop = int(val.split('-')[1])+1
        return range (start, stop)

    result = []

    if '-' in val and ',' not in val:
        result = get_range_hypen(val)
        
    elif ',' in val:
        items = val.split(',')
        for item in items:
            if '-' in item:
                result.extend(get_range_hypen(item))
            else:
                result.append(int(item))
    else:
        result.append(int(val))

    return result


def get_root_dbfile(ntppaths, datatype, label):
    dbpath = os.path.join(atr.datpath, 'db', label)
    dbname = datatype+'.db'
    dbfile = check_and_join(dbpath, dbname)

    db = open(dbfile, 'w')
    path_name_size_list = []
    for ntppath in ntppaths:
        cmd = 'rfdir %s' % ntppath
        output = proc_cmd(cmd)
        name_size_list = parse_ls(output)
        for name, size in name_size_list:
            filename = os.path.join(ntppath, name)
            path_name_size_list.append((filename, size))

    json.dump(path_name_size_list, db)
    db.close()
    
    check_update_status(dbfile)
    return dbfile 


def get_rootfiles_from_db(dbfile):
    db = open(dbfile)
    jobs = json.load(db)
    db.close()
    rootfiles = []
    for k, v in jobs.items():
        rootfile = v['rootfile']
        rfs = [r[0].encode('ascii') for r in rootfile]
        rootfiles.extend(rfs)

    return rootfiles


def get_rootfiles_sizes_from_paths(paths, lscmd):
    rootfiles_sizes = []
    for p in paths:
        cmd = '%s %s' % (lscmd, p)
        output = proc_cmd(cmd)

        if 'lcg-ls' in lscmd:
            f = LcgFile(data=output.split('\n'))
            name_size_list = f.name_size_list

        elif 'rfdir' in lscmd:
            name_size_list = parse_ls(output, ignore='.test')

        else:
            raise NameError(lscmd)

        for name, size in name_size_list:
            filename = os.path.join(p, name)
            rootfiles_sizes.append((filename, size))

    return rootfiles_sizes


def get_events_set(txtfile):
    f = UserFile(txtfile)
    run_evts = []
    n = 0
    for l in f.data:
        n += 1
        run_evts.append(tuple(l.split()))
    evts_set = set(run_evts)
    duplicate = len(run_evts) - len(evts_set)
    if duplicate > 0:
        sys.stdout.write('    There are %s duplicated events in %s. \n'
                         % (duplicate, txtfile))
    return set(run_evts)


def get_file_groups_by_size(files_name_size_list, size_max='2G', nfile_max=500):
    if size_max == '2G':
        size_max = 2*1073741824
        
    groups = []

    group = []    
    size_sum = 0

    for name_size in files_name_size_list:
        name = name_size[0]
        size = name_size[1]
        
        name = name.encode('ascii')
        if size_sum < size_max and len(group) < nfile_max:
            group.append(name)
            #group.append(name_size)
            size_sum += float(size)
        else:
            groups.append(group)
            group = []
            size_sum = 0
            group.append(name)
            #group.append(name_size)            
            size_sum += float(size)

        if name_size == files_name_size_list[-1]:
            groups.append(group)

    return groups

def get_selfile(datatype, label, job_opt=False):
    if datatype == 'dt' and label in ['5ifbv2.3.1', '5ifbv2.3.2']:
        jobs = range(87)
    elif datatype == 'dt' and label in ['5ifbv2.6.1']:
        jobs = range(80)
    elif datatype == 'data' and label in [
            'Run2012B-PromptReco-v1.5.4.1',
            'Run2012B-PromptReco-v1.5.5.1',
            ]:
        fpath = os.path.join(atr.datpath, 'ntp/data/Run2012B-PromptReco-v1_5')
        if label == 'Run2012B-PromptReco-v1.5.4.1':
            fname = 'SingleCandB0KstMuMu_v4.root'
        elif label == 'Run2012B-PromptReco-v1.5.5.1':
            fname = 'SingleCandB0KstMuMu_v5.root'
        else:
            raise NameError(label)
        
        ffile = os.path.join(fpath, fname)
        return ffile

    elif datatype == 'data' and label in [
            'Run2011-Merged-v1.1.1.1',
            ]:
        fpath = os.path.join(atr.datpath, 'ntp/data/Run2011-Merged-v1_1')
        fname = 'SingleCandB0KstMuMu_v1.root'
        
        ffile = os.path.join(fpath, fname)
        return ffile
    else:
        raise NameError(datatype, label)

    if job_opt:
        filepath = atr.selpath(label)
    else:
        filepath = os.path.join(atr.datpath, 'sel')

    f = []
    for job in jobs:
        name = 'tree_%s.root' %(job+1)
        f.append(os.path.join(filepath, datatype, label, name))
        
    return f


def get_skmfile(datatype, label):
    if datatype == 'dt':
        name = '*.root'
    else:
        raise NameError(datatype)
    filepath = atr.skmpath(label)
    f = os.path.join(filepath, datatype, label, name)
    return f

def get_y_title(h, unit=''):
    number = ( h.GetXaxis().GetXmax() - h.GetXaxis().GetXmin() )/h.GetXaxis().GetNbins()
    if unit != '':
        title = 'Events / %s (%s) ' % (number, unit)
    else:
        title = 'Events / %s ' % number

    return title


def group_files_by_size(path_name_size_list, size_max='2G'):
    if size_max == '2G':
        size_max = 2*1073741824
        
    groups = []

    group = []    
    size_sum = 0

    for name_size in path_name_size_list:
        name = name_size[0]
        size = name_size[1]
        
        #name = name.encode('ascii')
        if size_sum < size_max:
            #group.append(name)
            group.append(name_size)
            size_sum += float(size)
        else:
            groups.append(group)
            group = []
            size_sum = 0
            #group.append(name)
            group.append(name_size)            
            size_sum += float(size)

        if name_size == path_name_size_list[-1]:
            groups.append(group)

    return groups

def hists_draw(hist_list, legend=None, xtitle=''):
    hist_list.sort()
    hist_list.reverse()

    if legend is None:
        from ROOT import TLegend
        legend = ROOT.TLegend()
    
    i = 0
    for hist, val in  hist_list:
        i += 1
        if i == 1:
            hist.SetLineWidth(3)
            hist.SetLineColor(i)
            hist.GetXaxis().SetTitle(xtitle)
            hist.Draw()
        else:
            hist.SetLineWidth(3)
            hist.SetLineColor(i)
            hist.Draw("same")

        legend.AddEntry(hist, hist.GetName())

    legend.SetFillColor(0)
    legend.SetBorderSize(1)
    legend.Draw()

def init_logging(args):
    #loglevel = get_options(args, 'loglevel')
    loglevel = 'INFO'
    if option_exists(args, '-loglevel'):
        get_option(args, '-loglevel')
        
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=numeric_level)
    return logging


def label_to_name(label):
    name = label.replace('/', '_')
    name = name.replace('.', '_')
    return name 


def list_to_file(li, f, listname='l', prefix=''):
    fo = UserFile()
    fo.append('%s = [\n' %listname)
    num_files = len(li)
    for i in range(num_files):
        rf = li[i]
        if i != num_files-1:
            fo.append('\'%s%s\', \n' % (prefix, rf))
        else:
            fo.append('\'%s%s\'] \n' % (prefix, rf))
    fo.output(f)


def merge_root_files(srcdir, dstdir, job=None, size_max=2*9.8e8, nfile_max=500):

    dbname = '.files.db'
    dbfile = check_and_join(dstdir, dbname)
    

    files_name_size_list = get_files_name_size_list_from_path(srcdir)
    set_obj_db(files_name_size_list, dbfile)

    # if force_clean_db:
    #     os.remove(dbfile)
        
    # if os.access(dbfile, os.F_OK):
    #     sys.stdout.write('using existing db file...\n')
    #     files_name_size_list = get_obj_db(dbfile)
    # else:
    #     sys.stdout.write('creating db file...\n')
    #     files_name_size_list = get_files_name_size_list_from_path(srcdir)
    #     set_obj_db(files_name_size_list, dbfile)

    size = [int(s) for n, s in files_name_size_list]
    totalsize = sum(size)
    sys.stdout.write('total of %s files, size %s.\n' % (len(files_name_size_list),
                                                      convert_bytes(totalsize)))

    com_name = files_name_size_list[0][0].split('_')[0] 
    rootfile_groups = get_file_groups_by_size(files_name_size_list, 
                                              size_max=size_max, nfile_max=nfile_max)

    ngroups = len(rootfile_groups)
    dst_files = get_files_name_size_list_from_path(dstdir)
    last_target_number = len(dst_files)

    if job is None:
        jobs = range(1, ngroups+1)
    else: 
        jobs = get_range_from_str(job)

    for i in jobs: 
        sourcefiles = ' '.join(rootfile_groups[i-1])
        targetname = '%s_merged_%s.root' % (com_name, i+last_target_number)
        targetfile = os.path.join(dstdir, targetname)
        cmd = 'hadd %s %s' %(targetfile, sourcefiles)
        sys.stdout.write('merging %s of %s ...\n' %(i, ngroups))
        sys.stdout.flush()
        output = proc_cmd(cmd, procdir=srcdir)
        print output

    sys.stdout.write('Merged files: \n')
    cmd = 'ls -lh %s' % dstdir
    output = proc_cmd(cmd)
    print output 


def numbers_to_string(numbers):
    strnums = map(str, numbers)
    string = ','.join(strnums)
    return string


def option_exists(args, key): 
    exists = False 
    if key in args:
        exists = True
    return exists


def output_set(s, fname):
    f = UserFile()
    for l in s:
        f.append('%s %s \n' %l)
    f.output(fname, verbose=1)
    
    
def parse_ls(content, ignore=None, poz=None):
    name_size_list = []
    content = content.strip()
    lines = content.split('\n')
    for l in lines:
        if ignore is not None and ignore in l:
                continue
        if poz is not None: 
            size = l.split()[poz[0]]
            name = l.split()[poz[1]]

        else:
            try:
                size = l.split()[4]
                name = l.split()[-1]
            except IndexError:
                print l
                raise
            
        name_size_list.append((name, size))

    return name_size_list

    
def proc_cmd(cmd, test=False, verbose=1, procdir=None, shell=False):
    if test:
        sys.stdout.write(cmd+'\n')
        return 

    # No need to add quote in the cmd! 
    cwd = os.getcwd()
    if procdir is not None:
        os.chdir(procdir)
    
    args = cmd.split()

    if isinstance(cmd, list) or shell:
        args = cmd 

    process = subprocess.Popen(args, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=shell)
    stdout = process.communicate()[0]

    if 'error' in stdout:
        sys.stdout.write(stdout)
        
    if procdir is not None:
        os.chdir(cwd)

    return stdout

def root_chain(rootfile, obj, debug=False):
    from ROOT import TChain
    chain = TChain(obj)
    num = 0
    if isinstance(rootfile, str):
        num = chain.Add(rootfile)
    if isinstance(rootfile, list):
        for rf in rootfile:
            rf = str(rf)
            num += chain.Add(rf)

    if debug:
        sys.stdout.write('\ndebug: %s root file(s) added:\n%s\n'
                         % (num, rootfile))
    if num < 1:
        raise NameError('No rootfile added!', rootfile)
    
    return chain


    
def save_fit_results(pars, txtfile, err_type='SYMM', verbose=0):
    f = UserFile()
    if err_type == 'SYMM':
        f.append('Name\t|| Value\t|| Error\n')
        for par in pars:
            line = '%s\t| %s\t| %s\n' %(
                par.GetName(), par.getVal(), par.getError())
            f.append(line)

    elif err_type == 'ASYM':
        f.append('Name\t|| Value\t|| Error\t|| Low\t || High\n')
        for par in pars:
            line = '%s\t| %s\t| %s\t| %s\t| %s\n' %(
                par.GetName(), par.getVal(), par.getError(),
                par.getAsymErrorLo(), par.getAsymErrorHi())
            f.append(line)

    else:
        raise NameError(err_type)
    
    f.output(txtfile, verbose=verbose)


def scram_build(cmssw):
    procdir = os.path.join(atr.afbpath, 'rel', cmssw, 'src')
    cmd = 'scram b'
    output = proc_cmd(cmd, procdir=procdir)
    

def set_file(path, label, name, suffix, mode='', test=False):
    path_ = os.path.join(path, label)

    if test:
        name += '_test'
    name_ = name + suffix

    file_ = check_and_join(path_, name_, mode=mode)
    return file_


def set_figname(args):
    args_ = list(args)
    if '-b' in args_:
        args_.remove('-b')
    if '-t' in args_:
        args_.remove('-t')
    if '-f' in args_:
        args_.remove('-f')

    name = '_'.join(args_)
    name = name.replace('/', '_')
    name = name.replace('.', '_')
    return name 


def set_obj_db(obj, dbfile):
    db = open(dbfile, 'w')
    json.dump(obj, db)
    db.close()
    check_update_status(dbfile)
    

# def set_logfile(path_, name, job=None, test=False):
#     if test:
#         name += '_test'

#     if job is not None:
#         name += '_%s' % job

#     name_ = name + '.txt'
#     file_ = check_and_join(path_, name_)
#     return file_


def set_root_style(stat=1110, grid=0, PadTopMargin=0.08, PadLeftMargin=0.10,
                   TitleOffsetY=1.2):
    # must be used in the beginning
    from ROOT import gROOT, gStyle
    gROOT.SetBatch(1)
    gROOT.Reset()

    gStyle.SetOptTitle(1) 
    gStyle.SetTitleFillColor(0) 
    gStyle.SetTitleBorderSize(0)
    
    # gStyle.SetTitleFont(22)
    # gStyle.SetStatFont(22)
    # gStyle.SetStatColor(1)
    # gStyle.SetLabelFont(22,"X")
    # gStyle.SetLabelFont(22,"Y")
 
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetCanvasDefX(0)
    gStyle.SetCanvasDefY(0)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameBorderSize(1)
    gStyle.SetFrameFillColor(0)
    gStyle.SetFrameFillStyle(0)
    gStyle.SetFrameLineColor(1)
    gStyle.SetFrameLineStyle(1)
    gStyle.SetFrameLineWidth(1)

    gStyle.SetPadTopMargin(PadTopMargin) 
    gStyle.SetPadLeftMargin(PadLeftMargin) 
    gStyle.SetPadRightMargin(0.05) 

    gStyle.SetLabelSize(0.03, "XYZ") 
    gStyle.SetTitleSize(0.04, "XYZ") 
    gStyle.SetTitleOffset(TitleOffsetY, "Y") 

    gStyle.SetPadBorderMode(0) 
    gStyle.SetPadColor(0) 
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadGridX(grid)
    gStyle.SetPadGridY(grid)

    gStyle.SetOptStat(stat)
    gStyle.SetStatColor(0)
    gStyle.SetStatBorderSize(1)

    
 
def set_selfile(datatype, label, path, name, job=None, test=False):
    print 'Please use the atr.sel.setrootfile!'
    raise NameError()

    # if job is not None:
    #     name += '_%s' % job
    # if test:
    #     name += '_test'
    # name_ = name + '.root'

    # if label in ['5ifbv2.3']:
    #     path_ = os.path.join(path, 'sel', datatype, label)
    # else:
    #     path_ = os.path.join(path, datatype, label)
        
    # file_ = check_and_join(path_, name_)
    # if '/castor' in file_:
    #     file_ = 'rfio:%s' %file_

    # return file_


def set_skmfile(datatype, label, name, job=None, test=False):
    if job is not None:
        name += '_%s' % job
    if test:
        name += '_test'
        
    filename = name + '.root'
    skmpath = atr.skmpath(label)
    filepath = os.path.join(skmpath, datatype, label)
    skimfile = check_and_join(filepath, filename)
    if '/castor' in skimfile:
        skimfile = 'rfio:%s' %skimfile
    return skimfile


def stager_cmd(cmdname, rootfile, test=False, verbose=False):
    cmd = '%s -M %s' % (cmdname, rootfile)
    stdout = proc_cmd(cmd, test)

    if stdout is not None and verbose:
        sys.stdout.write(stdout)

    return stdout 


def update_bashfile_cmd(bashfile, cmd, pre='setafb 0', test=False):
    bash_content =  '''#!/bin/sh
date
hostname

.  ~/.bashrc
%s

%s

date
''' % (pre, cmd)

    f = UserFile()
    f.append(bash_content) 
    f.output(bashfile)
    os.chmod(bashfile, 0755)
    return bashfile

   
def update_jobs_db(stage, datatype, label, test=False):
    if stage == 'srm':
        paths = atr.srmpath(label)
        lscmd = 'lcg-ls -l'
    elif stage == 'ntp':
        paths = atr.ntppath(label, datatype)
        lscmd = 'rfdir'
    elif stage == 'skm':
        skmpath = atr.skmpath(label)
        paths = [os.path.join(skmpath, datatype, label)]
        lscmd = 'rfdir'
    else:
        raise NameError(stage)

    dbpath = os.path.join(atr.datpath, stage, 'db', label)
    dbname = 'jobs_%s.db' % datatype
    if test:
        dbname += '.test'
    dbfile = check_and_join(dbpath, dbname)

    old_rootfiles = get_rootfiles_from_db(dbfile)
    rootfiles_sizes  = get_rootfiles_sizes_from_paths(paths, lscmd) 
    path_name_size_list = []
    for rf, sz in rootfiles_sizes:
        if rf in old_rootfiles:
            continue
        path_name_size_list.append((rf, sz))

    rootfiles = group_files_by_size(path_name_size_list)

    db = open(dbfile)
    jobs = json.load(db)
    db.close()

    num_old_jobs = len(jobs)
    job = num_old_jobs 
    for rootfile in rootfiles:
        job += 1
        jobs[job] = {}
        jobs[job]['rootfile'] = rootfile
        jobs[job]['status'] = 'created'

    db = open(dbfile, 'w')    
    json.dump(jobs, db)
    db.close()

    check_update_status(dbfile)
    sys.stdout.write('updated %s jobs.\n' % (len(jobs) - num_old_jobs))
    return dbfile 

    
