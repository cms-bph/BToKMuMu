"""
Providing Tools for file handling.

"""

import os
import sys
import filecmp
import shutil
from ConfigParser import RawConfigParser


__author__ = "Xin Shi <Xin.Shi@cern.ch>"
__copyright__ = "Copyright (c) Xin Shi"


# ------------------------------------------------------------
# Functions 
# ------------------------------------------------------------
def make_tmpfile(f):
    path, name = os.path.split(f)
    tmpname = '.tmp_' + name
    tmpfile = os.path.join(path, tmpname)
    return tmpfile


def check_update_status(f, verbose=1):
    tmpfile = make_tmpfile(f)
    if not os.access(tmpfile, os.F_OK):
        message = 'created %s ...\n' %f
        shutil.copy2(f, tmpfile)
    elif filecmp.cmp(f, tmpfile, shallow=False):
        message = 'up-to-date: %s\n' % f
    else:
        message = 'updated %s ...\n' %f
        #os.remove(tmpfile)

    if verbose > 0 :
        sys.stdout.write(message)

    return message


def check_and_copy(source_file, dest_file, verbose=1):
    message = ''
    local_source_file = None
    if '@' in source_file and ':' in source_file:
        local_source_file = make_tmpfile(dest_file)
        try:
            p = subprocess.Popen(['scp', source_file, local_source_file])
            os.waitpid(p.pid, 0)
            source_file = local_source_file
        except IOError:
            sys.stdout.write('Skipping %s \n' %source_file)
            return 

    if os.access(dest_file, os.F_OK) :
        if filecmp.cmp(source_file, dest_file, shallow=False):
            message =  'up-to-date: %s' % dest_file
        else:
            message = 'Updating %s ...' %dest_file
    else:
        message = 'Writing %s ...' %dest_file
    try:
        shutil.copy2(source_file, dest_file)
    except IOError:
        sys.stdout.write('Skipping %s \n' %source_file)
        return 
            
    if local_source_file != None:
        os.remove(local_source_file)

    if verbose > 0:
        sys.stdout.write(message+'\n')
        
    return message


def check_and_join(filepath, filename, mode=''):
    if '/castor' in filepath:
        cmd = 'rfdir %s' % filepath
        stdout = proc_cmd(cmd)
        if 'No such file or directory' in stdout:
            sys.stdout.write('creating dir (chmod 775) %s ...' % filepath)
            cmd = 'rfmkdir -p %s' % filepath
            proc_cmd(cmd)
            cmd = 'rfchmod 775 %s' % filepath
            proc_cmd(cmd)
            sys.stdout.write(' OK.\n')
    elif 'root://' in filepath:
        sys.stdout.write('Skip checking root:// ...\n')
        
    else: 
        if not os.access(filepath, os.F_OK):
            sys.stdout.write('creating dir %s ...' % filepath)
            os.makedirs(filepath)
            sys.stdout.write(' OK.\n')
        
    file_ = os.path.join(filepath, filename)
    if os.access(file_, os.F_OK) :
        tmpfile = make_tmpfile(file_)
        shutil.copy2(file_, tmpfile)
        if mode == 'w':
            os.remove(file_)

    return file_

def check_file(f):
    print f
    sys.exit()



# ------------------------------------------------------------
# Classes 
# ------------------------------------------------------------


class UserFile(object):
    '''Class to handle file  '''
    def __init__(self, filename=None, data=None):
        self.data = []
        if data != None:
            self.data = data
            
        if filename:
            self.input(filename)
            self.file = filename

    def append(self, content):
        self.data.append(content)

    def backup(self, prefix, verbose=0):
        path, name =  os.path.split(self.file)
        backup_file  = os.path.join(path, prefix + name)

        if 'default' in prefix and os.path.exists(backup_file):
            if verbose > 0: 
                sys.stdout.write('\nDefault file exits! : %s\n' % backup_file)
            return

        self.output(backup_file , verbose)

    def extend(self, content):
        self.data.extend(content)

    def find(self, pat):
        for line in self.data:
            if pat in line:
                return True
        return False

    def input(self, filename, verbose=0):
        fi = open(filename, 'r')
        for line in fi:
            self.data.append(line)
        fi.close()

    def input_data(self, data):
        self.data = data 

    def insert(self, index, newline):
        line_num = 0
        for line in self.data:
            line_num = line_num + 1
            if index in line:
                self.data.insert(line_num-1, newline)
                return 
        
    def output(self, f=None, verbose=1):
        fo = sys.stdout
        if f != None:
            filepath, filename = os.path.split(f)
            f = check_and_join(filepath, filename) 
            fo = open(f ,'w')
            
        for line in self.data:
            fo.write(line)
        fo.close()
        
        if f != None:
            message = check_update_status(f, verbose=verbose)
            return message


    def replace(self, old, new):
        line_num = 0
        for line in self.data:
            line_num = line_num + 1
            if old in line:
                self.data[line_num-1] = line.replace(old, new)

    def restore(self, prefix, verbose=0):
        path, name  =  os.path.split(self.file)
        backup_file =  os.path.join(path, prefix + name)
        f = UserFile(backup_file)
        self.data = f.data
            

    def set_key_value(self, key, value):
        for line in self.data:
            if key in line:
                old = line
                new = '%s = %s \n' % (key, value)
                self.replace(old, new)

 
class LogFile(UserFile):
    "Handle log file"

    def __init__(self, filename=None):
        self.processed = 'N/A'
        self.skimmed = 'N/A'
        self.selected = 'N/A'
        self.duration = 'N/A'
        
        try:
            UserFile.__init__(self, filename)
            self.parse()
        except IOError:
            pass
        
    def parse(self):
        "parse log file"
        line_no = -1
        found_stream_event = False
        start_lumi_info = False
        for line in self.data:
            line_no += 1
            line = line.strip()
            if 'processed' in line:
                self.processed = line.split(' ')[1]
                
            if 'skimmed' in line:
                self.skimmed = line.split(' ')[1]

            if 'selected' in line:
                self.selected = line.split(' ')[1]

            if 'duration' in line:
                self.duration = line.replace('duration ', '')


class LcgFile(UserFile):
    "Handle Lcg output file"

    def __init__(self, filename=None, data=None):
        UserFile.__init__(self, filename, data)
        self.parse()
        
    def parse(self):
        "parse Lcg file"
        line_no = -1
        self.name_size_list = [] 
        for line in self.data:
            line_no += 1
            if '/pnfs' in line:
                line = line.split() 
                size = line[4]
                name = line[6]
                self.name_size_list.append((name, size))


class CfgFile(UserFile, RawConfigParser):
    "handle Cfg file"

    def __init__(self, filename=None):
        super(CfgFile,self).__init__(filename)
        RawConfigParser.__init__(self)


    def output(self, f=None): 
        for section in self._sections:
            self.append("[%s]\n" % section)
            for (key, value) in self._sections[section].items():
                if key == "__name__":
                    continue
                if (value is not None) or (self._optcre == self.OPTCRE):
                    key = " = ".join((key, str(value).replace('\n', '\n\t')))
                self.append("%s\n" % (key))

        super(CfgFile,self).output(f)
