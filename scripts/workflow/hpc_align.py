import subprocess
import os
import sys
import datetime
from configparser import ConfigParser

def check_exist(cmd, thing):
    try:
        subprocess.check_output('%s %s' % (cmd, thing), shell=True)
    except subprocess.CalledProcessError:
        print("Error: did not find %s in path." % thing)
        sys.exit(0)

def log_error(cmd, exec_output, exec_error):
        with open(LOG_FILE, 'a') as f:
                f.write('time: %s\ncmd: %s\noutput: %s\nexec error:%s\n' % (str(datetime.datetime.now()), cmd, exec_output, exec_error))
                
def log_final(no_error, argv):
    log_output = os.path.join(SCRIPT_DIR, 'log_align_analyze_sort.txt')
    with open(log_output, 'a') as f:
        f.write('%s %s %s %s\n' % (no_error, argv[0], argv[1], str(datetime.datetime.now())))

if len(sys.argv) != 3:
    print('Usage: python', sys.argv[0], 'config_file.txt','read_file.txt')
    sys.exit(0)

if len(sys.argv) != 3:
    print('Usage: python', sys.argv[0], 'config_file.txt','read_file.txt')
    sys.exit(0)

#--------------------------------------------------------------
# read config file
#--------------------------------------------------------------
config = ConfigParser()
config.readfp(open(sys.argv[1]))

READS_DIR = config.get('config', 'READS_DIR')
LOG_FILE = config.get('config', 'LOG_FILE')

ref = config.get('config', 'REF')
OUTPUT_DIR = config.get('config', 'OUTPUT_DIR')

#--------------------------------------------------------------

SCRIPT_DIR = os.getcwd()
read_file = open(sys.argv[2])

check_exist('which', 'bwa')
check_exist('which', 'samtools')
check_exist('ls', ref)


if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

for line in read_file:
    read1 = os.path.join(READS_DIR, line.strip() + '_1.fastq')
    read2 = os.path.join(READS_DIR, line.strip() + '_2.fastq')
    check_exist('ls', read1)
    check_exist('ls', read2)
    
    name = read1.split('/')[-1].split('_R1')[0]
    out_sam = os.path.join(OUTPUT_DIR, name+'.sam')
    out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_f2_q20.sam')
    

    output = 'None'

    # 01_alignment      
    if os.path.exists(out_sam):
        print('Alignment might have been done already.  Skip bwa.')
    else:
        cmd = 'bwa mem %s %s %s' % (ref,read1,read2)
        try:
            output = subprocess.check_call(cmd, shell=True, stdout=open(out_sam, 'w'))
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    # 02_filter_by_samtools
    if os.path.exists(out_filtered_sam):
        print('Alignment might have been filtered already.  Skip samtools.')
    else:
        print("Filter bwa's output")
        cmd = 'samtools view -f 2 -q 20 %s' % out_sam
        try:
            ouptut = subprocess.check_call(cmd, shell=True, stdout=open(out_filtered_sam, 'w'))
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    
    print ("Finished %s. " %(line))
