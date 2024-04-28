import subprocess
import os
import sys
import datetime
from configparser import ConfigParser
import time
import threading

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

# maxthreads = 4
# sema = threading.Semaphore(value=maxthreads)
def process(*args, **kwargs):
    # sema.acquire()
    config_file = kwargs['config_file']
    read_ID = kwargs['read_ID']
    print("\nConfig file:", config_file)
    print("Processing read", read_ID)
    #--------------------------------------------------------------
    # read config file
    #--------------------------------------------------------------
    SCRIPT_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))+'/'
    config = ConfigParser()
    # config.readfp(open(SCRIPT_DIR+'defaults.ini'))
    config.read_file(open(SCRIPT_DIR+'defaults.ini'))
    default_alignment_quality = config.get('defaults', 'alignment_quality')
    default_chloroplast = config.get('defaults', 'chloroplast')
    default_mitochondria = config.get('defaults', 'mitochondria')

    # config.readfp(open(sys.argv[1]))
    # config.readfp(open(config_file))
    config.read_file(open(config_file))

    READS_DIR = config.get('config', 'READS_DIR')
    LOG_FILE = config.get('config', 'LOG_FILE')

    ref = config.get('config', 'REF')
    OUTPUT_DIR = config.get('config', 'OUTPUT_DIR')
    is_pair_read = int(config.get('config', 'PE'))

    # quality for SAM filter
    try:
        alignment_quality = config.get('config', 'alignment_quality')
    except:
        alignment_quality = default_alignment_quality

    try: 
        chloroplast = config.get('config', 'chloroplast')
    except:
        chloroplast = default_chloroplast

    try:
        mitochondria = config.get('config', 'mitochondria')
    except:
        mitochondria = default_mitochondria

    #--------------------------------------------------------------

    # SCRIPT_DIR = os.getcwd()
    # read_file = open(sys.argv[2])
    # read_file = open(input_read_file)

    check_exist('which', 'bwa-meme')
    check_exist('which', 'samtools')
    check_exist('ls', ref)
    filter_cp_mt = os.path.join(SCRIPT_DIR, 'filter_samfiles_cp_mt.py')
    check_exist('ls', filter_cp_mt)

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    SAM_DIR = os.path.dirname(os.path.realpath(ref))+'/'
    start_time = time.time()
    # for line in read_file:
    if is_pair_read == 1:
        read1 = os.path.join(READS_DIR, read_ID + '_R1.fastq')
        read2 = os.path.join(READS_DIR, read_ID + '_R2.fastq')
        check_exist('ls', read1)
        check_exist('ls', read2)
        name = read1.split('/')[-1].split('_R1')[0]
    else:
        read = os.path.join(READS_DIR, read_ID + '.fastq')
        check_exist('ls', read)
        name = read.split('/')[-1].split('.')[0]

    # out_sam = os.path.join(OUTPUT_DIR, name+'.sam')
    out_sam = os.path.join(SAM_DIR, name+'.sam')
    # out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_f2_q'+alignment_quality+'.sam')
    out_filtered_sam = os.path.join(SAM_DIR, name+'_F0x900_F0x04_q'+alignment_quality+'.sam')
    output = 'None'

    # if '_MT' in name:
        # is_pair_read = 1

    if is_pair_read == 1:
        bwacmd = 'bwa-meme_mode2 mem -7 -t 2 %s %s %s' % (ref,read1,read2)
    else:
        bwacmd = 'bwa-meme_mode2 mem -7 -t 2 %s %s' % (ref,read)
        # out_fastq = os.path.join(SAM_DIR, name+'.fastq')
        # if os.path.exists(out_fastq):
            # print('cat fastq might have been done already.  Skip cat.')
        # else:
            # cmd = 'cat %s %s' % (read1,read2)
            # try:
                # output = subprocess.check_call(cmd, shell=True, stdout=open(out_fastq, 'w'))
                # # output.wait()
            # except:
                # no_error = False
                # log_error(cmd, output, sys.exc_info())
        # bwacmd = 'bwa-meme mem -7 %s %s' % (ref,out_fastq)

    # 01_alignment      
    if os.path.exists(out_sam):
        print('Alignment might have been done already.  Skip bwa-meme.')
    else:
        # cmd = 'bwa-meme mem %s %s %s' % (ref,read1,read2)
        cmd = bwacmd
        try:
            output = subprocess.check_call(cmd, shell=True, stdout=open(out_sam, 'w'))
            # output.wait()
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    alignment_time = time.time()
    # print("Alignment time for ", line.strip(), ": ", alignment_time-start_time)
    print("Alignment time for ", read_ID, ": ", alignment_time-start_time)

    # 02_filter_by_samtools
    if os.path.exists(out_filtered_sam):
        print('Alignment might have been filtered already.  Skip samtools.')
    else:
        print("Filter bwa-meme's output")
        # cmd = 'samtools view -f 2 -q %s %s' % (alignment_quality , out_sam)
        # cmd = 'samtools view -f 2 -F 0x900 -q %s %s' % (alignment_quality , out_sam)
        cmd = 'samtools view -h -F 0x904 -q %s %s' % (alignment_quality , out_sam)
        try:
            output = subprocess.check_call(cmd, shell=True, stdout=open(out_filtered_sam, 'w'))
            # output.wait()
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    # select reads that mapped to chloroplast and mitochondria
    print('Filter alignments for chloroplast and mitochondrial genomes.')
    # cmd = 'python filter_samfiles_cp_mt.py %s %s %s %s' %(out_filtered_sam, OUTPUT_DIR, chloroplast, mitochondria)
    cmd = 'python %s %s %s %s %s' %(filter_cp_mt, out_filtered_sam, OUTPUT_DIR, chloroplast, mitochondria)
    try:
        output = subprocess.check_call(cmd, shell=True)
        # output.wait()
    except:
        no_error = False
        log_error(cmd, output, sys.exc_info())            

    filter_time = time.time()
    # print("Filter time for ", line.strip(), ": ", filter_time-alignment_time)
    print("Filter time for ", read_ID, ": ", filter_time-alignment_time)

    print ("Finished %s. " %(line))
    # sema.release()


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python', sys.argv[0], 'config_file.txt','read_ID')
        sys.exit(0)
    
    kw = {'config_file': sys.argv[1], 'read_ID': sys.argv[2]}
    main(kw)
