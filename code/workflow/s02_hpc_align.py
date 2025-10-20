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
    default_rm_alter_align = config.get('defaults', 'rm_alter_align')
    # config.readfp(open(sys.argv[1]))
    # config.readfp(open(config_file))
    config.read_file(open(config_file))

    READS_DIR = config.get('config', 'READS_DIR')
    LOG_FILE = config.get('config', 'LOG_FILE')

    ref = config.get('config', 'REF')
    OUTPUT_DIR = config.get('config', 'OUTPUT_DIR')
    is_pair_read = int(config.get('config', 'PE'))

    try:
        rm_alter_align = int(config.get('config', 'rm_alter_align'))
    except:
        rm_alter_align = int(default_rm_alter_align)

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

    if rm_alter_align == 1:
        filter_cp_mt = os.path.join(SCRIPT_DIR, 'filter_samfiles_rm_numt.py')
    else:
        filter_cp_mt = os.path.join(SCRIPT_DIR, 'filter_samfiles_cp_mt.py')
    check_exist('ls', filter_cp_mt)

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    SAM_DIR = os.path.dirname(os.path.realpath(ref))+'/'

    start_time = time.time()

    # for line in read_file:
    if is_pair_read == 1:
        read1 = os.path.join(READS_DIR, read_ID + '_1.fastq.gz')
        read2 = os.path.join(READS_DIR, read_ID + '_2.fastq.gz')
        # check_exist('ls', read1)
        # check_exist('ls', read2)
        name = read1.split('/')[-1].split('_1')[0]
        sam_flag = 'f2_F0x900'
    else:
        read = os.path.join(READS_DIR, read_ID + '.fastq.gz')
        # check_exist('ls', read)
        name = read.split('/')[-1].split('.')[0]
        sam_flag = 'F0x904'

    out_sam = os.path.join(OUTPUT_DIR, name+'.sam')
    out_filtered_sam = os.path.join(OUTPUT_DIR, name+'_'+sam_flag+'_q'+alignment_quality+'.sam')
    out_bam = os.path.join(OUTPUT_DIR, name+'.bam')
    out_markdup = os.path.join(OUTPUT_DIR, name+'.markdup')
    out_sorted_bam = os.path.join(OUTPUT_DIR, name+'_sorted.bam')
    out_sorted_bai = os.path.join(OUTPUT_DIR, name+'_sorted.bam.bai')
    if mitochondria != 'None':
        mt_out = os.path.join(OUTPUT_DIR,'mitochondria')
        out_csv_filtered_sam = os.path.join(mt_out, name+'_'+sam_flag+'_q'+alignment_quality+'.sam')
        out_threshold_file = os.path.join(mt_out, name+'_threshold.txt')
    if chloroplast != 'None':
        cp_out = os.path.join(OUTPUT_DIR,'chloroplast')
        out_csv_filtered_sam = os.path.join(cp_out, name+'_'+sam_flag+'_q'+alignment_quality+'.sam')
        out_threshold_file = os.path.join(cp_out, name+'_threshold.txt')
    # out_sam = os.path.join(SAM_DIR, name+'.sam')
    # out_filtered_sam = os.path.join(SAM_DIR, name+'_f2_F0x900_q'+alignment_quality+'.sam')
    # out_bam = os.path.join(SAM_DIR, name+'.bam')
    # out_markdup = os.path.join(SAM_DIR, name+'.markdup')
    # out_sorted_bam = os.path.join(SAM_DIR, name+'_sorted.bam')
    # out_threshold_file = os.path.join(SAM_DIR, name+'_threshold.txt')
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
    if os.path.exists(out_csv_filtered_sam) or os.path.exists(out_filtered_sam):
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

        # Check if the output file exists, and if it does, skip the operation
        for out_file, cmd in [(out_bam, 'samtools view -@ 2 -hb %s > %s' % (out_sam, out_bam)),
                              (out_markdup, 'sambamba markdup -t 2 -r -p %s %s' % (out_bam, out_markdup)),
                              (out_sorted_bam, 'sambamba sort -t 2 -p -o %s %s' % (out_sorted_bam, out_markdup))]:
            if os.path.exists(out_file):
                print(f'Skip {out_file}.')
            else:
                try:
                    subprocess.check_call(cmd, shell=True)
                except subprocess.CalledProcessError as e:
                    no_error = False
                    log_error(cmd, e.output, e)

    alignment_time = time.time()
    # print("Alignment time for ", line.strip(), ": ", alignment_time-start_time)
    print("Alignment time for ", read_ID, ": ", alignment_time-start_time)

    # 02_cal_cov_by_samtools
    if mitochondria != 'None':
        if os.path.exists(out_threshold_file):
            print('Threshold might have been calculated already.  Skip Calculate.')
        else:
            print("Calculating percentage threshold")
            if os.path.exists(out_sorted_bam):
                # Calculate the average coverage for the nuclear genome and mitogenome.
                avg_cov_genome = subprocess.check_output(f'samtools depth -a {out_sorted_bam} | grep -v "^{mitochondria}" | awk \'BEGIN {{count=0; sum=0}} {{sum+=$3; count++}} END {{if (count > 0) print sum/count; else print 0}}\'', shell=True)
                avg_cov_mt = subprocess.check_output(f'samtools depth -a {out_sorted_bam} | grep "^{mitochondria}" | awk \'BEGIN {{count=0; sum=0}} {{sum+=$3; count++}} END {{if (count > 0) print sum/count; else print 0}}\'', shell=True)

                # Convert to float
                avg_cov_genome = float(avg_cov_genome.strip())
                avg_cov_mt = float(avg_cov_mt.strip())

                # percentage_threshold = 1 / (float(avg_cov_mt) / float(avg_cov_genome) + 1)
                if avg_cov_genome != 0 and avg_cov_mt != 0:
                    percentage_threshold = 1 / (avg_cov_mt / avg_cov_genome + 1)
                else:
                    percentage_threshold = 0

                # Write the values to the output file
                with open(out_threshold_file, 'w') as f:
                    f.write(f"percentage_threshold: {percentage_threshold}\n")
                    f.write(f"avg_cov_genome: {avg_cov_genome}\n")
                    f.write(f"avg_cov_mt: {avg_cov_mt}\n")

                print("Calculate percentage threshold done!")

    # if chloroplast != 'None':
        # to be continue

    # 03_filter_by_samtools
    if os.path.exists(out_csv_filtered_sam) or os.path.exists(out_filtered_sam):
        print('Alignment might have been filtered already.  Skip samtools.')
    else:
        print("Filter bwa-meme's output")
        if is_pair_read == 1:
            cmd = 'samtools view -@ 2 -h -f 2 -F 0x900 -q %s %s' % (alignment_quality , out_sorted_bam)
        else:
            cmd = 'samtools view -@ 2 -h -F 0x904 -q %s %s' % (alignment_quality , out_sorted_bam)
        try:
            output = subprocess.check_call(cmd, shell=True, stdout=open(out_filtered_sam, 'w'))
            # output.wait()
        except:
            no_error = False
            log_error(cmd, output, sys.exc_info())

    if os.path.exists(out_csv_filtered_sam):
        print('Filter alignments have been filtered already.  Skip Filter alignments.')
    else:
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

        # rm all other sam files
        if os.path.exists(out_sorted_bai):
            try:
                subprocess.check_call('rm %s && rm %s && rm %s && rm %s && rm %s' % (out_sam, out_bam, out_markdup, out_sorted_bam, out_sorted_bai), shell=True)
                # subprocess.check_call('rm %s && rm %s && rm %s' % (out_sam, out_bam, out_markdup), shell=True)
            except subprocess.CalledProcessError as e:
                no_error = False
                log_error(cmd, e.output, e)

    # # If 'rm_alter_align' is set to 0, then after deleting `out_csv_filtered_sam`, the filtering process can be repeated.
    # if rm_alter_align == 1:
    if os.path.exists(out_csv_filtered_sam):
        if os.path.exists(out_filtered_sam):
            try:
                subprocess.check_call('rm %s' % out_filtered_sam, shell=True)
            except subprocess.CalledProcessError as e:
                no_error = False
                log_error(cmd, e.output, e)

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
