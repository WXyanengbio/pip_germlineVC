from __future__ import barry_as_FLUFL

__all__  =  ['bwa_dir' , 'samtools_dir' ,'ref_fa_file' , 'ref_index_name','total_ref_fa_file' ,'exome_target_bed' ,'read1' , 'read2' , 'out_file' , 'num_threads' , 'logger_bwa_process' , 'logger_bwa_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'


import os
import logging
import time
import sys
import multiprocessing
import shlex
import subprocess


# samtools faidx to get target reference fasta based on the bed
def samtools_faidx(bed_dict,logger_bwa_process, logger_bwa_errors):
    samtools_dir = bed_dict[0]
    total_ref_fa_file = bed_dict[1]
    bed = bed_dict[2]
    tmp_target_ref = bed_dict[3]
    command = '{0} faidx {1} {2} >> {3}'.format(samtools_dir, total_ref_fa_file, bed, tmp_target_ref)
    stdout, stderr = stdout_err(command)
    logger_bwa_process.info(stdout)
    logger_bwa_errors.info(stderr)

# use samtools faidx to get target reference fasta based on the bed file
def filter_ref_fa_by_bed(samtools_dir, ref_fa_file, ref_index_name, 
                         exome_target_bed,total_ref_fa_file, 
                         logger_bwa_process, 
                    logger_bwa_errors):
    if not os.path.isfile(exome_target_bed):
        logger_bwa_errors.error('%s does not exist!', exome_target_bed)
        print("Error: cannot find target bed file!")
        exit()
    if not os.path.isfile(ref_fa_file):
        logger_bwa_errors.error('%s does not exist!', ref_fa_file)
        print("Error: cannot find reference genome file!")
        exit()
    if not os.path.isfile(ref_index_name + '.refSeq.fa'):
        tmp_ref_dir = os.path.dirname(ref_index_name) + '/' + 'tmp_target_ref.fa'
        bed= open(exome_target_bed, 'U')
        bed_dicts = []
        for line in bed:
            if line.startswith('chr'):
                line= line.strip().split('\t')
                outname= line[0] + ':' + line[1]+ '-' + str(int(line[2])+1)
                #print('>' + outname)
                bed_dicts.append([samtools_dir, total_ref_fa_file, outname, tmp_ref_dir])
            else:
                continue
        bed.close()
        pool = multiprocessing.Pool(processes=1)
        pool.map(samtools_faidx, bed_dicts)
        print("--" * 20)
        pool.close()   # 关闭pool, 则不会有新的进程添加进去
        pool.join()    # 必须在join之前close, 然后join等待pool中所有的线程执行完毕
        fasta= open(tmp_ref_dir, 'U')
        ref_dir = open(ref_index_name + '.refSeq.fa','w')
        for line in fasta:
            line= line.strip()
            if line == '':
                continue
            if line.startswith('>'):
                seqname= line
                seqname1= seqname.replace(':', '_')
                seqname1= seqname1.replace('-', '_')
                ref_dir.write(seqname1+'\n')
            else:
                ref_dir.write(line +'\n')
        ref_dir.close()
        fasta.close()
        os.system('rm -rf {0}'.format(tmp_ref_dir))
    return ref_index_name + '.refSeq.fa'

def align_reads_bwa(bwa_dir, samtools_dir,
                    ref_fa_file, ref_index_name, exome_target_bed,total_ref_fa_file,
                    read1, read2, 
                    out_file, num_threads, logger_bwa_process, 
                    logger_bwa_errors): 
    if not os.path.isfile(read1):
        logger_bwa_errors.error('%s does not exist!', read1)
        print("Error: cannot find NGS read file!")
        return 1
    ref_fa_file_bed = filter_ref_fa_by_bed(samtools_dir, ref_fa_file, ref_index_name, exome_target_bed, total_ref_fa_file, logger_bwa_process, 
                                           logger_bwa_errors)
    index_file_extensions = ['.pac', '.amb','.ann','.bwt', '.sa']
    genome_indexed = True
    for extension in index_file_extensions:
        if not os.path.isfile(ref_index_name + extension):
            genome_indexed = False
            break
    if not genome_indexed:
        bwa_index_command = '{0} index -p {1} {2}'.format(
            bwa_dir, ref_index_name, ref_fa_file_bed)
        (status, output) = subprocess.getstatusoutput(bwa_index_command)
        logger_bwa_process.info(output)
        print('BWA genome index files have been built.')
    else:
        print('BWA genome index files exist.')

    bwa_align_command = '{0} mem {1} {2} {3} -t {4} > {5}'.format(bwa_dir, ref_index_name, read1, read2, num_threads, out_file)
    logger_bwa_process.info(bwa_align_command)
    #os.system(bwa_align_command)

    (status, output) = subprocess.getstatusoutput(bwa_align_command)
    logger_bwa_process.info(output)

    print('BWA alignment has been completed.')
    logger_bwa_process.info('BWA alignment has been completed.')
    return 0

def align_reads_bwa_based_all(bwa_dir, total_ref_fa_file,
                              read1, read2, 
                              out_file, num_threads, 
                              logger_bwa_process, logger_bwa_errors):
    if not os.path.isfile(read1):
        logger_bwa_errors.error('%s does not exist!', read1)
        print("Error: cannot find NGS read file!")
        return 1
    if not os.path.isfile(total_ref_fa_file):
        logger_bwa_errors.error('%s does not exist!', total_ref_fa_file)
        print("Error: cannot find reference fasta!")
        return 1
    ref_index_name = total_ref_fa_file.rstrip('.fasta')
    index_file_extensions = ['.pac', '.amb','.ann','.bwt', '.sa']
    genome_indexed = True
    for extension in index_file_extensions:
        if not os.path.isfile(ref_index_name + extension):
            genome_indexed = False
            break
    if not genome_indexed:
        bwa_index_command = '{0} index -p {1} {2}'.format(
            bwa_dir, ref_index_name, ref_fa_file_bed)
        logger_bwa_process.info(bwa_index_command)
        (status, output) = subprocess.getstatusoutput(bwa_index_command)
        logger_bwa_process.info(output)
        print('BWA genome index files have been built.')
    else:
        print('BWA genome index files exist.')

    bwa_align_command = '{0} mem -t {1} {2} {3} {4} > {5}'.format(
        bwa_dir, num_threads, ref_index_name, read1, read2, out_file)
    (status, output) = subprocess.getstatusoutput(bwa_align_command)
    logger_bwa_process.info(output)
    logger_bwa_errors.info(stderr)
    print('BWA alignment has been completed.')
    logger_bwa_process.info('BWA alignment has been completed.')
    return 0