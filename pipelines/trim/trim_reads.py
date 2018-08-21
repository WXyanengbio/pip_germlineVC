from __future__ import barry_as_FLUFL

__all__  =  ['read1' , 'read2' , 'trimmed1' ,'un_trimmed1', 'trimmed2' ,'un_trimmed2', 'min_read_len' , 'common_seq1' , 'common_seq2' , 'stats_file' , 'logger_trim_process' , 'logger_trim_errors']
__version__  =  '1.0'
__author__  =  'Maggie Ruimin Sun'

import logging
import os
import re
import sys
import time
import gzip
import itertools
import shlex
import subprocess

#put the info output to the log
def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr
#----------------------------------------------------------
def read_fq(file_name, logger_trim_process, logger_trim_errors):
    if not os.path.isfile(file_name):
        logger_trim_errors.error("%s does not exist!\n", file_name)
        print(file_name + ' does not exist!')
    if re.search('.gz$', file_name):
        fastq = gzip.open(file_name, 'r')
    else:
        fastq = open(file_name)

    with fastq as f:
        while True:
            l1 = str(f.readline(), 'utf-8')
            if not l1:
                break
            l2 = str(f.readline(), 'utf-8')
            l3 = str(f.readline(), 'utf-8')
            l4 = str(f.readline(), 'utf-8')
            yield [l1, l2, l3, l4]


def trim_read1(r1, common_seq2, mt_barcode):
    l1 = len(r1[1].strip())
    r1_end3 = r1[1].strip()[(l1 - 23):]
    trim_seq = common_seq2 + mt_barcode
    for i in range(23):
        if r1_end3[i:] == trim_seq[0:(23 - i)]:
            break
    pos_trim_r1 = l1 - 23 + i
    return pos_trim_r1, [r1[0], r1[1].strip()[0:pos_trim_r1] + '\n', r1[2], r1[3].strip()[0:pos_trim_r1] + '\n']


def trim_read2(r2, common_seq1):
    l2 = len(r2[1].strip())
    r2_end3 = r2[1].strip()[(l2 - 21):]
    for i in range(21):
        if r2_end3[i:] == common_seq1[0:(21 - i)]:
            break
    pos_trim_r2 = l2 - 21 + i
    return pos_trim_r2, [r2[0], r2[1].strip()[0:pos_trim_r2] + '\n', r2[2], r2[3].strip()[0:pos_trim_r2] + '\n']


def trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len, common_seq1,
                    common_seq2, stats_file, logger_trim_process, logger_trim_errors):
    time_start = time.time()
    num_total_reads = 0
    num_short_reads = 0
    num_error_reads1 = 0
    num_error_reads2 = 0
    fout1 = open(trimmed1, 'w')
    fout2 = open(trimmed2, 'w')

    for r1, r2 in zip(read_fq(read1, logger_trim_process, logger_trim_errors),
                      read_fq(read2, logger_trim_process, logger_trim_errors)):
        num_total_reads += 1
        if r1[0][0] != '@' or r2[0][0] != '@':
            num_error_reads1 += 1
            logger_trim_process.info("Error read pair: \n\t"
                                     + '\t'.join(r1) + '\n\t' + '\t'.join(r2) + '\n')
        else:
            start_common = r2[1].find(common_seq2)
            if start_common < 12:
                num_error_reads2 += 1
                #logger_trim_process.info("Error barcode/common seqs:" + str(start_common) + "\n\t"
                #                         + '\t'.join(r1) + '\n\t' + '\t'.join(r2) + '\n')
            else:
                umi = r2[1][(start_common - 12):start_common]
                qua = r2[3][(start_common - 12):start_common]
                r2[1] = r2[1][(start_common + 11):]
                r2[3] = r2[3][(start_common + 11):]
                pos_trim_r1, r1 = trim_read1(r1, common_seq2, umi)
                pos_trim_r2, r2 = trim_read2(r2, common_seq1)
                if pos_trim_r1 < min_read_len or pos_trim_r2 < min_read_len:
                    num_short_reads += 1
                   # logger_trim_process.info("Short read pair: \n\t"
                    #                         + '\t'.join(r1) + '\n\t' + '\t'.join(r2) + '\n')
                else:
                    #umi = umi + ';' + qua
                    h1 = r1[0].split(' ')[0] + '_' + umi + ' ' + r1[0].split(' ')[1]
                    h2 = r2[0].split(' ')[0] + '_' + umi + ' ' + r2[0].split(' ')[1]
                    fout1.write(h1 + r1[1] + r1[2] + r1[3])
                    fout2.write(h2 + r2[1] + r2[2] + r2[3])
    fout1.close()
    fout2.close()
    print('Total number of reads == ' + str(num_total_reads))
    print('Number of short reads (<{0}bp) == {1}'.format(min_read_len, num_short_reads))
    print('Number of error reads == ' + str(num_error_reads1+num_error_reads2))

    stats_out = open(stats_file, 'w')
    stats_out.write('Total number of reads == ' + str(num_total_reads) + '\n')
    stats_out.write('Number of short reads (either read_length <{0}bp) == {1}\n'.format(
        min_read_len, num_short_reads))
    stats_out.write('Number of unproper read pairs (containing incorrect headers) == ' + str(num_error_reads1) + '\n')
    stats_out.write('Number of read pairs without correct common sequences/MTs == ' + str(num_error_reads2) + '\n')
    stats_out.write('The time of trimming is %s minutes.' % str((time.time() - time_start) / 60))

def trim_read_pairs_by_trimmomatic(trimmomatic_dir,
                                   read1, read2, 
                                   trimmed1, un_trimmed1,
                                   trimmed2, un_trimmed2,
                                   min_read_len,
                                   stats_file, logger_trim_process,
                                   logger_trim_errors):
    if not os.path.isfile(read1):
        logger_trim_errors.error('%s does not exist!', read1)
        print("Error: cannot find NGS read file!")
        exit()
    if not os.path.isfile(trimmomatic_dir):
        logger_trim_errors.error('%s does not exist!', trimmomatic_dir)
        print("Error: cannot find trimmomatic.jar!")
        exit()
    #trimmomatic_log = os.path.dirname(trimmed1) + '/' +'trimmomatic.log'
    command = 'java -jar {0} PE -threads 1 -phred33 -summary {1} {2} {3} {4} {5} {6} {7} ILLUMINACLIP:{8}:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:{9} '.format(
        trimmomatic_dir, stats_file, read1, read2, trimmed1, un_trimmed1, trimmed2, un_trimmed2, os.path.dirname(trimmomatic_dir) + '/adapters/TruSeq3-PE.fa', min_read_len)
    stdout, stderr = stdout_err(command)
    logger_trim_process.info(stdout)
    logger_trim_errors.info(stderr)
    #os.system(command)