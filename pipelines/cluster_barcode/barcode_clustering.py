from __future__ import barry_as_FLUFL

__all__  =  ['read1', 'read2', 'read1_out', 'read2_out', 'consolidate_read1', 'consolidate_read2','min_qual', 'min_freq', 'logger_umi_process' , 'logger_umi_errors']
__version__  =  '1.0'
__author__  =  'Maggie Ruimin Sun'

import os
import re
import sys
import gzip
import time
import logging
import subprocess
from difflib import SequenceMatcher

def read_fq(filename):
    if re.search('.gz$', filename):
        fastq = gzip.open(filename, 'rb')
    else:
        fastq = open(filename)
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]

def umi_tagging(read1, read2, read1_out, read2_out, logger_umi_process, logger_umi_errors):
    
    # Sort FASTQ files based on molecular barcodes.
    cmd = 'cat '+read1+' | paste - - - - | sort -k 3,3 -k 1,1 | tr "\t" "\n" > '+read1_out 
    logger_umi_process.info(cmd)
    subprocess.check_call(cmd, shell=True, env=os.environ.copy())
    cmd = 'cat '+read2+' | paste - - - - | sort -k 3,3 -k 1,1 | tr "\t" "\n" > '+read2_out
    logger_umi_process.info(cmd)
    subprocess.check_call(cmd, shell=True, env=os.environ.copy())

def edit_dist(s1, s2):
    matcher = SequenceMatcher(None, s1, s2)
    dist = len(s1) + len(s2)
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == 'equal':
            dist -= (i2 + j2 - i1 -j1)
        elif tag == 'replace':
            dist -= (i2-i1)
    return dist

# Based on the molecular barcode tagged in the previous step, wrap up 
# all reads having the most similar molecular barcodes (edit distance < 1).
# Returned values: the molecular ID corresponding to the barcode sequences ;
#                   the sample ID, which is that of the very first reads;
#                   the bins of reads of the similar barcodes.
def read_bins(fastq_file):
    num_read = 0
    fastq = open(fastq_file)
    reads = fastq.readlines() 
    total_reads = len(reads) / 4
    fastq.close()
    header_out = open(fastq_file+'_ids.txt', 'w')

    i = 0
    cur_molecular_id = ''
    while i < total_reads:
        ri = 4*i
        header = reads[ri]
        seq = reads[ri+1]
        qual = reads[ri+3]
        num_read += 1
        sample_id, info, barcodes = header.strip().split(' ')
        molecular_id = barcodes[0:12]
        umi_qual = barcodes[13:]
        qual = [str(ord(x)-33) for x in qual.strip()]
        read = [header.rstrip(), seq.rstrip(),'+',qual]
        header_out.write(sample_id+'\t'+molecular_id+'\n')
        i += 1
        if molecular_id == cur_molecular_id:
        #if edit_dist(molecular_id, cur_molecular_id) < 1:
            bin_reads.append(read)
        else:
            if cur_molecular_id != '':
                yield cur_molecular_id, cur_sample_id, cur_umi_qual, bin_reads
            cur_molecular_id = molecular_id
            cur_sample_id = sample_id
            cur_umi_qual = umi_qual
            bin_reads = [read]
    header_out.close()
    yield cur_molecular_id, cur_sample_id, cur_umi_qual, bin_reads

def consolidate_position(bases, quals, min_qual, min_freq):
    num = {}
    qual = {}
    num['A'] = num['C'] = num['G'] = num['T'] = num['N'] = 0
    qual['A'] = qual['C'] = qual['G'] = qual['T'] = qual['N'] = 0
    for bb, qq in zip(bases, quals):
        qq = int(qq)
        if qq > min_qual:
            num[bb] += 1
        if qq > qual[bb]:
            qual[bb] = qq
    most_common_base = max(num.keys(), key=(lambda key:num[key]))
    freq = float(num[most_common_base]) / len(bases)
    if freq > min_freq:
        return True, most_common_base, qual[most_common_base]
    else:
        return False, 'N', 0

def consolidate(fastq_file, consolidate_fq_file, min_qual, min_freq, 
    logger_umi_process, logger_umi_errors):
    logger_umi_process.info('Consolidating reads in %s', fastq_file)
    consolidate = open(consolidate_fq_file, 'w')
    bins = read_bins(fastq_file)
    num_input_reads = 0
    num_consolidate_reads = 0
    num_success = 0
    num_bases = 0

    for cur_molecular_id, cur_sample_id, cur_umi_qual, reads in bins:
        num_input_reads += len(reads)
        if len(reads) < 5:
            continue
        num_consolidate_reads += 1
        read_bases = zip(*[list(read[1]) for read in reads])
        read_quals = zip(*[list(read[3]) for read in reads])
        consolidation_success, cons_seq, cons_qual = zip(*[
            consolidate_position( bases, quals, min_qual, min_freq) 
            for bases, quals in zip(read_bases, read_quals)
        ])
        num_success += sum(consolidation_success)
        num_bases += len(consolidation_success)
        num_Ns = len(consolidation_success) - sum(consolidation_success)
        if num_Ns > 20:
            continue
        consolidate.write('%s;%s;%s\n' % (cur_sample_id, cur_molecular_id, cur_umi_qual))
        consolidate.write(''.join(cons_seq)+'\n+\n')
        consolidate.write(''.join([chr(q+33) for q in cons_qual])+'\n')
    consolidate.close()
    print('Number of input reads == '+str(num_input_reads))
    print('Number of consolidated reads == '+str(num_consolidate_reads))
    logger_umi_process.info("Analyzed %d input reads", num_input_reads)
    logger_umi_process.info("Wrote %d consolidated reads to %s", 
        num_consolidate_reads, consolidate_fq_file)
    logger_umi_process.info("Successfully consolidated %d bases out of %d (%.2f%%)",
        num_success, num_bases, 100*float(num_success)/num_bases)