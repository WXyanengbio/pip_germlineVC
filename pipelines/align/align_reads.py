from __future__ import barry_as_FLUFL

__all__  =  ['bwa_dir' , 'ref_fa_file' , 'ref_index_name' , 'read1' , 'read2' , 'out_file' , 'num_threads' , 'logger_bwa_process' , 'logger_bwa_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'


import os
import logging
import time
import sys


def align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, read1, read2, 
                    out_file, num_threads, logger_bwa_process, 
                    logger_bwa_errors): 
    if not os.path.isfile(read1):
        logger_bwa_errors.error('%s does not exist!', read1)
        print("Error: cannot find NGS read file!")
        return 1

    index_file_extensions = ['.pac', '.amb','.ann','.bwt', '.sa']
    genome_indexed = True
    for extension in index_file_extensions:
        if not os.path.isfile(ref_index_name + extension):
            genome_indexed = False
            break
    if not genome_indexed:
        bwa_index_command = '{0} index -p {1} {2}'.format(
            bwa_dir, ref_index_name, ref_fa_file)
        logger_bwa_process.info(bwa_index_command)
        os.system(bwa_index_command)
        print('BWA genome index files have been built.')
    else:
        print('BWA genome index files exist.')

    bwa_align_command = '{0} mem -t {1} {2} {3} {4} > {5}'.format(
        bwa_dir, num_threads, ref_index_name, read1, read2, out_file)
    print(bwa_align_command)
    os.system(bwa_align_command)
    print('BWA alignment has been completed.')
    logger_bwa_process.info('BWA alignment has been completed.')
    return 0