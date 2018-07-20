from __future__ import barry_as_FLUFL

__all__  =  ['benchmark_dir','truth_vcf', 'filter_vcf', 'confident_region_bed', 'benchmarking_dir', 'total_ref_fa_file', 'Exon_Interval', 'logger_benchmark_process', 'logger_benchmark_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import os
import re
import sys
import time
import gzip
import itertools



def hap_py(benchmark_dir,truth_vcf, filter_vcf, confident_region_bed, benchmarking_dir, total_ref_fa_file, Exon_Interval, logger_benchmark_process, logger_benchmark_errors):
    Exon_bed = benchmarking_dir + '/' + os.path.basename(Exon_Interval).rstrip(".list") + ".bed"
    if not os.path.isfile(Exon_bed):
        logger_benchmark_errors.error("%s does not exist!\n", Exon_bed)
        print(Exon_bed + ' does not exist!')
        cmd1 = 'grep ^chr {0}.list | cut -f 1-3 > {1}'.format(Exon_Interval, Exon_bed)
        logger_benchmark_process.info('Build the target bed to benchmark.')
        os.system(cmd1)
    if not os.path.isfile(filter_vcf):
        logger_benchmark_errors.error("%s does not exist!\n", filter_vcf)
        print(filter_vcf + ' does not exist!')
    else:
        benchmarout = benchmarking_dir + '/' + os.path.basename(filter_vcf).rstrip(".vcf")
        cmd2 = '{0} {1} {2} -f {3} -r {4} -T {5} -o {6}'.format(benchmark_dir, truth_vcf, filter_vcf,confident_region_bed,total_ref_fa_file ,Exon_bed, benchmarout)
        logger_benchmark_process.info('hap.py benchmark.')
        os.system(cmd2)
    

