from __future__ import barry_as_FLUFL

__all__  =  ['benchmark_dir','truth_vcf', 'filter_vcf', 'confident_region_bed', 'benchmarking_dir', 'total_ref_fa_file', 'exon_interval', 'num_threads','logger_benchmark_process', 'logger_benchmark_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import os
import sys


def hap_py(benchmark_dir,truth_vcf, filter_vcf, confident_region_bed, benchmarking_dir, total_ref_fa_file, exon_interval, num_threads, logger_benchmark_process, logger_benchmark_errors):
    if exon_interval != 'all':
        exon_bed = benchmarking_dir + '/' + os.path.basename(exon_interval).rstrip(".list") + ".bed"
        if not os.path.isfile(exon_bed):
            logger_benchmark_errors.error("%s does not exist!\n", exon_bed)
            print(exon_bed + ' does not exist!')
            cmd1 = 'grep ^chr {0} | cut -f 1-3 > {1}'.format(exon_interval, exon_bed)
            logger_benchmark_process.info('Build the target bed to benchmark.')
            os.system(cmd1)
    if not os.path.isfile(filter_vcf):
        logger_benchmark_errors.error("%s does not exist!\n", filter_vcf)
        print(filter_vcf + ' does not exist!')
    else:
        benchmarout = benchmarking_dir + '/' + os.path.basename(filter_vcf).rstrip(".vcf")
        if exon_interval != 'all':
            cmd2 = '{0} {1} {2} -f {3} -r {4} -T {5} -o {6} --threads {7}'.format(benchmark_dir, truth_vcf, filter_vcf,confident_region_bed,total_ref_fa_file ,exon_bed, benchmarout,num_threads)
        else:
            cmd2 = '{0} {1} {2} -f {3} -r {4} -o {5} --threads {6}'.format(benchmark_dir, truth_vcf, filter_vcf,confident_region_bed,total_ref_fa_file , benchmarout,num_threads)
        logger_benchmark_process.info('hap.py benchmark.')
        os.system(cmd2)

