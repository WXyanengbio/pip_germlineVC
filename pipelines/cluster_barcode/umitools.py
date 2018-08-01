from __future__ import barry_as_FLUFL

__all__  = ['samtools_dir', 'umitools_dir', 'filtered_sam' ,'filtered_bam' , 'sorted_bam', 'umitool_stats' , 'umis_sam', 'edit_dist', 'logger_umi_process', 'logger_umi_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'


import os
import sys


def umitool(samtools_dir, umitools_dir, filtered_sam ,filtered_bam , sorted_bam, umitool_stats , umis_sam, edit_dist, logger_umi_process, logger_umi_errors):
    cmd1 = samtools_dir + ' view -bS ' + filtered_sam + ' > ' + filtered_bam
    logger_umi_process.info('Samtools transform sam to bam.')
    os.system(cmd1)
    cmd2 = samtools_dir + ' sort ' + filtered_bam + ' > ' + sorted_bam
    logger_umi_process.info('Samtools sort bam.')
    os.system(cmd2)
    cmd3 = samtools_dir + ' index ' + sorted_bam
    logger_umi_process.info('Samtools build index of bam.')
    os.system(cmd3)
    cmd4 = 'python3.6 ' + umitools_dir + ' dedup -I ' + sorted_bam + ' --output-stats=' + umitool_stats + ' -S ' + filtered_bam + ' --edit-distance-threshold ' + str(edit_dist) + ' --paired True'
    logger_umi_process.info('UMIs-tools cluster bam.')
    os.system(cmd4)
    cmd5 = samtools_dir + ' view -h ' + filtered_bam + ' > ' + umis_sam
    logger_umi_process.info('Samtools transform umis_bam to umis_sam.')
    os.system(cmd5)
