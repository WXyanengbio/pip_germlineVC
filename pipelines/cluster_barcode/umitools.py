from __future__ import barry_as_FLUFL

__all__  = ['samtools_dir', 'umitools_dir', 'filtered_sam' ,'filtered_bam' , 'sorted_bam', 'umitool_stats' , 'umis_sam', 'edit_dist', 'logger_umi_process', 'logger_umi_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'


import os
import sys
import shlex
import subprocess

#put the info output to the log
def stdout_err(command):
    command_pope = shlex.split(command)
    print(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr

def umitool(samtools_dir, umitools_dir, filtered_sam ,filtered_bam , sorted_bam, umitool_stats , umis_sam, edit_dist, logger_umi_process, logger_umi_errors):
    command1 = samtools_dir + ' view -bS ' + filtered_sam + ' > ' + filtered_bam
    logger_umi_process.info('Samtools transform sam to bam.')
    #(status, output) = subprocess.getstatusoutput(command1)
    os.system(command1)
    command2 = samtools_dir + ' sort ' + filtered_bam + ' > ' + sorted_bam
    logger_umi_process.info('Samtools sort bam.')
    #output2 = os.popen(command2)
    os.system(command2)
    #print(output2)
   # logger_umi_process.info(output2)
    command3 = samtools_dir + ' index ' + sorted_bam
    logger_umi_process.info('Samtools build index of bam.')
    os.system(command3)
    command4 = 'python3.6 {0}  dedup -I {1} --output-stats={2} -S {3} --edit-distance-threshold {4} --paired True'.format(umitools_dir,
                sorted_bam,umitool_stats,filtered_bam, edit_dist)
    logger_umi_process.info('UMIs-tools cluster bam.')
    (status, output) = subprocess.getstatusoutput(command4)
    logger_umi_process.info(output)
    command5 = samtools_dir + ' view -h ' + filtered_bam + ' > ' + umis_sam
    logger_umi_process.info('Samtools transform umis_bam to umis_sam.')
    os.system(command5)
