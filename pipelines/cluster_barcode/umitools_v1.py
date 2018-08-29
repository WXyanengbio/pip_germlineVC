from __future__ import barry_as_FLUFL

__all__  = ['samtools_dir', 'umitools_dir', 'filtered_sam' ,'filtered_bam' , 'sorted_bam', 'umitool_stats' , 'umis_sam', 'edit_dist', 'logger_umi_process', 'logger_umi_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'


import os
import sys
import shlex
import subprocess
sys.path.append("..")
from pipelines.log.log_v1 import store_cluster_logs
#put the info output to the log
def stdout_err(command):
    command_pope = shlex.split(command)
    print(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr

def umitool(samtools_dir, umitools_dir, filtered_sam ,filtered_bam , sorted_bam, umitool_stats , umis_sam, edit_dist, logger_umi_process, logger_umi_errors):
    command1 = samtools_dir + ' view -bS ' + filtered_sam + ' -o ' + filtered_bam
    store_cluster_logs(logger_umi_process,'null','Samtools transform sam to bam.')
    (status, output) = subprocess.getstatusoutput(command1)
    command2 = samtools_dir + ' sort ' + filtered_bam + ' -o ' + sorted_bam
    store_cluster_logs(logger_umi_process,'null','Samtools sort bam.')
    (status, output) = subprocess.getstatusoutput(command2)
    command3 = samtools_dir + ' index ' + sorted_bam
    store_cluster_logs(logger_umi_process,'null','Samtools build index of bam.')
    os.system(command3)
    command4 = 'python3.6 {0}  dedup -I {1} --output-stats={2} -S {3} --edit-distance-threshold {4} --paired True'.format(umitools_dir,
                sorted_bam,umitool_stats,filtered_bam, edit_dist)
    store_cluster_logs(logger_umi_process,'null','UMIs-tools cluster bam.')
    (status, output) = subprocess.getstatusoutput(command4)
    store_cluster_logs(logger_umi_process,'null',output)
    command5 = samtools_dir + ' view -h ' + filtered_bam + ' -o ' + umis_sam
    store_cluster_logs(logger_umi_process,'null','Samtools transform umis_bam to umis_sam.')
    (status, output) = subprocess.getstatusoutput(command5)
