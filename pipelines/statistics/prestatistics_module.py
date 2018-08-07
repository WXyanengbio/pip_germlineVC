__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
import re
import sys
import time
import itertools
from itertools import groupby
import shlex
import subprocess

#put the info output to the log
def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr

# get the range of the positions of the N bases in the reads 
def split_n_bases(nbases):
    if nbases == 'NULL':
        return ["NULL" ,"NULL"]
    else:
        bases = nbases.split(",")
       # print(bases)
        bases_posi =[]
        for posi in range(0,len(bases)):
            posi = bases[posi]
            #print(posi)
            if "-" in str(posi):
                for i in range(int(posi.split("-")[0]), int(posi.split("-")[1])+1):
                    bases_posi.append(i)
            else:
                bases_posi.append(int(posi))
    #print(bases_posi)
        l = bases_posi
        a = list()
        b = []
        result = []
        function = lambda x: x[1] - x[0]
        for k, g in groupby(enumerate(l), function):
            g = list(g)
            b.append(k)
            a.append(g)
        for c in range(len(b)):
            if [v for i, v in a[c]][0] != [v for i, v in a[c]][-1]:
                result.append("%d-%d" % ([v for i, v in a[c]][0], [v for i, v in a[c]][-1]))
            if [v for i, v in a[c]][0] == [v for i, v in a[c]][-1]:
                result.append([v for i, v in a[c]][0])
    #print(result)
        if "-" not in str(result[len(result)-1]):
            result_last = result[len(result)-1]
        else:
            result_last = result[len(result)-1].split('-')[1]
        return ["[" + str(result[0]) + "]" , str(result_last)]

# get informations from the qc result of the reads
def getinfo(fastqc,logger_statistics_process, logger_statistics_errors):
    qc_read = fastqc.split(".fastq")[0] + '_fastqc.zip'
    command1 = 'unzip' + ' -o ' + qc_read  + ' -d ' + os.path.dirname(qc_read)
    #os.system(command1)
    stdout, stderr = stdout_err(command1)
    logger_statistics_process.info(stdout)
    logger_statistics_errors.info(stderr)
    print('{0} has been completed.'.format(qc_read))
    qc_data = qc_read.rstrip(".zip") + '/' + 'fastqc_data.txt'
    qcdata = open(qc_data,"r")
    modules = 0
    # 定义正则表达式
    value = re.compile(r'\d+')
    perbasesequencequalit_posi = []
    perbasesequencequalit = []
    persequencequalityscores = []
    persequencequalityreads = []
    per_basen1 = []
    perbase_ncontent = []
    f = qcdata.readlines()
    for lines in range(0 ,len(f)):
        line = f[lines]
        line = line.strip()
        if '>>END_MODULE' in line:
            modules = modules + 1
        if 'Total Sequences' in line:
            raw_reads = re.findall('\d+',line)
        if 'Sequence length' in line:
            #mode = re.compile(r'\d+-\d+')
            #seq_length = mode.findall(line)
            seq_length = line.split('\t')[1]
            if '-' in seq_length:
                seq_length_min,seq_length_max = seq_length.split('-')
            else:
                seq_length_min = str(seq_length)
                seq_length_max = str(seq_length)
            #print(seq_length)
        if '%GC' in line:
            gc = re.findall('\d+',line)
        if modules == 1 and value.match(line[0]):
            perbasesequencequalit_posi.append(line.split('\t')[0])
            perbasesequencequalit.append(float(line.split('\t')[1]))
        if modules == 3 and value.match(line[0]):
            persequencequalityscores.append(line.split('\t')[0])
            persequencequalityreads.append(float(line.split('\t')[1]))
        if modules == 6 and value.match(line[0]):
            per_basen1.append(line.split('\t')[0])
            perbase_ncontent.append(line.split('\t')[1])
    #--Per base sequence quality
    perbasesequencequalit_mean = sum(perbasesequencequalit)/len(perbasesequencequalit)
    lowqualit_bases= []
    for i in range(0,len(perbasesequencequalit)):
        if float(perbasesequencequalit[i]) < perbasesequencequalit_mean-2:
            lowqualit_bases.append(perbasesequencequalit_posi[i])
    if len(lowqualit_bases) == 0:
        lowqualit_bases.append('NULL')
    #--per sequence quality
    persequencequalityreads1 = []
    for i in range(0,len(persequencequalityreads)):
        persequencequalityreads1.append(sum(persequencequalityreads[i:])/int(raw_reads[0]))
    for i in range(0,len(persequencequalityscores)):
        if persequencequalityscores[i] == '20':
            q20 = persequencequalityreads1[i]
        if persequencequalityscores[i] == '30':
            q30 = persequencequalityreads1[i]
    if '20' not in persequencequalityscores:
        q20 = 1.0000
    #--N content
    nbases = []
    for i in range(0,len(perbase_ncontent)):
        if float(perbase_ncontent[i]) > 0:
            nbases.append(per_basen1[i])
    if len(nbases) == 0:
        nbases.append('NULL')
    #---
    #print(raw_reads)
    #print(seq_length)
    #print(perbasesequencequalit_mean)
    #print([raw_reads[0], seq_length[0], GC[0], str(round(Perbasesequencequalit_mean,3)), ','.join(lowqualitBases), str('%.3f%%' % (Q20 * 100)), str('%.3f%%' % (Q30 * 100)), ','.join(nbases)]
    return [raw_reads[0], seq_length_min, seq_length_max, gc[0], str(round(perbasesequencequalit_mean,3)), 
            split_n_bases(','.join(lowqualit_bases))[0], str('%.3f%%' % (q20 * 100)), str('%.3f%%' % (q30 * 100)), split_n_bases(','.join(nbases))[0], split_n_bases(','.join(nbases))[1]] 

#--QC reads by fastQC
def qc_raw_reads(fastQC_dir, out_dir, sample, module, read1, read2,logger_statistics_process, logger_statistics_errors):
    command1 = '{0} {1} {2} -o {3}'.format(fastQC_dir, read1, read2, out_dir)
    stdout, stderr = stdout_err(command1)
    logger_statistics_process.info(stdout)
    logger_statistics_errors.info(stderr)
    #os.system(command1)
    print('QC-{0} has been completed.\n'.format(module))
    qc_statistics = out_dir + '/' + sample + '.' + module +'.statistics.txt'
    qc_result1 =getinfo(out_dir + '/' + os.path.basename(read1),logger_statistics_process, logger_statistics_errors)
    qc_result2 =getinfo(out_dir + '/' + os.path.basename(read2),logger_statistics_process, logger_statistics_errors)
    fout = open(qc_statistics, 'w')
    fout.write('\t'.join(['SampleID','Sequence direction','raw reads', 'min length','max length', 'GC content', 'mean of Per base qualit', 'low qualit Bases position', 'Q20', 'Q30', 'N_bases position']) + '\n')
    fout.write('\t'.join([sample, read1.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result1[0:9])]) + '\n')
    fout.write('\t'.join([sample, read2.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result2[0:9])]) + '\n')
    fout.close()
    return qc_result1, qc_result2

#-get the depth and coverage of the mapping result
def statistics_depth_coverage(samtools_dir, sam_bam, out_dir,sample, module, exome_target_bed, logger_statistics_process, logger_statistics_errors):
    if not os.path.isfile(sam_bam):
        logger_statistics_errors.error("%s does not exist!\n", sam_bam)
        print(sam_bam + ' does not exist!')
    if re.search('.sam$', sam_bam):
        bam = sam_bam.rstrip('.sam') + '.bam'
        command1 = samtools_dir + ' view -bS ' + sam_bam + ' > ' + bam
        os.system(command1)
        logger_statistics_process.info('{0} has been tranformed to bam.'.format(sam_bam))
        sorted_bam = sam_bam.rstrip('.sam') + '_sorted.bam'
    elif re.search('.bam$', sam_bam):
        bam = sam_bam
        sorted_bam = sam_bam.rstrip('.bam') + '_sorted.bam'
    if not os.path.isfile(sorted_bam):
        print(sorted_bam + ' does not exist!')
        command2 = samtools_dir + ' sort ' + bam + ' > ' + sorted_bam
        os.system(command2)
    sorted_bam_index  = sorted_bam +  '.bai'
    if not os.path.isfile(sorted_bam_index):
        print(sorted_bam_index + ' does not exist!')
        command3 = samtools_dir + ' index ' + sorted_bam
        os.system(command3)
    #-numbers of reads in target region
    num_reads_in_target_region = out_dir + '/' + sample +'_'+ module + '_numbersReadsInTargetRegion.txt'
    command4 = samtools_dir + ' idxstats ' + sorted_bam + ' > ' + num_reads_in_target_region
    os.system(command4)
    #-coverage of reads in target region
    coverage_in_target_region = out_dir + '/' + sample +'_'+ module + '_covergerInTargetRegion.txt'
    command5 ='{0} mpileup {1} | perl -alne \'{2}\' > {3}'.format(samtools_dir, sorted_bam,
        '{$pos{$F[0]}++;$depth{$F[0]}+=$F[3]} END{print "$_\t$pos{$_}\t$depth{$_}" foreach sort keys %pos}',coverage_in_target_region)
    stdout, stderr = stdout_err(command5)
    logger_statistics_process.info(stdout)
    logger_statistics_errors.info(stderr)
    #-statistics and plot of  the depth and coverage in target region
    statistics_plot = out_dir + '/' + sample +'_'+ module + '_depth_coverageInTargetRegion'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    command6 = 'Rscript ' + scriptdir + '/statistics_depth_coverage.R' + ' -p ' + num_reads_in_target_region + ' -s ' + coverage_in_target_region + ' -r '+ exome_target_bed + ' -o ' + statistics_plot
    stdout, stderr = stdout_err(command6)
    logger_statistics_process.info(stdout)
    logger_statistics_errors.info(stderr)
    #depth of the bases in target region
    bases_depth_in_target_region = out_dir + '/' + sample +'_'+ module + '_basesDepthInTargetRegion.txt'
    command7 ='{0} mpileup {1} | perl -alne \'{2}\' > {3}'.format(samtools_dir, sorted_bam,
        '{$depth{$F[3]}++}END{print "$_\t$depth{$_}" foreach sort{$a <=> $b}keys %depth}',bases_depth_in_target_region)
    stdout, stderr = stdout_err(command7)
    logger_statistics_process.info(stdout)
    logger_statistics_errors.info(stderr)
    #-statistics and plot of  the depth and coverage in target region
    statistics_plot1 = out_dir + '/' + sample +'_'+ module + '_basesDepthInTargetRegion'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    command8 = 'Rscript ' + scriptdir + '/statistics_bases_depth.R' + ' -p ' + bases_depth_in_target_region + ' -o ' + statistics_plot1
    stdout, stderr = stdout_err(command8)
    logger_statistics_process.info(stdout)
    logger_statistics_errors.info(stderr)
    return sorted_bam
    
#--get the mapping result 
def statistics_sam_bam(samtools_dir, sam_bam, out_dir, sample, module, logger_statistics_process, logger_statistics_errors):
    #-statistics of
    if not os.path.isfile(sam_bam):
        logger_statistics_errors.error("%s does not exist!\n", sam_bam)
        print(sam_bam + ' does not exist!')
    align_statistics = out_dir + '/' + sample +'_'+ module +  '_statistics.txt'
    command = samtools_dir + ' stats ' + sam_bam +' | grep ^SN | cut -f 2-3  > ' + align_statistics
    os.system(command)
    return align_statistics

#--statistics of the time cost by the pipeline
def statistics_time(out_dir, sample, process, logger_statistics_process, logger_statistics_errors):
    if not os.path.isfile(process):
        logger_statistics_errors.error("%s does not exist!\n", process)
        print(process + ' does not exist!')
    time_statistics = out_dir + '/' + sample +'_'+ 'time_Cost'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    command = 'Rscript ' + scriptdir + '/statistics_time.R' + ' -p ' + process + ' -o ' + time_statistics
    os.system(command)

#--merge the statistics of sam built by modules in pipelines
def merge_statistics_sam_bam(logger_statistics_process, logger_statistics_errors, out_dir, sample, names,*args):
    for arg in args:
        if not os.path.isfile(arg):
            logger_statistics_errors.error("%s does not exist!\n", arg)
            print(arg + ' does not exist!')
    statisticsfiles = ','.join(args)
    mergestatistics = out_dir + '/' + sample +'_merge_sam_bam_statisticsfile.txt'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    command = 'Rscript ' + scriptdir + '/statistics_merge_sam_bam_statisticsfile.R' + ' -p ' + statisticsfiles + ' -g '+ names +' -o ' + mergestatistics
    os.system(command)

