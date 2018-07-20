__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
import re
import sys
import time
import gzip
import itertools
from itertools import groupby
import argparse

def split_N_bases(Nbases):
    bases = Nbases.split(",")
    print(bases)
    bases_posi =[]
    for posi in range(0,len(bases)):
        posi = bases[posi]
        print(posi)
        if "-" in str(posi):
            for i in range(int(posi.split("-")[0]), int(posi.split("-")[1])+1):
                bases_posi.append(i)
        else:
            bases_posi.append(int(posi))
    print(bases_posi)
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
    print(result)
    if "-" not in str(result[len(result)-1]):
        result_last = result[len(result)-1]
    else:
        result_last = result[len(result)-1].split('-')[1]
    return result_last

def getinfo(fastqc):
    qc_read = fastqc.split(".fastq")[0] + '_fastqc.zip'
    commond1 = 'unzip' + ' -o ' + qc_read  + ' -d ' + os.path.dirname(qc_read)
    os.system(commond1)
    print('{0} has been completed.'.format(qc_read))
    qc_data = qc_read.rstrip(".zip") + '/' + 'fastqc_data.txt'
    qcdata = open(qc_data,"r")
    modules = 0
    # 定义正则表达式
    value = re.compile(r'\d+')
    Perbasesequencequalit_posi = []
    Perbasesequencequalit = []
    Persequencequalityscores = []
    Persequencequalityreads = []
    per_basen1 = []
    PerbaseNcontent = []
    f = qcdata.readlines()
    for lines in range(0 ,len(f)):
        line = f[lines]
        line = line.strip()
        if '>>END_MODULE' in line:
            modules = modules + 1
        if 'Total Sequences' in line:
            raw_reads = re.findall('\d+',line)
        if 'Sequence length' in line:
            mode = re.compile(r'\d+-\d+')
            seq_length = mode.findall(line)
        if '%GC' in line:
            GC = re.findall('\d+',line)
        if modules == 1 and value.match(line[0]):
            Perbasesequencequalit_posi.append(line.split('\t')[0])
            Perbasesequencequalit.append(float(line.split('\t')[1]))
        if modules == 3 and value.match(line[0]):
            Persequencequalityscores.append(line.split('\t')[0])
            Persequencequalityreads.append(float(line.split('\t')[1]))
        if modules == 6 and value.match(line[0]):
            per_basen1.append(line.split('\t')[0])
            PerbaseNcontent.append(line.split('\t')[1])
    #--Per base sequence quality
    Perbasesequencequalit_mean = sum(Perbasesequencequalit)/len(Perbasesequencequalit)
    lowqualitBases= []
    for i in range(0,len(Perbasesequencequalit)):
        if float(Perbasesequencequalit[i]) < Perbasesequencequalit_mean-2:
            lowqualitBases.append(Perbasesequencequalit_posi[i])
    if len(lowqualitBases) == 0:
        lowqualitBases.append('NULL')
    #--Per sequence quality
    Persequencequalityreads1 = []
    for i in range(0,len(Persequencequalityreads)):
        Persequencequalityreads1.append(sum(Persequencequalityreads[i:])/int(raw_reads[0]))
    for i in range(0,len(Persequencequalityscores)):
        if Persequencequalityscores[i] == '20':
            Q20 = Persequencequalityreads1[i]
        if Persequencequalityscores[i] == '30':
            Q30 = Persequencequalityreads1[i]
    if '20' not in Persequencequalityscores:
        Q20 = 1.0000
    #--N content
    Nbases = []
    for i in range(0,len(PerbaseNcontent)):
        if float(PerbaseNcontent[i]) > 0:
            Nbases.append(per_basen1[i])
    if len(Nbases) == 0:
        Nbases.append('NULL')
    #---
    #print([raw_reads[0], seq_length[0], GC[0], str(round(Perbasesequencequalit_mean,3)), ','.join(lowqualitBases), str('%.3f%%' % (Q20 * 100)), str('%.3f%%' % (Q30 * 100)), ','.join(Nbases)]
    return [raw_reads[0], seq_length[0], GC[0], str(round(Perbasesequencequalit_mean,3)), ','.join(lowqualitBases), str('%.3f%%' % (Q20 * 100)), str('%.3f%%' % (Q30 * 100)), ','.join(Nbases)]

def qc_raw_reads(fastQC_dir, out_dir, sample, module, read1, read2,logger_static_process, logger_static_errors):
    commond1 = '{0} {1} {2} -o {3}'.format(fastQC_dir, read1, read2, out_dir)
    os.system(commond1)
    print('QC-{0} has been completed.'.format(module))
    qc_static = out_dir + '/' + sample + '.' + module +'.static.txt'
    qc_result1 =getinfo(out_dir + '/' + os.path.basename(read1))
    qc_result2 =getinfo(out_dir + '/' + os.path.basename(read2))
    fout = open(qc_static, 'w')
    fout.write('\t'.join(['SampleID','Sequence direction','raw reads', 'seq length', 'GC content', 'mean of Per base qualit', 'lowqualitBases', 'Q20', 'Q30', 'N_bases']) + '\n')
    fout.write('\t'.join([sample, read1.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result1)]) + '\n')
    fout.write('\t'.join([sample, read2.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result2)]) + '\n')
    fout.close()
    return qc_result1, qc_result2

def static_depth_coverage(samtools_dir, sam_bam, out_dir,sample, module, logger_static_process, logger_static_errors):
    if not os.path.isfile(sam_bam):
        logger_static_errors.error("%s does not exist!\n", sam_bam)
        print(sam_bam + ' does not exist!')
    if re.search('.sam$', sam_bam):
        bam = sam_bam.rstrip('.sam') + '.bam'
        commond1 = samtools_dir + ' view -bS ' + sam_bam + ' > ' + bam
        os.system(commond1)
        logger_static_process.info('{0} has been tranformed to bam.'.format(sam_bam))
        sorted_bam = sam_bam.rstrip('.sam') + '_sorted.bam'
    elif re.search('.bam$', sam_bam):
        bam = sam_bam
        sorted_bam = sam_bam.rstrip('.bam') + '_sorted.bam'
    if not os.path.isfile(sorted_bam):
        print(sorted_bam + ' does not exist!')
        commond2 = samtools_dir + ' sort ' + bam + ' > ' + sorted_bam
        os.system(commond2)
    sorted_bam_index  = sorted_bam +  '.bai'
    if not os.path.isfile(sorted_bam_index):
        print(sorted_bam_index + ' does not exist!')
        commond3 = samtools_dir + ' index ' + sorted_bam
        os.system(commond3)
    #-numbers of reads in target region
    numReadsInTargetRegion = out_dir + '/' + sample +'_'+ module + '_numbersReadsInTargetRegion.txt'
    commond4 = samtools_dir + ' idxstats ' + sorted_bam + ' > ' + numReadsInTargetRegion
    os.system(commond4)
    #-coverage of reads in target region
    coverageInTargetRegion = out_dir + '/' + sample +'_'+ module + '_covergerInTargetRegion.txt'
    commond5 ='{0} mpileup {1} | perl -alne \'{2}\' > {3}'.format(samtools_dir, sorted_bam,
        '{$pos{$F[0]}++;$depth{$F[0]}+=$F[3]} END{print "$_\t$pos{$_}\t$depth{$_}" foreach sort keys %pos}',coverageInTargetRegion)
    os.system(commond5)
    #-static and plot of  the depth and coverage in target region
    static_plot = out_dir + '/' + sample +'_'+ module + '_depth_coverageInTargetRegion'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond6 = 'Rscript ' + scriptdir + '/static_depth_coverage.R' + ' -p ' + numReadsInTargetRegion + ' -s ' + coverageInTargetRegion + ' -o ' + static_plot
    os.system(commond6)
    #depth of the bases in target region
    basesDepthInTargetRegion = out_dir + '/' + sample +'_'+ module + '_basesDepthInTargetRegion.txt'
    commond7 ='{0} mpileup {1} | perl -alne \'{2}\' > {3}'.format(samtools_dir, sorted_bam,
        '{$depth{$F[3]}++}END{print "$_\t$depth{$_}" foreach sort{$a <=> $b}keys %depth}',basesDepthInTargetRegion)
    os.system(commond7)
    #-static and plot of  the depth and coverage in target region
    static_plot1 = out_dir + '/' + sample +'_'+ module + '_basesDepthInTargetRegion'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond8 = 'Rscript ' + scriptdir + '/static_bases_depth.R' + ' -p ' + basesDepthInTargetRegion + ' -o ' + static_plot1
    os.system(commond8)
    return sorted_bam
    
def static_sam_bam(samtools_dir, sam_bam, out_dir, sample, module, logger_static_process, logger_static_errors):
    #-static of
    if not os.path.isfile(sam_bam):
        logger_static_errors.error("%s does not exist!\n", sam_bam)
        print(sam_bam + ' does not exist!')
    align_static = out_dir + '/' + sample +'_'+ module +  '_Static.txt'
    commond = samtools_dir + ' stats ' + sam_bam +' | grep ^SN | cut -f 2-3  > ' + align_static
    os.system(commond)
    return align_static

def static_time(out_dir, sample, process, logger_static_process, logger_static_errors):
    if not os.path.isfile(process):
        logger_static_errors.error("%s does not exist!\n", process)
        print(process + ' does not exist!')
    time_static = out_dir + '/' + sample +'_'+ 'time_Cost'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond = 'Rscript ' + scriptdir + '/static_time.R' + ' -p ' + process + ' -o ' + time_static
    os.system(commond)

def merge_static_sam_bam(logger_static_process, logger_static_errors, out_dir, sample, names,*args):
    for arg in args:
        if not os.path.isfile(arg):
            logger_static_errors.error("%s does not exist!\n", arg)
            print(arg + ' does not exist!')
    staticfiles = ','.join(args)
    mergestatic = out_dir + '/' + sample +'_merge_sam_bam_staticfile.txt'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond = 'Rscript ' + scriptdir + '/static_merge_sam_bam_staticfile.R' + ' -p ' + staticfiles + ' -g '+ names +' -o ' + mergestatic
    os.system(commond)

