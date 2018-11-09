__version__ = '1.0'
__author__ = 'Wang Xian'

import logging
import os
import re
import sys
import time
import itertools
from itertools import groupby
import shlex
import subprocess
sys.path.append("..")
from pipelines.log.log_v1 import store_statistics_logs


# put the info output to the log
def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr


# get the range of the positions of bases in the reads
def split_n_bases(nbases):
    if nbases == 'NULL':
        return ["NULL", "NULL"]
    else:
        bases = nbases.split(",")
        bases_posi = []
        for posi in range(0, len(bases)):
            posi = bases[posi]
            if "-" in str(posi):
                for i in range(int(posi.split("-")[0]), int(posi.split("-")[1])+1):
                    bases_posi.append(i)
            else:
                bases_posi.append(int(posi))
        l1 = bases_posi
        a = list()
        b = []
        result = []
        function1 = lambda x: x[1] - x[0]
        for k, g in groupby(enumerate(l1), function1):
            g = list(g)
            b.append(k)
            a.append(g)
        for c in range(len(b)):
            if [v for i, v in a[c]][0] != [v for i, v in a[c]][-1]:
                result.append("%d-%d" % ([v for i, v in a[c]][0], [v for i, v in a[c]][-1]))
            if [v for i, v in a[c]][0] == [v for i, v in a[c]][-1]:
                result.append([v for i, v in a[c]][0])
        if "-" not in str(result[len(result)-1]):
            result_last = result[len(result)-1]
        else:
            result_last = result[len(result)-1].split('-')[1]
        return ["[" + str(result[0]) + "]", str(result_last)]


# get informations from the qc result of the reads
def getinfo(fastqc, logger_statistics_process, logger_statistics_errors):
    qc_read = fastqc.split(".fastq")[0] + '_fastqc.zip'
    if not os.path.isfile(fastqc.split(".fastq")[0] + '_fastqc'):
        command1 = 'unzip' + ' -o ' + qc_read + ' -d ' + os.path.dirname(qc_read)
        stdout, stderr = stdout_err(command1)
        store_statistics_logs(logger_statistics_process, 'null', stdout)
        store_statistics_logs('null', logger_statistics_errors, stderr)
    # print('{0} has been completed.'.format(qc_read))
    store_statistics_logs(logger_statistics_process, 'null', '{0} has been completed.\n'.format(qc_read))
    qc_data = qc_read.rstrip(".zip") + '/' + 'fastqc_data.txt'
    qcdata = open(qc_data, "r")
    modules = 0
    #
    value = re.compile(r'\d+')
    perbasesequencequalit_posi = []
    perbasesequencequalit = []
    persequencequalityscores = []
    persequencequalityreads = []
    per_basen1 = []
    perbase_ncontent = []
    f = qcdata.readlines()
    for lines in range(0, len(f)):
        line = f[lines]
        line = line.strip()
        if '>>END_MODULE' in line:
            modules = modules + 1
        if 'Total Sequences' in line:
            raw_reads = re.findall('\d+', line)
        if 'Sequence length' in line:
            seq_length = line.split('\t')[1]
            if '-' in seq_length:
                seq_length_min, seq_length_max = seq_length.split('-')
            else:
                seq_length_min = str(seq_length)
                seq_length_max = str(seq_length)
        if '%GC' in line:
            gc = re.findall('\d+', line)
        if modules == 1 and value.match(line[0]):
            perbasesequencequalit_posi.append(line.split('\t')[0])
            perbasesequencequalit.append(float(line.split('\t')[1]))
        if modules == 3 and value.match(line[0]):
            persequencequalityscores.append(line.split('\t')[0])
            persequencequalityreads.append(float(line.split('\t')[1]))
        if modules == 6 and value.match(line[0]):
            per_basen1.append(line.split('\t')[0])
            perbase_ncontent.append(line.split('\t')[1])

    # --mean of Per base sequence quality
    all_perbasesequencequalit = 0
    for i in range(0, len(perbasesequencequalit)):
        if '-' not in perbasesequencequalit_posi[i]:
            qual = float(perbasesequencequalit[i])
            all_perbasesequencequalit += qual
        else:
            num1, num2 = perbasesequencequalit_posi[i].split('-')
            qual = float(perbasesequencequalit[i])
            all_perbasesequencequalit += qual*(int(num2)-int(num1)+1)
    perbasesequencequalit_mean = all_perbasesequencequalit/int(seq_length_max)

    # low qual base : qual < 25
    lowqualit_bases = []
    for i in range(0, len(perbasesequencequalit)):
        if float(perbasesequencequalit[i]) < 25:
            lowqualit_bases.append(perbasesequencequalit_posi[i])
    if len(lowqualit_bases) == 0:
        lowqualit_bases.append('NULL')
    lowqualit_bases_region = split_n_bases(','.join(lowqualit_bases))[0]

    # percent of Q20 and Q30
    persequencequalityreads1 = []
    for i in range(0, len(persequencequalityreads)):
        persequencequalityreads1.append(sum(persequencequalityreads[i:])/int(raw_reads[0]))
    for i in range(0, len(persequencequalityscores)):
        if persequencequalityscores[i] == '20':
            q20 = persequencequalityreads1[i]
        if persequencequalityscores[i] == '30':
            q30 = persequencequalityreads1[i]
    if '20' not in persequencequalityscores:
        q20 = 1.0000

    # N content
    nbases = []
    for i in range(0, len(perbase_ncontent)):
        if float(perbase_ncontent[i]) > 0:
            nbases.append(per_basen1[i])
    if len(nbases) == 0:
        nbases.append('NULL')
    nbase_region = split_n_bases(','.join(nbases))

    return [raw_reads[0], seq_length_min, seq_length_max, gc[0], str(round(perbasesequencequalit_mean, 3)),
            lowqualit_bases_region, str('%.3f%%' % (q20 * 100)),
            str('%.3f%%' % (q30 * 100)), nbase_region[0]]


# --QC reads by fastQC
def qc_raw_reads(fastQC_dir, out_dir, sample, module, read1, read2,
                 logger_statistics_process, logger_statistics_errors):
    qc_read = out_dir + '/' + os.path.basename(read1).split(".fastq")[0] + '_fastqc.zip'
    if not os.path.isfile(qc_read):
        command1 = '{0} {1} {2} -o {3}'.format(fastQC_dir, read1, read2, out_dir)
        stdout, stderr = stdout_err(command1)
        store_statistics_logs(logger_statistics_process, 'null', stdout)
        store_statistics_logs('null', logger_statistics_errors, stderr)
        store_statistics_logs(logger_statistics_process, 'null', 'QC-{0} has been completed.\n'.format(module))
    else:
        store_statistics_logs(logger_statistics_process, 'null', 'QC-{0} exists.\n'.format(module))
    qc_statistics = out_dir + '/' + sample + '.' + module + '.statistics.txt'
    qc_result1 = getinfo(out_dir + '/' + os.path.basename(read1), logger_statistics_process, logger_statistics_errors)
    qc_result2 = getinfo(out_dir + '/' + os.path.basename(read2), logger_statistics_process, logger_statistics_errors)
    fout = open(qc_statistics, 'w')
    fout.write('\t'.join(['SampleID', 'Sequence direction', 'raw reads', 'min length', 'max length',
                          'GC content', 'mean of Per base qualit', 'low qualit Bases position',
                          'Q20', 'Q30', 'N_bases position']) + '\n')
    fout.write('\t'.join([sample, read1.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result1[0:9])]) + '\n')
    fout.write('\t'.join([sample, read2.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result2[0:9])]) + '\n')
    fout.close()
    return qc_result1, qc_result2


# -get the depth and coverage of the mapping result
def statistics_depth_coverage(samtools_dir, sam_bam, out_dir, sample, module, exome_target, exome_target_bed,
                              logger_statistics_process, logger_statistics_errors, renew):
    # get the path
    scriptdir = os.path.dirname(os.path.abspath(__file__))

    if not os.path.isfile(sam_bam):
        store_statistics_logs('null', logger_statistics_errors, sam_bam + " does not exist!\n")
    sorted_bam = sam_bam.rstrip('.sam') + '_sorted.bam'
    bam = sam_bam.rstrip('.sam') + '.bam'
    if not os.path.isfile(sorted_bam) or renew is 'T':
        # print(sorted_bam + ' does not exist!')
        if not os.path.isfile(bam):
            bam = sam_bam.rstrip('.sam') + '.bam'
            command1 = samtools_dir + ' view -bS ' + sam_bam + ' -o ' + bam
            stdout, stderr = stdout_err(command1)
            store_statistics_logs(logger_statistics_process, 'null', stdout)
            store_statistics_logs('null', logger_statistics_errors, stderr)
            store_statistics_logs(logger_statistics_process, 'null',
                                  '{0} has been tranformed to bam.\n'.format(sam_bam))
            command2 = samtools_dir + ' sort ' + bam + ' -o ' + sorted_bam
            stdout, stderr = stdout_err(command2)
            store_statistics_logs(logger_statistics_process, 'null', stdout)
            store_statistics_logs('null', logger_statistics_errors, stderr)
            os.system('rm -rf {0}'.format(bam))
        #else:
        #    command2 = samtools_dir + ' sort ' + bam + ' -o ' + sorted_bam
        #    stdout, stderr = stdout_err(command2)
        #    store_statistics_logs(logger_statistics_process, 'null', stdout)
        #    store_statistics_logs('null', logger_statistics_errors, stderr)
    sorted_bam_index = sorted_bam + '.bai'
    if not os.path.isfile(sorted_bam_index) or renew is 'T':
        # print(sorted_bam_index + ' does not exist!')
        command3 = samtools_dir + ' index ' + sorted_bam
        stdout, stderr = stdout_err(command3)
        store_statistics_logs(logger_statistics_process, 'null', stdout)
        store_statistics_logs('null', logger_statistics_errors, stderr)
    # -numbers of reads in target region
    num_reads_in_target_region = out_dir + '/' + sample + '_' + module + '_numbersReadsInTargetRegion.txt'
    if not os.path.isfile(num_reads_in_target_region) or renew is 'T':
        command4 = samtools_dir + ' idxstats ' + sorted_bam + ' -o ' + num_reads_in_target_region
        stdout, stderr = stdout_err(command4)
        store_statistics_logs(logger_statistics_process, 'null', stdout)
        store_statistics_logs('null', logger_statistics_errors, stderr)
    # -coverage of reads in target region
    coverage_in_target_region = out_dir + '/' + sample + '_' + module + '_covergerInTargetRegion.txt'
    if not os.path.isfile(coverage_in_target_region) or renew is 'T':
        command5 = '{0} mpileup {1} | perl -alne \'{2}\' > {3}'.format(
            samtools_dir, sorted_bam,
            '{$pos{$F[0]}++;$depth{$F[0]}+=$F[3]} END{print "$_\t$pos{$_}\t$depth{$_}" foreach sort keys %pos}',
            coverage_in_target_region)
        os.system(command5)
    # -
    # -statistics and plot of  the depth and coverage in target region
    statistics_plot = out_dir + '/' + sample + '_' + module + '_depth_coverageInTargetRegion'
    if not os.path.isfile(statistics_plot + '.pdf') or renew is 'T':
        # scriptdir = os.path.dirname(os.path.abspath(__file__))
        command6 = 'Rscript ' + scriptdir + '/statistics_depth_coverage.R' + ' -p ' \
                   + num_reads_in_target_region + ' -s ' + coverage_in_target_region\
                   + ' -r ' + exome_target_bed + ' -o ' + statistics_plot
        stdout, stderr = stdout_err(command6)
        store_statistics_logs(logger_statistics_process, 'null', stdout)
        store_statistics_logs('null', logger_statistics_errors, stderr)
    # - statistics of the depth of the baes in region
    if module == 'Fliter':
        bases_depth_in_region = out_dir + '/' + sample + '_' + module + '_basesDepthInRegion.txt'
        if not os.path.isfile(bases_depth_in_region) or renew is 'T':
            command7 = '{0} depth {1} > {2}'.format(samtools_dir, sorted_bam, bases_depth_in_region)
            os.system(command7)
        # statistics of the depth of exon region
        statistics_plot1 = out_dir + '/' + sample + '_' + module + '_basesDepthInTargetExon'
        command7 = ''.join(['Rscript ',  scriptdir, '/statistics_bases_coverage_on_target_exon.R', ' -p ',
                           bases_depth_in_region, ' -s ', exome_target, ' -t True ',
                           ' -o ', statistics_plot1])
        os.system(command7)
    # depth of the bases in target region
    bases_depth_in_target_region = out_dir + '/' + sample + '_' + module + '_basesDepthInTargetRegion.txt'
    if not os.path.isfile(bases_depth_in_target_region) or renew is 'T':
        command7 = '{0} mpileup {1} | perl -alne \'{2}\' > {3}'.format(
            samtools_dir, sorted_bam, '{$depth{$F[3]}++}END{print "$_\t$depth{$_}" foreach sort{$a <=> $b}keys %depth}',
            bases_depth_in_target_region)
        os.system(command7)
    # -statistics and plot of the depth and coverage in target region
    statistics_plot2 = out_dir + '/' + sample + '_' + module + '_basesDepthInTargetRegion'
    command8 = 'Rscript ' + scriptdir + '/statistics_bases_depth.R' + ' -p ' \
               + bases_depth_in_target_region + ' -o ' + statistics_plot2
    stdout, stderr = stdout_err(command8)
    store_statistics_logs(logger_statistics_process, 'null', stdout)
    store_statistics_logs('null', logger_statistics_errors, stderr)
    return sorted_bam



# -get the depth and coverage of MT
def statistics_mtdepth_coverage(germline_vc_dir, out_dir, sample, exome_target,
                                logger_statistics_process, logger_statistics_errors):
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    mtdp = germline_vc_dir + '/' + sample + ".smCounter.all.txt"
    if os.path.isfile(mtdp):
        statistics_plot1 = out_dir + '/' + sample + '_bases_MTDepthInTargetExon'
        command7 = ''.join(['Rscript ',  scriptdir, '/statistics_basesMT_coverage_on_target_exon.R', ' -p ',
                           mtdp, ' -s ', exome_target, ' -t True ',
                           ' -o ', statistics_plot1])
        os.system(command7)

# --get the mapping result
def statistics_sam_bam(samtools_dir, sam_bam, out_dir, sample, module,
                       logger_statistics_process, logger_statistics_errors, renew):
    # -statistics of
    if not os.path.isfile(sam_bam):
        store_statistics_logs('null', logger_statistics_errors, sam_bam + " does not exist!\n")
        # print(sam_bam + ' does not exist!')
    align_statistics = out_dir + '/' + sample + '_' + module + '_statistics.txt'
    if not os.path.isfile(align_statistics) or renew is 'T':
        command = samtools_dir + ' stats ' + sam_bam + ' | grep ^SN | cut -f 2-3  > ' + align_statistics
        os.system(command)
    return align_statistics


# --statistics of the time cost by the pipeline
def statistics_time(out_dir, sample, process, logger_statistics_process, logger_statistics_errors):
    if not os.path.isfile(process):
        store_statistics_logs('null', logger_statistics_errors, process + " does not exist!\n")
        # print(process + ' does not exist!')
    time_statistics = out_dir + '/' + sample + '_' + 'time_Cost'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    command = 'Rscript ' + scriptdir + '/statistics_time.R' + ' -p ' + process + ' -o ' + time_statistics
    os.system(command)


# --merge the statistics of sam built by modules in pipelines
def merge_statistics_sam_bam(logger_statistics_process, logger_statistics_errors, out_dir, sample, names, *args):
    for arg in args:
        if not os.path.isfile(arg):
            store_statistics_logs('null', logger_statistics_errors, arg+" does not exist!\n")
            # print(arg + ' does not exist!')
    statisticsfiles = ','.join(args)
    mergestatistics = out_dir + '/' + sample + '_merge_sam_bam_statisticsfile.txt'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    command = 'Rscript ' + scriptdir + '/statistics_merge_sam_bam_statisticsfile.R' + \
              ' -p ' + statisticsfiles + ' -g ' + names + ' -o ' + mergestatistics
    os.system(command)


def merge_statistics(sample_preinfo, exome_target, qc, trim_qc,
                     trim_statis, filter_statis, primer_statis, umi_statis, align_base, mt_depth):
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    command = 'Rscript {0}/statistics_preinfo_v1.R -q {1} -t {2} -c {3} -f {4} -p {5} ' \
              '-u {6}  -a {7} -m {8} -e {9} -o {10}'.format(scriptdir, qc, trim_qc, trim_statis, filter_statis,
                                                   primer_statis, umi_statis, align_base, mt_depth, exome_target, sample_preinfo)
    os.system(command)
