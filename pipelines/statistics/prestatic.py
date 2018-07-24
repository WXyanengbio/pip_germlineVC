__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
import re
import sys
import time
import gzip
import itertools
import argparse

#sys.path.append(os.path.split(os.path.realpath(__file__))[0])
#import the logging functions
#from pipelines.log.log import store_static_logs


def script_information():
    print ("\nApplication: pipelines of QIAseq Targeted DNA Panel\n")
    print ("=====================================================================")
    print ("Required environment: python \ fastQC \ samtools \ GATK")

parser = argparse.ArgumentParser(usage = "\n\npython %(prog)s --source --sample_name --tailname --output \
                                          --samtools_dir --fastQC_dir --gatk_dir")
parser.add_argument("--source", help = "Path to input reads in FASTA format", type = str)
parser.add_argument("--sample_name", help = "the sample name of raw reads", type = str)
parser.add_argument("--tailname", help = "the tailname of sample raw reads", type =str)
parser.add_argument("--output", help = "Path of output file", type = str)
#parser.add_argument("--ref_index_name", help = "the path of ref index--if there isn't a ref index, it will make a index in the path of ref fasta by bwa", type = str)
#parser.add_argument("--ref_fa_file", help = "the path of ref fasta", type = str)
#parser.add_argument("--total_ref_fa_file", help = "the path of ref total ref fa", type = str)
#parser.add_argument("--total_ref_fa_dict", help = "the path of ref total ref fa dict", type = str)
parser.add_argument("--samtools_dir", help = "the install path of samtools", type = str)
#parser.add_argument("--primers_file", help = "Load all primer sequences in the panel", type = str)
parser.add_argument("--fastQC_dir", help = "the install path of fastQC", type = str)
parser.add_argument("--gatk_dir", help = "the install path of GATK4", type = str)
#parser.add_argument("--exome_target_bed", help = "the bed file of exome intervals", type = str) 
parser.add_argument("-v", '-version', action = 'version', version =' %(prog)s 1.0')
parser.add_argument("--test", help = "the subprocess of the script", type = int, default = 1)

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    script_information()
    args = parser.parse_args()
if "-v" in sys.argv[1:]:
    args = parser.parse_args()
try:
    args = parser.parse_args()
except SystemExit:
    script_information()
    parser.print_help()
    exit()

if len(sys.argv) == 1:
    script_information()
    parser.print_help()
    exit()

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger

def store_static_logs(log_dir):
    formatter_pipeline_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_pipeline_errors = logging.Formatter(
        "%(asctime)s;%(levelname)s;                                             %(message)s")
    logger_static_process = setup_logger('Running Messages of static',
                                       log_dir + '/static_process.log',
                                       formatter_pipeline_process)
    logger_static_errors = setup_logger('Errors & Warnings of pipeline',
                                      log_dir + '/static_errors.log',
                                      formatter_pipeline_errors)
    return  logger_static_process, logger_static_errors

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
    #--N content
    Nbases = []
    for i in range(0,len(PerbaseNcontent)):
        if float(PerbaseNcontent[i]) > 0:
            Nbases.append(per_basen1[i])
    if len(Nbases) == 0:
        Nbases.append('NULL')
    #---
    print([raw_reads[0], seq_length[0], GC[0], str(round(Perbasesequencequalit_mean,3)), ','.join(lowqualitBases), str('%.3f%%' % (Q20 * 100)), str('%.3f%%' % (Q30 * 100)), ','.join(Nbases)]
    return [raw_reads[0], seq_length[0], GC[0], str(round(Perbasesequencequalit_mean,3)), ','.join(lowqualitBases), str('%.3f%%' % (Q20 * 100)), str('%.3f%%' % (Q30 * 100)), ','.join(Nbases)]

def qc_raw_reads(fastQC_dir, out_dir, sample, module, read1, read2,logger_trim_process, logger_static_errors):
    commond1 = '{0} {1} {2} -o {3}'.format(fastQC_dir, read1, read2, out_dir)
    os.system(commond1)
    print('QC has been completed.')
    qc_static = out_dir + '/' + sample + '.' + module +'.static.txt'
    qc_result1 =getinfo(out_dir + '/' + os.path.filename(read1))
    qc_result2 =getinfo(out_dir + '/' + os.path.filename(read2))
    fout = open(qc_static, 'w')
    fout.write('\t'.join(['SampleID','Sequence direction','raw reads', 'seq length', 'GC content', 'mean of Per base qualit', 'lowqualitBases', 'Q20', 'Q30', 'N_bases']) + '\n')
    fout.write('\t'.join([sample, read1.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result1)]) + '\n')
    fout.write('\t'.join([sample, read2.lstrip(sample + '_').rstrip(".fastq.gz"), '\t'.join(qc_result2)]) + '\n')
    fout.close()
    return qc_result1,qc_result2
    
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
    commond2 = samtools_dir + ' sort ' + bam + ' > ' + sorted_bam
    os.system(commond2)
    commond3 = samtools_dir + ' index ' + sorted_bam
    os.system(commond3)
    #-numbers of reads in target region
    numReadsInTargetRegion = out_dir + '/' + module + '_numbersReadsInTargetRegion.txt'
    commond4 = samtools_dir + ' idxstats ' + sorted_bam + ' > ' + numReadsInTargetRegion
    os.system(commond4)
    #-coverage of reads in target region
    coverageInTargetRegion = out_dir + '/' + module + '_covergerInTargetRegion.txt'
    commond5 = samtools_dir + ' mpileup ' + sorted_bam +' | perl -alne '{$pos{$F[0]}++;$depth{$F[0]}+=$F[3]} END{print "$_\t$pos{$_}\t$depth{$_}" foreach sort keys %pos}' >' + coverageInTargetRegion
    os.system(commond5)
    #-static and plot of  the depth and coverage in target region
    static_plot = out_dir + '/' + module + '_depth_coverageInTargetRegion'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond6 = 'Rscript ' + scriptdir + '/static_depth_coverage.R' + ' -p ' + numReadsInTargetRegion + ' -s ' + coverageInTargetRegion + ' -o ' + static_plot
    os.system(commond6)
    #depth of the bases in target region
    basesDepthInTargetRegion = out_dir + '/' + module + '_basesDepthInTargetRegion.txt'
    commod7 = samtools_dir + ' mpileup ' + sorted_bam +' | perl -alne '{$depth{$F[3]}++}END{print "$_\t$depth{$_}" foreach sort{$a <=> $b}keys %depth}' >' + basesDepthInTargetRegion
    os.system(commond7)
    #-static and plot of  the depth and coverage in target region
    static_plot1 = out_dir + '/' + module + '_basesDepthInTargetRegion'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond8 = 'Rscript ' + scriptdir + '/static_bases_depth.R' + ' -p ' + basesDepthInTargetRegion + ' -o ' + static_plot1
    os.system(commond8)
    
def static_sam_bam(samtools_dir, sam_bam, out_dir,sample, module, logger_static_process, logger_static_errors):
    #-static of
    if not os.path.isfile(sam_bam):
        logger_static_errors.error("%s does not exist!\n", sam_bam)
        print(sam_bam + ' does not exist!')
    align_static = out_dir + '/' + module +  '_Static.txt'
    commond = samtools_dir + ' stats ' + sam_bam +' | grep ^SN | cut -f 2-3  > ' align_static
    os.system(commond)

def static_time(out_dir, process, logger_static_process, logger_static_errors):
    if not os.path.isfile(sam_bam):
        logger_static_errors.error("%s does not exist!\n", sam_bam)
        print(sam_bam + ' does not exist!')
    time_static = out_dir + '/' + 'time_Cost'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond = 'Rscript ' + scriptdir + '/static_time.R' + ' -p ' + process + ' -o ' + time_static
    os.system(commond)

def merge_static_sam_bam(logger_static_process, logger_static_errors,*args):
    for arg in args:
        if not os.path.isfile(arg):
            logger_static_errors.error("%s does not exist!\n", arg)
            print(arg + ' does not exist!')
    staticfiles = ','.join(args)
    mergestatic = out_dir + '/' + '_merge_sam_bam_staticfile.txt'
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    commond = 'Rscript ' + scriptdir + '/static_merge_sam_bam_staticfile.R' + ' -p ' + staticfiles + ' -o ' + mergestatic
    os.system(commond)
    
def main():
    #time cost
    time_start1 = time.time()
    #---check the outputdir
    out_dir = args.output
    #----pipeline log file
    log_dir = out_dir + '/' + 'log/'
    logger_static_process, logger_static_errors = store_static_logs(log_dir)
    #---input
    source = args.source
    sample = args.sample_name
    tailname = args.tailname
    fastQC_dir = args.fastQC_dir
    #trim
    #undetermined_dir = source + '/'+ 'undetermined'
    #stats_file = undetermined_dir + '/' + sample + '_basic_stats.txt'
    #bwa---align
    #num_threads = args.num_threads
    #ref_fa_file = args.ref_fa_file
    #ref_index_name = args.ref_index_name
    #post-align
    samtools_dir = args.samtools_dir
    #primers_file = args.primers_file
    #--clustering
   #--variant calling
    #memorySize = args.memorySize
    #memorySize = '-Xmx' + memorySize + ' ' + '-Djava.io.tmpdir=./'
    #gatk_dir = args.gatk_dir
    samtools_dir = args.samtools_dir
    #total_ref_fa_dict = args.total_ref_fa_dict
    #total_ref_fa_file = args.total_ref_fa_file
    #read_length = args.read_length
    #exome_target_bed = args.exome_target_bed
    test_level = args.test
    ##########################################################################################
    #---QC
    ##########################################################################################
    #static_dir
    static_dir = out_dir + '/'+ 'static'
    if not os.path.exists(static_dir):
        os.makedirs(static_dir)
    sample = sample + '_' + tailname
    read1 = source + '/' + sample + '_R1_001.fastq.gz'
    read2 = source + '/' + sample  + '_R2_001.fastq.gz'
    module = 'qc'
    qc_result1 , qc_result2 = qc_raw_reads(fastQC_dir, static_dir, sample, module, read1, read2)
    
    logger_static_process.info("QC of reads is completed after %.2f min.", (time.time()-time_start1)/60)
    
    if test_level == 1:
        print("Test QC module!")
        exit()
    ##########################################################################################
    #---trim
    #########################################################################################
    
if __name__ == '__main__':
    main()