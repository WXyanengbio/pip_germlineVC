__author__ = 'Wang Xian'
__version__ = 'v1.0'

import os
import re
import sys
import time
import gzip
import itertools
import argparse

time_start = time.time()
def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def store_logs(log_dir):
    formatter_g_variantCalling_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_g_variantCalling_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_g_variantCalling_process = setup_logger('Germline variant calling Messages',
                                         log_dir + '/Germline_variant_calling_process.log',
                                         formatter_g_variantCalling_process)
    logger_g_variantCalling_errors = setup_logger('Errors & Warnings of Germline_variant_calling',
                                        log_dir + '/Germline_variant_calling_errors.log',
                                        formatter_fg_variantCalling_errors)
    return logger_g_variantCalling_process, logger_g_variantCalling_errors


def script_information():
    print ("\nApplication: Germline variant calling of QIAseq Targeted DNA Panel\n")
    print ("=====================================================================")
    print ("Required environment：python \  samtools \ GATK4")

def sam_to_bem(gatk_dir, samtools_dir,
                            sam, sample, 
                            output, memorySize, 
                            exome_target_bed, known_sites, 
                            read_length, 
                            logger_g_variantCalling_process, logger_g_variantCalling_errors):
    # sam to bam By SortSam in GATK
    bam = output + '/' + sample + '.bam'
    command_count = '{0} --java-options "{1}" SortSam -SO coordinate -I {2} -O {3}'.format(
        gatk_dir, memorySize, sam, bam)
    logger_g_variantCalling_process.info('GATK sort sam to bam.')
    os.system(command_count)
    #statistics of coverages
    cov_file =  output + '/' + sample + '.cov.txt'
    command_count_1 = '{0} CollectHsMetrics -BI {1} -TI {2} -I {3} -O {4}'.format(
        gatk_dir, exome_target_bed, exome_target_bed, bam, cov_file)
    logger_g_variantCalling_process.info('GATK count the coverage.')
    os.system(command_count_1)
    #build index of bam By samtools
    command_count1 = '{0} index {1}'.format(samtools_dir, bam)
    logger_g_variantCalling_process.info('Samtools build the index of bam.')
    os.system(command_count1)
    #MarkDuplicates
    mark_bam = output + '/' + sample + '_marked.bam'
    bam_metrics = output + '/' + sample + '.metrics'
    command_count2 ='{0} --java-options "{1}" MarkDuplicates -I {2} -O {3} -M {4}'.format(
        gatk_dir, memorySize, bam, mark_bam, bam_metrics)
    logger_g_variantCalling_process.info('GATK marks the duplicates.')
    os.system(command_count2)
    #build the index of the marked_fixed bam
    command_count4 = '{0} index {1}'.format(samtools_dir, mark_bam)
    logger_g_variantCalling_process.info('Samtools build the index of  marked bam.')
    os.system(command_count4)
    #---------------------------------------
    #Base(Quality Score) Recalibration
    #BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model
    recal_data_table =  output + '/' + sample + '.recal_data.table'
    command_count5 ='{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} -L {4} -ip {5} --known-sites {6} -O {7}'.format(
        gatk_dir, memorySize, ref_fa_file, mark_bam, exome_target_bed, read_length, known_sites, recal_data_table)
    logger_g_variantCalling_process.info('GATK generate Base Quality Score Recalibration (BQSR) model.')
    os.system(command_count5)
    #GatherBQSRReports--Generate report of Base Quality Score Recalibration (BQSR) model
    BQSRReports = output + '/' + sample + '.BQSR.report'
    command_count6 ='{0} --java-options "{1}" GatherBQSRReports -I {2} -O {3}'.format(
        gatk_dir, memorySize,  recal_data_table, BQSRReports)
    logger_g_variantCalling_process.info('GATK generate report of Base Quality Score Recalibration (BQSR) model.')
    os.system(command_count6)
    #ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model
    bgsr_bam =  output + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
    command_count7 ='{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4}  -L {5} -ip {6} -O {7}'.format(
        gatk_dir, memorySize,  ref_fa_file, mark_bam, recal_data_table, exome_target_bed, read_length, bgsr_bam)
    logger_g_variantCalling_process.info('GATK apply Base Quality Score Recalibration (BQSR) model.')
    os.system(command_count7)
    #time used for translating the format of the sam by samtools and GATK
    logger_g_variantCalling_process.info('Time cost at  translating the format of the sam by samtools and GATK == ' + str((time.time() - time_start) / 60) + 'min')
    print('Compeleted translating the format of the sam by samtools and GATK .')

def germline_variant_calling(gatk_dir, marked_BQSR_bam,
                            sample, output, 
                            memorySize, ref_fa_file,
                            exome_target_bed, ERC,
                            read_filter, read_length,
                            reduce_logs, create_output_variant_index,
                            logger_g_variantCalling_process, logger_g_variantCalling_errors)
    command_count ='{0} --java-options "{1}" HaplotypeCaller  -R {2} -I {3} -L {4} -ip {5}'.format(
        gatk_dir, memorySize, ref_fa_file, marked_BQSR_bam, exome_target_bed, read_length)
    logger_g_variantCalling_process.info('Begin to confirm the options parameters of running HaplotypeCaller.')
    if ERC == NULL:
        logger_g_variantCalling_process.info('Runs HaplotypeCaller in default mode on a single input BAM file containing sequence data!')
        vcf =  output + '/' + sample + '_variants.vcf'
        command_count + = vcf
    else:
        if ERC == 'GVCF' 
            logger_g_variantCalling_process.info('Runs HaplotypeCaller in GVCF mode!')
            vcf =  output + '/' + sample + '_variants.g.vcf'
            command_count + = vcf + '-ERC {0}'.format(ERC)
        else:
            logger_g_variantCalling_process.info('{0} is not a HaplotypeCaller model. Please check the input parameter of ERC!'.format(ERC))

    if read_filter == NULL:
        logger_g_variantCalling_process.info('Run HaplotypeCaller without a read filter!')
    else:
        logger_g_variantCalling_process.info('Run HaplotypeCaller with a read filter that deals with some problems in bam!')
        command_count + = '--read-filter {0}'.format(read_filter)
    if reduce_logs == 'no':
        logger_g_variantCalling_process.info('Run HaplotypeCaller without reduceing the amount of chatter in the logs!')
    else:
        logger_g_variantCalling_process.info('Run HaplotypeCaller with reduceing the amount of chatter in the logs!')
        command_count + = '--QUIET'
    if create_output_variant_index == 'FALSE':
        logger_g_variantCalling_process.info('Run HaplotypeCaller with turnning off automatic variant index creation!')
        command_count + = '--create-output-variant-index {0}'.format(create_output_variant_index)
    else:
        logger_g_variantCalling_process.info('Run HaplotypeCaller with turnning on automatic variant index creation!')
        command_count + = '--create-output-variant-index {0}'.format(create_output_variant_index)
    logger_g_variantCalling_process.info('Begin to do germline variant calling.')
    os.system(command_count)
    logger_g_variantCalling_process.info('Time cost at germline variant calling == ' + str((time.time() - time_start) / 60) + 'min')
    print('Compeleted germline variant calling by GATK.')

def main():
    parser = argparse.ArgumentParser(usage = "\n\npython3 %(prog)s --source --sample_name  --output  --memorySize --gatk_dir  --ref_fa_file --exome_target_bed --ERC --samtools_dir --read_filter --reduce_logs --create_output_variant_index ")
    parser.add_argument("--source", help = "Path to input reads in FASTA format", type = str)
    parser.add_argument("--sample_name", help = "the sample name of raw reads", type = str)
    parser.add_argument("--output", help = "Path of output file", type = str)
    parser.add_argument("--memorySize ", help = "the cutoff of Java memory", type = str, default = '4G')
    parser.add_argument("--gatk_dir", help = "the install path of GATK4", type = str)
    parser.add_argument("--ref_fa_file", help = "the path of ref fasta", type = str)
    #parser.add_argument("--exome_target", help = "the list of exome intervals", type = str, default = NULL)
    parser.add_argument("--known_sites", help = "the list of --known-sites , sep=',' ", type = str)
    parser.add_argument("--read_length", help = "the length of reads ", type = int)
    parser.add_argument("--exome_target_bed", help = "the bed file of exome intervals", type = str) 
    parser.add_argument("--ERC", help = "switch to running HaplotypeCaller in GVCF mode", type = str, default = NULL)
    parser.add_argument("--samtools_dir", help = "the install path of samtools", type = str)
    parser.add_argument("--read_filter", help = "add a read filter that deals with some problems", type = str, default = NULL)
    parser.add_argument("--reduce_logs", help = "reduce the amount of chatter in the logs: --QUIET : yes", type = str, default = 'no' , choice = ['no', 'yes'])
    parser.add_argument("--create_output_variant_index", help = "turn off automatic variant index creation", type = str, default = 'FALSE' , choice = ['FALSE', 'TRUE'])
    parser.add_argument("-v", '-version', action = 'version', version =' %(prog)s 1.0')


#---input
    source = args.source
    sample = args.sample_name
    output = args.output
    memorySize = '-Xmx' + args.memorySize + ' ' + '-Djava.io.tmpdir=./'
    gatk_dir = args.gatk_dir
    samtools_dir = args.samtools_dir
    ref_fa_file = args.ref_fa_file
    known_sites = args.known_sites
    read_length = args.read_length
    exome_target_bed = args.exome_target_bed
    ERC = args.ERC
    read_filter = args.read_filter
    reduce_logs = args.reduce_logs
    create_output_variant_index = args.create_output_variant_index
 #---modify the known-sites
    known_sites=known_sites.replace(',',' --known-sites ')
 #---logging
    germline_VC_dir = out_dir + '/'+ 'germline_VC'
    os.makedirs(germline_VC_dir)

    log_dir = germline_VC_dir + 'log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    logger_g_variantCalling_process, logger_g_variantCalling_errors = store_logs(log_dir)
#---
    sam = source + '/' + sample + '_filtered.sam'
    sam_to_bem(gatk_dir, samtools_dir,
                            sam, sample, output, memorySize,
                            exome_target_bed, 
                            known_sites, read_length,
                            logger_g_variantCalling_process, logger_g_variantCalling_errors)
    marked_BQSR_bam = output + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
    germline_variant_calling(gatk_dir, marked_BQSR_bam,
                            sample, output, 
                            memorySize, ref_fa_file, 
                            exome_target_bed, ERC,
                            read_filter, read_length,
                            reduce_logs, create_output_variant_index,
                            logger_g_variantCalling_process, logger_g_variantCalling_errors)

if __name__ == '__main__':
    main()    