from __future__ import barry_as_FLUFL

__all__  =  [ 'sample' , 'output ' , 'memorySize' , 'gatk_dir ' , 'vready_sam', 'marked_BQSR_bam',
              'ref_fa_file' , 'exome_target_bed' , 'ERC' , 'samtools_dir' ,
              'known_sites', 'read_length', 
              'logger_g_variantCalling_process', 'logger_g_variantCalling_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import os
import re
import sys
import time
import gzip
import itertools

time_start = time.time()

def sam_to_bem(gatk_dir, samtools_dir,
               vready_sam, sample, 
               output, memorySize, 
               exome_target_bed, known_sites, 
               read_length, 
               logger_g_variantCalling_process, logger_g_variantCalling_errors):
    # sam to bam By SortSam in GATK
    bam = output + '/' + sample + '.bam'
    command_count = '{0} --java-options "{1}" SortSam -SO coordinate -I {2} -O {3}'.format(
        gatk_dir, memorySize, sam, bam)
    
    time_start1= time.time()
    os.system(command_count)
    logger_g_variantCalling_process.info("ATK sort sam to bam---cost %.2f min.", (time.time()-time_start1)/60)
    #
    #statistics of coverages
    cov_file =  output + '/' + sample + '.cov.txt'
    command_count_1 = '{0} CollectHsMetrics -BI {1} -TI {2} -I {3} -O {4}'.format(
        gatk_dir, exome_target_bed, exome_target_bed, bam, cov_file)
    time_start1 = time.time()
    os.system(command_count_1)
    logger_g_variantCalling_process.info('GATK count the coverage----cost %.2f min.', (time.time()-time_start1)/60)
    #build index of bam By samtools
    command_count1 = '{0} index {1}'.format(samtools_dir, bam)
    time_start1 = time.time()
    os.system(command_count1)
    logger_g_variantCalling_process.info('Samtools build the index of bam----cost %.2f min.', (time.time()-time_start1)/60)
    #MarkDuplicates
    mark_bam = output + '/' + sample + '_marked.bam'
    bam_metrics = output + '/' + sample + '.metrics'
    command_count2 ='{0} --java-options "{1}" MarkDuplicates -I {2} -O {3} -M {4}'.format(
        gatk_dir, memorySize, bam, mark_bam, bam_metrics)
    time_start1= time.time()
    os.system(command_count2)
    logger_g_variantCalling_process.info("GATK marks the duplicates----cost %.2f min.", (time.time()-time_start1)/60)
    #build the index of the marked_fixed bam
    command_count4 = '{0} index {1}'.format(samtools_dir, mark_bam)
    time_start1= time.time()
    os.system(command_count4)
    logger_g_variantCalling_process.info("Samtools build the index of  marked bam----cost %.2f min.", (time.time()-time_start1)/60)
    #---------------------------------------
    #Base(Quality Score) Recalibration
    #BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model
    recal_data_table =  output + '/' + sample + '.recal_data.table'
    command_count5 ='{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} -L {4} -ip {5} --known-sites {6} -O {7}'.format(
        gatk_dir, memorySize, ref_fa_file, mark_bam, exome_target_bed, read_length, known_sites, recal_data_table)
    time_start1= time.time()
    os.system(command_count5)
    logger_g_variantCalling_process.info("GATK generate Base Quality Score Recalibration (BQSR) model----cost %.2f min.", (time.time()-time_start1)/60)
    #GatherBQSRReports--Generate report of Base Quality Score Recalibration (BQSR) model
    BQSRReports = output + '/' + sample + '.BQSR.report'
    command_count6 ='{0} --java-options "{1}" GatherBQSRReports -I {2} -O {3}'.format(
        gatk_dir, memorySize,  recal_data_table, BQSRReports)
    time_start1= time.time()
    os.system(command_count6)
    logger_g_variantCalling_process.info("ATK generate report of Base Quality Score Recalibration (BQSR) model----cost %.2f min.", (time.time()-time_start1)/60)
    #ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model
    bgsr_bam =  output + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
    command_count7 ='{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4}  -L {5} -ip {6} -O {7}'.format(
        gatk_dir, memorySize, ref_fa_file, mark_bam, recal_data_table, exome_target_bed, read_length, bgsr_bam)
    time_start1= time.time()
    os.system(command_count7)
    logger_g_variantCalling_process.info("GATK apply Base Quality Score Recalibration (BQSR) model----cost %.2f min.", (time.time()-time_start1)/60)
    #time used for translating the format of the sam by samtools and GATK
    logger_g_variantCalling_process.info('Compeleted translating the format of the sam by samtools and GATK.')

def germline_variant_calling(gatk_dir, marked_BQSR_bam,
                            sample, output, 
                            memorySize, ref_fa_file,
                            exome_target_bed, ERC,
                            read_filter, read_length,
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
    logger_g_variantCalling_process.info('Begin to do germline variant calling.')
    time_start1 = time.time()
    os.system(command_count)
    logger_g_variantCalling_process.info('Compeleted germline variant calling by GATK.')
    logger_g_variantCalling_process.info('Time cost at germline variant calling == ' + str((time.time() - time_start1) / 60) + 'min')
