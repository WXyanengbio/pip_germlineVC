from __future__ import barry_as_FLUFL

__all__  =  [ 'sample' , 'output ' , 'memorySize' , 'gatk_dir ' , 'vready_sam', 'marked_BQSR_bam',
              'ref_fa_file' ,'ref_fa_dict', 'exome_target_bed' , 'ERC' , 'samtools_dir' ,
              'known_sites', 'read_filter', 'snp_filter' , 'indel_filter',
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
               exome_target_bed, 
               ref_fa_file,ref_fa_dict,
               known_sites, 
               logger_g_variantCalling_process, logger_g_variantCalling_errors):
    # sam to bam By SortSam in GATK
    sortedsam = output + '/' + sample + '_sorted.sam'
    command_count = '{0} --java-options "{1}" SortSam -SO coordinate -I {2} -O {3}'.format(
        gatk_dir, memorySize, vready_sam, sortedsam)
    time_start1= time.time()
    os.system(command_count)
    logger_g_variantCalling_process.info("GATK sort sam--cost %.2f min.", (time.time()-time_start1)/60)
    # check the indexs of reference geonome fasta
    if not os.path.exists(ref_fa_dict):
        command_count_r1 = '{0} --java-options "{1}" CreateSequenceDictionary -R {2} -O {3}'.format(
        gatk_dir, memorySize, ref_fa_file, ref_fa_dict)
        time_start1 = time.time()
        os.system(command_count_r1)
        logger_g_variantCalling_process.info('GATK build the dict of genome----cost %.2f min.', (time.time()-time_start1)/60)
    #--VCF.idx
    genome_idx = ref_fa_file + '.idx'
    if not os.path.exists(genome_idx):
        command_count_r2 = samtools_dir +' faidx '+ ref_fa_file
        time_start1 = time.time()
        os.system(command_count_r2)
        logger_g_variantCalling_process.info('Samtools build the fai of genome----cost %.2f min.', (time.time()-time_start1)/60)
    
    #target bed to exon intervel list
    Exon_Interval = output + '/' + 'target_interval.list'
    command_count_r3 = '{0} --java-options "{1}" BedToIntervalList -I {2} -O {3} -SD {4}'.format(
        gatk_dir, memorySize, exome_target_bed, Exon_Interval, ref_fa_dict)
    time_start1 = time.time()
    os.system(command_count_r3)
    logger_g_variantCalling_process.info('GATK target bed to exon intervel list----cost %.2f min.', (time.time()-time_start1)/60)
    #--sam to bem
    bam = output + '/' + sample + '_sorted.bam'
    command_count_r4 = samtools_dir+' view -bS '+ sortedsam +' > '+ bam
    time_start1 = time.time()
    os.system(command_count_r4)
    logger_g_variantCalling_process.info('Samtools transform sam to bam----cost %.2f min.', (time.time()-time_start1)/60)
    #statistics of coverages
    cov_file =  output + '/' + sample + '.cov.txt'
    command_count_1 = '{0} CollectHsMetrics -BI {1} -TI {2} -I {3} -O {4}'.format(
        gatk_dir, Exon_Interval, Exon_Interval, bam, cov_file)
    time_start1 = time.time()
    os.system(command_count_1)
    logger_g_variantCalling_process.info('GATK count the coverage----cost %.2f min.', (time.time()-time_start1)/60)
    #build index of bam By samtools
    command_count1 = '{0} index {1}'.format(samtools_dir, bam)
    time_start1 = time.time()
    os.system(command_count1)
    logger_g_variantCalling_process.info('Samtools build the index of bam----cost %.2f min.', (time.time()-time_start1)/60)
    #MarkDuplicates
    mark_bam = output + '/' + sample + '_sorted.MarkDuplicates.bam'
    bam_metrics = output + '/' + sample + '_sorted.MarkDuplicates.metrics'
    command_count2 ='{0} --java-options "{1}" MarkDuplicates -I {2} -O {3} -M {4} --REMOVE_SEQUENCING_DUPLICATES false'.format(
        gatk_dir, memorySize, bam, mark_bam, bam_metrics)
    time_start1= time.time()
    os.system(command_count2)
    logger_g_variantCalling_process.info("GATK marks the duplicates----cost %.2f min.", (time.time()-time_start1)/60)
    #add row 'RG' in head
    mark_RG_bam = output + '/' + sample + '_sorted.MarkDuplicates.RG.bam'
    command_count3 = '{0} --java-options "{1}" AddOrReplaceReadGroups -I {2} -O {3} -LB lib1 -PL illumina -PU unit1 -SM {4}'.format(
        gatk_dir, memorySize, mark_bam, mark_RG_bam, bam_metrics)
    os.system(command_count3)
    #build the index of the marked_fixed bam
    command_count4 = '{0} index {1}'.format(samtools_dir, mark_RG_bam)
    time_start1= time.time()
    os.system(command_count4)
    logger_g_variantCalling_process.info("Samtools build the index of  marked bam----cost %.2f min.", (time.time()-time_start1)/60)
    #---------------------------------------
    #Base(Quality Score) Recalibration
    #BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model
    recal_data_table =  output + '/' + sample + '.recal_data.table'
    command_count5 ='{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} -L {4} --known-sites {5} -O {6}'.format(
        gatk_dir, memorySize, ref_fa_file, mark_RG_bam, Exon_Interval, known_sites, recal_data_table)
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
    command_count7 ='{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -L {5} -O {6}'.format(
        gatk_dir, memorySize, ref_fa_file, mark_RG_bam, recal_data_table, Exon_Interval, bgsr_bam)
    time_start1= time.time()
    os.system(command_count7)
    logger_g_variantCalling_process.info("GATK apply Base Quality Score Recalibration (BQSR) model----cost %.2f min.", (time.time()-time_start1)/60)
    #time used for translating the format of the sam by samtools and GATK
    logger_g_variantCalling_process.info('Compeleted translating the format of the sam by samtools and GATK.')

def germline_variant_calling(gatk_dir, marked_BQSR_bam,
                            sample, output, 
                            memorySize, ref_fa_file,
                            Exon_Interval, ERC,
                            read_filter,
                            snp_filter,indel_filter,
                            logger_g_variantCalling_process, logger_g_variantCalling_errors):
    command_count ='{0} --java-options "{1}" HaplotypeCaller -R {2} -I {3} -L {4}'.format(
        gatk_dir, memorySize, ref_fa_file, marked_BQSR_bam, Exon_Interval)
    logger_g_variantCalling_process.info('Begin to confirm the options parameters of running HaplotypeCaller.')
    vcf1 =  output + '/' + sample + '.raw_variants.vcf'
    if ERC == 'no':
        logger_g_variantCalling_process.info('Runs HaplotypeCaller in default mode on a single input BAM file containing sequence data!')
        vcf =  output + '/' + sample + '.raw_variants.vcf'
        command_count = command_count + ' ' + vcf
    else:
        if ERC == 'GVCF':
            logger_g_variantCalling_process.info('Runs HaplotypeCaller in GVCF mode!')
            vcf =  output + '/' + sample + '.raw_variants.g.vcf'
            command_count = command_count + ' -O {0} -ERC {1}'.format(vcf, ERC)
        else:
            logger_g_variantCalling_process.info('{0} is not a HaplotypeCaller model. Please check the input parameter of ERC!'.format(ERC))

    if read_filter == 'no':
        logger_g_variantCalling_process.info('Run HaplotypeCaller without a read filter!')
    else:
        logger_g_variantCalling_process.info('Run HaplotypeCaller with a read filter that deals with some problems in bam!')
        command_count = command_count + ' ' + '--read-filter {0}'.format(read_filter)
    logger_g_variantCalling_process.info('Begin to do germline variant calling.')
    time_start1 = time.time()
    os.system(command_count)
    logger_g_variantCalling_process.info('Compeleted HaplotypeCaller by GATK.')
    logger_g_variantCalling_process.info('Time cost at HaplotypeCaller == ' + str((time.time() - time_start1) / 60) + 'min')
    #-Joint-Call Cohort
    if ERC == 'GVCF':
        command_count1 ='{0} --java-options "{1}" GenotypeGVCFs -R {2} --variant {3} -O {4}'.format(
        gatk_dir, memorySize, ref_fa_file, vcf,vcf1)
        time_start1 = time.time()
        os.system(command_count1)
        logger_g_variantCalling_process.info('Compeleted GenotypeGVCFs by GATK.')
        logger_g_variantCalling_process.info('Time cost at GenotypeGVCFs == ' + str((time.time() - time_start1) / 60) + 'min')
    #-SNP
    snp_vcf = output + '/' + sample + '.raw_variants_SNP.vcf'
    command_count2 ='{0} --java-options "{1}" SelectVariants -R {2} --variant {3} -O {4} --select-type-to-include SNP'.format(
        gatk_dir, memorySize, ref_fa_file, vcf1, snp_vcf)
    os.system(command_count2)
    logger_g_variantCalling_process.info('Compeleted SelectVariants SNP by GATK.')
    logger_g_variantCalling_process.info('Time cost at SelectVariants SNP == ' + str((time.time() - time_start1) / 60) + 'min')
    #filter variant in SNP
    filter_snp_vcfs = output + '/' + sample + '.filter_SNP.vcf'
    command_count2fs ='{0} --java-options "{1}" VariantFiltration -R {2} --variant {3} -O {4} --filter-expression "{5}" --filter-name "my_snp_filter"'.format(
        gatk_dir, memorySize, ref_fa_file, snp_vcf,filter_snp_vcfs, snp_filter)
    os.system(command_count2fs)
    logger_g_variantCalling_process.info('Compeleted filter SNP by GATK.')
    logger_g_variantCalling_process.info('Time cost at filter SNP == ' + str((time.time() - time_start1) / 60) + 'min')
    #indel
    indel_vcf = output + '/' + sample + '.raw_variants_indel.vcf'
    command_count3 ='{0} --java-options "{1}" SelectVariants -R {2} --variant {3} -O {4} --select-type-to-include INDEL'.format(
        gatk_dir, memorySize, ref_fa_file, vcf1, indel_vcf)
    os.system(command_count3)
    logger_g_variantCalling_process.info('Compeleted SelectVariants indel by GATK.')
    logger_g_variantCalling_process.info('Time cost at SelectVariants indel == ' + str((time.time() - time_start1) / 60) + 'min')
    #filter variant in indel
    filter_indel_vcfi = output + '/' + sample + '.filter_indel.vcf'
    command_count2fi ='{0} --java-options "{1}" VariantFiltration -R {2} --variant {3} -O {4} --filter-expression "{5}" --filter-name "my_indel_filter"'.format(
        gatk_dir, memorySize, ref_fa_file, indel_vcf, filter_indel_vcfi, indel_filter)
    os.system(command_count2fi)
    logger_g_variantCalling_process.info('Compeleted filter SNP by GATK.')
    logger_g_variantCalling_process.info('Time cost at filter SNP == ' + str((time.time() - time_start1) / 60) + 'min')
