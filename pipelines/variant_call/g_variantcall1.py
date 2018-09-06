from __future__ import barry_as_FLUFL

__all__  =  [ 'sample', 'output ' , 'memory_size' , 'gatk_dir ' , 'vready_sam', 'marked_bqsr_bam',
              'ref_fa_file' ,'ref_fa_dict', 'exome_target_bed' , 'erc' , 'samtools_dir' ,
              'known_sites', 'snp_filter' , 'indel_filter',
              'strelka2_dir','bgzip','tabix','total_ref_chrom_fa_file',
              'logger_g_variantcalling_process', 'logger_g_variantcalling_errors','bqsr']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import os
import sys
import time
import shlex
import subprocess

#put the info output to the log
def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr


time_start = time.time()

def check_variant_exit(tmpvcf):
    headlines = 0
    f1 = open(tmpvcf)
    lens= len(f1.readlines())
    #print(lens)
    f1.close()
    f1 = open(tmpvcf)
    for line in f1.readlines():
        #print(line)
        if line.startswith('#'):
            headlines = headlines +1
    f1.close()
    #print(headlines)
    if headlines == lens:
        check_value = 0
    else:
        check_value = 1
    return check_value

def sam_to_bam(gatk_dir, samtools_dir,
               vready_sam, sample,
               output, memory_size,
               exome_target_bed, 
               ref_fa_file,ref_fa_dict,
               known_sites, 
               logger_g_variantcalling_process, logger_g_variantcalling_errors, bqsr):
    # sam to bam By SortSam in GATK
    sortedsam = output + '/' + sample + '_sorted.sam'
    command_count = '{0} --java-options "{1}" SortSam -SO coordinate -I {2} -O {3} --showHidden true'.format(
        gatk_dir, memory_size, vready_sam, sortedsam)
    time_start1= time.time()
    stdout, stderr = stdout_err(command_count)
    logger_g_variantcalling_process.info(stdout)
    logger_g_variantcalling_errors.info(stderr)
    logger_g_variantcalling_process.info("--{0}--GATK sort sam--cost {1} min.".format(sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # check the indexs of reference geonome fasta
    if not os.path.exists(ref_fa_dict):
        command_count_r1 = '{0} --java-options "{1}" CreateSequenceDictionary -R {2} -O {3} --showHidden true'.format(
        gatk_dir, memory_size, ref_fa_file, ref_fa_dict)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count_r1)
        logger_g_variantcalling_process.info(stdout)
        logger_g_variantcalling_errors.info(stderr)
        logger_g_variantcalling_process.info('--{0}--GATK build the dict of genome----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #--VCF.idx
    genome_idx = ref_fa_file + '.idx'
    if not os.path.exists(genome_idx):
        command_count_r2 = samtools_dir +' faidx '+ ref_fa_file
        time_start1 = time.time()
        os.system(command_count_r2)
        logger_g_variantcalling_process.info('--{0}--Samtools build the fai of genome----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    
    #target bed to exon intervel list
    if os.path.basename(exome_target_bed) != 'all':
        exon_interval = output + '/' + 'target_interval.list'
        command_count_r3 = '{0} --java-options "{1}" BedToIntervalList -I {2} -O {3} -SD {4} --showHidden true'.format(
            gatk_dir, memory_size, exome_target_bed, exon_interval, ref_fa_dict)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count_r3)
        logger_g_variantcalling_process.info(stdout)
        logger_g_variantcalling_errors.info(stderr)
        logger_g_variantcalling_process.info('--{0}--GATK target bed to exon intervel list----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #--sam to bem
    bam = output + '/' + sample + '_sorted.bam'
    command_count_r4 = samtools_dir+' view -bS '+ sortedsam +' > '+ bam
    time_start1 = time.time()
    os.system(command_count_r4)
    
    logger_g_variantcalling_process.info('--{0}--Samtools transform sam to bam----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #statistics of coverages
    if os.path.basename(exome_target_bed) != 'all':
        cov_file =  output + '/' + sample + '.cov.txt'
        command_count_1 = '{0} CollectHsMetrics -BI {1} -TI {2} -I {3} -O {4} --showHidden true'.format(
            gatk_dir, exon_interval, exon_interval, bam, cov_file)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count_1)
        logger_g_variantcalling_process.info(stdout)
        logger_g_variantcalling_errors.info(stderr)
        logger_g_variantcalling_process.info('--{0}--GATK count the coverage----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #build index of bam By samtools
    command_count1 = '{0} index {1}'.format(samtools_dir, bam)
    time_start1 = time.time()
    os.system(command_count1)
    logger_g_variantcalling_process.info('--{0}--Samtools build the index of bam----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #MarkDuplicates
    mark_bam = output + '/' + sample + '_sorted.MarkDuplicates.bam'
    bam_metrics = output + '/' + sample + '_sorted.MarkDuplicates.metrics'
    command_count2 ='{0} --java-options "{1}" MarkDuplicates -I {2} -O {3} -M {4} --REMOVE_SEQUENCING_DUPLICATES false --showHidden true'.format(
        gatk_dir, memory_size, bam, mark_bam, bam_metrics)
    time_start1= time.time()
    stdout, stderr = stdout_err(command_count2)
    logger_g_variantcalling_process.info(stdout)
    logger_g_variantcalling_errors.info(stderr)
    logger_g_variantcalling_process.info("--{0}--GATK marks the duplicates----cost {1} min.".format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #add row 'RG' in head
    mark_rg_bam = output + '/' + sample + '_sorted.MarkDuplicates.RG.bam'
    command_count3 = '{0} --java-options "{1}" AddOrReplaceReadGroups -I {2} -O {3} -LB lib1 -PL illumina -PU unit1 \
                     -SM {4} --showHidden true'.format(gatk_dir, memory_size, mark_bam, mark_rg_bam, bam_metrics)
    stdout, stderr = stdout_err(command_count3)
    logger_g_variantcalling_process.info(stdout)
    logger_g_variantcalling_errors.info(stderr)
    #build the index of the marked_fixed bam
    command_count4 = '{0} index {1}'.format(samtools_dir, mark_rg_bam)
    time_start1= time.time()
    os.system(command_count4)
    logger_g_variantcalling_process.info("--{0}--Samtools build the index of  marked bam----cost {1} min.".format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #---------------------------------------
    # Base(Quality Score) Recalibration
    print(bqsr)
    if bqsr == 'yes':
        original_recal_data_table = output + '/' + sample + '_original.recal_data.csv'
        original_bgsr_bam = output + '/' + sample + '_sorted.MarkDuplicates.originalBQSR.bam'
        recal_data_table = output + '/' + sample + '_bqsr.recal_data.csv'
        bgsr_bam = output + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
        if os.path.basename(exome_target_bed) != 'all':
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--first
            command_count5_1 = '{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} -L {4} --known-sites {5} -O {6} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, exon_interval, known_sites, original_recal_data_table)
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--first
            command_count7_1 = '{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -L {5} -O {6} --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, original_recal_data_table, exon_interval,
                original_bgsr_bam)
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--second
            command_count5_2 = '{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} -L {4} --known-sites {5} -O {6} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, exon_interval, known_sites, recal_data_table)
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--second
            command_count7_2 = '{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -L {5} -O {6} --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, recal_data_table, exon_interval, bgsr_bam)
        else:
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model
            command_count5_1 = '{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} --known-sites {4} -O {5} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, known_sites, original_recal_data_table)
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--first
            command_count7_1 = '{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -O {5} --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, original_recal_data_table, original_bgsr_bam)
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model
            command_count5_2 = '{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} --known-sites {4} -O {5} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, known_sites, recal_data_table)
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--second
            command_count7_2 = '{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -O {5} --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, recal_data_table, bgsr_bam)
        stdout, stderr = stdout_err(command_count5_1)
        stdout, stderr = stdout_err(command_count7_1)
        stdout, stderr = stdout_err(command_count5_2)
        stdout, stderr = stdout_err(command_count7_2)
        base_reports = output + '/' + sample + '.BQSR.before_VS_after.pdf'
        base_reports_csv = output + '/' + sample + '.BQSR.before_VS_after.csv'
        command_count6 = '{0} --java-options "{1}" AnalyzeCovariates --before-report-file {2} --after-report-file {3} --plots-report-file {4} --intermediate-csv-file {5}'.format(
            gatk_dir, memory_size, original_recal_data_table, recal_data_table, base_reports, base_reports_csv)
        stdout, stderr = stdout_err(command_count6)
    else:
        bgsr_bam = 'null'
    return mark_rg_bam, bgsr_bam
    #time used for translating the format of the sam by samtools and GATK
    logger_g_variantcalling_process.info('Compeleted translating the format of the sam by samtools and GATK.')

def germline_variant_calling(gatk_dir, marked_BQSR_bam,
                            sample, output, 
                            memory_size, ref_fa_file,
                            exon_interval, erc,
                            snp_filter,indel_filter,
                            logger_g_variantcalling_process, logger_g_variantcalling_errors):
    if exon_interval != 'all':
         command_count ='{0} --java-options "{1}" HaplotypeCaller -R {2} -I {3} -L {4} --minimum-mapping-quality 20 --showHidden true'.format(
            gatk_dir, memory_size, ref_fa_file, marked_BQSR_bam, exon_interval)
    else:
         command_count ='{0} --java-options "{1}" HaplotypeCaller -R {2} -I {3} --minimum-mapping-quality 20 --showHidden true'.format(
            gatk_dir, memory_size, ref_fa_file, marked_BQSR_bam)
    logger_g_variantcalling_process.info('Begin to confirm the options parameters of running HaplotypeCaller.')
    vcf1 =  output + '/' + sample + '.raw_variants.vcf'
    if erc == 'no':
        logger_g_variantcalling_process.info('Runs HaplotypeCaller in default mode on a single input BAM file containing sequence data!')
        vcf =  output + '/' + sample + '.raw_variants.vcf'
        command_count = command_count + ' -O ' + vcf
    else:
        if erc == 'GVCF':
            logger_g_variantcalling_process.info('Runs HaplotypeCaller in GVCF mode!')
            vcf =  output + '/' + sample + '.raw_variants.g.vcf'
            command_count = command_count + ' -O {0} -ERC {1}'.format(vcf, erc)
        else:
            logger_g_variantcalling_process.info('{0} is not a HaplotypeCaller model. Please check the input parameter of ERC!'.format(ERC))
    logger_g_variantcalling_process.info('Begin to do germline variant calling.')
    time_start1 = time.time()
    stdout, stderr = stdout_err(command_count)
    logger_g_variantcalling_process.info(stdout)
    logger_g_variantcalling_errors.info(stderr)
    logger_g_variantcalling_process.info('Compeleted HaplotypeCaller by GATK.')
    logger_g_variantcalling_process.info('--{0}--Time cost at HaplotypeCaller ----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #-Joint-Call Cohort
    if erc == 'GVCF':
        command_count1 ='{0} --java-options "{1}" GenotypeGVCFs -R {2} --variant {3} -O {4}'.format(
        gatk_dir, memory_size, ref_fa_file, vcf, vcf1)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count1)
        logger_g_variantcalling_process.info(stdout)
        logger_g_variantcalling_errors.info(stderr)
        logger_g_variantcalling_process.info('Compeleted GenotypeGVCFs by GATK.')
        logger_g_variantcalling_process.info('--{0}--Time cost at GenotypeGVCFs----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #-SNP
    snp_vcf = output + '/' + sample + '.raw_variants_SNP.vcf'
    command_count2 ='{0} --java-options "{1}" SelectVariants -R {2} --variant {3} -O {4} --select-type-to-include SNP --showHidden true'.format(
        gatk_dir, memory_size, ref_fa_file, vcf1, snp_vcf)
    stdout, stderr = stdout_err(command_count2)
    logger_g_variantcalling_process.info(stdout)
    logger_g_variantcalling_errors.info(stderr)
    logger_g_variantcalling_process.info('Compeleted SelectVariants SNP by GATK.')
    logger_g_variantcalling_process.info('--{0}--Time cost at SelectVariants SNP----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))
    #indel
    indel_vcf = output + '/' + sample + '.raw_variants_indel.vcf'
    command_count3 ='{0} --java-options "{1}" SelectVariants -R {2} --variant {3} -O {4} --select-type-to-include INDEL --showHidden true'.format(
        gatk_dir, memory_size, ref_fa_file, vcf1, indel_vcf)
    stdout, stderr = stdout_err(command_count3)
    logger_g_variantcalling_process.info(stdout)
    logger_g_variantcalling_errors.info(stderr)
    logger_g_variantcalling_process.info('Compeleted SelectVariants indel by GATK.')
    logger_g_variantcalling_process.info('--{0}--Time cost at SelectVariants indel----cost {1} min.'.format(sample,str('%.3f' % ((time.time()-time_start1)/60))))

def strelka2_call(strelka2_dir,
                  bgzip,
                  tabix,
                  total_ref_chrom_fa_file,
                  outputdir,
                 sample, bam,
                 exome_target_bed):
    os.system('cp {0} {1}'.format(exome_target_bed, outputdir))
    command1= bgzip + ' '+ outputdir+'/'+os.path.basename(exome_target_bed)
    print(command1)
    os.system(command1)
    command2= tabix + ' '+ outputdir+'/'+ os.path.basename(exome_target_bed) +'.gz'
    print(command2)
    os.system(command2)
    command3='python {0} --bam {1} --callRegions {2} --exome --referenceFasta {3} --runDir {4}'.format(
        strelka2_dir, bam, outputdir+'/'+os.path.basename(exome_target_bed)+'.gz', total_ref_chrom_fa_file, outputdir)
    stdout, stderr = stdout_err(command3)
    #print(stderr)
    run_py = outputdir+'/'+ 'runWorkflow.py'
    command4='python {0} -m local -j 1'.format(run_py)
    stdout, stderr = stdout_err(command4)