from __future__ import barry_as_FLUFL

__all__ = ['sample', 'output ', 'memory_size',  'vready_sam', 'marked_bam',
           'gatk_dir ', 'samtools_dir', 'bcftools_dir', 'varsan2_dir',
           'ref_fa_file', 'ref_fa_dict', 'exome_target_bed', 'exon_interval', 'erc', 
           'known_sites', 'snp_filter', 'indel_filter',
           'strelka2_dir', 'bgzip', 'tabix', 'total_ref_chrom_fa_file',
           'logger_g_variantcalling_process', 'logger_g_variantcalling_errors', 'bqsr',
           'smcounter', 'mtdepth', 'rpb', 'ncpu', 'minbq', 'minmq', 'hplen', 'mismatchthr',
           'mtdrop', 'maxmt', 'primerdist', 'threshold', 
           'bedtandemrepeats', 'bedrepeatmaskersubset', 'bedtools_dir', 'logfile', 'renew']
__version__ = '1.0'
__author__ = 'Wang Xian'

import os
import sys
import time
import shlex
import subprocess
sys.path.append("..")
from pipelines.log.log_v1 import store_germline_vc_logs


def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr


def check_variant_exit(tmpvcf):
    headlines = 0
    f1 = open(tmpvcf)
    lens = len(f1.readlines())
    # print(lens)
    f1.close()
    f1 = open(tmpvcf)
    for line in f1.readlines():
        # print(line)
        if line.startswith('#'):
            headlines = headlines +1
    f1.close()
    # print(headlines)
    if headlines == lens:
        check_value = 0
    else:
        check_value = 1
    return check_value


def sam_to_bam(gatk_dir, samtools_dir, vready_sam, sample, output, memory_size, exome_target_bed,
               ref_fa_file, ref_fa_dict, known_sites,
               logger_g_variantcalling_process, logger_g_variantcalling_errors, bqsr, renew):
    # check the indexs of reference geonome fasta
    if not os.path.exists(ref_fa_dict):
        command_count_r1 = '{0} --java-options "{1}" CreateSequenceDictionary -R {2} -O {3} --showHidden true'.format(
            gatk_dir, memory_size, ref_fa_file, ref_fa_dict)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count_r1)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
        store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                               '--{0}--GATK build the dict of genome----cost {1} min.'.format(
                                   sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # --VCF.idx
    genome_idx = ref_fa_file + '.idx'
    if not os.path.exists(genome_idx) or renew == 'T':
        command_count_r2 = samtools_dir + ' faidx ' + ref_fa_file
        time_start1 = time.time()
        os.system(command_count_r2)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                               '--{0}--Samtools build the fai of genome----cost {1} min.'.format(
                                   sample, str('%.3f' % ((time.time()-time_start1)/60))))
    
    # target bed to exon intervel list
    if os.path.basename(exome_target_bed) != 'all':
        exon_interval = output + '/' + 'target_interval.list'
        if not os.path.exists(exon_interval) or renew == 'T':
            command_count_r3 = '{0} --java-options "{1}" BedToIntervalList -I {2} -O {3} -SD {4} --showHidden true'.format(
              gatk_dir, memory_size, exome_target_bed, exon_interval, ref_fa_dict)
            time_start1 = time.time()
            stdout, stderr = stdout_err(command_count_r3)
            store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
            store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
            store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                               '--{0}--GATK target bed to exon intervel list----cost {1} min.'.format(
                                   sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # --sam to bam
    bam = output + '/' + sample + '.bam'
    if not os.path.exists(bam) or renew == 'T':
        command_count_r4 = samtools_dir +' view -bS ' + vready_sam + ' > ' + bam
        time_start1 = time.time()
        os.system(command_count_r4)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           '--{0}--Samtools transform sam to bam----cost {1} min.'.format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # bam to sortedbam By SortSam in GATK
    sortedbam = output + '/' + sample + '_sorted.bam'
    if not os.path.exists(sortedbam) or renew == 'T':
        command_count_r4 = samtools_dir +' sort ' + bam + ' > ' + sortedbam
        time_start1 = time.time()
        os.system(command_count_r4)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           '--{0}--Samtools sort bam----cost {1} min.'.format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # statistics of coverages
    if os.path.basename(exome_target_bed) != 'all':
        cov_file = output + '/' + sample + '.cov.txt'
        command_count_1 = '{0} CollectHsMetrics -BI {1} -TI {2} -I {3} -O {4} --showHidden true'.format(
            gatk_dir, exon_interval, exon_interval, sortedbam, cov_file)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count_1)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
        store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                               '--{0}--GATK count the coverage----cost {1} min.'.format(
                                   sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # build index of bam By samtools
    bamindex = output + '/' + sample + '_sorted.bam.bai'
    if not os.path.exists(bamindex) or renew == 'T':
        command_count1 = '{0} index {1}'.format(samtools_dir, sortedbam)
        time_start1 = time.time()
        os.system(command_count1)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           '--{0}--Samtools build the index of bam----cost {1} min.'.format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # MarkDuplicates
    mark_bam = output + '/' + sample + '_sorted.MarkDuplicates.bam'
    bam_metrics = output + '/' + sample + '_sorted.MarkDuplicates.metrics'
    if not os.path.exists(mark_bam) or renew == 'T':
        command_count2 = ' '.join([gatk_dir,
                               '--java-options "{0}" MarkDuplicates'.format(memory_size),
                               '-I', sortedbam, '-O', mark_bam, '-M', bam_metrics,
                               '--REMOVE_SEQUENCING_DUPLICATES false'])
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count2)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
        store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           "--{0}--GATK marks the duplicates----cost {1} min.".format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # add row 'RG' in head
    mark_rg_bam = output + '/' + sample + '_sorted.MarkDuplicates.RG.bam'
    if not os.path.exists(mark_rg_bam) or renew == 'T':
        command_count3 = ' '.join([gatk_dir,
                               '--java-options "{0}" AddOrReplaceReadGroups'.format(memory_size),
                               '-I', mark_bam, '-O', mark_rg_bam, '-LB lib1 -PL illumina -PU unit1',
                               '-SM', bam_metrics,
                               '--showHidden true'])
        stdout, stderr = stdout_err(command_count3)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
        store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
    # build the index of the marked_fixed bam
        command_count4 = '{0} index {1}'.format(samtools_dir, mark_rg_bam)
        time_start1 = time.time()
        os.system(command_count4)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           "--{0}--Samtools build the index of  marked bam----cost {1} min.".format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # Base(Quality Score) Recalibration
    # print(bqsr)
    if bqsr == 'yes':
        original_recal_data_table = output + '/' + sample + '_original.recal_data.csv'
        original_bgsr_bam = output + '/' + sample + '_sorted.MarkDuplicates.originalBQSR.bam'
        recal_data_table = output + '/' + sample + '_bqsr.recal_data.csv'
        bgsr_bam = output + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
        if os.path.basename(exome_target_bed) != 'all':
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--first
            command_count5_1 = ' '.join([gatk_dir,
                                         '--java-options "{0}" BaseRecalibrator'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '-L', exon_interval, '--known-sites', known_sites,
                                         '-O', original_recal_data_table, '--bqsr-baq-gap-open-penalty 40',
                                         '--showHidden true'])
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--first
            command_count7_1 = ' '.join([gatk_dir,
                                         '--java-options "{0}" ApplyBQSR'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '-bqsr', original_recal_data_table,
                                         '-L', exon_interval,
                                         '-O', original_bgsr_bam,
                                         '--showHidden true'])
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--second
            command_count5_2 = ' '.join([gatk_dir,
                                         '--java-options "{0}" BaseRecalibrator'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '-L', exon_interval, '--known-sites', known_sites,
                                         '-O', recal_data_table, '--bqsr-baq-gap-open-penalty 40',
                                         '--showHidden true'])
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--second
            command_count7_2 = ' '.join([gatk_dir,
                                         '--java-options "{0}" ApplyBQSR'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '-bqsr', recal_data_table,
                                         '-L', exon_interval,
                                         '-O', bgsr_bam,
                                         '--showHidden true'])
        else:
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--first
            command_count5_1 = ' '.join([gatk_dir,
                                         '--java-options "{0}" BaseRecalibrator'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '--known-sites', known_sites,
                                         '-O', original_recal_data_table, '--bqsr-baq-gap-open-penalty 40',
                                         '--showHidden true'])
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--first
            command_count7_1 = ' '.join([gatk_dir,
                                         '--java-options "{0}" ApplyBQSR'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '-bqsr', original_recal_data_table,
                                         '-O', original_bgsr_bam,
                                         '--showHidden true'])
            # BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--second
            command_count5_2 = ' '.join([gatk_dir,
                                         '--java-options "{0}" BaseRecalibrator'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '--known-sites', known_sites,
                                         '-O', recal_data_table, '--bqsr-baq-gap-open-penalty 40',
                                         '--showHidden true'])
            # ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--second
            command_count7_2 = ' '.join([gatk_dir,
                                         '--java-options "{0}" ApplyBQSR'.format(memory_size),
                                         '-R', ref_fa_file,
                                         '-I', mark_rg_bam, '-bqsr', recal_data_table,
                                         '-O', bgsr_bam,
                                         '--showHidden true'])
        stdout, stderr = stdout_err(command_count5_1)
        stdout, stderr = stdout_err(command_count7_1)
        stdout, stderr = stdout_err(command_count5_2)
        stdout, stderr = stdout_err(command_count7_2)
        base_reports = output + '/' + sample + '.BQSR.before_VS_after.pdf'
        base_reports_csv = output + '/' + sample + '.BQSR.before_VS_after.csv'
        command_count6 = ' '.join([gatk_dir,
                                   '--java-options "{0}" AnalyzeCovariates'.format(memory_size),
                                   '--before-report-file', original_recal_data_table,
                                   '--after-report-file', recal_data_table,
                                   '--plots-report-file', base_reports,
                                   '--intermediate-csv-file', base_reports_csv])
        stdout, stderr = stdout_err(command_count6)
    else:
        bgsr_bam = 'null'
    store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           'Compeleted translating the format of the sam by samtools and GATK.')
    return mark_rg_bam, bgsr_bam


def germline_variant_calling(gatk_dir, marked_bam,
                             sample, output,
                             memory_size, ref_fa_file,
                             exon_interval, erc,
                             snp_filter, indel_filter,
                             logger_g_variantcalling_process, logger_g_variantcalling_errors):
    if exon_interval != 'all':
        command_count = ' '.join([gatk_dir,
                                  '--java-options "{0}" HaplotypeCaller'.format(memory_size),
                                  '-R', ref_fa_file,
                                  '-I', marked_bam,
                                  '-L', exon_interval,
                                  '--minimum-mapping-quality 20 --showHidden true'])
    else:
        command_count = ' '.join([gatk_dir,
                                  '--java-options "{0}" HaplotypeCaller'.format(memory_size),
                                  '-R', ref_fa_file,
                                  '-I', marked_bam,
                                  '--minimum-mapping-quality 20 --showHidden true'])
    store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           'Begin to confirm the options parameters of running HaplotypeCaller.')
    vcf1 = output + '/' + sample + '.raw_variants.vcf'
    if erc == 'no':
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                               'Runs HaplotypeCaller in default mode on a single BAM file containing sequence data!')
        vcf = output + '/' + sample + '.raw_variants.vcf'
        command_count = command_count + ' -O ' + vcf
    else:
        if erc == 'GVCF':
            store_germline_vc_logs(logger_g_variantcalling_process, 'null', 'Runs HaplotypeCaller in GVCF mode!')
            vcf = output + '/' + sample + '.raw_variants.g.vcf'
            command_count = command_count + ' -O {0} -ERC {1}'.format(vcf, erc)
        else:
            store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                                   ERC + ' is not a HaplotypeCaller model. Please check the input parameter of ERC!')
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', 'Begin to do germline variant calling.')
    time_start1 = time.time()
    stdout, stderr = stdout_err(command_count)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
    store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', 'Compeleted HaplotypeCaller by GATK.')
    store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           '--{0}--Time cost at HaplotypeCaller ----cost {1} min.'.format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # -Joint-Call Cohort
    if erc == 'GVCF':
        command_count1 = '{0} --java-options "{1}" GenotypeGVCFs -R {2} --variant {3} -O {4}'.format(
            gatk_dir, memory_size, ref_fa_file, vcf, vcf1)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count1)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
        store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', 'Compeleted GenotypeGVCFs by GATK.')
        store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                               '--{0}--Time cost at GenotypeGVCFs----cost {1} min.'.format(
                                   sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # -SNP
    snp_vcf = output + '/' + sample + '.raw_variants_SNP.vcf'
    command_count2 = ' '.join([gatk_dir,
                               '--java-options "{0}" SelectVariants'.format(memory_size),
                               '-R', ref_fa_file,
                               '--variant', vcf1,
                               '-O', snp_vcf,
                               '--select-type-to-include SNP',
                               '--showHidden true'])
    stdout, stderr = stdout_err(command_count2)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
    store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', 'Compeleted SelectVariants SNP by GATK.')
    store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           '--{0}--Time cost at SelectVariants SNP----cost {1} min.'.format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))
    # indel
    indel_vcf = output + '/' + sample + '.raw_variants_indel.vcf'
    command_count3 = ' '.join([gatk_dir,
                               '--java-options "{0}" SelectVariants'.format(memory_size),
                               '-R', ref_fa_file,
                               '--variant', vcf1,
                               '-O', indel_vcf,
                               '--select-type-to-include INDEL',
                               '--showHidden true'])
    stdout, stderr = stdout_err(command_count3)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
    store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', 'Compeleted SelectVariants indel by GATK.')
    store_germline_vc_logs(logger_g_variantcalling_process, 'null',
                           '--{0}--Time cost at SelectVariants indel----cost {1} min.'.format(
                               sample, str('%.3f' % ((time.time()-time_start1)/60))))


def strelka2_call(strelka2_dir, bgzip, tabix, total_ref_chrom_fa_file, output, sample, marked_bam,
                  exome_target_bed, logger_g_variantcalling_process, logger_g_variantcalling_errors, renew):
    os.system('cp {0} {1}'.format(exome_target_bed, output))
    if not os.path.exists(output+'/'+os.path.basename(exome_target_bed)+'.gz') or renew == 'T':
        command1 = bgzip + ' ' + output+'/'+os.path.basename(exome_target_bed)
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', command1)
        os.system(command1)
    if not os.path.exists(output+'/'+os.path.basename(exome_target_bed)+'.gz.tbi') or renew == 'T':
        command2 = tabix + ' ' + output+'/' + os.path.basename(exome_target_bed) + '.gz'
        store_germline_vc_logs(logger_g_variantcalling_process, 'null', command2)
        os.system(command2)
    command3 = 'python {0} --bam {1} --callRegions {2} --exome --referenceFasta {3} --runDir {4}'.format(
        strelka2_dir, marked_bam, output+'/'+os.path.basename(exome_target_bed)+'.gz', total_ref_chrom_fa_file, output)
    stdout, stderr = stdout_err(command3)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
    store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)
    run_py = output+'/' + 'runWorkflow.py'
    command4 = 'python {0} -m local -j 1'.format(run_py)
    stdout, stderr = stdout_err(command4)
    store_germline_vc_logs(logger_g_variantcalling_process, 'null', stdout)
    store_germline_vc_logs('null', logger_g_variantcalling_errors, stderr)


def samtools_call(samtools_dir, bcftools_dir, marked_bam, sample, output,
                  ref_fa_file, logger_g_variantcalling_process, logger_g_variantcalling_errors):
    bcf = output + '/' + sample + '_raw.bcf'
    command1 = bcftools_dir + ' mpileup -Ob -o ' + bcf + ' -f ' + ref_fa_file + ' ' + marked_bam
    os.system(command1)
    vcf = output + '/' + sample + '.raw_samtools.vcf'
    command2 = bcftools_dir + ' call -vmO v -o ' + vcf + ' ' + bcf
    os.system(command2)


def varsan2_call(samtools_dir, varsan2_dir, ref_fa_file, output,
                 sample, marked_bam, logger_g_variantcalling_process, logger_g_variantcalling_errors):
    pileup = output + '/' + sample + '.mpileup'
    command1 = samtools_dir + ' mpileup -B -f ' + ref_fa_file + ' ' + marked_bam + ' > ' + pileup
    os.system(command1)
    snp = output + '/' + sample + '.raw_varscan2.snp.vcf'
    indel = output + '/' + sample + '.raw_varscan2.indel.vcf'
    command2 = 'java -jar ' + varsan2_dir + ' mpileup2snp ' + pileup + ' > ' + snp
    os.system(command2)
    # stdout, stderr = stdout_err(command2)
    command3 = 'java -jar ' + varsan2_dir + ' mpileup2indel ' + pileup + ' > ' + indel
    os.system(command3)

def smcounter_call(smcounter, output, marked_bam, exome_target_bed, mtdepth, rpb, ncpu, minbq, minmq, hplen, mismatchthr,
              mtdrop, maxmt, primerdist, threshold, ref_fa_file,
              bedtandemrepeats, bedrepeatmaskersubset, bedtools_dir, logfile,
              logger_g_variantcalling_process, logger_g_variantcalling_errors):
    cmd = 'python2.7 ' + smcounter + \
          ' --outPrefix ' + output + \
          ' --bamFile ' + marked_bam + \
          ' --bedTarget ' + exome_target_bed + \
          ' --mtDepth ' + str(mtdepth) + \
          ' --rpb ' + str(rpb) + \
          ' --nCPU ' + str(ncpu) + \
          ' --minBQ ' + str(minbq) + \
          ' --minMQ ' + str(minmq) + \
          ' --hpLen ' + str(hplen) + \
          ' --mismatchThr ' + str(mismatchthr) + \
          ' --mtDrop ' + str(mtdrop) + \
          ' --maxMT ' + str(maxmt) + \
          ' --primerDist ' + str(primerdist) + \
          ' --threshold ' + str(threshold) + \
          ' --refGenome ' + ref_fa_file + \
          ' --bedTandemRepeats ' + bedtandemrepeats + \
          ' --bedRepeatMaskerSubset ' + bedrepeatmaskersubset + \
          ' --bedtoolsPath ' + bedtools_dir + \
          ' --runPath ./ ' + \
          ' --logFile ' + logfile
    #if not os.path.exists(output + '.smCounter.cut.txt'):
    os.system(cmd)