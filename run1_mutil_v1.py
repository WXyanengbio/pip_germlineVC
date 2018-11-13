__version__ = '1.0'
__author__ = 'Wang Xian'

import logging
import os
import re
import sys
import time
import argparse
import multiprocessing
import psutil
import yaml

# sys.path.append(os.path.split(os.path.realpath(__file__))[0])
# import the logging functions
from pipelines.log.log_v1 import store_pipeline_logs, store_trim_logs, store_filter_logs, store_align_logs, \
    store_cluster_logs, store_reformat_logs, store_germline_vc_logs, store_annotation_logs,\
    store_statistics_logs, store_benchmark_logs
# import the trim function
from pipelines.trim.trim_reads_v1 import trim_read_pairs
# import the align function
from pipelines.align.align_reads_v1 import align_reads_bwa, align_reads_bwa_based_all
# import the post-align functions
from pipelines.post_align.post_alignment_v1 import filter_alignment_samtools, identify_gs_primers
# import the barcode clusering function
from pipelines.cluster_barcode.umitools_v1 import umitool
# import the reformat sam function
from pipelines.reformat.reformat_sam import reformat_sam
# import the variant calling funcitons
from pipelines.variant_call.g_variantcall1_v1 import sam_to_bam, germline_variant_calling, strelka2_call, \
    samtools_call, varsan2_call, smcounter_call
# import the annotation variant funcitons
from pipelines.variant_call.annotation_gatk_hc_v1 import annotationmain, read_vcf_filter
# import the statistics functions
from pipelines.statistics.prestatistics_module_v1 import qc_raw_reads, statistics_depth_coverage, statistics_sam_bam,\
    statistics_time, merge_statistics_sam_bam, merge_statistics, statistics_mtdepth_coverage
# import the benchmaking funciton
# from pipelines.benchmark.hap_benchmark import hap_py


def script_information():
    print("\nApplication: multiprocessing of pipelines of QIAseq Targeted DNA Panel\n")
    print("=====================================================================")
    print("Required environment: python \ JAVA\ R \ bwa \ samtools \ GATK \ UMI-tools \ hay.py ")


parser = argparse.ArgumentParser(usage="\n\npython3.6 %(prog)s --yaml")

parser.add_argument("--yaml", 
                    type=str,
                    default='null',
                    help="The yaml file of the parameters")
parser.add_argument("-v", 
                    '-version', 
                    action='version',
                    version=' %(prog)s 1.0')

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


def main_run_germline_variant_calling(path_sampleID_sub):
    time_start1 = time.time()
    source = path_sampleID_sub[0].split('\t')[0]
    sample = path_sampleID_sub[0].split('\t')[1]
    if len(path_sampleID_sub[0].split('\t')) > 2:
        tailname = path_sampleID_sub[0].split('\t')[2]
        sample = sample + '_' + tailname
    # parameters
    (output, fastqc_dir, primers_file, exome_target_bed, min_read_len, common_seq1, common_seq2, num_threads, edit_dist,
     min_mapq, max_soft_clip, max_dist, memory_size, snp_filter, indel_filter, ref_ens, bwa_dir, samtools_dir,
     umitools_dir, gatk_dir, ref_index_name, ref_fa_file, total_ref_fa_file, total_ref_fa_dict, known_sites,
     erc, db_cosmic, db_clinvar, db_g1000, test_level, exome_target, calling, tabix, bgzip, bcftools_dir,
     varsan2_dir, strelka2_dir, total_ref_chrom_fa_file, datasets_dir,smcounter, mtdepth, rpb, ncpu, minbq, minmq, hplen, mismatchthr,
     mtdrop, maxmt, primerdist, bedtandemrepeats, bedrepeatmaskersubset, bedtools_dir, renew) = path_sampleID_sub[1:55]
    # check the output
    out_dir = output + '/' + sample
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # pipeline log file
    log_dir = out_dir + '/' + 'log'
    if not os.path.isdir(log_dir):
        try:
            os.makedirs(log_dir)
        except OSError as e:
            if e.errno == 17:
                logger_pipeline_process = log_dir
                logger_pipeline_errors = log_dir
    logger_pipeline_process = log_dir
    logger_pipeline_errors = log_dir
    # time cost
    time_start = time.time()
    module = "QC"

    read1 = source + '/' + sample + '_R1_001.fastq.gz'
    read2 = source + '/' + sample + '_R2_001.fastq.gz'

    #if tools in ['all', 'qc']:
    if 'qc' in tools or 'all' in tools:
        print("Test QC module!\n")
    # qc_dir
        qc_dir = out_dir + '/' + 'QC'
        if not os.path.exists(qc_dir):
            os.makedirs(qc_dir)

        logger_statistics_process = log_dir
        logger_statistics_errors = log_dir
        qc_result1, qc_result2 = qc_raw_reads(fastqc_dir, qc_dir,
                                              sample, module, read1, read2,
                                              logger_statistics_process, logger_statistics_errors)
    # check the quality of the raw reads
        if float(qc_result1[7].strip('%')) > 70 and float(qc_result2[7].strip('%')) > 70:
            print("The ratio of read1 and read2 with Q30 quality are both higher than 70%.")
        else:
            exit("The ratio of read1 and read2 with Q30 quality are both lower than 80%!!!!!!!")
        store_pipeline_logs(logger_pipeline_process, 'null',
                            "--{0}--QC of reads is completed after {1} min.\n".format(
                                sample, str('%.3f' % ((time.time()-time_start)/60))))
        store_statistics_logs(logger_statistics_process, 'null',
                              "--{0}--QC of reads is completed after {1} min.\n".format(
                                  sample, str('%.3f' % ((time.time()-time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    # trim
    ##########################################################################################
    # time cost
    time_start = time.time()
    # undetermined_dir
    undetermined_dir = out_dir + '/' + 'undetermined'
    trimmed1 = undetermined_dir + '/' + sample + '_R1_undetermined.fastq'
    trimmed2 = undetermined_dir + '/' + sample + '_R2_undetermined.fastq'
    stats_file = undetermined_dir + '/' + sample + '_basic_stats.txt'
    # if tools in ['all', 'trim']:
    if 'trim' in tools or 'all' in tools:
        print("please check the QC subprocess result--the min read length!")
        print("The cutoff of the min read length is the default: {0}".format(min_read_len))
        print("Test trim module!\n\n\n")
        # mkdir undetermined_dir
        if not os.path.exists(undetermined_dir):
            os.makedirs(undetermined_dir)
        
        logger_trim_process = log_dir
        logger_trim_errors = log_dir
        
        trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len,
                        common_seq1, common_seq2, stats_file, logger_trim_process,
                        logger_trim_errors)
        store_pipeline_logs(logger_pipeline_process, 'null',
                            "--{0}--Trim of reads is completed after {1} min.\n".format(
                                sample, str('%.3f' % ((time.time()-time_start)/60))))
        store_trim_logs(logger_trim_process, 'null',
                        "--{0}--Trimming of reads is completed after {1} min.\n".format(
                            sample, str('%.3f' % ((time.time()-time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    # align
    ##########################################################################################
    # time cost
    time_start = time.time()
    # aligned_dir
    aligned_dir = out_dir + '/' + 'aligned'
    trim_read1 = out_dir + '/' + 'undetermined' + '/' + sample + '_R1_undetermined.fastq'
    trim_read2 = out_dir + '/' + 'undetermined' + '/' + sample + '_R2_undetermined.fastq'

    out_file = aligned_dir + '/' + sample + '_aligned.sam'
    # if tools in ['all', 'align']:
    if 'align' in tools or 'all' in tools:
        print("please check the Trim subprocess result--undetermined.fastq!")
        print("Test align module!\n")
        if not os.path.exists(aligned_dir):
            os.makedirs(aligned_dir)
        logger_bwa_process = log_dir
        logger_bwa_errors = log_dir
        align_reads_bwa(bwa_dir, samtools_dir, ref_fa_file, ref_index_name,
                        exome_target_bed, total_ref_fa_file, trim_read1, trim_read2,
                        out_file, num_threads, logger_bwa_process, logger_bwa_errors, renew)
        store_align_logs(logger_bwa_process, 'null',
                         "--{0}--Alignment of reads is completed after {1} min.".format(
                             sample, str('%.3f' % ((time.time()-time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            "--{0}--Align of reads is completed after {1} min.".format(
                                sample, str('%.3f' % ((time.time()-time_start)/60))))
        print("--" * 20 + '\n\n')

    # #####################################annotationmain####################################################
    # post_align
    # #########################################################################################
    # time cost
    time_start = time.time()
    # post_aligned_dir
    filtered_dir = out_dir + '/' + 'filtered'
    # out_file from the align
    alignment_sam = out_file
    out_file1 = filtered_dir + '/' + sample + '_tmp.sam'
    stats_file = filtered_dir + '/' + sample + '_align_stats.txt'
    primer_stats_file = filtered_dir + '/' + sample + '_primer_stats.csv'
    out_file2 = filtered_dir + '/' + sample + '_filtered.sam'

    # if tools in ['all', 'post_align']:
    if 'post_align' in tools or 'all' in tools:
        print("please check the Algin subprocess result--aligned.sam!")
        print("Test post align module!\n")
        if not os.path.exists(filtered_dir):
            os.makedirs(filtered_dir)
        logger_filter_process = log_dir
        logger_filter_errors = log_dir
        filter_alignment_samtools(samtools_dir, alignment_sam, min_mapq,
                                  max_soft_clip, out_file1, stats_file,
                                  logger_filter_process, logger_filter_errors)
        identify_gs_primers(samtools_dir, out_file1, primers_file, max_dist, out_file2,
                            stats_file, primer_stats_file, logger_filter_process,
                            logger_filter_errors)
        store_filter_logs(logger_filter_process, 'null',
                          "--{0}--Post Alignment of reads is completed after {1} min.".format(
                              sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            "--{0}--Post_align of reads is completed after {1} min.".format(
                                sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    # barcode clustering
    ##########################################################################################
    # time cost
    time_start = time.time()
    # clustering_dir
    clustered_dir = out_dir + '/' + 'clustered'
    filtered_sam = out_dir + '/' + 'filtered' + '/' + sample + '_filtered.sam'
    filtered_bam = clustered_dir + '/' + sample + '_filtered.bam'
    sorted_bam = clustered_dir + '/' + sample + '_filtered_sorted.bam'
    umitool_stats = clustered_dir + '/' + sample + '_deduplicated'
    #umitool_stats = clustered_dir + '/' + sample + '_group.tsv'
    umis_sam = clustered_dir + '/' + sample + '_umis.sam'
    #if tools in ['all', 'cluster']:
    if 'cluster' in tools or 'all' in tools:
        print("please check the post algin subprocess result--filtered.sam!")
        print("Test cluster module!\n")
        if not os.path.exists(clustered_dir):
            os.makedirs(clustered_dir)
        logger_umi_process = log_dir
        logger_umi_errors = log_dir
        umitool(samtools_dir, umitools_dir, filtered_sam, filtered_bam, sorted_bam,
                umitool_stats, umis_sam, edit_dist, logger_umi_process, logger_umi_errors)
        store_cluster_logs(logger_umi_process, 'null',
                           "--{0}--UMIs tools clustering of reads is completed after {1} min.".format(
                               sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            "--{0}--Cluster of reads is completed after {1} min.".format(
                                sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    # reformat
    ##########################################################################################
    # time cost
    time_start = time.time()
    # reformated_dir
    reformated_dir = out_dir + '/' + 'reformated'
    alignment_sam = out_dir + '/' + 'clustered' + '/' + sample + '_umis.sam'
    output_sam = reformated_dir + '/' + sample + '_vcready.sam'
    # if tools in ['all', 'reformat']:
    if 'reformat' in tools or 'all' in tools:
        print("please check the cluster subprocess result--umis.sam!")
        print("Test reformat module!\n")
        if not os.path.exists(reformated_dir):
            os.makedirs(reformated_dir)
        logger_reformat_process = log_dir
        logger_reformat_errors = log_dir
        reformat_sam(alignment_sam, output_sam, logger_reformat_process, logger_reformat_errors)
        store_reformat_logs(logger_reformat_process, 'null',
                            '--{0}--Finish reformating alignment SAM file is completed after {1} min.'.format(
                                sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            '--{0}--Reformat alignment SAM file is completed after {1} min.'.format(
                                sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    # #########################################################################################
    # Germline variant calling
    # #########################################################################################
    # time cost
    time_start = time.time()
    # germline_vc_dir
    germline_vc_dir = out_dir + '/' + 'germline_vc'
    # modify the known-sites
    known_sites = known_sites.replace(',', ' --known-sites ' + datasets_dir + '/')
    known_sites = datasets_dir + '/' + known_sites
    vready_sam = out_dir + '/' + 'reformated' + '/' + sample + '_vcready.sam'
    # marked_bqsr_bam = germline_vc_dir + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
    if os.path.basename(exome_target_bed) != 'all':
        exon_interval = germline_vc_dir + '/' + 'target_interval.list'
    else:
        exon_interval = 'all'
    # if tools in ['all', 'variant_call']:
    if 'variant_call' in tools or 'all' in tools:
        print("please check the reformat subprocess result--vcready.sam!")
        print("Test variant_call module!\n")
        if not os.path.exists(germline_vc_dir):
            os.makedirs(germline_vc_dir)
        logger_germline_vc_process = log_dir
        logger_germline_vc_errors = log_dir

        bqsr = 'n'
        bam_to_variant, bqsr_bam_to_variant = sam_to_bam(gatk_dir, samtools_dir, vready_sam, sample,
                                                         germline_vc_dir, memory_size, exome_target_bed,
                                                         total_ref_fa_file, total_ref_fa_dict, known_sites,
                                                         logger_germline_vc_process, logger_germline_vc_errors, bqsr, renew)
        callings = calling.split(',')
        # if calling == 'GATK':
        if 'GATK' in callings:
            germline_variant_calling(gatk_dir, bam_to_variant, sample, germline_vc_dir, memory_size,
                                     total_ref_fa_file, exon_interval, erc, snp_filter, indel_filter,
                                     logger_germline_vc_process, logger_germline_vc_errors)
        # elif calling == 'strelka2':
        if 'strelka2' in callings:
            strelka2_call(strelka2_dir, bgzip, tabix, total_ref_chrom_fa_file, germline_vc_dir, sample, bam_to_variant,
                          exome_target_bed, logger_germline_vc_process, logger_germline_vc_errors, renew)
        # elif calling == 'samtools':
        if 'samtools' in callings:
            samtools_call(samtools_dir, bcftools_dir, bam_to_variant, sample, germline_vc_dir,
                          total_ref_fa_file, logger_germline_vc_process, logger_germline_vc_errors)
        # elif calling == 'varscan2':
        if 'varscan2' in callings:
            varsan2_call(samtools_dir, varsan2_dir, total_ref_fa_file, germline_vc_dir, sample,
                         bam_to_variant, logger_germline_vc_process, logger_germline_vc_errors)
        if 'smcounter' in callings:
            bam_to_variant = bam_to_variant.rstrip(".MarkDuplicates.RG.bam") + '.bam'
            threshold = 0
            logfile = germline_vc_dir + '/smcountlog'
            smcounter_call(smcounter, germline_vc_dir + '/' + sample, bam_to_variant, exome_target_bed, mtdepth,
                      rpb, ncpu, minbq, minmq, hplen, mismatchthr,
                      mtdrop, maxmt, primerdist, threshold, total_ref_fa_file,
                      bedtandemrepeats, bedrepeatmaskersubset, bedtools_dir, logfile,
                      logger_germline_vc_process, logger_germline_vc_errors, renew)
        store_germline_vc_logs(logger_germline_vc_process, 'null',
                               '--{0}--Germline variant calling is completed after {1} min.'.format(
                                   sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            '--{0}--variant_calling is completed after {1} min.'.format(
                                sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    # #########################################################################################
    # Annotation variant calling
    # #########################################################################################
    # time cost
    time_start = time.time()
    # Annotation dir
    annotation_dir = out_dir + '/' + 'annotation'
    raw_vcf = germline_vc_dir + '/' + sample + '.raw_variants.vcf'
    snp_vcf = germline_vc_dir + '/' + sample + '.raw_variants_SNP.vcf'
    indel_vcf = germline_vc_dir + '/' + sample + '.raw_variants_indel.vcf'
    # annotation
    # if tools in ['all', 'annotation']:
    if 'annotation' in tools or 'all' in tools:
        print("please check the variant_call subprocess result--VCF!")
        print("Test annotation module!\n")
    # Annotation dir
        if not os.path.exists(annotation_dir):
            os.makedirs(annotation_dir)
        logger_annotation_process = log_dir
        logger_annotation_errors = log_dir
        snp_limit, indel_limit = read_vcf_filter(snp_filter, indel_filter)
        callings = calling.split(',')
        for callingsub in callings:
            annotationmain(db_cosmic, db_clinvar, db_g1000, ref_ens, raw_vcf, sample, snp_limit, indel_limit,
                       annotation_dir, logger_annotation_process, logger_annotation_errors, callingsub)
        if 'GATK' in callings:
            callingsub = 'GATK'
            annotationmain(db_cosmic, db_clinvar, db_g1000, ref_ens, snp_vcf, sample, snp_limit, indel_limit,
                           annotation_dir, logger_annotation_process, logger_annotation_errors, callingsub)
            annotationmain(db_cosmic, db_clinvar, db_g1000, ref_ens, indel_vcf, sample, snp_limit, indel_limit,
                           annotation_dir, logger_annotation_process, logger_annotation_errors, callingsub)

        store_annotation_logs(logger_annotation_process, 'null',
                              '--{0}--Finish annotation variant  is completed after {1} min.'.format(
                                  sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            '--{0}--Annotation variant is completed after {1} min.'.format(
                                sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    # statistics of the variant calling pipeline
    ##########################################################################################
    # time cost
    time_start = time.time()
    # statistics dir
    statistics_dir = out_dir + '/' + 'statistics'
    if not os.path.exists(statistics_dir):
        os.makedirs(statistics_dir)
    # statistics the clean reads
    # statistics trim dir
    statistics_trim_dir = statistics_dir + '/' + 'trim_QC'
    if not os.path.exists(statistics_trim_dir):
        os.makedirs(statistics_trim_dir)
    # if tools in ['all', 'statis']:
    if 'statis' in tools or 'all' in tools:
        print("please check the others subprocess results!")
        print("Test statistics module!\n")
        #if tools == 'statis':
        logger_statistics_process = log_dir
        logger_statistics_errors = log_dir
    # statistics dir
        if not os.path.exists(statistics_dir):
            os.makedirs(statistics_dir)
        module = "Trim"
        trim_result1, trim_result2 = qc_raw_reads(fastqc_dir, statistics_trim_dir, sample, module,  trimmed1, trimmed2,
                                                  logger_statistics_process, logger_statistics_errors)
    # statistics the align
        module1 = "Align"
        align_sorted_bam = statistics_depth_coverage(samtools_dir, out_file, statistics_dir, sample, module1,
                                                     exome_target, exome_target_bed, logger_statistics_process,
                                                     logger_statistics_errors, renew)
        align_statistics = statistics_sam_bam(samtools_dir, align_sorted_bam, statistics_dir, sample, module1,
                                              logger_statistics_process,
                                              logger_statistics_errors, renew)
    # statistics the filter
    # cluster module would build the filter sorted bam, but it has been changed UMIs-tools
        module2 = "Fliter"
        filtered_sorted_bam = statistics_depth_coverage(samtools_dir, filtered_sam, statistics_dir, sample, module2,
                                                        exome_target, exome_target_bed, logger_statistics_process,
                                                        logger_statistics_errors, renew)
        fliter_statistics = statistics_sam_bam(samtools_dir, filtered_sorted_bam, statistics_dir, sample, module2,
                                               logger_statistics_process, logger_statistics_errors, renew)
    # statistics the umi-tools
        module3 = "Cluster_reformat"
        cr_sorted_bam = statistics_depth_coverage(samtools_dir, vready_sam, statistics_dir, sample, module3,
                                                  exome_target, exome_target_bed, logger_statistics_process,
                                                  logger_statistics_errors, renew)
        cr_statistics = statistics_sam_bam(samtools_dir, cr_sorted_bam, statistics_dir, sample, module3,
                                           logger_statistics_process, logger_statistics_errors, renew)
    # staistics the bases MT depth
        statistics_mtdepth_coverage(germline_vc_dir, statistics_dir, sample, exome_target,
                                    logger_statistics_process, logger_statistics_errors)
    # merge the sorted bam
        merge_statistics_sam_bam(logger_statistics_process, logger_statistics_errors, statistics_dir, sample,
                                 ','.join([module1, module2, module3]), align_statistics,
                                 fliter_statistics, cr_statistics)
        store_statistics_logs(logger_statistics_process, 'null',
                              '--{0}--Statistics is completed after {1} min.'.format(
                                  sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            '--{0}--Statistics is completed after {1} min.'.format(
                                sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',
                            '--{0}--Pipeline : {0} is completed after {1} min.'.format(
                                sample, ('%.2f' % ((time.time() - time_start1)/60))))
    # statistics the time cost
        process_log = out_dir + '/' + 'log' + '/' + 'process.log'
        statistics_time(statistics_dir, sample, process_log, logger_statistics_process, logger_statistics_errors)
        print("--" * 20 + '\n\n')
    return 'Pipeline : {0} is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start1)/60)))


cores = multiprocessing.cpu_count()
data = psutil.virtual_memory()
total = data.total
free = data.available


if __name__ == '__main__':
    # input
    yaml_file = args.yaml
    if yaml_file != 'null':
        f = open(yaml_file)
        file_yaml = yaml.load(f)
    if 'file_list' in file_yaml.keys():
        file_list = file_yaml['file_list']
    # check the outputdir
    if 'output' in file_yaml.keys():
        out_dir = file_yaml['output']
    if 'threads' in file_yaml.keys():
        threads = file_yaml['threads']
    if 'fastqc_dir' in file_yaml.keys():
        fastqc_dir = file_yaml['fastqc_dir']
    if 'datasets_dir' in file_yaml.keys():
        datasets_dir = file_yaml['datasets_dir']
    # trim
    if 'min_read_len' in file_yaml.keys():
        min_read_len = file_yaml['min_read_len']
    if 'common_seq1' in file_yaml.keys():
        common_seq1 = file_yaml['common_seq1']
    if 'common_seq2' in file_yaml.keys():
        common_seq2 = file_yaml['common_seq2']
    # bwa-align
    if 'num_threads' in file_yaml.keys():
        num_threads = file_yaml['num_threads']
    if 'bwa_dir' in file_yaml.keys():
        bwa_dir = file_yaml['bwa_dir']
    if 'ref_fa_file' in file_yaml.keys():
        ref_fa_file = datasets_dir + '/' + file_yaml['ref_fa_file']
    if 'ref_index_name' in file_yaml.keys():
        ref_index_name = datasets_dir + '/' + file_yaml['ref_index_name']
    # post-align
    if 'samtools_dir' in file_yaml.keys():
        samtools_dir = file_yaml['samtools_dir']
    if 'min_mapq' in file_yaml.keys():
        min_mapq = file_yaml['min_mapq']
    if 'max_soft_clip' in file_yaml.keys():
        max_soft_clip = file_yaml['max_soft_clip']
    if 'max_dist' in file_yaml.keys():
        max_dist = file_yaml['max_dist']
    if 'primers_file' in file_yaml.keys():
        primers_file = file_yaml['primers_file']
    # clustering
    if 'umitools_dir' in file_yaml.keys():
        umitools_dir = file_yaml['umitools_dir']
    if 'edit_dist' in file_yaml.keys():
        edit_dist = file_yaml['edit_dist']
    # variant calling
    if 'memory_size' in file_yaml.keys():
        memory_size = file_yaml['memory_size']
        memory_size = '-Xmx' + str(memory_size) + 'G ' + '-Djava.io.tmpdir=./'
    if 'gatk_dir' in file_yaml.keys():
        gatk_dir = file_yaml['gatk_dir']
    # strelka2
    if 'strelka2_dir' in file_yaml.keys():
        strelka2_dir = file_yaml['strelka2_dir']
    if 'bgzip' in file_yaml.keys():
        bgzip = file_yaml['bgzip']
    if 'tabix' in file_yaml.keys():
        tabix = file_yaml['tabix']
    if 'total_ref_chrom_fa_file' in file_yaml.keys():
        total_ref_chrom_fa_file = datasets_dir + '/' + file_yaml['total_ref_chrom_fa_file']
    # samtools
    if 'bcftools_dir' in file_yaml.keys():
        bcftools_dir = file_yaml['bcftools_dir']
    # varscan2
    if 'varsan2_dir' in file_yaml.keys():
        varsan2_dir = file_yaml['varsan2_dir']
    #smcounter
    smcounter = file_yaml['smcounter_path']
    mtdepth = file_yaml['mt_depth']
    rpb = file_yaml['rpb']
    ncpu = file_yaml['n_cpu']
    minbq = file_yaml['min_bq']
    minmq = file_yaml['min_mq']
    hplen = file_yaml['hp_len']
    mismatchthr = file_yaml['mismatch']
    mtdrop = file_yaml['mt_drop']
    maxmt = file_yaml['max_mt']
    primerdist = file_yaml['primerdist']
    bedtandemrepeats = file_yaml['bed_tandem_repeats']
    bedrepeatmaskersubset = file_yaml['bed_repeat_masker_subset']
    bedtools_dir = file_yaml['bedtools_dir']
    # ref database for variant calling
    if 'total_ref_fa_dict' in file_yaml.keys():
        total_ref_fa_dict = datasets_dir + '/' + file_yaml['total_ref_fa_dict']

    if 'total_ref_fa_file' in file_yaml.keys():
        total_ref_fa_file = datasets_dir + '/' + file_yaml['total_ref_fa_file']
    if 'known_sites' in file_yaml.keys():
        known_sites = file_yaml['known_sites']
    # read_length = args.read_length
    if 'exome_target' in file_yaml.keys():
        exome_target = datasets_dir + '/' + file_yaml['exome_target']
    if 'exome_target_bed' in file_yaml.keys():
        exome_target_bed = datasets_dir + '/' + file_yaml['exome_target_bed']
    # GATK Ha
    if 'erc' in file_yaml.keys():
        erc = file_yaml['erc']
    else:
        erc = args.erc
    # variant filters
    if 'snp_filter' in file_yaml.keys():
        snp_filter = file_yaml['snp_filter']
    if 'indel_filter' in file_yaml.keys():
        indel_filter = file_yaml['indel_filter']
    # annotation
    if 'ref_ens' in file_yaml.keys():
        ref_ens = datasets_dir + '/' + file_yaml['ref_ens']
    if 'db_cosmic' in file_yaml.keys():
        db_cosmic = datasets_dir + '/' + file_yaml['db_cosmic']
    if 'db_clinvar' in file_yaml.keys():
        db_clinvar = datasets_dir + '/' + file_yaml['db_clinvar']
    if 'db_g1000' in file_yaml.keys():
        db_g1000 = datasets_dir + '/' + file_yaml['db_g1000']
    # subprocess of the pipeline
    if 'tools' in file_yaml.keys():
        tools = file_yaml['tools']
    # the variant calling model
    if 'calling' in file_yaml.keys():
        calling = file_yaml['calling']
    # check the cpus and memery
    time_start_run = time.time()
    num_thread_of_process = int(cores/threads)
    if num_thread_of_process*threads > 0.7*cores:
        num_thread_of_process = num_thread_of_process - 1
    memory_size_of_process = int(free*0.6/(1024*1024*1024*num_thread_of_process))
    # check the outputdir
    # variant calling
    memory_size = memory_size_of_process
    if memory_size > 4:
        memory_size = 4
        memory_size = '-Xmx' + str(memory_size) + 'G '
    else:
        if threads > 1:
            memory_size = '-Xmx' + str(memory_size) + 'G '
        else:
            memory_size = '-Xmx4G '
    # if renew the sam to variant calling
    renew = file_yaml['renew']
    
    file_list1 = open(file_list, 'r')
    path_sampleID = []
    for line in file_list1.readlines():
        path_sampleID.append([line.strip(), out_dir, fastqc_dir, primers_file, exome_target_bed,
                              min_read_len, common_seq1, common_seq2, num_threads, edit_dist, min_mapq,
                              max_soft_clip, max_dist, memory_size, snp_filter, indel_filter, ref_ens,
                              bwa_dir, samtools_dir, umitools_dir, gatk_dir, ref_index_name,
                              ref_fa_file, total_ref_fa_file, total_ref_fa_dict, known_sites, erc, db_cosmic,
                              db_clinvar, db_g1000, tools, exome_target, calling, tabix, bgzip,
                              bcftools_dir, varsan2_dir,
                              strelka2_dir, total_ref_chrom_fa_file, datasets_dir, smcounter,
                              mtdepth, rpb, ncpu, minbq, minmq, hplen, mismatchthr, mtdrop, 
                              maxmt, primerdist, bedtandemrepeats, bedrepeatmaskersubset, bedtools_dir, renew])
    
    pool = multiprocessing.Pool(processes=threads)
    results = pool.map(main_run_germline_variant_calling, path_sampleID)
    print("--" * 20)
    pool.close()
    pool.join()
    
    # if tools in ['all', 'statis', 'summary']:
    if 'all' in tools or 'summary' in tools:
        print("please check the others subprocess results!")
        print("Test statistics module!\n")
        # -merge the statis info
        sample_preinfo = out_dir + '/pre_summary.txt'
        qc = out_dir + '/QC.statistics.list'
        trim_qc = out_dir + '/' + 'QC_Trim.statistics.txt'
        trim_statis = out_dir + '/trim.basic_stats.txt'
        filter_statis = out_dir + '/align.basic_stats.txt'
        primer_statis = out_dir + '/primer.basic_stats.txt'
        umi_statis = out_dir + '/umis.basic_stats.txt'
        align_base = out_dir + '/basesDepthInRegion.stats.txt'
        mt_depth = out_dir + '/bases_MT_DepthInRegion.stats.txt'
        
        for line in range(0, len(path_sampleID)):
            parameters = path_sampleID[line][0]
            sample = parameters.split('\t')[1]
            if line is 0:
                os.system('ls {0}/{1}/QC/*QC.statistics.txt > {2}'.format(out_dir, sample, qc))
                os.system('ls {0}/{1}/statistics/trim_QC/*Trim.statistics.txt > {2}'.format(out_dir, sample, trim_qc))
                os.system('ls {0}/{1}/undetermined/*_basic_stats.txt > {2}'.format(out_dir, sample, trim_statis))
                os.system('ls {0}/{1}/filtered/*_align_stats.txt > {2}'.format(out_dir, sample, filter_statis))
                os.system('ls {0}/{1}/filtered/*_primer_stats.csv > {2}'.format(out_dir, sample, primer_statis))
                os.system('ls {0}/{1}/clustered/*_deduplicated_per_umi.tsv > {2}'.format(out_dir, sample, umi_statis))
                os.system('ls {0}/{1}/statistics/*_Fliter_basesDepthInRegion.txt > {2}'.format(out_dir, sample, align_base))
                os.system('ls {0}/{1}/statistics/*MTDepthInTargetExon.txt > {2}'.format(out_dir, sample, mt_depth))
            else:
                os.system('ls {0}/{1}/QC/*QC.statistics.txt >> {2}'.format(out_dir, sample, qc))
                os.system('ls {0}/{1}/statistics/trim_QC/*Trim.statistics.txt >> {2}'.format(out_dir, sample, trim_qc))
                os.system('ls {0}/{1}/undetermined/*_basic_stats.txt >> {2}'.format(out_dir, sample, trim_statis))
                os.system('ls {0}/{1}/filtered/*_align_stats.txt >> {2}'.format(out_dir, sample, filter_statis))
                os.system('ls {0}/{1}/filtered/*_primer_stats.csv >> {2}'.format(out_dir, sample, primer_statis))
                os.system('ls {0}/{1}/clustered/*_deduplicated_per_umi.tsv >> {2}'.format(out_dir, sample, umi_statis))
                os.system('ls {0}/{1}/statistics/*_Fliter_basesDepthInRegion.txt >> {2}'.format(out_dir, sample, align_base))
                os.system('ls {0}/{1}/statistics/*MTDepthInTargetExon.txt >> {2}'.format(out_dir, sample, mt_depth))
        
        merge_statistics(sample_preinfo, exome_target, qc, trim_qc, trim_statis,
                         filter_statis, primer_statis, umi_statis, align_base, mt_depth)

    print("All process done.")
    print("Return results: ")
    # pipeline log file
    log_dir_run = out_dir + '/' + 'log/'
    if not os.path.exists(log_dir_run):
        os.makedirs(log_dir_run)
    logger_pipeline_process_run = log_dir_run
    logger_pipeline_error_run = log_dir_run
    for i in results:
        print(i)
        store_pipeline_logs(logger_pipeline_process_run, 'null', i+'\n')
    store_pipeline_logs(logger_pipeline_process_run, 'null',
                        'All samples are completed after {0} min.\n'.format(
                            ('%.2f' % ((time.time() - time_start_run)/60))))
