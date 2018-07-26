__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
import re
import sys
import time
import argparse
import multiprocessing
import psutil

#sys.path.append(os.path.split(os.path.realpath(__file__))[0])
#import the logging functions
from pipelines.log.log import store_pipeline_logs , store_trim_logs, store_filter_logs, store_align_logs , store_cluster_logs , \
                              store_reformat_logs, store_germline_vc_logs, store_annotation_logs, store_statistics_logs, store_benchmark_logs
#import the trim function
from pipelines.trim.trim_reads import trim_read_pairs
#import the align function
from pipelines.align.align_reads import align_reads_bwa
#import the post-align functions
from pipelines.post_align.post_alignment import filter_alignment_samtools , identify_gs_primers
#import the barcode clusering function
from pipelines.cluster_barcode.umitools import umitool
#import the reformat sam function
from pipelines.reformat.reformat_sam import reformat_sam
#import the variant calling funcitons
from pipelines.variant_call.g_variantcall1 import sam_to_bem , germline_variant_calling
#import the annotation variant funcitons
from pipelines.variant_call.annotation_gatk_hc import annotationmain
#import the statistics functions
from pipelines.statistics.prestatistics_module import qc_raw_reads, statistics_depth_coverage, statistics_sam_bam, statistics_time, merge_statistics_sam_bam
#import the benchmaking funciton
from pipelines.benchmark.hap_benchmark import hap_py

def script_information():
    print ("\nApplication: multiprocessing of pipelines of QIAseq Targeted DNA Panel\n")
    print ("=====================================================================")
    print ("Required environment: python \ JAVA\ R \ bwa \ samtools \ GATK \ UMI-tools \ hay.py ")

parser = argparse.ArgumentParser(usage = "\n\npython3.6 %(prog)s --file_list --threads --output")

parser.add_argument("--file_list", 
                    help = "Path to file contian tow colown: one is the path of raw reads, two is the sampleID of raw reads", 
                    type = str)
parser.add_argument("--threads", 
                    help = "the threads of multiprocessing ", 
                    type = int)
#parser.add_argument("--source", 
#                    help = "Path to input reads in FASTA format", 
#                    type = str)
#parser.add_argument("--sample_name", 
#                    help = "the sample name of raw reads", 
#                    type = str)
#parser.add_argument("--tailname", 
#                    type =str,
#                    default = 'null',
#                    help = "the tailname of sample raw read")
parser.add_argument("--datasets_dir", 
                    type =str,
                    default = '/home/dell/Works/Projects/Datasets',
                    help = "the tailname of sample raw read")
parser.add_argument("--common_seq1",
                    type = str, 
                    default = 'CAAAACGCAATACTGTACATT',
                    help = "the common seq1 of QIAseq Targeted DNA Panel")
parser.add_argument("--common_seq2",
                    type = str, 
                    default = 'ATTGGAGTCCT',
                    help = "the common seq2 of QIAseq Targeted DNA Panel")
parser.add_argument("--output", 
                    type = str,
                    default = 'null',
                    help = "Path of output file")
parser.add_argument("--fastqc_dir", 
                    type = str,
                    default= 'fastqc',
                    help = "the install path of fastQC")
parser.add_argument("--bwa_dir",
                    type = str,
                    default= 'bwa',
                    help = "the install path of bwa"
                    )
parser.add_argument("--samtools_dir", 
                    type = str,
                    default= 'samtools',
                    help = "the install path of samtools"
                    )
parser.add_argument("--gatk_dir", 
                    type = str,
                    default= 'gatk',
                    help = "the install path of GATK4"
                    )
parser.add_argument("--umitools_dir", 
                    type = str,
                    default= '/home/dell/.local/bin/umi_tools',
                    help = "the install path of umitools"
                    )
parser.add_argument("--benchmark_dir", 
                    type = str, 
                    default = '/home/dell/Works/Softwares/bioapps/benchmarking-tools/tools/hap.py_build/bin/hap.py',
                    help = "the install path of benchmark"
                    )
parser.add_argument("--ref_index_name", 
                    type = str,
                    default= 'genome/target_breast/target_breast_BRCA',
                    help = "the path of ref index--if there isn't a ref index, it will make a index in the path of ref fasta by bwa"
                    )
parser.add_argument("--ref_fa_file",
                    type = str,
                    default= 'genome/target_breast/target_breast_BRCA.refSeq.fa',
                    help = "the path of ref fasta"
                    )
parser.add_argument("--total_ref_fa_file", 
                    type = str,
                    default= 'genome/ucsc.hg19.fasta',
                    help = "the path of ref total ref fa"
                    )
parser.add_argument("--total_ref_fa_dict", 
                    type = str,
                    default= 'genome/ucsc.hg19.dict',
                    help = "the path of ref total ref fa dict"
                    )
parser.add_argument("--truth_vcf", 
                    type = str,
                    default= 'truthDB/samll_variant/NA12878.vcf.gz',
                    help = "the path of truth VCF of variant for benchmarking"
                    )
parser.add_argument("--confident_region_bed", 
                    type = str,
                    default= 'truthDB/confidentregion/ConfidentRegions.bed',
                    help = "the path of confident region bed file of the truth VCF of variant for benchmarking"
                    )
parser.add_argument("--min_read_len", 
                    type = int, 
                    default = 40,
                    help = "the cutoff of the min read length"
                   )
parser.add_argument("--num_threads",
                    type = int, 
                    default = 4,
                    help = "the number of threads to align"
                    )
parser.add_argument("--min_mapq",
                    type = int, 
                    default = 17,
                    help = "the parameter of filter alignment_sam"
                    )
parser.add_argument("--max_soft_clip",
                    type = int,
                    default = 10,
                    help = "the parameter of filter alignment_sam"
                    )
parser.add_argument("--max_dist",
                    type = int, 
                    default = 2,
                    help = "the parameter of filter alignment_sam"
                    )
parser.add_argument("--primers_file",
                    type = str,
                    default = 'DHS-001Z_primers_target.csv',
                    help = "Load all primer sequences in the panel")
parser.add_argument("--edit_dist",
                    type = int, 
                    default = 2,
                    help = "the parameter of edit distance between barcodes"
                    )
parser.add_argument("--memory_size", 
                    type = str, 
                    default = '4',
                    help = "the cutoff of Java memory"
                    )
parser.add_argument("--known_sites", 
                    type = str,
                    default ='known_sites/1000G_phase1.snps.high_confidence.hg19.sites.vcf,known_sites/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf,known_sites/dbsnp_138.hg19.vcf',
                    help = "the list of --known-sites , sep by: , "
                   )
                    
parser.add_argument("--exome_target_bed",
                    type = str,
                    default = 'genome/target_breast/target_breast_BRCA.bed',
                    help = "the bed file of exome intervals") 
parser.add_argument("--erc", 
                    type = str, 
                    default = 'GVCF',
                    help = "switch to running HaplotypeCaller in GVCF mode"
                    )
parser.add_argument("--read_filter", 
                    type = str, 
                    default = 'no',
                    help = "add a read filter that deals with some problems"
                    )
parser.add_argument("--snp_filter", 
                    type = str, 
                    default = 'QD < 5.0 || FS > 60.0 || MQ < 50.0 || SOR > 3.0 || MQRankSum < -2.5 || ReadPosRankSum < -3.0',
                    help = "add parameters for filtering SNPs"
                    )
parser.add_argument("--indel_filter",
                    type = str, 
                    default = 'QD < 5.0 || FS > 200 || ReadPosRankSum < -3.0 || SOR > 10.0',
                    help = "add parameters for filtering Indels"
                    )
parser.add_argument("--db_cosmic", 
                    type = str,
                    default ='annonDB/COSMIC_variant.csv',
                    help = "add cosmic databases of variants")
parser.add_argument("--db_clinvar", 
                    type = str,
                    default ='annonDB/breast_cancer_variant_clinvar.csv',
                    help = "add clinvar databases of variants")
parser.add_argument("--db_g1000", 
                    type = str,
                    default ='annonDB/breast_cancer_variant_1000genomes.csv',
                    help = "add g1000 databases of variants")
parser.add_argument("--anno_geneID", 
                    type = str,
                    default ='annonDB/geneid_cancer_qiagen.csv',
                    help = "add annotation gene ID of variants")
parser.add_argument("-v", 
                    '-version', 
                    action = 'version', 
                    version =' %(prog)s 1.0')
parser.add_argument("--test",
                    type = int,
                    default = 9,
                    help = "the subprocess of the script")

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
    #---
    output = path_sampleID_sub[1]
    fastqc_dir = path_sampleID_sub[2]
    primers_file = source + '/' + path_sampleID_sub[3]
    exome_target_bed = path_sampleID_sub[4]
    min_read_len = path_sampleID_sub[5]
    common_seq1 = path_sampleID_sub[6]
    common_seq2 = path_sampleID_sub[7]
    num_threads = path_sampleID_sub[8]
    edit_dist = path_sampleID_sub[9]
    min_mapq = path_sampleID_sub[10]
    max_soft_clip = path_sampleID_sub[11]
    max_dist = path_sampleID_sub[12]
    memory_size = path_sampleID_sub[13]
    snp_filter = path_sampleID_sub[14]
    indel_filter = path_sampleID_sub[15]
    ref_ens = path_sampleID_sub[16]
    bwa_dir = path_sampleID_sub[17]
    samtools_dir = path_sampleID_sub[18]
    umitools_dir = path_sampleID_sub[19]
    gatk_dir = path_sampleID_sub[20]
    benchmark_dir = path_sampleID_sub[21]
    ref_index_name = path_sampleID_sub[22]
    ref_fa_file = path_sampleID_sub[23]
    total_ref_fa_file = path_sampleID_sub[24]
    total_ref_fa_dict = path_sampleID_sub[25]
    known_sites = path_sampleID_sub[26]
    erc = path_sampleID_sub[27]
    db_cosmic = path_sampleID_sub[28]
    db_clinvar = path_sampleID_sub[29]
    db_g1000 = path_sampleID_sub[30]
    truth_vcf = path_sampleID_sub[31]
    confident_region_bed = path_sampleID_sub[32]
    test_level = path_sampleID_sub[33]
    #---check the output
    out_dir1 = output  + '/' + sample
    if not os.path.exists(out_dir1):
        os.makedirs(out_dir1)
    #----pipeline log file
    log_dir = out_dir1 + '/' + 'log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    logger_pipeline_process, logger_pipeline_errors = store_pipeline_logs(log_dir)
    ##########################################################################################
    #---QC
    ##########################################################################################
    #time cost
    time_start = time.time()
    #qc_dir
    module = "QC"
    qc_dir = out_dir1 + '/'+ 'QC'
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)

    read1 = source + '/' + sample + '_R1_001.fastq.gz'
    read2 = source + '/' + sample  + '_R2_001.fastq.gz'
    #---
    logger_statistics_process, logger_statistics_errors = store_statistics_logs(log_dir)
    qc_result1, qc_result2 = qc_raw_reads(fastqc_dir, qc_dir, 
                 sample, module, 
                 read1, read2,
                 logger_statistics_process, logger_statistics_errors)
    #--check the quality of the raw reads
    if float(qc_result1[6].strip('%')) > 70 and  float(qc_result2[6].strip('%')) > 70:
        print("The ratio of read1 and read2 with Q30 quality are both higher than 70%.")
    else:
        exit("The ratio of read1 and read2 with Q30 quality are both lower than 80%!!!!!!!")
    #--statistics the N base in raw reads and set the cutoff of the min read length
    if max(int(qc_result1[8]), int(qc_result2[8])) < min_read_len:
        print("The cutoff of the min read length is the default: {0}".format(min_read_len))
    else:
        min_read_len = max(int(qc_result1[8]), int(qc_result2[8]))
        print("The cutoff of the min read length is based on the N base in the reads: {0}".format(min_read_len))

    logger_statistics_process.info("QC of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("QC of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level >= 0:
        print("Test QC module!\n\n\n")
    else:
        exit()
    ##########################################################################################
    #---trim
    ##########################################################################################
    #time cost
    time_start = time.time()
    #undetermined_dir
    undetermined_dir = out_dir1 + '/'+ 'undetermined'
    if not os.path.exists(undetermined_dir):
        os.makedirs(undetermined_dir)

    trimmed1 = undetermined_dir + '/' + sample + '_R1_undetermined.fastq'
    trimmed2 = undetermined_dir + '/' + sample + '_R2_undetermined.fastq'
    stats_file = undetermined_dir + '/' + sample + '_basic_stats.txt'

    logger_trim_process, logger_trim_errors = store_trim_logs(log_dir)
    trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len,
                           common_seq1, common_seq2, stats_file, logger_trim_process,
                           logger_trim_errors)
    
    logger_trim_process.info("Trimming of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("Trim of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level >= 1:
        print("Test trim module!\n\n\n")
    else:
        exit()
    ##########################################################################################
    #---align
    ##########################################################################################
    #time cost
    time_start = time.time()
    #aligned_dir
    aligned_dir = out_dir1 + '/'+ 'aligned'
    if not os.path.exists(aligned_dir):
        os.makedirs(aligned_dir)

    trim_read1 = trimmed1
    trim_read2 = trimmed2
    out_file = aligned_dir + '/' + sample + '_aligned.sam'

    logger_bwa_process, logger_bwa_errors = store_align_logs(log_dir)
    
    returncode = align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, exome_target_bed, trim_read1, trim_read2, 
                                                out_file, num_threads, logger_bwa_process, logger_bwa_errors)
    
    logger_bwa_process.info("Alignment of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("Align of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level >= 2:
        print("Test align module!\n\n\n")
    else:
        exit()

    ######################################annotationmain####################################################
    #---post_align
    ##########################################################################################
    #time cost
    time_start = time.time()
    #post_aligned_dir
    filtered_dir = out_dir1 + '/'+ 'filtered'
    if not os.path.exists(filtered_dir):
        os.makedirs(filtered_dir)
    
    logger_filter_process, logger_filter_errors = store_filter_logs(log_dir)
    #out_file from the align
    alignment_sam = out_file
    #primers_file = primers_file
    #min_mapq=17
    #max_soft_clip=10
    out_file1 = filtered_dir + '/' + sample + '_tmp.sam'
    stats_file = filtered_dir + '/' + sample+ '_align_stats.txt'
    primer_stats_file = filtered_dir + '/' + sample + '_primer_stats.csv'
    #max_dist = 2
    out_file2 = filtered_dir + '/' + sample + '_filtered.sam'
    filter_alignment_samtools(samtools_dir, alignment_sam, min_mapq,
                              max_soft_clip, out_file1, stats_file,
                              logger_filter_process, logger_filter_errors)
    identify_gs_primers(samtools_dir, out_file1, primers_file, max_dist, out_file2,
                        stats_file, primer_stats_file, logger_filter_process,
                        logger_filter_errors)
    
    logger_bwa_process.info("Post Alignment of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("Post_align of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level >= 3:
        print("Test psot align module!\n\n\n")
    else:
        exit()
    ##########################################################################################
    #---barcode clustering
    ##########################################################################################
    #time cost
    time_start = time.time()
    #clustering_dir
    clustered_dir = out_dir1 + '/'+ 'clustered'
    if not os.path.exists(clustered_dir):
        os.makedirs(clustered_dir)

    logger_umi_process, logger_umi_errors = store_cluster_logs(log_dir)
    
    filtered_sam = out_file2
    filtered_bam = clustered_dir + '/' + sample + '_filtered.bam'
    sorted_bam = clustered_dir + '/' + sample + '_filtered_sorted.bam'
    umitool_stats = clustered_dir + '/' + sample + '_deduplicated'
    umis_sam = clustered_dir + '/' + sample + '_umis.sam'
    
    umitool(samtools_dir, umitools_dir, filtered_sam ,filtered_bam , sorted_bam, umitool_stats , umis_sam, edit_dist, logger_umi_process, logger_umi_errors)
    logger_umi_process.info("UMIs tools clustering of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("Cluster of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level >= 4:
        print("Test barcode clustering module!\n\n\n")
    else:
        exit()
    ##########################################################################################
    #---reformat 
    ##########################################################################################
    #time cost
    time_start = time.time()
    #reformated_dir
    reformated_dir = out_dir1 + '/'+ 'reformated'
    if not os.path.exists(reformated_dir):
        os.makedirs(reformated_dir)

    logger_reformat_process, logger_reformat_errors = store_cluster_logs(log_dir)
    
    alignment_sam = umis_sam
    output_sam = reformated_dir + '/' + sample + '_vcready.sam'
    reformat_sam(alignment_sam, output_sam, logger_reformat_process, logger_reformat_errors)
    
    logger_reformat_process.info('Finish reformating alignment SAM file is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Reformat alignment SAM file is completed after %.2f min.',(time.time() - time_start)/60)
    
    if test_level >= 5:
        print("Test reformat sam module!\n\n\n")
    else:
        exit()
    ##########################################################################################
    #---Germline variant calling
    ##########################################################################################
    #time cost
    time_start = time.time()
    #germline_vc_dir
    germline_vc_dir = out_dir1 + '/'+ 'germline_vc'
    if not os.path.exists(germline_vc_dir):
        os.makedirs(germline_vc_dir)
    
    logger_germline_vc_process, logger_germline_vc_errors = store_germline_vc_logs(log_dir)
    
    #---modify the known-sites
    known_sites = known_sites.replace(',' , ' --known-sites ' + args.datasets_dir + '/'  )
    known_sites = args.datasets_dir + '/' + known_sites
    vready_sam = output_sam
    sam_to_bem(gatk_dir, samtools_dir,
               vready_sam, sample,
               germline_vc_dir, memory_size,
               exome_target_bed, 
               total_ref_fa_file, total_ref_fa_dict,
               known_sites,
               logger_germline_vc_process, logger_germline_vc_errors)

    marked_bqsr_bam = germline_vc_dir + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
    exon_interval = germline_vc_dir + '/' + 'target_interval.list'
    germline_variant_calling(gatk_dir, marked_bqsr_bam,
                             sample, germline_vc_dir, 
                             memory_size, total_ref_fa_file, 
                             exon_interval, erc,
                             read_filter,
                             snp_filter,indel_filter,
                             logger_germline_vc_process, logger_germline_vc_errors)
    logger_germline_vc_process.info('Germline variant calling is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('variant_calling is completed after %.2f min.',(time.time() - time_start)/60)

    if test_level >= 6:
        print("Test variant calling module!\n\n\n")
    else:
        exit()
    
    ##########################################################################################
    #---Annotation variant calling
    ##########################################################################################
    #time cost
    time_start = time.time()
    #Annotation dir
    annotation_dir = out_dir1 + '/'+ 'annotation'
    if not os.path.exists(annotation_dir):
        os.makedirs(annotation_dir)
    
    logger_annotation_process, logger_annotation_errors = store_annotation_logs(log_dir)
    
    snp_vcf = germline_vc_dir + '/'  + sample + '.raw_variants_SNP.vcf'
    filter_snp = germline_vc_dir + '/'  + sample + '.filter_SNP.vcf'
    indel_vcf = germline_vc_dir + '/'  + sample + '.raw_variants_indel.vcf'
    filter_indel = germline_vc_dir + '/'  + sample + '.filter_indel.vcf'
    #annotation
    annotationmain(db_cosmic, db_clinvar, db_g1000, 
                   ref_ens,
                   snp_vcf, sample,
                   annotation_dir, logger_annotation_process, logger_annotation_errors)
    annotationmain(db_cosmic, db_clinvar, db_g1000, 
                   ref_ens,
                   filter_snp, sample,
                   annotation_dir, logger_annotation_process, logger_annotation_errors)
    annotationmain(db_cosmic, db_clinvar, db_g1000, 
                   ref_ens,
                   indel_vcf, sample,
                   annotation_dir, logger_annotation_process, logger_annotation_errors)
    annotationmain(db_cosmic, db_clinvar, db_g1000, 
                   ref_ens,
                   filter_indel, sample,
                   annotation_dir, logger_annotation_process, logger_annotation_errors)
    logger_annotation_process.info('Finish annotation variant  is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Annotation variant is completed after %.2f min.',(time.time() - time_start)/60)
    if test_level >= 7:
        print("Test annotation variant module!\n\n\n")
    else:
        exit()
    ##########################################################################################
    #---benchmarking 
    ##########################################################################################
    #time cost
    time_start = time.time()
    #benchmarking_dir
    benchmarking_dir = out_dir1 + '/'+ 'benchmarking'
    if not os.path.exists(benchmarking_dir):
        os.makedirs(benchmarking_dir)

    logger_benchmark_process, logger_benchmark_errors = store_benchmark_logs(log_dir)
    
    hap_py(benchmark_dir,truth_vcf, filter_snp, confident_region_bed, benchmarking_dir, total_ref_fa_file, exon_interval, num_threads, logger_benchmark_process, logger_benchmark_errors)
    hap_py(benchmark_dir,truth_vcf, filter_indel, confident_region_bed, benchmarking_dir, total_ref_fa_file, exon_interval,num_threads, logger_benchmark_process, logger_benchmark_errors)
    
    logger_benchmark_process.info('benchmarking is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Benchmark is completed after %.2f min.',(time.time() - time_start)/60)
    
    if test_level >= 8:
        print("Test benchmark module!")
    else:
        exit()
    ##########################################################################################
    #---statistics of the variant calling pipeline
    ##########################################################################################
    #time cost
    time_start = time.time()
    #statistics dir
    statistics_dir = out_dir1 + '/'+ 'statistics'
    if not os.path.exists(statistics_dir):
        os.makedirs(statistics_dir)
    #-statistics the clean reads
    #statistics trim dir
    statistics_trim_dir = statistics_dir + '/'+ 'trim_QC'
    if not os.path.exists(statistics_trim_dir):
        os.makedirs(statistics_trim_dir)
    module = "Trim"
    trim_result1, trim_result2 = qc_raw_reads(fastqc_dir, statistics_trim_dir, 
                 sample, module, 
                 trimmed1, trimmed2,
                 logger_statistics_process, logger_statistics_errors)
    #--statistics the align
    module1 = "Align"
    align_sorted_bam = statistics_depth_coverage(samtools_dir, out_file, statistics_dir, sample, module1,exome_target_bed, logger_statistics_process, logger_statistics_errors)
    align_statistics = statistics_sam_bam(samtools_dir, sorted_bam, statistics_dir,sample, module1, logger_statistics_process, logger_statistics_errors)
    #--statistics the filter
    #----cluster module would build the filter sorted bam, but it has been changed UMIs-tools
    module2 = "Fliter"
    filtered_sorted_bam = statistics_depth_coverage(samtools_dir, filtered_sam, statistics_dir, sample, module2,exome_target_bed, logger_statistics_process, logger_statistics_errors)
    fliter_statistics = statistics_sam_bam(samtools_dir, filtered_sorted_bam, statistics_dir, sample, module2, logger_statistics_process, logger_statistics_errors)
    # statistics the umi-tools
    module3 = "Cluster_reformat"
    cr_sorted_bam = statistics_depth_coverage(samtools_dir, vready_sam, statistics_dir, sample, module3, exome_target_bed,logger_statistics_process, logger_statistics_errors)
    cr_statistics = statistics_sam_bam(samtools_dir, cr_sorted_bam, statistics_dir, sample, module3, logger_statistics_process, logger_statistics_errors)
    #-merge the sorted bam
    merge_statistics_sam_bam(logger_statistics_process, logger_statistics_errors, statistics_dir, sample, ','.join([module1,module2,module3]),align_statistics,fliter_statistics,cr_statistics)
    
    logger_statistics_process.info('Statistics is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Statistics is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Pipeline : {0} is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start1)/60))))
    #--statistics the time cost
    process_log = out_dir + '/'+ 'log' + '/' + 'process.log'
    statistics_time(statistics_dir, sample, process_log, logger_statistics_process, logger_statistics_errors)
    #---
    if test_level == 9:
        print("Test statistics module!\n\n\n")
    return '{0} is complete'.format(sample)

#--
cores = multiprocessing.cpu_count()
data = psutil.virtual_memory()
total = data.total #总内存,单位为byte
free = data.available #可以内存

if __name__ == '__main__':
    #---check the cpus and memery
    threads = args.threads
    num_thread_of_process = int(cores/threads)
    if num_thread_of_process*threads > 0.75*cores:
        num_thread_of_process = num_thread_of_process -1
    memory_size_of_process = int(free*0.75/(1024*1024*1024*num_thread_of_process))
    #---check the outputdir
    out_dir = args.output
    if args.output is 'null':
        out_dir = "output"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #----pipeline log file
    #log_dir = out_dir + '/' + 'log/'
    #if not os.path.exists(log_dir):
    #    os.makedirs(log_dir)

    #logger_pipeline_process, logger_pipeline_errors = store_pipeline_logs(log_dir)
    #QC
    fastqc_dir = args.fastqc_dir
    #trim
    min_read_len = args.min_read_len
    common_seq1 = args.common_seq1
    common_seq2 = args.common_seq2
    #bwa---align
    num_threads = num_thread_of_process
    bwa_dir = args.bwa_dir
    ref_fa_file = args.datasets_dir + '/' + args.ref_fa_file
    ref_index_name = args.datasets_dir + '/' +args.ref_index_name
    #post-align
    samtools_dir = args.samtools_dir
    min_mapq = args.min_mapq
    max_soft_clip = args.max_soft_clip
    max_dist = args.max_dist
    primers_file = args.primers_file
    #--clustering
    umitools_dir = args.umitools_dir
    edit_dist = args.edit_dist
    #--variant calling
    memory_size = memory_size_of_process
    if memory_size > 4:
        memory_size = 4
        memory_size = '-Xmx' + str(memory_size) + 'G '
    print(memory_size)
    gatk_dir = args.gatk_dir
    samtools_dir = args.samtools_dir
    total_ref_fa_dict = args.datasets_dir + '/' + args.total_ref_fa_dict
    total_ref_fa_file = args.datasets_dir + '/' + args.total_ref_fa_file
    known_sites = args.known_sites
    #read_length = args.read_length
    exome_target_bed = args.datasets_dir + '/' + args.exome_target_bed
    erc = args.erc
    read_filter = args.read_filter
    snp_filter = args.snp_filter
    indel_filter = args.indel_filter
    #--annotation
    ref_ens = args.datasets_dir + '/' + args.anno_geneID
    db_cosmic = args.datasets_dir + '/' + args.db_cosmic
    db_clinvar = args.datasets_dir + '/' + args.db_clinvar
    db_g1000 =  args.datasets_dir + '/' + args.db_g1000
    #--benchmark
    benchmark_dir = args.benchmark_dir
    confident_region_bed = args.datasets_dir + '/' + args.confident_region_bed
    truth_vcf = args.datasets_dir + '/' + args.truth_vcf
    #---
    test_level = args.test
    #-----
    #--read the sample list 
    file_list = open(args.file_list,'r')
    path_sampleID = []
    for line in file_list.readlines():
        path_sampleID.append([line.strip(), out_dir, fastqc_dir, primers_file, exome_target_bed,
                                      min_read_len, common_seq1, common_seq2, num_threads, edit_dist, min_mapq,
                                      max_soft_clip, max_dist, memory_size, snp_filter, indel_filter, ref_ens,
                                      bwa_dir, samtools_dir, umitools_dir, gatk_dir, benchmark_dir, ref_index_name,
                                      ref_fa_file, total_ref_fa_file, total_ref_fa_dict, known_sites, erc,db_cosmic, 
                                      db_clinvar, db_g1000, truth_vcf, confident_region_bed, test_level])
    
    pool = multiprocessing.Pool(processes=threads)
    results = pool.map(main_run_germline_variant_calling, path_sampleID)
    print("--" * 20)
    pool.close()   # 关闭pool, 则不会有新的进程添加进去
    pool.join()    # 必须在join之前close, 然后join等待pool中所有的线程执行完毕
    print("All process done.")

    print("Return results: ")
    for i in results:
        print(i)   # 获得进程的执行结果