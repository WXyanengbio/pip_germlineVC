__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
import re
import sys
import time
import gzip
import itertools
from difflib import SequenceMatcher
import argparse

#sys.path.append(os.path.split(os.path.realpath(__file__))[0])
#import the logging functions
from pipelines.log.log import store_pipeline_logs , store_trim_logs, store_filter_logs, store_align_logs , store_cluster_logs , \
                              store_reformat_logs, store_germline_VC_logs, store_annotation_logs, store_static_logs, store_benchmark_logs
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
from pipelines.variant_call.g_variantCall1 import sam_to_bem , germline_variant_calling
#import the annotation variant funcitons
from pipelines.variant_call.annotation_gatk_HC import annotationmain
#import the static functions
from pipelines.static.prestatic_module import qc_raw_reads, static_depth_coverage, static_sam_bam, static_time, merge_static_sam_bam, split_N_bases
#import the benchmaking funciton
from pipelines.benchmark.hap_benchmark import hap_py
def script_information():
    print ("\nApplication: pipelines of QIAseq Targeted DNA Panel\n")
    print ("=====================================================================")
    print ("Required environment: python \ bwa \ samtools \ GATK")

parser = argparse.ArgumentParser(usage = "\n\npython %(prog)s --source --sample_name --tailname --primers_file --exome_target_bed --output \
--fastQC_dir --bwa_dir --samtools_dir --umitools_dir --gatk_dir --benchmark_dir --ref_index_name --ref_fa_file --total_ref_fa_file \
--total_ref_fa_dict --known_sites --ERC --db_cosmic --db_clinvar --db_g1000 --anno_geneID --truth_vcf --confident_region_bed")

parser.add_argument("--source", 
                    help = "Path to input reads in FASTA format", 
                    type = str)
parser.add_argument("--sample_name", 
                    help = "the sample name of raw reads", 
                    type = str)
parser.add_argument("--tailname", 
                    help = "the tailname of sample raw reads", 
                    type =str)
parser.add_argument("--common_seq1",
                    type = str, 
                    default = 'CAAAACGCAATACTGTACATT',
                    help = "the common seq1 of QIAseq Targeted DNA Panel")
parser.add_argument("--common_seq2",
                    type = str, 
                    default = 'ATTGGAGTCCT',
                    help = "the common seq2 of QIAseq Targeted DNA Panel")
parser.add_argument("--output", 
                    help = "Path of output file", 
                    type = str)
parser.add_argument("--fastQC_dir", 
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
                    default= '/home/dell/Works/Projects/Datasets/genome/target_breast/target_breast',
                    help = "the path of ref index--if there isn't a ref index, it will make a index in the path of ref fasta by bwa"
                    )
parser.add_argument("--ref_fa_file",
                    type = str,
                    default= '/home/dell/Works/Projects/Datasets/genome/target_breast/target_breast.refSeq.fa',
                    help = "the path of ref fasta"
                    )
parser.add_argument("--total_ref_fa_file", 
                    type = str,
                    default= '/home/dell/Works/Projects/Datasets/genome/ucsc.hg19.fasta',
                    help = "the path of ref total ref fa"
                    )
parser.add_argument("--total_ref_fa_dict", 
                    type = str,
                    default= '/home/dell/Works/Projects/Datasets/genome/ucsc.hg19.dict',
                    help = "the path of ref total ref fa dict"
                    )
parser.add_argument("--truth_vcf", 
                    type = str,
                    default= '/home/dell/Works/Projects/Datasets/truthDB/samll_variant/NA12878.vcf.gz',
                    help = "the path of truth VCF of variant for benchmarking"
                    )
parser.add_argument("--confident_region_bed", 
                    type = str,
                    default= '/home/dell/Works/Projects/Datasets/truthDB/confidentregion/ConfidentRegions.bed',
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
                    help = "Load all primer sequences in the panel", 
                    type = str)
parser.add_argument("--edit_dist",
                    type = int, 
                    default = 2,
                    help = "the parameter of edit distance between barcodes"
                    )
parser.add_argument("--memorySize", 
                    type = str, 
                    default = '4G',
                    help = "the cutoff of Java memory"
                    )
parser.add_argument("--known_sites", 
                    type = str,
                    default ='/home/dell/Works/Projects/Datasets/db/known_sites/1000G_phase1.snps.high_confidence.hg19.sites.vcf,/home/dell/Works/Projects/Datasets/db/known_sites/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf,/home/dell/Works/Projects/Datasets/db/known_sites/dbsnp_138.hg19.vcf',
                    help = "the list of --known-sites , sep by: , "
                   )
                    
parser.add_argument("--exome_target_bed",
                    help = "the bed file of exome intervals", 
                    type = str) 
parser.add_argument("--ERC", 
                    type = str, 
                    default = 'no',
                    help = "switch to running HaplotypeCaller in GVCF mode"
                    )
parser.add_argument("--read_filter", 
                    type = str, 
                    default = 'no',
                    help = "add a read filter that deals with some problems"
                    )
parser.add_argument("--snp_filter", 
                    type = str, 
                    default = 'QD < 9.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0',
                    help = "add parameters for filtering SNPs"
                    )
parser.add_argument("--indel_filter",
                    type = str, 
                    default = 'QD < 9.0 || FS > 200 || ReadPosRankSum < -20.0 || SOR > 10.0',
                    help = "add parameters for filtering Indels"
                    )
parser.add_argument("--db_cosmic", 
                    type = str,
                    default ='/home/dell/Works/Projects/Datasets/db/annonDB/COSMIC_variant.csv',
                    help = "add cosmic databases of variants")
parser.add_argument("--db_clinvar", 
                    type = str,
                    default ='/home/dell/Works/Projects/Datasets/db/annonDB/breast_cancer_variant_clinvar.csv',
                    help = "add clinvar databases of variants")
parser.add_argument("--db_g1000", 
                    type = str,
                    default ='/home/dell/Works/Projects/Datasets/db/annonDB/breast_cancer_variant_1000genomes.csv',
                    help = "add g1000 databases of variants")
parser.add_argument("--anno_geneID", 
                    type = str,
                    default ='/home/dell/Works/Projects/Datasets/db/annonDB/geneid_cancer_qiagen.csv',
                    help = "add annotation gene ID of variants")
parser.add_argument("-v", 
                    '-version', 
                    action = 'version', 
                    version =' %(prog)s 1.0')
parser.add_argument("--test",
                    type = int,
                    default = 0,
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

def main():
    #time cost
    time_start1 = time.time()
    #---check the outputdir
    out_dir = args.output
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #----pipeline log file
    log_dir = out_dir + '/' + 'log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logger_pipeline_process, logger_pipeline_errors = store_pipeline_logs(log_dir)

    #---input
    source = args.source
    sample = args.sample_name
    tailname = args.tailname
    #QC
    fastQC_dir = args.fastQC_dir
    #trim
    min_read_len = args.min_read_len
    common_seq1 = args.common_seq1
    common_seq2 = args.common_seq2
    #bwa---align
    num_threads = args.num_threads
    bwa_dir = args.bwa_dir
    ref_fa_file = args.ref_fa_file
    ref_index_name = args.ref_index_name
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
    memorySize = args.memorySize
    memorySize = '-Xmx' + memorySize + ' ' + '-Djava.io.tmpdir=./'
    gatk_dir = args.gatk_dir
    samtools_dir = args.samtools_dir
    total_ref_fa_dict = args.total_ref_fa_dict
    total_ref_fa_file = args.total_ref_fa_file
    known_sites = args.known_sites
    #read_length = args.read_length
    exome_target_bed = args.exome_target_bed
    ERC = args.ERC
    read_filter = args.read_filter
    snp_filter = args.snp_filter
    indel_filter = args.indel_filter
    #--annotation
    ref_ens = args.anno_geneID
    db_cosmic = args.db_cosmic
    db_clinvar = args.db_clinvar
    db_g1000 = args.db_g1000
    #--benchmark
    benchmark_dir = args.benchmark_dir
    confident_region_bed = args.confident_region_bed
    truth_vcf = args.truth_vcf
    #---
    test_level = args.test
    ##########################################################################################
    #---QC
    ##########################################################################################
    #time cost
    time_start = time.time()
    #qc_dir
    module = "QC"
    qc_dir = out_dir + '/'+ 'QC'
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    sample = sample + '_' + tailname
    read1 = source + '/' + sample + '_R1_001.fastq.gz'
    read2 = source + '/' + sample  + '_R2_001.fastq.gz'
    #---
    logger_static_process, logger_static_errors = store_static_logs(log_dir)
    qc_result1, qc_result2 = qc_raw_reads(fastQC_dir, qc_dir, 
                 sample, module, 
                 read1, read2,
                 logger_static_process, logger_static_errors)
    #--check the quality of the raw reads
    if float(qc_result1[6].strip('%')) > 80 and  float(qc_result2[6].strip('%')) > 80:
        print("The ratio of read1 and read2 with Q30 quality are both higher than 80%.")
    else:
        exit("The ratio of read1 and read2 with Q30 quality are both lower than 80%!!!!!!!")
    #--static the N base in raw reads and set the cutoff of the min read length
    if max(int(split_N_bases(qc_result1[7])), int(split_N_bases(qc_result2[7]))) < min_read_len:
        print("The cutoff of the min read length is the default: {0}".format(min_read_len))
    else:
        min_read_len = max(int(split_N_bases(qc_result1[7])), int(split_N_bases(qc_result2[7])))
        print("The cutoff of the min read length is based on the N base in the reads: {0}".format(min_read_len))

    logger_static_process.info("QC of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("QC of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level >= 0:
        print("Test QC module!")
    else:
        exit()
    ##########################################################################################
    #---trim
    ##########################################################################################
    #time cost
    time_start = time.time()
    #undetermined_dir
    undetermined_dir = out_dir + '/'+ 'undetermined'
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
        print("Test trim module!")
    else:
        exit()
    ##########################################################################################
    #---align
    ##########################################################################################
    #time cost
    time_start = time.time()
    #aligned_dir
    aligned_dir = out_dir + '/'+ 'aligned'
    if not os.path.exists(aligned_dir):
        os.makedirs(aligned_dir)

    trim_read1 = trimmed1
    trim_read2 = trimmed2
    out_file = aligned_dir + '/' + sample + '_aligned.sam'

    logger_bwa_process, logger_bwa_errors = store_align_logs(log_dir)
    
    returncode = align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, trim_read1, trim_read2, 
                                                out_file, num_threads, logger_bwa_process, logger_bwa_errors)
    
    logger_bwa_process.info("Alignment of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("Align of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level >= 2:
        print("Test align module!")
    else:
        exit()

    ######################################annotationmain####################################################
    #---post_align
    ##########################################################################################
    #time cost
    time_start = time.time()
    #post_aligned_dir
    filtered_dir = out_dir + '/'+ 'filtered'
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
        print("Test psot align module!")
    else:
        exit()
    ##########################################################################################
    #---barcode clustering
    ##########################################################################################
    #time cost
    time_start = time.time()
    #clustering_dir
    clustered_dir = out_dir + '/'+ 'clustered'
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
        print("Test barcode clustering module!")
    else:
        exit()
    ##########################################################################################
    #---reformat 
    ##########################################################################################
    #time cost
    time_start = time.time()
    #reformated_dir
    reformated_dir = out_dir + '/'+ 'reformated'
    if not os.path.exists(reformated_dir):
        os.makedirs(reformated_dir)

    logger_reformat_process, logger_reformat_errors = store_cluster_logs(log_dir)
    
    alignment_sam = umis_sam
    output_sam = reformated_dir + '/' + sample + '_vcready.sam'
    reformat_sam(alignment_sam, output_sam, logger_reformat_process, logger_reformat_errors)
    
    logger_reformat_process.info('Finish reformating alignment SAM file is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Reformat alignment SAM file is completed after %.2f min.',(time.time() - time_start)/60)
    
    if test_level >= 5:
        print("Test reformat sam module!")
    else:
        exit()
    ##########################################################################################
    #---Germline variant calling
    ##########################################################################################
    #time cost
    time_start = time.time()
    #germline_VC_dir
    germline_VC_dir = out_dir + '/'+ 'germline_VC'
    if not os.path.exists(germline_VC_dir):
        os.makedirs(germline_VC_dir)
    
    logger_germline_VC_process, logger_germline_VC_errors = store_germline_VC_logs(log_dir)
    
    #---modify the known-sites
    known_sites = known_sites.replace(',' , ' --known-sites ')
    
    vready_sam = output_sam
    sam_to_bem(gatk_dir, samtools_dir,
               vready_sam, sample,
               germline_VC_dir, memorySize,
               exome_target_bed, 
               total_ref_fa_file, total_ref_fa_dict,
               known_sites,
               logger_germline_VC_process, logger_germline_VC_errors)

    marked_BQSR_bam = germline_VC_dir + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
    Exon_Interval = germline_VC_dir + '/' + 'target_interval.list'
    germline_variant_calling(gatk_dir, marked_BQSR_bam,
                             sample, germline_VC_dir, 
                             memorySize, total_ref_fa_file, 
                             Exon_Interval, ERC,
                             read_filter,
                             snp_filter,indel_filter,
                             logger_germline_VC_process, logger_germline_VC_errors)
    logger_germline_VC_process.info('Germline variant calling is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('variant_calling is completed after %.2f min.',(time.time() - time_start)/60)

    if test_level >= 6:
        print("Test variant calling module!")
    else:
        exit()
    
    ##########################################################################################
    #---Annotation variant calling
    ##########################################################################################
    #time cost
    time_start = time.time()
    #Annotation dir
    annotation_dir = out_dir + '/'+ 'annotation'
    if not os.path.exists(annotation_dir):
        os.makedirs(annotation_dir)
    
    logger_annotation_process, logger_annotation_errors = store_annotation_logs(log_dir)
    
    snp_vcf = germline_VC_dir + '/'  + sample + '.raw_variants_SNP.vcf'
    filter_snp = germline_VC_dir + '/'  + sample + '.filter_SNP.vcf'
    indel_vcf = germline_VC_dir + '/'  + sample + '.raw_variants_indel.vcf'
    filter_indel = germline_VC_dir + '/'  + sample + '.filter_indel.vcf'
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
        print("Test annotation variant module!")
    else:
        exit()
    ##########################################################################################
    #---benchmarking 
    ##########################################################################################
    #time cost
    time_start = time.time()
    #benchmarking_dir
    benchmarking_dir = out_dir + '/'+ 'benchmarking'
    if not os.path.exists(benchmarking_dir):
        os.makedirs(benchmarking_dir)

    logger_benchmark_process, logger_benchmark_errors = store_benchmark_logs(log_dir)
    
    hap_py(benchmark_dir,truth_vcf, filter_snp, confident_region_bed, benchmarking_dir, total_ref_fa_file, Exon_Interval, logger_benchmark_process, logger_benchmark_errors)
    hap_py(benchmark_dir,truth_vcf, filter_indel, confident_region_bed, benchmarking_dir, total_ref_fa_file, Exon_Interval, logger_benchmark_process, logger_benchmark_errors)
    
    logger_benchmark_process.info('benchmarking is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Benchmark is completed after %.2f min.',(time.time() - time_start)/60)
    
    if test_level >= 8:
        print("Test benchmark module!")
    else:
        exit()
    ##########################################################################################
    #---Static of the variant calling pipeline
    ##########################################################################################
    #time cost
    time_start = time.time()
    #static dir
    static_dir = out_dir + '/'+ 'static'
    if not os.path.exists(static_dir):
        os.makedirs(static_dir)
    #-static the clean reads
    #static trim dir
    static_trim_dir = static_dir + '/'+ 'trim_QC'
    if not os.path.exists(static_trim_dir):
        os.makedirs(static_trim_dir)
    module = "Trim"
    trim_result1, trim_result2 = qc_raw_reads(fastQC_dir, static_trim_dir, 
                 sample, module, 
                 trimmed1, trimmed2,
                 logger_static_process, logger_static_errors)
    #--static the align
    module1 = "Align"
    align_sorted_bam = static_depth_coverage(samtools_dir, out_file, static_dir, sample, module1, logger_static_process, logger_static_errors)
    align_static = static_sam_bam(samtools_dir, sorted_bam, static_dir,sample, module1, logger_static_process, logger_static_errors)
    #--static the filter
    #----cluster module would build the filter sorted bam, but it has been changed UMIs-tools
    module2 = "Fliter"
    filtered_sorted_bam = static_depth_coverage(samtools_dir, filtered_sam, static_dir, sample, module2, logger_static_process, logger_static_errors)
    fliter_static = static_sam_bam(samtools_dir, filtered_sorted_bam, static_dir, sample, module2, logger_static_process, logger_static_errors)
    # static the umi-tools
    module3 = "Cluster_reformat"
    cr_sorted_bam = static_depth_coverage(samtools_dir, vready_sam, static_dir, sample, module3, logger_static_process, logger_static_errors)
    cr_static = static_sam_bam(samtools_dir, cr_sorted_bam, static_dir, sample, module3, logger_static_process, logger_static_errors)
    #-merge the sorted bam
    merge_static_sam_bam(logger_static_process, logger_static_errors, static_dir, sample, ','.join([module1,module2,module3]),align_static,fliter_static,cr_static)
    
    logger_static_process.info('static is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Static is completed after %.2f min.',(time.time() - time_start)/60)
    logger_pipeline_process.info('Pipeline : {0} is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start1)/60))))
    #--static the time cost
    process_log = out_dir + '/'+ 'log' + '/' + 'process.log'
    static_time(static_dir, sample, process_log, logger_static_process, logger_static_errors)
    #---
    if test_level == 9:
        print("Test Static module!")
if __name__ == '__main__':
    main()