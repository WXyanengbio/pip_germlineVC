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
                              store_reformat_logs, store_germline_VC_logs, store_annotation_logs
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

def script_information():
    print ("\nApplication: pipelines of QIAseq Targeted DNA Panel\n")
    print ("=====================================================================")
    print ("Required environment: python \ bwa \ samtools \ GATK")

parser = argparse.ArgumentParser(usage = "\n\npython %(prog)s --source --sample_name --tailname --primers_file --exome_target_bed --output \
                                          --bwa_dir --samtools_dir --umitools_dir --gatk_dir --ref_index_name --ref_fa_file --total_ref_fa_file \
                                          --total_ref_fa_dict --known_sites --ERC --db_cosmic --db_clinvar --db_g1000 --anno_geneID")
parser.add_argument("--source", help = "Path to input reads in FASTA format", type = str)
parser.add_argument("--sample_name", help = "the sample name of raw reads", type = str)
parser.add_argument("--tailname", help = "the tailname of sample raw reads", type =str)
parser.add_argument("--common_seq1", help = "the common seq1 of QIAseq Targeted DNA Panel", type = str, default = 'CAAAACGCAATACTGTACATT')
parser.add_argument("--common_seq2", help = "the common seq2 of QIAseq Targeted DNA Panel", type = str, default = 'ATTGGAGTCCT')
parser.add_argument("--output", help = "Path of output file", type = str)
parser.add_argument("--min_read_len", help = "the cutoff of the min read length", type = int, default = 40)
parser.add_argument("--bwa_dir", help = "the install path of bwa", type = str)
parser.add_argument("--ref_index_name", help = "the path of ref index--if there isn't a ref index, it will make a index in the path of ref fasta by bwa", type = str)
parser.add_argument("--ref_fa_file", help = "the path of ref fasta", type = str)
parser.add_argument("--total_ref_fa_file", help = "the path of ref total ref fa", type = str)
parser.add_argument("--total_ref_fa_dict", help = "the path of ref total ref fa dict", type = str)
parser.add_argument("--num_threads", help = "the number of threads to align", type = int, default = 4)
parser.add_argument("--samtools_dir", help = "the install path of samtools", type = str)
parser.add_argument("--min_mapq", help = "the parameter of filter alignment_sam", type = int, default = 17)
parser.add_argument("--max_soft_clip", help = "the parameter of filter alignment_sam", type = int, default = 10)
parser.add_argument("--max_dist", help = "the parameter of filter alignment_sam", type = int, default = 2)
parser.add_argument("--primers_file", help = "Load all primer sequences in the panel", type = str)
parser.add_argument("--umitools_dir", help = "the install path of umitools", type = str)
parser.add_argument("--edit_dist", help = "the parameter of edit distance between barcodes", type = int, default = 2)
parser.add_argument("--memorySize", help = "the cutoff of Java memory", type = str, default = '4G')
parser.add_argument("--gatk_dir", help = "the install path of GATK4", type = str)
parser.add_argument("--known_sites", help = "the list of --known-sites , sep=',' ", type = str)
parser.add_argument("--exome_target_bed", help = "the bed file of exome intervals", type = str) 
parser.add_argument("--ERC", help = "switch to running HaplotypeCaller in GVCF mode", type = str, default = 'no')
parser.add_argument("--read_filter", help = "add a read filter that deals with some problems", type = str, default = 'no')
parser.add_argument("--snp_filter", help = "add parameters for filtering SNPs", type = str, default = 'DP < 100 || QD < 2.0 || FS > 20 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0')
parser.add_argument("--indel_filter", help = "add parameters for filtering Indels", type = str, default = 'DP < 100 || QD < 2.0 || FS > 20 || ReadPosRankSum < -20.0')
parser.add_argument("--db_cosmic", help = "add cosmic databases of variants", type = str)
parser.add_argument("--db_clinvar", help = "add clinvar databases of variants", type = str)
parser.add_argument("--db_g1000", help = "add g1000 databases of variants", type = str)
parser.add_argument("--anno_geneID", help = "add annotation gene ID of variants", type = str)
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
    #---
    test_level = args.test
    ##########################################################################################
    #---trim
    ##########################################################################################
    #time cost
    time_start = time.time()
    #undetermined_dir
    undetermined_dir = out_dir + '/'+ 'undetermined'
    if not os.path.exists(undetermined_dir):
        os.makedirs(undetermined_dir)
    
    sample = sample + '_' + tailname
    read1 = source + '/' + sample + '_R1_001.fastq.gz'
    read2 = source + '/' + sample  + '_R2_001.fastq.gz'
    trimmed1 = undetermined_dir + '/' + sample + '_R1_undetermined.fastq'
    trimmed2 = undetermined_dir + '/' + sample + '_R2_undetermined.fastq'
    stats_file = undetermined_dir + '/' + sample + '_basic_stats.txt'

    logger_trim_process, logger_trim_errors = store_trim_logs(log_dir)
    trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len,
                           common_seq1, common_seq2, stats_file, logger_trim_process,
                           logger_trim_errors)
    
    logger_trim_process.info("Trimming of reads is completed after %.2f min.", (time.time()-time_start)/60)
    logger_pipeline_process.info("Trimming of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level == 1:
        print("Test trim module!")
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
    logger_pipeline_process.info("Alignment of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level == 2:
        print("Test align module!")
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
    logger_pipeline_process.info("Post Alignment of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level == 3:
        print("Test psot align module!")
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
    logger_pipeline_process.info("UMIs tools clustering of reads is completed after %.2f min.", (time.time()-time_start)/60)
    
    if test_level == 4:
        print("Test barcode clustering module!")
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
    
    logger_reformat_process.info('Finish reformating alignment SAM file after {0} min.'.format((time.time() - time_start) / 60))
    logger_pipeline_process.info('Finish reformating alignment SAM file after {0} min.'.format((time.time() - time_start) / 60))
    
    if test_level == 5:
        print("Test reformat sam module!")
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
    logger_germline_VC_process.info('Finish germline variant calling after {0} min.'.format((time.time() - time_start) / 60))
    logger_pipeline_process.info('Finish germline variant callin after {0} min.'.format((time.time() - time_start) / 60))

    if test_level == 6:
        print("Test variant calling module!")
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
    logger_annotation_process.info('Finish annotation variant after {0} min.'.format((time.time() - time_start) / 60))
    logger_pipeline_process.info('Sample: {0}  has been processed after {1} min.'.format( sample , (time.time() - time_start1) / 60))

    if test_level == 7:
        print("Test annotation variant module!")
        exit()
if __name__ == '__main__':
    main()