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
import yaml

#sys.path.append(os.path.split(os.path.realpath(__file__))[0])
#import the logging functions
from pipelines.log.log_v1 import store_pipeline_logs,store_trim_logs, store_filter_logs, store_align_logs , store_cluster_logs , \
                              store_reformat_logs, store_germline_vc_logs, store_annotation_logs, store_statistics_logs, store_benchmark_logs
#import the trim function
from pipelines.trim.trim_reads_v1 import trim_read_pairs
#import the align function
from pipelines.align.align_reads_v1 import align_reads_bwa
#import the post-align functions
from pipelines.post_align.post_alignment_v1 import filter_alignment_samtools, identify_gs_primers
#import the barcode clusering function
from pipelines.cluster_barcode.umitools_v1 import umitool
#import the reformat sam function
from pipelines.reformat.reformat_sam import reformat_sam
#import the variant calling funcitons
from pipelines.variant_call.g_variantcall1_v1 import sam_to_bam, germline_variant_calling, strelka2_call, samtools_call, varsan2_call
#import the annotation variant funcitons
from pipelines.variant_call.annotation_gatk_hc_v1 import annotationmain,read_vcf_filter
#import the statistics functions
from pipelines.statistics.prestatistics_module_v1 import qc_raw_reads, statistics_depth_coverage, statistics_sam_bam, statistics_time, merge_statistics_sam_bam, merge_statistics
#import the benchmaking funciton
#from pipelines.benchmark.hap_benchmark import hap_py

def script_information():
    print ("\nApplication: multiprocessing of pipelines of QIAseq Targeted DNA Panel\n")
    print ("=====================================================================")
    print ("Required environment: python \ JAVA\ R \ bwa \ samtools \ GATK \ UMI-tools \ hay.py ")

parser = argparse.ArgumentParser(usage = "\n\npython3.6 %(prog)s --yaml --threads")

parser.add_argument("--yaml", 
                    type =str,
                    default = 'null',
                    help = "The yaml file of the parameters")
parser.add_argument("--file_list", 
                    help = "Path to file contian tow colown: one is the path of raw reads, two is the sampleID of raw reads", 
                    type = str)
parser.add_argument("--threads", 
                    help = "the threads of multiprocessing ", 
                    type = int)
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
parser.add_argument("--exome_target",
                    type = str,
                    default = 'genome/target_breast/BRCA.csv',
                    help = "the posi of exon to statistics the depth of exon bases") 
parser.add_argument("--exome_target_bed",
                    type = str,
                    default = 'genome/target_breast/target_breast_BRCA.bed',
                    help = "the bed file of exome intervals") 
parser.add_argument("--erc", 
                    type = str, 
                    default = 'GVCF',
                    help = "switch to running HaplotypeCaller in GVCF mode"
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
parser.add_argument("--tools",
                    type = str,
                    default = 'all',
                    help = "the subprocess of the script, such as qc, trim, align, post_align, cluster,reformat, variant_call, annotation, statis")

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
    primers_file = path_sampleID_sub[3]
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
    ref_index_name = path_sampleID_sub[21]
    ref_fa_file = path_sampleID_sub[22]
    total_ref_fa_file = path_sampleID_sub[23]
    total_ref_fa_dict = path_sampleID_sub[24]
    known_sites = path_sampleID_sub[25]
    erc = path_sampleID_sub[26]
    db_cosmic = path_sampleID_sub[27]
    db_clinvar = path_sampleID_sub[28]
    db_g1000 = path_sampleID_sub[29]
    test_level = path_sampleID_sub[30]
    exome_target = path_sampleID_sub[31]
    calling = path_sampleID_sub[32]
    tabix = path_sampleID_sub[33]
    bgzip = path_sampleID_sub[34]
    bcftools_dir = path_sampleID_sub[35]
    varsan2_dir = path_sampleID_sub[36]
    strelka2_dir = path_sampleID_sub[37]
    total_ref_chrom_fa_file = path_sampleID_sub[38]
    #---check the output
    out_dir = output  + '/' + sample
    #out_dir1 = output
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #----pipeline log file
    log_dir = out_dir + '/' + 'log'
    if not os.path.isdir(log_dir):
        try:
            os.makedirs(log_dir)
        except OSError as e:
            if e.errno == 17:  # errno.EEXIST
                logger_pipeline_process = log_dir
                logger_pipeline_errors = log_dir
    logger_pipeline_process = log_dir
    logger_pipeline_errors = log_dir
     #time cost
    time_start = time.time()
    module = "QC"
    #---
    read1 = source + '/' + sample + '_R1_001.fastq.gz'
    read2 = source + '/' + sample  + '_R2_001.fastq.gz'
    #---
    if tools in ['all','qc']:
        print("Test QC module!\n")
    #qc_dir
        qc_dir = out_dir + '/'+ 'QC'
        if not os.path.exists(qc_dir):
            os.makedirs(qc_dir)
    #---
        logger_statistics_process =log_dir
        logger_statistics_errors =log_dir
        qc_result1, qc_result2 = qc_raw_reads(fastqc_dir, qc_dir, 
                 sample, module, 
                 read1, read2,
                 logger_statistics_process, logger_statistics_errors)
    #--check the quality of the raw reads
        if float(qc_result1[7].strip('%')) > 70 and  float(qc_result2[7].strip('%')) > 70:
            print("The ratio of read1 and read2 with Q30 quality are both higher than 70%.")
        else:
            exit("The ratio of read1 and read2 with Q30 quality are both lower than 80%!!!!!!!")
    #--statistics the N base in raw reads and set the cutoff of the min read length
        if max(int(qc_result1[9]), int(qc_result2[9])) < min_read_len:
            print("The cutoff of the min read length is the default: {0}".format(min_read_len))
        elif max(int(qc_result1[9]), int(qc_result2[9])) > 1.25*min_read_len:
            #min_read_len = max(int(qc_result1[9]), int(qc_result2[9]))
            print("The sites in raw reads have N base: 1-{0}".format(min_read_len))
            print("The cutoff of the min read length is based on the N base in the reads: {0}".format(min_read_len))
        store_pipeline_logs(logger_pipeline_process,'null',"--{0}--QC of reads is completed after {1} min.\n".format(sample,str('%.3f' % ((time.time()-time_start)/60))))
        store_statistics_logs(logger_statistics_process,'null',"--{0}--QC of reads is completed after {1} min.\n".format(sample,str('%.3f' % ((time.time()-time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    #---trim
    ##########################################################################################
    #time cost
    time_start = time.time()
    #undetermined_dir
    undetermined_dir = out_dir + '/'+ 'undetermined'
    trimmed1 = undetermined_dir + '/' + sample + '_R1_undetermined.fastq'
    trimmed2 = undetermined_dir + '/' + sample + '_R2_undetermined.fastq'
    stats_file = undetermined_dir + '/' + sample + '_basic_stats.txt'
    if tools in ['all','trim']:
        print("please check the QC subprocess result--the min read length!")
        print("The cutoff of the min read length is the default: {0}".format(min_read_len))
        print("Test trim module!\n\n\n")
        # ,kdir undetermined_dir
        if not os.path.exists(undetermined_dir):
            os.makedirs(undetermined_dir)
        
        logger_trim_process = log_dir
        logger_trim_errors = log_dir
        
        trim_read_pairs(read1, read2, trimmed1, trimmed2, min_read_len,
                           common_seq1, common_seq2, stats_file, logger_trim_process,
                           logger_trim_errors)
        store_pipeline_logs(logger_pipeline_process,'null',"--{0}--Trim of reads is completed after {1} min.\n".format(sample,str('%.3f' % ((time.time()-time_start)/60))))
        store_trim_logs(logger_trim_process,'null',"--{0}--Trimming of reads is completed after {1} min.\n".format(sample,str('%.3f' % ((time.time()-time_start)/60))))
        #logger_pipeline_process.info("--{0}--Trim of reads is completed after {1} min.".format(sample,str('%.3f' % ((time.time()-time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    #---align
    ##########################################################################################
    #time cost
    time_start = time.time()
    #aligned_dir
    aligned_dir = out_dir + '/'+ 'aligned'
    trim_read1 = out_dir + '/'+ 'undetermined'+ '/' + sample + '_R1_undetermined.fastq'
    trim_read2 = out_dir + '/'+ 'undetermined'+ '/' + sample + '_R2_undetermined.fastq'

    out_file = aligned_dir + '/' + sample + '_aligned.sam'
    if tools in ['all','align']:
        print("please check the Trim subprocess result--undetermined.fastq!")
        print("Test align module!\n")
        if not os.path.exists(aligned_dir):
            os.makedirs(aligned_dir)
        logger_bwa_process = log_dir
        logger_bwa_errors = log_dir
        returncode = align_reads_bwa(bwa_dir, samtools_dir,ref_fa_file, ref_index_name, exome_target_bed, total_ref_fa_file, trim_read1, trim_read2, 
                                                out_file, num_threads, logger_bwa_process, logger_bwa_errors)
        store_align_logs(logger_bwa_process,'null',"--{0}--Alignment of reads is completed after {1} min.".format(sample,str('%.3f' % ((time.time()-time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',"--{0}--Align of reads is completed after {1} min.".format(sample,str('%.3f' % ((time.time()-time_start)/60))))
        print("--" * 20 + '\n\n')

    ######################################annotationmain####################################################
    #---post_align
    ##########################################################################################
    #time cost
    time_start = time.time()
    #post_aligned_dir
    filtered_dir = out_dir + '/'+ 'filtered'
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

    if tools in ['all','post_align']:
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
        store_filter_logs(logger_filter_process,'null',"--{0}--Post Alignment of reads is completed after {1} min.".format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',"--{0}--Post_align of reads is completed after {1} min.".format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    #---barcode clustering
    ##########################################################################################
    #time cost
    time_start = time.time()
    #clustering_dir
    clustered_dir = out_dir + '/'+ 'clustered'
    filtered_sam = out_dir + '/'+ 'filtered'+ '/' + sample + '_filtered.sam'
    filtered_bam = clustered_dir + '/' + sample + '_filtered.bam'
    sorted_bam = clustered_dir + '/' + sample + '_filtered_sorted.bam'
    umitool_stats = clustered_dir + '/' + sample + '_deduplicated'
    umis_sam = clustered_dir + '/' + sample + '_umis.sam'
    if tools in ['all','cluster']:
        print("please check the post algin subprocess result--filtered.sam!")
        print("Test cluster module!\n")
        if not os.path.exists(clustered_dir):
            os.makedirs(clustered_dir)
        logger_umi_process = log_dir
        logger_umi_errors = log_dir
        umitool(samtools_dir, umitools_dir, filtered_sam ,filtered_bam , sorted_bam, umitool_stats , umis_sam, edit_dist, logger_umi_process, logger_umi_errors)
        store_cluster_logs(logger_umi_process,'null',"--{0}--UMIs tools clustering of reads is completed after {1} min.".format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null',"--{0}--Cluster of reads is completed after {1} min.".format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    #---reformat 
    ##########################################################################################
    #time cost
    time_start = time.time()
    #reformated_dir
    reformated_dir = out_dir + '/'+ 'reformated'
    alignment_sam = out_dir + '/'+ 'clustered'+ '/' + sample + '_umis.sam'
    output_sam = reformated_dir + '/' + sample + '_vcready.sam'
    if tools in ['all','reformat']:
        print("please check the cluster subprocess result--umis.sam!")
        print("Test reformat module!\n")
        if not os.path.exists(reformated_dir):
            os.makedirs(reformated_dir)
        logger_reformat_process = log_dir
        logger_reformat_errors = log_dir
        reformat_sam(alignment_sam, output_sam, logger_reformat_process, logger_reformat_errors)
        store_reformat_logs(logger_reformat_process,'null','--{0}--Finish reformating alignment SAM file is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null','--{0}--Reformat alignment SAM file is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    #---Germline variant calling
    ##########################################################################################
    #time cost
    time_start = time.time()
    #germline_vc_dir
    germline_vc_dir = out_dir + '/'+ 'germline_vc'
    #---modify the known-sites
    known_sites = known_sites.replace(',' , ' --known-sites ' + args.datasets_dir + '/'  )
    known_sites = args.datasets_dir + '/' + known_sites
    vready_sam = out_dir + '/'+ 'reformated'+ '/' + sample + '_vcready.sam'
    marked_bqsr_bam = germline_vc_dir + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
    if os.path.basename(exome_target_bed) != 'all':
        exon_interval = germline_vc_dir + '/' + 'target_interval.list'
    else:
        exon_interval = 'all'
    if tools in ['all','variant_call']:
        print("please check the reformat subprocess result--vcready.sam!")
        print("Test variant_call module!\n")
        if not os.path.exists(germline_vc_dir):
            os.makedirs(germline_vc_dir)
        logger_germline_vc_process = log_dir
        logger_germline_vc_errors = log_dir

        bqsr = 'n'
        bam_to_variant, bqsr_bam_to_variant =sam_to_bam(gatk_dir, samtools_dir,
               vready_sam, sample,
               germline_vc_dir, memory_size,
               exome_target_bed, 
               total_ref_fa_file, total_ref_fa_dict,
               known_sites,
               logger_germline_vc_process, logger_germline_vc_errors,bqsr)
        if calling =='GATK':
            germline_variant_calling(gatk_dir, bam_to_variant,
                             sample, germline_vc_dir, 
                             memory_size, total_ref_fa_file, 
                             exon_interval, erc,
                             snp_filter,indel_filter,
                             logger_germline_vc_process, logger_germline_vc_errors)
        elif calling =='strelka2':
            strelka2_call(strelka2_dir,
                          bgzip,
                          tabix,
                          total_ref_fa_file,
                          germline_vc_dir,
                          sample, bam_to_variant,
                          exome_target_bed,logger_germline_vc_process,logger_germline_vc_errors)
        elif calling =='samtools':
            samtools_call(samtools_dir, 
                          bcftools_dir, 
                          bam_to_variant, 
                          sample, 
                          germline_vc_dir, 
                          total_ref_fa_file,logger_germline_vc_process,logger_germline_vc_errors)
        elif calling =='varscan2':
            varsan2_call(samtools_dir, 
                         varsan2_dir,
                         total_ref_fa_file,
                         germline_vc_dir,
                         sample, 
                         bam_to_variant,logger_germline_vc_process,logger_germline_vc_errors)
        store_germline_vc_logs(logger_germline_vc_process,'null','--{0}--Germline variant calling is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null', '--{0}--variant_calling is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    #---Annotation variant calling
    ##########################################################################################
    #time cost
    time_start = time.time()
    #Annotation dir
    annotation_dir = out_dir + '/'+ 'annotation'
    raw_vcf = germline_vc_dir + '/'  + sample + '.raw_variants.vcf'
    snp_vcf = germline_vc_dir + '/'  + sample + '.raw_variants_SNP.vcf'
    indel_vcf = germline_vc_dir + '/'  + sample + '.raw_variants_indel.vcf'
    #annotation
    if tools in ['all','annotation']:
        print("please check the variant_call subprocess result--VCF!")
        print("Test annotation module!\n")
    #Annotation dir
        if not os.path.exists(annotation_dir):
            os.makedirs(annotation_dir)
        logger_annotation_process = log_dir
        logger_annotation_errors = log_dir
        snp_limit, indel_limit= read_vcf_filter(snp_filter, indel_filter)
        annotationmain(db_cosmic, db_clinvar, db_g1000, 
                   ref_ens,
                   raw_vcf, sample,snp_limit, indel_limit,
                   annotation_dir, logger_annotation_process, logger_annotation_errors,calling)
        if 'GATK' in calling :
            calling1='GATK'
            annotationmain(db_cosmic, db_clinvar, db_g1000,
                   ref_ens,
                   snp_vcf, sample,snp_limit, indel_limit,
                   annotation_dir, logger_annotation_process, logger_annotation_errors,calling1)

            annotationmain(db_cosmic, db_clinvar, db_g1000,
                   ref_ens,
                   indel_vcf, sample,snp_limit, indel_limit,
                   annotation_dir, logger_annotation_process, logger_annotation_errors,calling1)

        store_annotation_logs(logger_annotation_process,'null', '--{0}--Finish annotation variant  is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null','--{0}--Annotation variant is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        print("--" * 20 + '\n\n')

    ##########################################################################################
    #---statistics of the variant calling pipeline
    ##########################################################################################
    #time cost
    time_start = time.time()
    #statistics dir
    statistics_dir = out_dir + '/'+ 'statistics'
    if not os.path.exists(statistics_dir):
        os.makedirs(statistics_dir)
    #-statistics the clean reads
    #statistics trim dir
    statistics_trim_dir = statistics_dir + '/'+ 'trim_QC'
    if not os.path.exists(statistics_trim_dir):
        os.makedirs(statistics_trim_dir)
    if tools in ['all','statis']:
        print("please check the others subprocess results!")
        print("Test statistics module!\n")
        if tools == 'statis':
            logger_statistics_process = log_dir
            logger_statistics_errors = log_dir
    #statistics dir
        if not os.path.exists(statistics_dir):
            os.makedirs(statistics_dir)
        module = "Trim"
        trim_result1, trim_result2 = qc_raw_reads(fastqc_dir, statistics_trim_dir, 
                 sample, module, 
                 trimmed1, trimmed2,
                 logger_statistics_process, logger_statistics_errors)
    #--statistics the align
        module1 = "Align"
        align_sorted_bam = statistics_depth_coverage(samtools_dir, out_file, statistics_dir, sample, module1, exome_target_bed,logger_statistics_process, logger_statistics_errors)
        align_statistics = statistics_sam_bam(samtools_dir, align_sorted_bam, statistics_dir,sample, module1, logger_statistics_process, logger_statistics_errors)
    #--statistics the filter
    #----cluster module would build the filter sorted bam, but it has been changed UMIs-tools
        module2 = "Fliter"
        filtered_sorted_bam = statistics_depth_coverage(samtools_dir, filtered_sam, statistics_dir, sample, module2,exome_target_bed, logger_statistics_process, logger_statistics_errors)
        fliter_statistics = statistics_sam_bam(samtools_dir, filtered_sorted_bam, statistics_dir, sample, module2, logger_statistics_process, logger_statistics_errors)
    # statistics the umi-tools
        module3 = "Cluster_reformat"
        cr_sorted_bam = statistics_depth_coverage(samtools_dir, vready_sam, statistics_dir, sample, module3,exome_target_bed, logger_statistics_process, logger_statistics_errors)
        cr_statistics = statistics_sam_bam(samtools_dir, cr_sorted_bam, statistics_dir, sample, module3, logger_statistics_process, logger_statistics_errors)
    #-merge the sorted bam
        merge_statistics_sam_bam(logger_statistics_process, logger_statistics_errors, statistics_dir, 
                                 sample, ','.join([module1,module2,module3]),align_statistics,fliter_statistics,cr_statistics)
        store_statistics_logs(logger_statistics_process, 'null','--{0}--Statistics is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null','--{0}--Statistics is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start)/60))))
        store_pipeline_logs(logger_pipeline_process, 'null','--{0}--Pipeline : {0} is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start1)/60))))
    #--statistics the time cost
        process_log = out_dir + '/'+ 'log' + '/' + 'process.log'
        statistics_time(statistics_dir, sample, process_log, logger_statistics_process, logger_statistics_errors)
        print("--" * 20 + '\n\n')
    return 'Pipeline : {0} is completed after {1} min.'.format(sample, ('%.2f' % ((time.time() - time_start1)/60)))

#--
cores = multiprocessing.cpu_count()
data = psutil.virtual_memory()
total = data.total #总内存,单位为byte
free = data.available #可以内存

if __name__ == '__main__':
    #---input
    yaml_file = args.yaml
    if yaml_file != 'null':
        f = open(yaml_file)
        file_yaml = yaml.load(f)
        #print(file_yaml)
    if yaml_file != 'null' and 'file_list' in file_yaml.keys():
        file_list = file_yaml['file_list']
    else:
        file_list = args.file_list
    #---check the outputdir
    if yaml_file != 'null' and 'output' in file_yaml.keys():
        out_dir = file_yaml['output']
    else:
        exit('No output path in the yaml')
    if yaml_file != 'null' and 'threads' in file_yaml.keys():
        threads = file_yaml['threads']
    else:
        threads = args.threads
    if yaml_file != 'null' and 'fastqc_dir' in file_yaml.keys():
        fastqc_dir = file_yaml['fastqc_dir']
    else:
        fastqc_dir = args.fastqc_dir
    #trim
    if yaml_file != 'null' and 'min_read_len' in file_yaml.keys():
        min_read_len = file_yaml['min_read_len']
    else:
        min_read_len = args.min_read_len
    if yaml_file != 'null' and 'common_seq1' in file_yaml.keys():
        common_seq1 = file_yaml['common_seq1']
    else:
        common_seq1 = args.common_seq1
    if yaml_file != 'null' and 'common_seq2' in file_yaml.keys():
        common_seq2 = file_yaml['common_seq2']
    else:
        common_seq2 = args.common_seq2
    #bwa---align
    if yaml_file != 'null' and 'num_threads' in file_yaml.keys():
        num_threads = file_yaml['num_threads']
    else:
        num_threads = args.num_threads
    if yaml_file != 'null' and 'bwa_dir' in file_yaml.keys():
        bwa_dir = file_yaml['bwa_dir']
    else:
        bwa_dir = args.bwa_dir
    if yaml_file != 'null' and 'ref_fa_file' in file_yaml.keys():
        ref_fa_file = args.datasets_dir + '/' + file_yaml['ref_fa_file']
    else:
        ref_fa_file = args.datasets_dir + '/' + args.ref_fa_file
    if yaml_file != 'null' and 'ref_index_name' in file_yaml.keys():
        ref_index_name = args.datasets_dir + '/' + file_yaml['ref_index_name']
    else:
        ref_index_name = args.datasets_dir + '/' +args.ref_index_name
    #post-align
    if yaml_file != 'null' and 'samtools_dir' in file_yaml.keys():
        samtools_dir = file_yaml['samtools_dir']
    else:
        samtools_dir = args.samtools_dir
    if yaml_file != 'null' and 'min_mapq' in file_yaml.keys():
        min_mapq = file_yaml['min_mapq']
    else:
        min_mapq = args.min_mapq
    if yaml_file != 'null' and 'max_soft_clip' in file_yaml.keys():
        max_soft_clip = file_yaml['max_soft_clip']
    else:
        max_soft_clip = args.max_soft_clip
    if yaml_file != 'null' and 'max_dist' in file_yaml.keys():
        max_dist = file_yaml['max_dist']
    else:
        max_dist = args.max_dist
    if yaml_file != 'null' and 'primers_file' in file_yaml.keys():
        primers_file = file_yaml['primers_file']
    else:
        primers_file = args.primers_file
    #--clustering
    if yaml_file != 'null' and 'umitools_dir' in file_yaml.keys():
        umitools_dir = file_yaml['umitools_dir']
    else:
        umitools_dir = args.umitools_dir
    if yaml_file != 'null' and 'edit_dist' in file_yaml.keys():
        edit_dist = file_yaml['edit_dist']
    else:
        edit_dist = args.edit_dist
    #--variant calling
    if yaml_file != 'null' and 'memory_size' in file_yaml.keys():
        memory_size = file_yaml['memory_size']
    else:
        memory_size = args.memory_size
    memory_size = '-Xmx' + str(memory_size) + 'G ' + '-Djava.io.tmpdir=./'
    if yaml_file != 'null' and 'gatk_dir' in file_yaml.keys():
        gatk_dir = file_yaml['gatk_dir']
    else:
        gatk_dir = args.gatk_dir
    #--strelka2
    if yaml_file != 'null' and 'strelka2_dir' in file_yaml.keys():
        strelka2_dir = file_yaml['strelka2_dir']
    if yaml_file != 'null' and 'bgzip' in file_yaml.keys():
        bgzip = file_yaml['bgzip']
    if yaml_file != 'null' and 'tabix' in file_yaml.keys():
        tabix = file_yaml['tabix']
    if yaml_file != 'null' and 'total_ref_chrom_fa_file' in file_yaml.keys():
        total_ref_chrom_fa_file = args.datasets_dir + '/' + file_yaml['total_ref_chrom_fa_file']
    #--samtools
    if yaml_file != 'null' and 'bcftools_dir' in file_yaml.keys():
        bcftools_dir = file_yaml['bcftools_dir']
    #varscan2
    if yaml_file != 'null' and 'varsan2_dir' in file_yaml.keys():
        varsan2_dir = file_yaml['varsan2_dir']
    #--ref database for variant calling
    if yaml_file != 'null' and 'total_ref_fa_dict' in file_yaml.keys():
        total_ref_fa_dict = args.datasets_dir + '/' + file_yaml['total_ref_fa_dict']
    else:
        total_ref_fa_dict = args.datasets_dir + '/' + args.total_ref_fa_dict
    if yaml_file != 'null' and 'total_ref_fa_file' in file_yaml.keys():
        total_ref_fa_file = args.datasets_dir + '/' + file_yaml['total_ref_fa_file']
    else:
        total_ref_fa_file = args.datasets_dir + '/' + args.total_ref_fa_file
    if yaml_file != 'null' and 'known_sites' in file_yaml.keys():
        known_sites = file_yaml['known_sites']
    else:
        known_sites = args.known_sites
    #read_length = args.read_length
    if yaml_file != 'null' and 'exome_target' in file_yaml.keys():
        exome_target = args.datasets_dir + '/' + file_yaml['exome_target']
    else:
        exome_target = args.datasets_dir + '/' + args.exome_target
    if yaml_file != 'null' and 'exome_target_bed' in file_yaml.keys():
        exome_target_bed = args.datasets_dir + '/' + file_yaml['exome_target_bed']
    else:
        exome_target_bed = args.datasets_dir + '/' + args.exome_target_bed
    #--GATK Ha
    if yaml_file != 'null' and 'erc' in file_yaml.keys():
        erc = file_yaml['erc']
    else:
        erc = args.erc
    #--variant filters
    if yaml_file != 'null' and 'snp_filter' in file_yaml.keys():
        snp_filter = file_yaml['snp_filter']
    else:
        snp_filter = args.snp_filter
    if yaml_file != 'null' and 'indel_filter' in file_yaml.keys():
        indel_filter = file_yaml['indel_filter']
    else:
        indel_filter = args.indel_filter
    #--annotation
    if yaml_file != 'null' and 'ref_ens' in file_yaml.keys():
        ref_ens = args.datasets_dir + '/' + file_yaml['ref_ens']
    else:
        ref_ens = args.datasets_dir + '/' + args.anno_geneID
    if yaml_file != 'null' and 'db_cosmic' in file_yaml.keys():
        db_cosmic = args.datasets_dir + '/' + file_yaml['db_cosmic']
    else:
        db_cosmic = args.datasets_dir + '/' + args.db_cosmic
    if yaml_file != 'null' and 'db_clinvar' in file_yaml.keys():
        db_clinvar = args.datasets_dir + '/' + file_yaml['db_clinvar']
    else:
        db_clinvar = args.datasets_dir + '/' + args.db_clinvar
    if yaml_file != 'null' and 'db_g1000' in file_yaml.keys():
        db_g1000 = args.datasets_dir + '/' + file_yaml['db_g1000']
    else:
        db_g1000 =  args.datasets_dir + '/' + args.db_g1000
    #--subprocess of the pipeline
    if yaml_file != 'null' and 'tools' in file_yaml.keys():
        tools = file_yaml['tools']
    else:
        tools = args.tools
    #--the variant calling model
    if yaml_file != 'null' and 'calling' in file_yaml.keys():
        calling = file_yaml['calling']
    #---check the cpus and memery
    time_start_run = time.time()
    #---
    #threads = args.threads
    
    num_thread_of_process = int(cores/threads)
    if num_thread_of_process*threads > 0.7*cores:
        num_thread_of_process = num_thread_of_process -1
    memory_size_of_process = int(free*0.6/(1024*1024*1024*num_thread_of_process))
    #---check the outputdir
    #--variant calling
    memory_size = memory_size_of_process
    if memory_size > 4:
        memory_size = 4
        memory_size = '-Xmx' + str(memory_size) + 'G '
    else:
        if threads > 1:
            memory_size = '-Xmx' + str(memory_size) + 'G '
        else:
            memory_size = '-Xmx4G '
    print(memory_size)

    #--check the variant calling model
    #if calling =='GATK':
    #    variant_call = gatk_dir
    #    tabix = 'null'
    #    bgzip = 'null'
    #elif calling == 'strelka2':
    #    variant_call = strelka2_dir
    #    total_ref_fa_file = total_ref_chrom_fa_file
    #elif
    #-----
    #--read the sample list 
    file_list1 = open(file_list,'r')
    path_sampleID = []
    for line in file_list1.readlines():
        path_sampleID.append([line.strip(), out_dir, fastqc_dir, primers_file, exome_target_bed,
                                      min_read_len, common_seq1, common_seq2, num_threads, edit_dist, min_mapq,
                                      max_soft_clip, max_dist, memory_size, snp_filter, indel_filter, ref_ens,
                                      bwa_dir, samtools_dir, umitools_dir, gatk_dir, ref_index_name,
                                      ref_fa_file, total_ref_fa_file, total_ref_fa_dict, known_sites, erc,db_cosmic, 
                                      db_clinvar, db_g1000, tools, exome_target,calling,tabix,bgzip,bcftools_dir,varsan2_dir,
                                      strelka2_dir,total_ref_chrom_fa_file
                                      ])
    
    pool = multiprocessing.Pool(processes=threads)
    results = pool.map(main_run_germline_variant_calling, path_sampleID)
    print("--" * 20)
    pool.close()   # 关闭pool, 则不会有新的进程添加进去
    pool.join()    # 必须在join之前close, 然后join等待pool中所有的线程执行完毕
    print("All process done.")
    print("Return results: ")
    #----pipeline log file
    log_dir_run = out_dir + '/' + 'log/'
    if not os.path.exists(log_dir_run):
        os.makedirs(log_dir_run)
    logger_pipeline_process_run = log_dir_run
    logger_pipeline_error_run = log_dir_run
    for i in results:
        print(i)   # 获得进程的执行结果
        store_pipeline_logs(logger_pipeline_process_run,'null',i+'\n')
    store_pipeline_logs(logger_pipeline_process_run,'null','All samples are completed after {0} min.\n'.format(('%.2f' % ((time.time() - time_start_run)/60))))
    #----
    if tools in ['all','statis','summary']:
        print("please check the others subprocess results!")
        print("Test statistics module!\n")
        #-merge the statis info 
        sample_preinfo = out_dir +'/pre_summary.txt'
        qc = out_dir +'/QC.statistics.list'
        os.system('ls {0}/*/QC/*QC.statistics.txt > {1}'.format(out_dir,qc ))
        trim_qc = out_dir +'/' + 'QC_Trim.statistics.txt'
        os.system('ls {0}/*/statistics/trim_QC/*Trim.statistics.txt > {1}'.format(out_dir,trim_qc))
        trim_statis = out_dir + '/trim.basic_stats.txt'
        os.system('ls {0}/*/undetermined/*_basic_stats.txt > {1}'.format(out_dir,trim_statis))
        filter_statis = out_dir + '/align.basic_stats.txt'
        os.system('ls {0}/*/filtered/*_align_stats.txt > {1}'.format(out_dir,filter_statis))
        primer_statis = out_dir + '/primer.basic_stats.txt'
        os.system('ls {0}/*/filtered/*_primer_stats.csv > {1}'.format(out_dir,primer_statis))
        umi_statis = out_dir + '/umis.basic_stats.txt'
        os.system('ls {0}/*/clustered/*_deduplicated_per_umi.tsv > {1}'.format(out_dir,umi_statis))
        align_base = out_dir + '/basesDepthInRegion.stats.txt'
        os.system('ls {0}/*/statistics/*_Align_basesDepthInRegion.txt > {1}'.format(out_dir,align_base))
        merge_statistics(sample_preinfo,exome_target,qc,trim_qc,trim_statis,filter_statis, primer_statis, umi_statis, align_base)
    #----
