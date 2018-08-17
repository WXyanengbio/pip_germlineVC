import os
import sys
import time
import shlex
import subprocess
import argparse
import yaml
import pysam
import logging

def setup_logger(name, log_file, formatter, level=logging.DEBUG):
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    
    return logger


def store_pipeline_logs(log_dir):
    formatter_pipeline_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_pipeline_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_pipeline_process = setup_logger('Running Messages of pipeline',
                                       log_dir + '/process.log',
                                       formatter_pipeline_process)
    logger_pipeline_errors = setup_logger('Errors & Warnings of pipeline',
                                      log_dir + '/errors.log',
                                      formatter_pipeline_errors)
    return logger_pipeline_process, logger_pipeline_errors

def script_information():
    print ("\nApplication: pipelines of benchmarking \n")
    print ("=====================================================================")
    print ("Required environment: \n\
            python\n \
            JAVA\n \
            R\n \
            bwa\n \
            samtools\n \
            GATK\n \
            bcftools\n \
            VarDict\n \
            Strelka2\n \
            VarScan2")
    print ("\n")
    print("To get help , type:\n")
    print("python3.6 run1.py -h")

parser = argparse.ArgumentParser(usage = "\n\npython3.6 %(prog)s --yaml")

parser.add_argument("--yaml", 
                    type =str,
                    default = 'null',
                    help = "The yaml file of the parameters")
parser.add_argument("--vready_bam", 
                    help = "Path to input test bam", 
                    type = str)
parser.add_argument("--sample_name", 
                    help = "the sample name of raw reads", 
                    type = str)
parser.add_argument("--tailname", 
                    type =str,
                    default = 'null',
                    help = "the tailname of sample raw read")
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
                    help = "the inos.path.dirnamestall path of fastQC")
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
                    default= 'genome/target_breast/target_breast.refSeq.fa',
                    help = "the path of ref fasta")

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
                    default = 'DHS-001Z_primers_target_BRCA.csv',
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
                    default = 'DP < 10 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0',
                    help = "add parameters for filtering SNPs"
                    )
parser.add_argument("--indel_filter",
                    type = str, 
                    default = 'DP < 10 || QD < 2.0 || FS > 200 || ReadPosRankSum < -20.0 || SOR > 10.0',
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
#-------------------------------
def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr

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

    outfile.close()
def sam_to_bam(gatk_dir, samtools_dir, vready_bam, sample, output, memory_size, exome_target_bed, 
               ref_fa_file,ref_fa_dict,
               known_sites, bqsr):
    #target bed to exon intervel list
    if os.path.basename(exome_target_bed) != 'all':
        exon_interval = output + '/' + 'target_interval.list'
        command_count_r3 = '{0} --java-options "{1}" BedToIntervalList -I {2} -O {3} -SD {4} --showHidden true'.format(
            gatk_dir, memory_size, exome_target_bed, exon_interval, ref_fa_dict)
        stdout, stderr = stdout_err(command_count_r3)
    bam = vready_bam
    #--sam to bem
    if '_sorted.bam' not in vready_bam:
        bam = output + '/' + sample + '_sorted.bam'
        command_count_r4 = samtools_dir+' sort '+ vready_bam +' > '+ bam
        os.system(command_count_r4)
    #build index of bam By samtools
    if not os.path.isfile(vready_bam +'.bai'):
        command_count1 = '{0} index {1}'.format(samtools_dir, vready_bam)
        os.system(command_count1)
    #MarkDuplicates
    mark_bam = output + '/' + sample + '_sorted.MarkDuplicates.bam'
    bam_metrics = output + '/' + sample + '_sorted.MarkDuplicates.metrics'
    command_count2 ='{0} --java-options "{1}" MarkDuplicates -I {2} -O {3} -M {4} --REMOVE_SEQUENCING_DUPLICATES false --showHidden true'.format(
        gatk_dir, memory_size, bam, mark_bam, bam_metrics)
    stdout, stderr = stdout_err(command_count2)
    #add row 'RG' in head
    #mark_rg_bam = output + '/' + sample + '_sorted.MarkDuplicates.RG.bam'
    #command_count3 = '{0} --java-options "{1}" AddOrReplaceReadGroups -I {2} -O {3} -LB lib1 -PL illumina -PU unit1 -SM {4} --showHidden true'.format(
    ##    gatk_dir, memory_size, mark_bam, mark_rg_bam, bam_metrics)
    #stdout, stderr = stdout_err(command_count3)
    #build the index of the marked_fixed bam
    command_count4 = '{0} index {1}'.format(samtools_dir, mark_bam)
    os.system(command_count4)
    #logger_g_variantcalling_process.info("Samtools build the index of  marked bam----cost %.2f min.", (time.time()-time_start1)/60)
    #---------------------------------------
    bam_to_variant = mark_bam
    #Base(Quality Score) Recalibration
    print(bqsr)
    if bqsr =='yes':
        original_recal_data_table =  output + '/' + sample + '_original.recal_data.csv'
        original_bgsr_bam =  output + '/' + sample + '_sorted.MarkDuplicates.originalBQSR.bam'
        recal_data_table =  output + '/' + sample + '_bqsr.recal_data.csv'
        bgsr_bam =  output + '/' + sample + '_sorted.MarkDuplicates.BQSR.bam'
        if os.path.basename(exome_target_bed) != 'all':
            #BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--first
            command_count5_1 ='{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} -L {4} --known-sites {5} -O {6} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, exon_interval, known_sites, original_recal_data_table)
            #ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--first
            command_count7_1 ='{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -L {5} -O {6} --showHidden true'.format(
              gatk_dir, memory_size, ref_fa_file, mark_bam, original_recal_data_table, exon_interval, original_bgsr_bam)
            #BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model--second
            command_count5_2 ='{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} -L {4} --known-sites {5} -O {6} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, exon_interval, known_sites, recal_data_table)
            #ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--second
            command_count7_2 ='{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -L {5} -O {6} --showHidden true'.format(
             gatk_dir, memory_size, ref_fa_file, mark_bam, recal_data_table, exon_interval, bgsr_bam)
        else:
            #BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model
            command_count5_1 ='{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} --known-sites {4} -O {5} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, known_sites, original_recal_data_table)
            #ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--first
            command_count7_1 ='{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -O {5} --showHidden true'.format(
              gatk_dir, memory_size, ref_fa_file, mark_bam, original_recal_data_table, original_bgsr_bam)
            #BaseRecalibrator---Generate Base Quality Score Recalibration (BQSR) model
            command_count5_2 ='{0} --java-options "{1}" BaseRecalibrator -R {2} -I {3} --known-sites {4} -O {5} --bqsr-baq-gap-open-penalty 40 --showHidden true'.format(
                gatk_dir, memory_size, ref_fa_file, mark_bam, known_sites, recal_data_table)
            #ApplyBQSR--Apply Base Quality Score Recalibration (BQSR) model--second
            command_count7_2 ='{0} --java-options "{1}" ApplyBQSR -R {2} -I {3} -bqsr {4} -O {5} --showHidden true'.format(
               gatk_dir, memory_size, ref_fa_file, mark_bam, recal_data_table, bgsr_bam)
        stdout, stderr = stdout_err(command_count5_1)
        stdout, stderr = stdout_err(command_count7_1)
        stdout, stderr = stdout_err(command_count5_2)
        stdout, stderr = stdout_err(command_count7_2)
        base_reports = output + '/' + sample + '.BQSR.before_VS_after.pdf'
        base_reports_csv = output + '/' + sample + '.BQSR.before_VS_after.csv'
        command_count6 ='{0} --java-options "{1}" AnalyzeCovariates --before-report-file {2} --after-report-file {3} --plots-report-file {4} --intermediate-csv-file {5}'.format(
            gatk_dir, memory_size, original_recal_data_table, recal_data_table, base_reports,base_reports_csv)
        stdout, stderr = stdout_err(command_count6)
    else:
        bgsr_bam= 'null'
    return mark_bam, bgsr_bam

def germline_variant_calling(gatk_dir, marked_BQSR_bam,
                             sample, output, 
                             memory_size, ref_fa_file,
                             exon_interval, erc,
                             snp_filter,indel_filter):
    #logger_g_variantcalling_process, #logger_g_variantcalling_errors):
    if exon_interval != 'all':
         command_count ='{0} --java-options "{1}" HaplotypeCaller -R {2} -I {3} -L {4} --minimum-mapping-quality 20 --showHidden true'.format(
            gatk_dir, memory_size, ref_fa_file, marked_BQSR_bam, exon_interval)
    else:
         command_count ='{0} --java-options "{1}" HaplotypeCaller -R {2} -I {3} --minimum-mapping-quality 20 --showHidden true'.format(
            gatk_dir, memory_size, ref_fa_file, marked_BQSR_bam)
    #logger_g_variantcalling_process.info('Begin to confirm the options parameters of running HaplotypeCaller.')
    vcf1 =  output + '/' + sample + '.raw_variants.vcf'
    if erc == 'no':
        #logger_g_variantcalling_process.info('Runs HaplotypeCaller in default mode on a single input BAM file containing sequence data!')
        vcf =  output + '/' + sample + '.raw_variants.vcf'
        command_count = command_count + ' -O ' + vcf
    else:
        if erc == 'GVCF':
            #logger_g_variantcalling_process.info('Runs HaplotypeCaller in GVCF mode!')
            vcf =  output + '/' + sample + '.raw_variants.g.vcf'
            command_count = command_count + ' -O {0} -ERC {1}'.format(vcf, erc)
    #logger_g_variantcalling_process.info('Begin to do germline variant calling.')
    stdout, stderr = stdout_err(command_count)
    #-Joint-Call Cohort
    if erc == 'GVCF':
        command_count1 ='{0} --java-options "{1}" GenotypeGVCFs -R {2} --variant {3} -O {4}'.format(
        gatk_dir, memory_size, ref_fa_file, vcf, vcf1)
        time_start1 = time.time()
        stdout, stderr = stdout_err(command_count1)
    #-SNP
    snp_vcf = output + '/' + sample + '.raw_variants_SNP.vcf'
    command_count2 ='{0} --java-options "{1}" SelectVariants -R {2} --variant {3} -O {4} --select-type-to-include SNP --showHidden true'.format(
        gatk_dir, memory_size, ref_fa_file, vcf1, snp_vcf)
    stdout, stderr = stdout_err(command_count2)
    #filter variant in SNP
    filter_snp_vcfs_tmp = output + '/' + sample + '.filter_SNP_tmp.vcf'
    command_count2fstmp ='{0} --java-options "{1}" VariantFiltration -R {2} --variant {3} -O {4} \
         --filter-expression "{5}" --filter-name "my_snp_filter" --showHidden true'.format(
        gatk_dir, memory_size, ref_fa_file, snp_vcf, filter_snp_vcfs_tmp, snp_filter)
    stdout, stderr = stdout_err(command_count2fstmp)
    #--check the raw SNPs number in tmp vcf
    if check_variant_exit(filter_snp_vcfs_tmp) == 0:
        print('There is no SNPs in {0}'.format(sample))
        #logger_g_variantcalling_process.info('There is no SNPs in {0}'.format(sample))
    elif check_variant_exit(filter_snp_vcfs_tmp) == 1:
        filter_snp_vcfs = output + '/' + sample + '.filter_SNPvardict_dir.vcf'
        command_count2fs ='{0} --java-options "{1}" SelectVariants --variant {2} -O {3} --exclude-filtered true --showHidden true'.format(
           gatk_dir, memory_size, filter_snp_vcfs_tmp,filter_snp_vcfs)
        stdout, stderr = stdout_err(command_count2fs)
    #indel
    indel_vcf = output + '/' + sample + '.raw_variants_indel.vcf'
    command_count3 ='{0} --java-options "{1}" SelectVariants -R {2} --variant {3} -O {4} --select-type-to-include INDEL --showHidden true'.format(
        gatk_dir, memory_size, ref_fa_file, vcf1, indel_vcf)
    stdout, stderr = stdout_err(command_count3)
    #filter variant in indel
    filter_indel_vcfi_tmp = output + '/' + sample + '.filter_indel_tmp.vcf'
    command_count2fitmp ='{0} --java-options "{1}" VariantFiltration -R {2} --variant {3} -O {4} \
      --filter-expression "{5}" --filter-name "my_indel_filter" --showHidden true'.format(
        gatk_dir, memory_size, ref_fa_file, indel_vcf, filter_indel_vcfi_tmp, indel_filter)
    stdout, stderr = stdout_err(command_count2fitmp)
    #--check the raw indel number in tmp vcf
    if check_variant_exit(filter_indel_vcfi_tmp) == 0:
        print('There is no SNPs in {0}'.format(sample))
        #logger_g_variantcalling_process.info('There is no indel in {0}'.format(sample))
    elif check_variant_exit(filter_indel_vcfi_tmp) == 1:
        filter_indel_vcfi = output + '/' + sample + '.filter_indel.vcf'
        command_count2fi ='{0} --java-options "{1}" SelectVariants --variant {2} -O {3} --exclude-filtered true --showHidden true'.format(
            gatk_dir, memory_size,filter_indel_vcfi_tmp, filter_indel_vcfi)
        stdout, stderr = stdout_err(command_count2fi)


def hap_py(benchmark_dir,truth_vcf, filter_vcf, confident_region_bed, 
           benchmarking_dir, total_ref_fa_file, exon_interval, num_threads):
    if exon_interval != 'all':
        exon_bed = benchmarking_dir + '/' + os.path.basename(exon_interval).rstrip(".list") + ".bed"
        if not os.path.isfile(exon_bed):
            #logger_benchmark_errors.error("%s does not exist!\n", exon_bed)
            print(exon_bed + ' does not exist!')
            cmd1 = 'grep ^chr {0} | cut -f 1-3 > {1}'.format(exon_interval, exon_bed)
            #logger_benchmark_process.info('Build the target bed to benchmark.')
            os.system(cmd1)
    if not os.path.isfile(filter_vcf):
        #logger_benchmark_errors.error("%s does not exist!\n", filter_vcf)
        print(filter_vcf + ' does not exist!')
    else:
        benchmarout = benchmarking_dir + '/' + os.path.basename(filter_vcf).rstrip(".vcf")
        if exon_interval != 'all':
            cmd2 = '{0} {1} {2} -f {3} -r {4} -T {5} -o {6} --threads {7}'.format(benchmark_dir, truth_vcf, 
                                                                                  filter_vcf,confident_region_bed,total_ref_fa_file ,
                                                                                  exon_bed, benchmarout,num_threads)
        else:
            cmd2 = '{0} {1} {2} -f {3} -r {4} -o {5} --threads {6}'.format(benchmark_dir, truth_vcf, 
                                                                           filter_vcf,confident_region_bed,total_ref_fa_file , benchmarout,num_threads)
        #logger_benchmark_process.info('hap.py benchmark.')
        os.system(cmd2)

def samtools_call(samtools_dir, bcftools_dir, bam, sample, outputdir, total_ref_fa_file):
    bcf = outputdir + '/'+ sample + '_raw.bcf'
    command1 =bcftools_dir +' mpileup -Ob -o '+ bcf+' -f '+ total_ref_fa_file + ' ' + bam
    os.system(command1)
    vcf = outputdir + '/'+ sample + '_raw.vcf'
    command2 = bcftools_dir + ' call -vmO v -o ' + vcf + ' ' + bcf
    os.system(command2)
    vchk = outputdir + '/'+ sample + '_raw.vchk '
    command3 = bcftools_dir + ' stats ' + vcf + ' > ' + vchk
    os.system(command3)
    command4 = os.path.dirname(bcftools_dir)+'/misc/plot-vcfstats '+ vchk + '-p '+ outputdir + '/plots/'
    os.system(command4)
    bcftools_filter_vcf=  outputdir + '/'+ sample + '_filter.vcf'
    command5 = bcftools_dir + ' filter -O v -o '+  bcftools_filter_vcf + ' -s LOWQUAL -e \'QUAL<30 || FMT/DP <5\' '+  vcf
    os.system(command4)
    #command4 = perl -ne 'print $_ if /DP4=(\d+),(\d+),(\d+),(\d+)/ && ($3+$4)>=10 && ($3+$4)/($1+$2+$3+$4)>=0.8' var.raw.vcf > snp_indel.final.vcf
    return vcf

def vardict_call(vardict_dir, total_ref_fa_file, outputdir,sample, bam,exome_target_bed, af):
    os.system('cp {0} {1}'.format(bam, outputdir))
    os.system('cp {0} {1}'.format(bam +'.bai', outputdir))
    os.chdir(outputdir)
    print(os.getcwd())
    bam1=outputdir+'/'+os.path.basename(bam)
    tmp=bam1.rstrip('.bam') +'_stdin.txt'
    command = '{0} -G {1} -N {2} -b {3} -f {4} -c 1 -S 2 -E 3 {5} > {6}'.format(
        vardict_dir,total_ref_fa_file, bam1.rstrip('.bam'), bam1, af, exome_target_bed ,tmp)
    print(command)
    os.system(command)
    os.system('cp {0} {1}'.format(tmp, outputdir+'/stdin.txt'))
    os.system('cd {0}'.format(outputdir))
    vcf=bam1.rstrip('.bam')+'_tmp.vcf'
    #---teststrandbias.R modified from VarDict
    scriptdir = '/home/dell/Works/Projects/pip_germlineVC/pipelines/benchmark'
    command1=  'Rscript ' + scriptdir+ '/teststrandbias.R' + ' | ' + os.path.dirname(vardict_dir)+'/var2vcf_valid.pl'+ ' -N '+ bam1.rstrip('.bam') + ' -f 0.01 >  ' + vcf
    print(command1)
    os.system(command1)
    #---remove the last line in the raw vcf file
    tmp_vcf = open(vcf)
    raw_vcf = bam1.rstrip('.bam')+'_raw.vcf'
    sam_out = open(raw_vcf, 'w')
    for line in tmp_vcf.readlines():
        if line.startswith('#') or line.startswith('chr'):
            sam_out.write(line)
    sam_out.close()
    tmp_vcf.close()
    os.system('rm -rf {0}'.format(tmp_vcf))
    
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
    #print(stderr)
    
def varsan2_call(samtools_dir, varsan2_dir,total_ref_fa_file,
                 outputdir,
                 sample, bam):
    pileup= outputdir+ '/' + sample+ '.mpileup'
    command1= samtools_dir + ' mpileup -B -f '+ total_ref_fa_file +' ' + bam + ' > ' + pileup
    os.system(command1)
    snp= outputdir+ '/' + sample+ '.snp.vcf'
    indel= outputdir+ '/' + sample+ '.indel.vcf'
    command2= 'java -jar ' + varsan2_dir + ' mpileup2snp ' + pileup + ' > ' + snp
    os.system(command2)
    #stdout, stderr = stdout_err(command2)
    command3= 'java -jar ' + varsan2_dir + ' mpileup2indel ' + pileup + ' > ' + indel
    os.system(command3)
    #print(stderr)
    
def get_target_reads_from_bam(target_region,infile,outfile):
    chrom, start, end = target_region
    for gtf in infile.fetch(chrom[3:],int(start),int(end)):
        outfile.write(gtf)
    
def reformat_sam(alignment_sam, output_sam):
    sam = open(alignment_sam)
    sam.readline()
    sam_out = open(output_sam, 'w')
    for row in sam:
        if row[0:3] == '@HD':
            sam_out.write(row)
            continue
        if row[0:3] == '@SQ':
            head, ref, ln, m5 = row.split('\t')
            sn, chrom = ref.split(':')
            if str.isdigit(chrom):
                #chrom = 'chr'+ chrom
                ref= sn+ ':' + 'chr'+ chrom
                sam_out.write('\t'.join([head, ref, ln]) + '\n')
            continue
        if row[0:3] == '@PG':
            sam_out.write(row)
            continue
        if row[0:3] == '@RG':
            sam_out.write(row)
            continue
        qname, flag, rname= row.strip().split()[0:3]
        #print(rname)
        chrom = 'chr'+rname
        sam_out.write('\t'.join([qname, flag, chrom]) + '\t' + '\t'.join(row.split()[3:]) + '\n')
    sam.close()
    sam_out.close()
    
def main():
    #time cost
    time_start1 = time.time()
    #---input
    yaml_file = args.yaml
    if yaml_file != 'null':
        f = open(yaml_file)
        file_yaml = yaml.load(f)
        #print(file_yaml)
    if yaml_file != 'null' and 'vready_bam' in file_yaml.keys():
        vready_bam = file_yaml['vready_bam']
    else:
        vready_bam = args.vready_bam
    if yaml_file != 'null' and 'sample_name' in file_yaml.keys():
        sample = file_yaml['sample_name']
    else:
        sample = args.sample_name
    if yaml_file != 'null' and 'tailname' in file_yaml.keys():
        tailname = file_yaml['tailname']
    else:
        tailname = args.tailname
    
    if tailname != 'null':
        sample = sample + '_' + tailname
    print(sample)
    #---check the outputdir
    if yaml_file != 'null' and 'output' in file_yaml.keys():
        out_dir = file_yaml['output']
        print(out_dir)
    elif args.output is 'null':
        out_dir = sample
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #----pipeline log file
    log_dir = out_dir + '/' + 'log/'
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    logger_pipeline_process, logger_pipeline_errors = store_pipeline_logs(log_dir)
    #QC
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
    if yaml_file != 'null' and 'bcftools_dir' in file_yaml.keys():
        bcftools_dir = file_yaml['bcftools_dir']
    if yaml_file != 'null' and 'vardict_dir' in file_yaml.keys():
        vardict_dir = file_yaml['vardict_dir']
    if yaml_file != 'null' and 'strelka2_dir' in file_yaml.keys():
        strelka2_dir = file_yaml['strelka2_dir']
    if yaml_file != 'null' and 'bgzip' in file_yaml.keys():
        bgzip = file_yaml['bgzip']
    if yaml_file != 'null' and 'varsan2_dir' in file_yaml.keys():
        varsan2_dir = file_yaml['varsan2_dir']
    if yaml_file != 'null' and 'tabix' in file_yaml.keys():
        tabix = file_yaml['tabix']
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
    if yaml_file != 'null' and 'total_ref_fa_dict' in file_yaml.keys():
        total_ref_fa_dict = args.datasets_dir + '/' + file_yaml['total_ref_fa_dict']
    else:
        total_ref_fa_dict = args.datasets_dir + '/' + args.total_ref_fa_dict
    if yaml_file != 'null' and 'total_ref_fa_file' in file_yaml.keys():
        total_ref_fa_file = args.datasets_dir + '/' + file_yaml['total_ref_fa_file']
    else:
        total_ref_fa_file = args.datasets_dir + '/' + args.total_ref_fa_file
    
    if yaml_file != 'null' and 'total_ref_chrom_fa_file' in file_yaml.keys():
        total_ref_chrom_fa_file = args.datasets_dir + '/' + file_yaml['total_ref_chrom_fa_file']
    
    if yaml_file != 'null' and 'known_sites' in file_yaml.keys():
        known_sites = file_yaml['known_sites']
    else:
        known_sites = args.known_sites
    #read_length = args.read_length
    if yaml_file != 'null' and 'exome_target_bed' in file_yaml.keys():
        exome_target_bed = args.datasets_dir + '/' + file_yaml['exome_target_bed']
    else:
        exome_target_bed = args.datasets_dir + '/' + args.exome_target_bed
    if yaml_file != 'null' and 'erc' in file_yaml.keys():
        erc = file_yaml['erc']
    else:
        erc = args.erc
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
    #--benchmark
    if yaml_file != 'null' and 'benchmark_dir' in file_yaml.keys():
        benchmark_dir = file_yaml['benchmark_dir']
    else:
        benchmark_dir = args.benchmark_dir
    if yaml_file != 'null' and 'confident_region_bed' in file_yaml.keys():
        confident_region_bed = file_yaml['confident_region_bed']
    else:
        confident_region_bed = args.datasets_dir + '/' + args.confident_region_bed
    if yaml_file != 'null' and 'truth_vcf' in file_yaml.keys():
        truth_vcf = file_yaml['truth_vcf']
    else:
        truth_vcf = args.datasets_dir + '/' + args.truth_vcf
    #---
    if yaml_file != 'null' and 'te\n \st' in file_yaml.keys():
        test_level = file_yaml['test']
    else:
        test_level = args.test
    #--------
    bam_dir = out_dir + '/'+ 'bam'
    if not os.path.exists(bam_dir):
        os.makedirs(bam_dir)
    #bqsr='no'
    print(vready_bam)
    target_vready_bam = bam_dir+'/'+os.path.basename(vready_bam)
    os.system('cp {0} {1}'.format(vready_bam,bam_dir))
    target_region = []
    for line in open(exome_target_bed,'r'):
        if line.startswith('chr'):
            chrom, start, end = line.strip().split('\t')
            target_region.append([chrom, start, end])
    infile=pysam.Samfile(vready_bam,'rb')
    outfile= pysam.Samfile(target_vready_bam,'wb',template=infile)
    for region in target_region:
        get_target_reads_from_bam(region, infile, outfile)
    #outfile.close()
    #--------
    target_vready_sam = target_vready_bam.rstrip('_sorted.bam') + '.sam'
    command1= samtools_dir + ' view -h ' + target_vready_bam + ' > ' + target_vready_sam
    os.system(command1)
    target_chrom_sam = target_vready_bam.rstrip('_sorted.bam') + '_vready.sam'
    reformat_sam(target_vready_sam, target_chrom_sam)
    #os.system('head -100 {0}'.format(target_chrom_sam))
    target_sort_bam = target_vready_bam.rstrip('_sorted.bam') + '_vready_sorted.bam'
    command2 = samtools_dir + ' sort ' + target_chrom_sam + ' > ' + target_sort_bam
    os.system(command2)
    #---modify the known-sites
    known_sites = known_sites.replace(',' , ' --known-sites ' + args.datasets_dir + '/'  )
    known_sites = args.datasets_dir + '/' + known_sites
    #----
    bqsr='yes'
    bam_to_variant, baser_bam_to_variant = sam_to_bam(gatk_dir, samtools_dir,
                              target_sort_bam, sample,
                              bam_dir, memory_size,
                              exome_target_bed, 
                              total_ref_fa_file, total_ref_fa_dict,
                              known_sites, bqsr)
    
    
    #GATK
      #time cost
    time_start = time.time()
    exon_interval = bam_dir + '/' + 'target_interval.list'
    germline_vc_dir= out_dir + '/'+ 'gatk'
    if not os.path.exists(germline_vc_dir):
        os.makedirs(germline_vc_dir)
    germline_variant_calling(gatk_dir, bam_to_variant,
                             sample, germline_vc_dir, 
                             memory_size, total_ref_fa_file, 
                             exon_interval, erc,
                             snp_filter,indel_filter)
    logger_pipeline_process.info('GATK is completed after %.2f min.', (time.time()-time_start)/60)
    
    
    #GATK-without BQSR
      #time cost
    time_start = time.time()
    germline_vc_dir1= out_dir + '/'+ 'gatk_with_bqsr'
    if not os.path.exists(germline_vc_dir1):
        os.makedirs(germline_vc_dir1)
    germline_variant_calling(gatk_dir, baser_bam_to_variant,
                             sample, germline_vc_dir1, 
                             memory_size, total_ref_fa_file, 
                             exon_interval, erc,
                             snp_filter,indel_filter)
    logger_pipeline_process.info('GATK_bqsr is completed after %.2f min.', (time.time()-time_start)/60)
    
    
    #samtools+bcftools
      #time cost
    time_start = time.time()
    samtools_dir1= out_dir + '/'+ 'samtools'
    if not os.path.exists(samtools_dir1):
        os.makedirs(samtools_dir1)
    samtools_call(samtools_dir, bcftools_dir, bam_to_variant, sample, samtools_dir1, total_ref_fa_file)
    logger_pipeline_process.info('Samtools is completed after %.2f min.', (time.time()-time_start)/60)
    
    
    #VarScan2+samtools
      #time cost
    time_start = time.time()
    varsan2_dir1= out_dir + '/'+ 'varsan2'
    if not os.path.exists(varsan2_dir1):
        os.makedirs(varsan2_dir1)
    varsan2_call(samtools_dir, varsan2_dir,total_ref_fa_file,
                 varsan2_dir1,
                 sample, bam_to_variant)
    logger_pipeline_process.info('VarScan2 is completed after %.2f min.', (time.time()-time_start)/60)
    
    
    #VarDict
      #time cost
    time_start = time.time()
    vardict_dir1= out_dir + '/'+ 'vardict'
    if not os.path.exists(vardict_dir1):
        os.makedirs(vardict_dir1)
    os.system('cp {0} {1}'.format(bam_to_variant,vardict_dir1))
    af=0.01
    bam= bam_to_variant
    vardict_call(vardict_dir, total_ref_fa_file, vardict_dir1, sample, bam, exome_target_bed, af)
    logger_pipeline_process.info('VarDict is completed after %.2f min.', (time.time()-time_start)/60)
    
    
    #Strelka2
      #time cost
    time_start = time.time()
    strelka2_call_dir1= out_dir + '/'+ 'strelka2'
    if not os.path.exists(strelka2_call_dir1):
        os.makedirs(strelka2_call_dir1)
    strelka2_call(strelka2_dir,
                  bgzip,
                  tabix,
                  total_ref_chrom_fa_file,
                  strelka2_call_dir1,
                 sample, bam_to_variant,
                 exome_target_bed)
    logger_pipeline_process.info('Strelka2 is completed after %.2f min.', (time.time()-time_start)/60)
    
    
if __name__ == '__main__':
    main()
