__version__  =  '1.0'
__author__  =  'Wang Xian'


import logging
import os
import re
import sys
import time
import argparse
import yaml
import pandas as pd

base_paired = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

def read_database(cosmic,clinvar,g1000):
    cos = open(cosmic, 'r')
    dict_cos = {}
    cos.readline()
    for db1 in cos.readlines():
        # if not db.startswith('Gene_name'):
        Gene_name, Accession_Number, Gene_CDS_length, HGNC_ID, Sample_name, ID_sample, ID_tumour, Primary_site, \
        Site_subtype1, Site_subtype2, Site_subtype3, Primary_histology, Histology_subtype1, Histology_subtype2, \
        Histology_subtype3, Genome_wide_screen, Mutation_ID, Mutation_CDS, Mutation_AA, Mutation_Description, \
        Mutation_zygosity, LOH, GRCh, Chr, Start, End, Mutation_strand, SNP, Resistance_Mutation, FATHMM_prediction, \
        FATHMM_score, Mutation_somatic_status, Pubmed_PMID, ID_STUDY, Sample_source, Tumor_origin, Age = db1.strip().split(',')
        if Mutation_CDS is 'NS':
            continue
        if 'del' in Mutation_CDS:
            Change = Mutation_CDS[Mutation_CDS.find('del'):]
        elif 'ins' in Mutation_CDS:
            Change = Mutation_CDS[Mutation_CDS.find('ins'):]
        elif '>' in Mutation_CDS:
            Change = Mutation_CDS[Mutation_CDS.find('>')-1:]
        key1 = [Chr, Start ,Change]
        value1 = [Mutation_ID, Mutation_Description, Accession_Number, Gene_name, Gene_CDS_length,
                  Mutation_zygosity, LOH, Mutation_strand, Mutation_CDS, Mutation_AA, FATHMM_prediction, FATHMM_score,
                  Mutation_somatic_status]
        dict_cos[','.join(key1)] = ','.join(value1)
    clin = open(clinvar, 'r')
    dict_clin = {}
    clin.readline()
    for db2 in clin.readlines():
        genename, chr, start, end, geneid, alleleid, rs, pos, ref, alt, af_esp, af_exac, af_tgp, \
        clndn, clnhgvs, clnsig = db2.strip().split(',')
        if len(ref) == len(alt):
            change = ref + '>' + alt
        elif len(ref) > len(alt) and len(alt) == 1:
            change = 'del' + ref[1:]
        elif len(ref) < len(alt) and len(ref) == 1:
            change = 'ins' + alt[1:]
        else:
            change = 'del' + ref + 'ins' + alt
        key2 = [chr, pos, change]
        value2 = [geneid, rs, clndn, clnhgvs, clnsig]
        dict_clin[','.join(key2)] = ','.join(value2)
    
    genomes1000 = open(g1000, 'r')
    dict_g1000 = {}
    genomes1000.readline()

    for db3 in genomes1000.readlines():
        genename1, chr1, start1, end1, variant_type, ref1, alt1, rs1, eas_af, eur_af, amr_af, sas_af, afr_af, hgvs1 = db3.strip().split(',')
        if variant_type == 'SNV':
            change1 = ref1 + '>' + alt1
        elif variant_type == 'insertion':
            change1 = 'ins' + alt1
        elif variant_type == 'deletion':
            change1 = 'del' + ref1
        else:
            pass
        key3 = [chr1, start1, change1]
        value3 = [genename1, rs1, eas_af, eur_af, amr_af, sas_af, afr_af]
        dict_g1000[','.join(key3)] = ','.join(value3)
    return dict_cos, dict_clin, dict_g1000



def get_info(info):
    vcf_parameters = ['AC','AF','AN','BaseQRankSum','ClippingRankSum','DP','DS','END','ExcessHet','FS','InbreedingCoeff','MLEAC','MLEAF','MQ','MQRankSum','QD','RAW_MQ','ReadPosRankSum','SOR']
    subparas = info.split(';')
    subparas_pre = list(map(lambda tag: tag[:tag.find('=')], subparas))
    paras = []
    for para in vcf_parameters:
        if para in subparas_pre:
            tag = subparas[subparas_pre.index(para)]
            paras.append(tag[tag.find('=')+1:])
        else:
            paras.append('-')
    return paras

#split mutation. From VCF built by GATK HaplotypeCaller
def split_variant(filter):
    if filter != 'my_snp_filter' and filter != 'my_indel_filter':
        ac, af, an, baseqranksum, clippingranksum, dp, ds, end, excesshet, fs, inbreedingcoeff, mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor = get_info(filter)
    num_mt = len(alt.split(','))
    #af = detail.split(':')[2]
    if num_mt is 1:
        return [[chrom, pos, ref, alt, filter, ac, af, an, baseqranksum, clippingranksum, dp, ds, end, excesshet, fs, inbreedingcoeff, mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor]]

def annotation(dict_cos,dict_clin,dict_g1000,variant_vcf,annotated_csv,stats_file,vcf_soft):
    key_list = []
    key = ''
    change1 = ''
    num_in_clinvar = 0
    num_in_cosmic = 0
    num_in_g1000 = 0
    num_unmatch = 0
    var = open(variant_vcf, 'r')
    output = open(annotated_csv, 'w')
    #output.write(
    #    'CHR,POS,REF,ALT,FILTER,AC,AF,AN,VF,BaseQRankSum,ClippingRankSum,DP,DS,END,ExcessHet,FS,InbreedingCoeff,MLEAC,MLEAF,MQ,MQRankSum,QD,RAW_MQ,ReadPosRankSum,SOR,Gene_ID,RS_ID,CLNDN,HGVS,CLNSIG,'
    #    'COSMIC_ID,Mutation_Description,Feature_ID,Gene_Name,Gene_CDS_Length,Mutation_Zygosity,LOH,Mutation_Strand,'
    #    'HGVS.c,HGVS.p,FATHMM_Prediction,FATHMM_Score,Mutation_Somatic_Status,Gene_Name1,RS_ID1,EAS_AF,EUR_AF,AMR_AF,'
    #    'SAS_AF,AFR_AF\n')
    if vcf_soft == 'gatk':
        output.write(
        'CHR,POS,REF,ALT,QUAL,FILTER,VF,GT,AD,DP,FS,MQ,MQRankSum,QD,ReadPosRankSum,SOR,Gene_ID,RS_ID,CLNDN,HGVS,CLNSIG,'
        'COSMIC_ID,Mutation_Description,Feature_ID,Gene_Name,Gene_CDS_Length,Mutation_Zygosity,LOH,Mutation_Strand,'
        'HGVS.c,HGVS.p,FATHMM_Prediction,FATHMM_Score,Mutation_Somatic_Status,Gene_Name1,RS_ID1,EAS_AF,EUR_AF,AMR_AF,'
        'SAS_AF,AFR_AF\n')
    if vcf_soft == 'truth':
        output.write(
        'CHR,POS,REF,ALT,QUAL,FILTER,DETAIL,Gene_ID,RS_ID,CLNDN,HGVS,CLNSIG,'
        'COSMIC_ID,Mutation_Description,Feature_ID,Gene_Name,Gene_CDS_Length,Mutation_Zygosity,LOH,Mutation_Strand,'
        'HGVS.c,HGVS.p,FATHMM_Prediction,FATHMM_Score,Mutation_Somatic_Status,Gene_Name1,RS_ID1,EAS_AF,EUR_AF,AMR_AF,'
        'SAS_AF,AFR_AF\n')
    for line in var:
        if not line.startswith('#'):
            chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
            if len(ref) == len(alt):
                change = ref + '>' + alt
                change1 = base_paired[ref] + '>' + base_paired[alt]
            elif len(ref) > len(alt) and len(alt) == 1:
                change = 'del' + ref[1:]
            elif len(ref) < len(alt) and len(ref) == 1:
                change = 'ins' + alt[1:]
            else:
                change = 'del' + ref + 'ins' + alt
            key = chrom[3:] + ',' + pos + ',' + change
            key1 = chrom[3:] + ',' + pos + ',' + change1
            print(vcf_soft)
            if vcf_soft == 'gatk':
                if filter != 'my_snp_filter' and filter != 'my_indel_filter':
                # if filter != "":
                    ac, af, an, baseqranksum, clippingranksum, dp, ds, end, excesshet, fs, inbreedingcoeff, mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor = get_info(info)
                    if len(detail.split(':'))>5:
                        gt, ad, dp= detail.split(':')[0:3]
                        pl = detail.split(':')[len(detail.split(':'))-1]
                        gq = ':'.join(detail.split(':')[3:len(detail.split(':'))])
                    else:
                        gt, ad, dp, gq, pl = detail.split(':')
                    vf = str(round(int(ad.split(',')[1])/int(dp),4))
                    ad = ad.replace(',' , '/')
                    chrom = chrom[3:]
                    value = [chrom,pos,ref,alt,qual,filter,vf,gt,ad,dp,fs,mq,mqranksum,qd,readposranksum,sor]
            if vcf_soft == 'truth':
                value = [chrom[3:],pos,ref,alt,qual,filter,detail]
            unmatch = 0
            # drop duplicate variant
            if key in key_list:
                continue
            if key in dict_clin:
                new = ','.join(value) + ',' + dict_clin[key] + ','
                num_in_clinvar += 1
            else:
                new = ','.join(value) + ',-,-,-,-,-,'
                unmatch += 1
            if key in dict_cos:
                new += dict_cos[key] + ','
                num_in_cosmic += 1
            elif key1 in dict_cos:
                new += dict_cos[key1] + ','
                num_in_cosmic += 1
            else:
                new += '-,-,-,-,-,-,-,-,-,-,-,-,-,'
                unmatch += 1
            if key in dict_g1000:
                new += dict_g1000[key] + '\n'
                num_in_g1000 += 1
            else:
                new += '-,-,-,-,-,-,-\n'
                unmatch += 1
            # 3 databases both unmatch
            if unmatch == 3:
                num_unmatch += 1
            else:
                output.write(new)
            key_list.append(key)
    output.close()
    print ('The sample has %s variants.\n' % len(key_list))
    print ('%s variants in COSMIC database\n' % num_in_cosmic)
    print ('%s variants in Clinvar database\n' % num_in_clinvar)
    print ('%s variants in G1000 database\n' % num_in_g1000)
    print ('%s variants unmatch in cosmic,clinvar and g1000\n.' % num_unmatch)
    stats_out = open(stats_file,'w')
    stats_out.write('#The sample has %s variants.\n' % len(key_list))
    stats_out.write('#Type\tvariants\n')
    stats_out.write('Variants\t%s\n' % len(key_list))
    stats_out.write('COSMIC\t%s\n' % num_in_cosmic)
    stats_out.write('Clinvar\t%s\n' % num_in_clinvar)
    stats_out.write('G1000\t%s\n' % num_in_g1000)
    stats_out.write('unmatch\t%s\n' % num_unmatch)
    stats_out.close()

#match genename,ENSG and ENST from ensembl.
def fill_table(annotated_csv, annotated_csv_add, non_rs, non_cos, ref_ens):
    g2n = {}
    g2t = {}
    n2g = {}
    n2t = {}
    f1 = open(ref_ens,'r')
    for line in f1.readlines():
        l1,l2,l3 = line.strip().split(',')
        g2n[l2] = l1
        g2t[l2] = l3
        n2g[l1] = l2
        n2t[l1] = l3
    df = pd.read_csv(annotated_csv)
    subframe = df[['Gene_Name','Gene_ID','Feature_ID','Gene_Name1','RS_ID','RS_ID1']] #n,g,t,n1
    #for name,id,transcript in subframe.iterrows():
    for num in range(0, len(subframe)):
        #subframe.iloc[i]['Gene_Name'], subframe.iloc[i]['Gene_ID'], subframe.iloc[i]['Feature_ID']
        if subframe.iloc[num]['Gene_Name'] is '-' and subframe.iloc[num]['Feature_ID'] is '-' and subframe.iloc[num]['Gene_ID'] is '-':
            subframe.iloc[num]['Gene_Name'] = subframe.iloc[num]['Gene_Name1']
            #print(subframe.iloc[num]['Gene_Name1'])
            #print(n2g)
            subframe.iloc[num]['Gene_ID'] = n2g[subframe.iloc[num]['Gene_Name1']]
            subframe.iloc[num]['Feature_ID'] = n2t[subframe.iloc[num]['Gene_Name1']]
        elif subframe.iloc[num]['Gene_Name'] is '-' and subframe.iloc[num]['Feature_ID'] is '-':
            subframe.iloc[num]['Gene_Name'] = g2n[subframe.iloc[num]['Gene_ID']]
            subframe.iloc[num]['Feature_ID'] = g2t[subframe.iloc[num]['Gene_ID']]
        elif subframe.iloc[num]['Gene_ID'] is '-':
            subframe.iloc[num]['Gene_ID'] = n2g[subframe.iloc[num]['Gene_Name']]
        if subframe.iloc[num]['RS_ID'] is '-' and subframe.iloc[num]['RS_ID1'] is not '-':
            subframe.iloc[num]['RS_ID'] = subframe.iloc[num]['RS_ID1']
    name = subframe['Gene_Name']
    ensg = subframe['Gene_ID']
    enst = subframe['Feature_ID']
    rs = subframe['RS_ID']
    df.drop(labels=['Gene_Name'], axis=1, inplace=True)
    df.drop(labels=['Gene_ID'], axis=1, inplace=True)
    df.drop(labels=['Feature_ID'], axis=1, inplace=True)
    df.drop(labels=['Gene_Name1'], axis=1, inplace=True)
    df.drop(labels=['RS_ID'], axis=1, inplace=True)
    df.drop(labels=['RS_ID1'], axis=1, inplace=True)
    df.insert(7, 'Gene_Name', name)
    df.insert(8, 'Gene_ID', ensg)
    df.insert(9, 'Feature_ID', enst)
    df.insert(7, 'RS_ID', rs)
    df.to_csv(annotated_csv_add, index=False, sep=',')
    df[(True^df['RS_ID'].isin(['-']))].to_csv(non_rs, index=False, sep=',')
    df[(True^df['COSMIC_ID'].isin(['-']))].to_csv(non_cos, index=False, sep=',')

#-annotation main
def annotationmain(cosmic, clinvar, g1000, 
                   ref_ens,
                   vcf, sample,
                   output,vcf_soft):
    if '.vcf' in vcf:
        annotated_csv = output + '/' + sample + '.annotated.csv'
        annotated_csv_add = output + '/' + sample + '.annotated_ensembl.csv'
        non_rs = output + '/' + sample + '_non_rs.csv'
        non_cos = output + '/' + sample + '_non_cos.csv'
        stats_file = output + '/' + sample + '_annotate_stats.txt'
    #-read the annotation database
    if not os.path.isfile(vcf):
        print(vcf + ' does not exist!')
    else:
        dict_cos, dict_clin, dict_g1000 = read_database(cosmic,clinvar,g1000)
        #--annotation
        annotation(dict_cos, dict_clin, dict_g1000, vcf, annotated_csv, stats_file,vcf_soft)
        #--add the annotation
        fill_table(annotated_csv, annotated_csv_add, non_rs, non_cos, ref_ens)




def script_information():
    print ("\nApplication: pipelines of QIAseq Targeted DNA Panel\n")
    print ("=====================================================================")
    print ("Required environment: python \ JAVA\ R \ bwa \ samtools \ GATK \ UMI-tools")
    print ("\n")
    print("To get help , type:\n")
    print("python3.6 run1.py -h")

parser = argparse.ArgumentParser(usage = "\n\npython3.6 %(prog)s --yaml")

parser.add_argument("--yaml", 
                    type =str,
                    default = 'null',
                    help = "The yaml file of the parameters")
parser.add_argument("--source", 
                    help = "Path to input reads in FASTA format", 
                    type = str)
parser.add_argument("--sample_name", 
                    help = "the sample name of raw reads", 
                    type = str)
parser.add_argument("--datasets_dir", 
                    type =str,
                    default = '/home/dell/Works/Projects/Datasets',
                    help = "the tailname of sample raw read")
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

def main():
    #time cost
    time_start1 = time.time()
    #---input
    yaml_file = args.yaml
    if yaml_file != 'null':
        f = open(yaml_file)
        file_yaml = yaml.load(f)
        #print(file_yaml)
    if yaml_file != 'null' and 'source' in file_yaml.keys():
        source = file_yaml['source']
    else:
        source = args.source
    if yaml_file != 'null' and 'sample_name' in file_yaml.keys():
        sample = file_yaml['sample_name']
    else:
        sample = args.sample_name
    #---check the outputdir
    if yaml_file != 'null' and 'output' in file_yaml.keys():
        out_dir = file_yaml['output']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
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
    #benchmark_dir = args.benchmark_dir
    #confident_region_bed = args.datasets_dir + '/' + args.confident_region_bed
    #truth_vcf = args.datasets_dir + '/' + args.truth_vcf
    #---
    if yaml_file != 'null' and 'tools' in file_yaml.keys():
        tools = file_yaml['tools']
    else:
        tools = args.tools
    if yaml_file != 'null' and 'vcf_soft' in file_yaml.keys():
        vcf_soft = file_yaml['vcf_soft']

    ##########################################################################################
    #---Annotation variant calling
    ##########################################################################################
    #time cost
    time_start = time.time()
    #Annotation dir
    annotation_dir = out_dir
    raw_vcf = source

    #annotation
    if tools in ['all','annotation']:
        print("please check the variant_call subprocess result--VCF!")
        print("Test annotation module!\n")
    #Annotation dir
        if not os.path.exists(annotation_dir):
            os.makedirs(annotation_dir)
        annotationmain(db_cosmic, db_clinvar, db_g1000, 
                   ref_ens,
                   raw_vcf, sample,
                   annotation_dir,vcf_soft)
        print("--" * 20 + '\n\n')

if __name__ == '__main__':
    main()