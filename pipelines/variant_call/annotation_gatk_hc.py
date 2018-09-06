#FILTER=<ID=LowQual,Description="Low quality">
#FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
#FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
#FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
#FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
#FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
#FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
#FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
#FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
#INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
#INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
#INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
#INFO= VF ,Description="Variant Frequency Allele bases/ (Allele bases + Ref Bases)">
#INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
#INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
#INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
#INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
#INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
#INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
#INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
#INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
#INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
#INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
#INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
#INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
#INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
#INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
#INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
#INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
# CLNDN--ClinVar`s preferred disease name for the concept specified by disease identifiers in CLNDISD.
# HGVS--Top-level(primary assembly,alt,or patch) HGVS expression.
# CLNSIG--Clinical significance for this single variant.
# Mutation_Description--Type of mutation at the amino acid level (substitution, deletion, insertion, complex, fusion, unknown etc.)
# Gene_CDS_Length--Length of the gene (base pair) counts.
# Mutation_Zygosity--Information on whether the mutation was reported to be homozygous , heterozygous or unknown within the sample.
# LOH--LOH Information on whether the gene was reported to have loss of heterozygosity in the sample: yes, no or unknown.
# Mutation_Strand--postive or negative.
# FATHMM_Prediction--Functional Analysis through Hidden Markov Models.
# FATHMM_Score--The scores are in the form of pvalues ranging from 0 to 1. Pvalues greater than 0.5 are pathogenic
# while less than 0.5 are benign. Pvalues close to 0 or 1 are the high confidence results which
# are more accurate. The results are annotated as 10 feature groups (separately for coding and
#  non coding variants) details of which can be found in the original FATHMM-MKL paper.
# Mutation_Somatic_Status--Information on whether the sample was reported to be Confirmed Somatic, Previously Reported or Variant of unknown origin.
from __future__ import barry_as_FLUFL

__all__  =  ['cosmic', 'clinvar', 'g1000' , 
             'ref_ens',
             'vcf', 
             'sample' , 'output' , 'snp_filter', 'indel_filter',
             'logger_annotation_process', 'logger_annotation_errors']
__version__  =  '1.0'
__author__  =  'Yang XueKai'

import os
import sys
import time
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
        genename, chr, start, end, geneid, alleleid, rs, pos, ref, alt, af_esp, af_exac, af_tgp, clndn, clnhgvs, clnsig = db2.strip().split(',')
        if rs is 'N':
            rs = '-'
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


#def get_info(tag):
#    return tag[tag.find('=')+1:]
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

#split mutation. such as, ref:CTTT  alt:C,CT
def split_variant(line):
    chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
    ac, af, an, baseqranksum, clippingranksum, dp, ds, end, excesshet, fs, inbreedingcoeff, mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor = get_info(info)
    num_mt = len(alt.split(','))
    if len(detail.split(':'))>5:
        gt, ad, dp= detail.split(':')[0:3]
        #pl = detail.split(':')[len(detail.split(':'))-1]
        #gq = ':'.join(detail.split(':')[3:len(detail.split(':'))])
    else:
        gt, ad, dp, gq, pl = detail.split(':')
    if num_mt is 1:
        dp0, dp1= ad.split(',')
        dp0 = float(dp0)
        dp1 = float(dp1)
        af=gt.replace('/',',')
        vf = str(round(dp1/(dp0+dp1),4))
        return [[chrom, pos, ref, alt, filter, qual, dp, af,vf, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor]]
    elif num_mt is 2:
        alt1, alt2 = alt.split(',')
        af1, af2 = af.split(',')
        af1 = '0,'+ af1
        af2 = '0,'+ af2
        dp0, dp1, dp2 = ad.split(',')
        dp0 = float(dp0)
        dp1 = float(dp1)
        dp2 = float(dp2)
        vf1 = str(round(dp1/(dp0+dp1+dp2),4))
        vf2 = str(round(dp2/(dp0+dp1+dp2),4))
        return [[chrom, pos, ref, alt1, filter, qual, dp, af1,vf1, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor],
                [chrom, pos, ref, alt2, filter, qual, dp, af2,vf2, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor]
                ]

def comp_filter(limists,value):
    result=[]
    count = 0
    for para in limists.keys():
        if para == 'DP':
            if value[6] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[6]):
                    count+=1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[6]):
                    count+=1
                else:
                    result.append('Low'+para)
        elif para == 'QUAL':
            if value[5] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[5]):
                    count+=1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[5]):
                    count+=1
                else:
                    result.append('Low'+para)
        elif para == 'QD':
            if value[14] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[14]):
                    count+=1
                else:
                    result.append('High'+para)
            elif limists[para][0] =='<':
                if float(limists[para][1]) <= float(value[14]):
                    count+=1
                else:
                    result.append('Low'+para)
        elif para == 'FS':
            if value[10] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[10]):
                    count+=1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[10]):
                    count+=1
                else:
                    result.append('Low'+para)
        elif para == 'MQ':
            if value[12] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[12]):
                    count+=1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[12]):
                    count+=1
                else:
                    result.append('Low'+para)
        elif para == 'SOR':
            if value[16] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[16]):
                    count+=1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[16]):
                    count+=1
                else:
                    result.append('Low'+para)
        elif para == 'MQRankSum':
            if value[13] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[13]):
                    count+=1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[13]):
                    count+=1
                else:
                    result.append('Low'+para)
        elif para == 'ReadPosRankSum':
            if value[15] == '-':
                result.append('No'+para)
            elif limists[para][0] =='>':
                if float(limists[para][1]) >= float(value[15]):
                    count+=1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[15]):
                    count+=1
                else:
                    result.append('Low'+para)
    if count == len(limists.keys()):
        return 'PASS'
    else:
        return '|'.join(result)

def read_vcf_filter(snp_filter, indel_filter):
    snp = snp_filter.split(' || ')
    snp_limit={}
    for limit in snp:
        para, char, value = limit.split(" ")
        snp_limit[para] = [char, value]

    indel = indel_filter.split(' || ')
    indel_limit={}
    for limit in indel:
        para, char, value = limit.split(" ")
        indel_limit[para] = [char, value]
    return snp_limit, indel_limit

def define_hgvs(chr, pos, ref, alt):
    chr_to_version = {'1': '10', '2': '11', '3': '11', '4': '11', '5': '9', '6': '11', '7': '13', '8': '10', '9': '11',
                      '10': '10', '11': '9', '12': '11', '13': '10', '14': '8', '15': '9', '16': '9', '17': '10',
                      '18': '9', '19': '9', '20': '10', '21': '8', '22': '10', 'X': '10'}
    version = chr_to_version[chr]
    if chr != 'X' and int(chr) < 10:
        chr = '0' + chr
    if len(ref) == 1 and len(alt) == 1:
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + pos + ref + '>' + alt
    if len(ref) > 1 and len(alt) == 1:
        deletion = ref[1:]
        if len(deletion) == 1:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+1) + 'del' + deletion
        else:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+1) + '_' + str(int(pos)+len(deletion)) + 'del' + deletion
    if len(ref) == 1 and len(alt) > 1:
        insertion = alt[1:]
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)) + '_' + str(int(pos)+1) + 'ins' + insertion
    if len(ref) > 1 and len(alt) > 1 and len(ref)>len(alt):
        deletion = ref[len(alt):]
        if len(deletion) == 1:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+len(alt)) + 'del' + deletion
        else:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+len(alt)) + '_' + str(int(pos)+len(alt)+len(deletion)-1) + 'del' + deletion
    return hgvs

def annotation(dict_cos,dict_clin,dict_g1000,variant_vcf,annotated_csv,stats_file, snp_filter, indel_filter,sample_name,logger_annotation_process):
    key_list = []
    key = ''
    change1 = ''
    num_in_clinvar = 0
    num_in_cosmic = 0
    num_in_g1000 = 0
    num_unmatch = 0
    var = open(variant_vcf, 'r')
    output = open(annotated_csv, 'w')
    output.write('\t'.join(['Sample','CHR','POS','REF','ALT','FILTER','QUAL','DP','GT','VF','Baseqranksum','FS','Inbreedingcoeff','MQ',
                            'MQRankSum','QD','ReadPosRankSum','SOR','Gene_ID','RS_ID','CLNDN','HGVS','CLNSIG',
        'COSMIC_ID','Mutation_Description','Feature_ID','Gene_Name','Gene_CDS_Length','Mutation_Zygosity','LOH','Mutation_Strand',
        'HGVS.c','HGVS.p','FATHMM_Prediction','FATHMM_Score','Mutation_Somatic_Status','Gene_Name1','RS_ID1','EAS_AF','EUR_AF','AMR_AF',
        'SAS_AF','AFR_AF\n']))
    for line in var:
        if not line.startswith('#'):
            for spl in split_variant(line):
                chrom, pos, ref, alt, filter, qual, dp, af,vf, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor = spl
                chrom = chrom[3:]
                value = [chrom, pos, ref, alt, filter, qual, dp, af, vf, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor]
                #print(value)
                if len(ref) == len(alt) and len(alt) == 1:
                    change = ref + '>' + alt
                    change1 = base_paired[ref] + '>' + base_paired[alt]
                    filter = comp_filter(snp_filter,value)
                elif len(ref) > len(alt) and len(alt) == 1:
                    change = 'del' + ref[1:]
                    filter = comp_filter(indel_filter,value)
                elif len(ref) < len(alt) and len(ref) == 1:
                    change = 'ins' + alt[1:]
                    filter = comp_filter(indel_filter,value)
                else:
                    change = 'del' + ref + 'ins' + alt
                    filter = comp_filter(indel_filter,value)
                key = chrom + ',' + pos + ',' + change
                key1 = chrom + ',' + pos + ',' + change1
                value = [sample_name,chrom, pos, ref, alt, filter, qual, dp, af, vf, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor]
                #print(value)
                unmatch = 0
                # drop duplicate variant
                if key in key_list:
                    continue
                if key in dict_clin:
                    print(dict_clin[key])
                    new = '\t'.join(value) + '\t' + dict_clin[key].replace(',' ,'\t') + '\t'
                    num_in_clinvar += 1
                else:
                    new = '\t'.join(value) + '\t'+ '-\t'*3+ define_hgvs(chrom, pos, ref, alt) + '\t-\t'
                    unmatch += 1
                if key in dict_cos:
                    new += dict_cos[key].replace(',' ,'\t') + '\t'
                    num_in_cosmic += 1
                elif key1 in dict_cos:
                    new += dict_cos[key1].replace(',' ,'\t') + '\t'
                    num_in_cosmic += 1
                else:
                    new += '-\t'*13
                    unmatch += 1
                if key in dict_g1000:
                    new += dict_g1000[key].replace(',' ,'\t') + '\n'
                    num_in_g1000 += 1
                else:
                    new += '-\t'*6+ '-\n'
                # 3 databases both unmatch
                if unmatch == 3:
                    num_unmatch += 1
                else:
                    output.write(new)
                key_list.append(key)
    output.close()
    logger_annotation_process.info('{0} has {1} variants.'.format(sample_name,key_list))
    logger_annotation_process.info('{0} has {1} variants in COSMIC database'.format(sample_name,num_in_cosmic))
    logger_annotation_process.info('{0} has {1} variants in Clinvar database'.format(sample_name,num_in_clinvar))
    logger_annotation_process.info('{0} has {1} variants in G1000 database'.format(sample_name, num_in_g1000))
    logger_annotation_process.info('{0} has {1} variants unmatch in cosmic and clinvar.'.format(sample_name,num_unmatch))
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
def fill_table(annotated_csv, annotated_csv_add, ref_ens):
    n2g = {}
    g2n = {}
    for line in open(ref_ens, 'r').readlines():
        name, ensg = line.strip().split(',')
        n2g[name] = ensg
        g2n[ensg] = name
    #df = pd.read_csv(annotated_csv)
    df = pd.read_table(annotated_csv)
    subframe = df[['Gene_Name','Gene_ID','Gene_Name1','RS_ID','RS_ID1']] #n,g,t,n1
    #for name,id,transcript in subframe.iterrows():
    for num in range(0, len(subframe)):
        if subframe.iloc[num, 0] is '-' and subframe.iloc[num, 1] is '-' and subframe.iloc[num, 2] is not '-':
            name = subframe.iloc[num, 2]
            ensg = n2g[name]
            subframe.iloc[num, 0] = name
            subframe.iloc[num, 1] = ensg
        elif subframe.iloc[num, 0] is '-' and subframe.iloc[num, 2] is '-' and subframe.iloc[num, 1] is not '-':
            ensg = subframe.iloc[num, 1]
            name = g2n[ensg]
            subframe.iloc[num, 0] = name
        elif subframe.iloc[num, 0] is not '-' and subframe.iloc[num, 1] is '-' and subframe.iloc[num, 2] is '-':
            name = subframe.iloc[num, 0]
            ensg = n2g[name]
            subframe.iloc[num, 1] = ensg
        elif subframe.iloc[num, 0] is '-':
            name = subframe.iloc[num, 2]
            subframe.iloc[num, 0] = name
        elif subframe.iloc[num, 1] is '-':
            name = subframe.iloc[num, 0]
            ensg = n2g[name]
            subframe.iloc[num, 1] = ensg
        if subframe.iloc[num, 3] is '-' and subframe.iloc[num, 4] is not '-':
            subframe.iloc[num, 3] = subframe.iloc[num, 4]
    name = subframe['Gene_Name']
    ensg = subframe['Gene_ID']
    rs = subframe['RS_ID']
    df.drop(labels=['Gene_Name'], axis=1, inplace=True)
    df.drop(labels=['Gene_ID'], axis=1, inplace=True)
    df.drop(labels=['Gene_Name1'], axis=1, inplace=True)
    df.drop(labels=['RS_ID'], axis=1, inplace=True)
    df.drop(labels=['RS_ID1'], axis=1, inplace=True)
    df.insert(10, 'Gene_Name', name)
    df.insert(11, 'Gene_ID', ensg)
    df.insert(10, 'RS_ID', rs)
    df.to_csv(annotated_csv, index=False, sep='\t')
    #df[(True^df['RS_ID'].isin(['-']))].to_csv(non_rs, index=False, sep='\t')
    #df[(True^df['COSMIC_ID'].isin(['-']))].to_csv(non_cos, index=False, sep='\t')

#-annotation main
def annotationmain(cosmic, clinvar, g1000, 
                   ref_ens,
                   vcf, sample,snp_filter, indel_filter,
                   output, logger_annotation_process, logger_annotation_errors):
    if 'raw_variants_SNP.vcf' in vcf:
        annotated_csv = output + '/' + sample + '.raw_SNP.annotated.txt'
        annotated_csv_add = output + '/' + sample + '.raw_SNP.annotated_ensembl.txt'
        non_rs = output + '/' + sample + '.raw_SNP_non_rs.txt'
        non_cos = output + '/' + sample + '.raw_SNP_non_cos.txt'
        stats_file = output + '/' + sample + '.raw_SNP_annotate_stats.txt'
    elif 'raw_variants_indel.vcf' in vcf:
        annotated_csv = output + '/' + sample + '.raw_indel.annotated.txt'
        annotated_csv_add = output + '/' + sample + '.raw_indel.annotated_ensembl.txt'
        non_rs = output + '/' + sample + '.raw_indel_non_rs.txt'
        non_cos = output + '/' + sample + '.raw_indel_non_cos.txt'
        stats_file = output + '/' + sample + '.raw_indel_annotate_stats.txt'
    elif 'filter_SNP.vcf' in vcf:
        annotated_csv = output + '/' + sample + '.filter_SNP.annotated.txt'
        annotated_csv_add = output + '/' + sample + '.filter_SNP.annotated_ensembl.txt'
        non_rs = output + '/' + sample + '.filter_SNP_non_rs.txt'
        non_cos = output + '/' + sample + '.filter_SNP_non_cos.txt'
        stats_file = output + '/' + sample + '.filter_SNP_annotate_stats.txt'
    elif 'filter_indel.vcf' in vcf:
        annotated_csv = output + '/' + sample + '.filter_indel.annotated.txt'
        annotated_csv_add = output + '/' + sample + '.filter_indel.annotated_ensembl.txt'
        non_rs = output + '/' + sample + '.filter_indel_non_rs.txt'
        non_cos = output + '/' + sample + '.filter_indel_non_cos.txt'
        stats_file = output + '/' + sample + '.filter_indel_annotate_stats.txt'
    elif 'raw_variants.vcf' in vcf:
        annotated_csv = output + '/' + sample + '.raw_variants.annotated.txt'
        annotated_csv_add = output + '/' + sample + '.raw_variants.annotated_ensembl.txt'
        non_rs = output + '/' + sample + '.raw_variants_non_rs.txt'
        non_cos = output + '/' + sample + '.raw_variants_non_cos.txt'
        stats_file = output + '/' + sample + '.raw_variants_annotate_stats.txt'
    #-read the annotation database
    if not os.path.isfile(vcf):
        logger_annotation_errors.error("%s does not exist!\n", vcf)
        print(vcf + ' does not exist!')
    else:
        dict_cos, dict_clin, dict_g1000 = read_database(cosmic,clinvar,g1000)
        #--annotation
        annotation(dict_cos, dict_clin, dict_g1000,vcf,annotated_csv,stats_file, snp_filter, indel_filter,sample,logger_annotation_process)
        #--add the annotation
        fill_table(annotated_csv, annotated_csv_add, ref_ens)

