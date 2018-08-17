# INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type: SNP or INDEL">
# INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
# INFO=<ID=MT,Number=1,Type=Integer,Description="Total MT depth">
# INFO=<ID=UMT,Number=1,Type=Integer,Description="Filtered MT depth">
# INFO=<ID=PI,Number=1,Type=Float,Description="Variant prediction index">
# INFO=<ID=THR,Number=1,Type=Integer,Description="Variant prediction index minimum threshold">
# INFO=<ID=VMT,Number=1,Type=Integer,Description="Variant MT depth">
# INFO=<ID=VMF,Number=1,Type=Float,Description="Variant MT fraction">
# INFO=<ID=VSM,Number=1,Type=Integer,Description="Variant strong MT depth">
# FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# FORMAT=<ID=AD,Number=.,Type=Integer,Description="Filtered allelic MT depths for the ref and alt alleles">
# FORMAT=<ID=VF,Number=1,Type=Float,Description="Variant MT fraction, same as VMF">
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
        genename, chr, start, end, geneid, alleleid, rs, pos, ref, alt, af_esp, af_exac, af_tgp, \
        clndn, clnhgvs, clnsig = db2.strip().split(',')
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
        print(alt)
        print(ad)
        alt1, alt2 = alt.split(',')
        af1, af2 = af.split(',')
        af1 = '0,'+ af1
        af2 = '0,'+ af2
        dp0, dp1, dp2 = ad.split(',')
        dp0 = float(dp0)
        dp1 = float(dp1)
        dp2 = float(dp2)
        vf1 = str(round(dp1/(dp0+dp1+dp2),4))
        vf2 = str(round(dp1/(dp0+dp1+dp2),4))
        return [[chrom, pos, ref, alt, filter, qual, dp, af1,vf1, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor],
                [chrom, pos, ref, alt, filter, qual, dp, af2,vf2, baseqranksum, fs, inbreedingcoeff, mq, mqranksum, qd, readposranksum, sor]
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


def annotation(dict_cos,dict_clin,dict_g1000,variant_vcf,annotated_csv,stats_file, snp_filter, indel_filter,sample_name):
    key_list = []
    key = ''
    change1 = ''
    num_in_clinvar = 0
    num_in_cosmic = 0
    num_in_g1000 = 0
    num_unmatch = 0
    var = open(variant_vcf, 'r')
    output = open(annotated_csv, 'w')
    output.write('\t'.join(['Sample','CHR','POS','REF','ALT','FILTER','QUAL','DP','AF','VF','Baseqranksum','FS','Inbreedingcoeff','MQ',
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
                if len(ref) == len(alt):
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
                    new = '\t'.join(value) + '\t' + dict_clin[key].replace(',' ,'\t') + '\t'
                    num_in_clinvar += 1
                else:
                    new = '\t'.join(value) + '\t'+ '-\t'*5
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
    print ('The sample has %s variants.' % len(key_list))
    print ('%s variants in COSMIC database' % num_in_cosmic)
    print ('%s variants in Clinvar database' % num_in_clinvar)
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
            print(subframe.iloc[num]['Gene_Name1'])
            print(n2g)
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

def main():
    time_start = time.time()
    (source, sample_name_file) = sys.argv[1:]
    cosmic= '/home/dell/Works/Projects/Datasets/annonDB/COSMIC_variant.csv'
    clinvar= '/home/dell/Works/Projects/Datasets/annonDB/breast_cancer_variant_clinvar.csv'
    g1000= '/home/dell/Works/Projects/Datasets/annonDB/breast_cancer_variant_1000genomes.csv'
    ref_ens= '/home/dell/Works/Projects/Datasets/annonDB/geneid_cancer_qiagen.csv'
    snp_filter='DP < 30 || QD < 2.0 || FS > 60.0 || MQ < 50.0 || SOR > 3.0 || MQRankSum < -2.5 || ReadPosRankSum < -3.0'
    indel_filter= 'DP < 30 || QD < 2.0 || FS > 200 || ReadPosRankSum < -3.0 || SOR > 10.0'
    #--
    dict_cos, dict_clin, dict_g1000 = read_database(cosmic,clinvar,g1000)
    #--
    snp_limit, indel_limit= read_vcf_filter(snp_filter, indel_filter)
    
    sample_name_list=open(sample_name_file,'r')
    for sample_name in sample_name_list:
        sample_name = sample_name.strip()
        variant_vcf = source + sample_name + '_L001.raw_variants.vcf'
        annotated_csv = source + sample_name + '.annotated.txt'
        annotated_csv_add = source + sample_name + '.annotated_add.txt'
        non_rs = source + sample_name + '.non_rs.csv'
        non_cos = source + sample_name + '.non_cos.csv'
        stats_file = source + sample_name + '.annotate_stats.txt'
        annotation(dict_cos,dict_clin,dict_g1000,variant_vcf,annotated_csv,stats_file, snp_limit, indel_limit,sample_name)
        #fill_table(annotated_csv,annotated_csv_add,non_rs,non_cos,ref_ens)
    print('The time of used annotation is %s minutes.' % str((time.time() - time_start) / 60))

if __name__ == '__main__':
    main()
