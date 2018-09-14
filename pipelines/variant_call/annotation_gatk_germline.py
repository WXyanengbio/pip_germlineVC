import os
import sys
import time
import pandas as pd
import numpy as np

# put the info output to the log
def stdout_err(command):
    command_pope = shlex.split(command)
    child = subprocess.Popen(command_pope, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = child.communicate()
    child.wait()
    return stdout, stderr


base_paired = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def read_database(cosmic, clinvar, g1000):
    cos = open(cosmic, 'r')
    dict_cos = {}
    cos.readline()
    for db1 in cos.readlines():
        # if not db.startswith('gene_name'):
        (gene_name, accession_number, gene_cds_length, hgnc_id, sample_name, id_sample, id_tumour, primary_site,
         site_subtype1, site_subtype2, site_subtype3, primary_histology, histology_subtype1, histology_subtype2,
         histology_subtype3, genome_wide_screen, mutation_id, mutation_cds, mutation_aa, mutation_description,
         mutation_zygosity, loh, grch, chr, start, end, mutation_strand, snp, resistance_mutation, fathmm_prediction,
         fathmm_score, mutation_somatic_status, pubmed_pmid, id_study,
         sample_source, tumor_origin, age) = db1.strip().split(',')
        if mutation_cds is 'ns':
            continue
        if 'del' in mutation_cds:
            change = mutation_cds[mutation_cds.find('del'):]
        elif 'ins' in mutation_cds:
            change = mutation_cds[mutation_cds.find('ins'):]
        elif '>' in mutation_cds:
            change = mutation_cds[mutation_cds.find('>')-1:]
        key1 = [chr, start, change]
        value1 = [mutation_id, mutation_description, accession_number, gene_name, gene_cds_length,
                  mutation_zygosity, loh, mutation_strand, mutation_cds, mutation_aa, fathmm_prediction, fathmm_score,
                  mutation_somatic_status]
        dict_cos[','.join(key1)] = ','.join(value1)
    clin = open(clinvar, 'r')
    dict_clin = {}
    clin.readline()
    for db2 in clin.readlines():
        (genename, chr, start, end, geneid, pos, ref, alt, alleleid,
         rs, af_esp, af_exac, af_tgp, clndn, clnhgvs, clnsig) = db2.strip().split(',')
        if rs is 'N':
            rs = '-'
        if len(ref) == len(alt):
            change = ref + '>' + alt
        elif (len(ref) > len(alt)) and len(alt) == 1:
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
        (genename1, chr1, start1, end1, variant_type, ref1, alt1,
         rs1, eas_af, eur_af, amr_af, sas_af, afr_af, hgvs1) = db3.strip().split(',')
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


# define HGVS by chrom, pos, ref, alt
def define_hgvs(chr, pos, ref, alt):
    chr_to_version = {'1': '10', '2': '11', '3': '11', '4': '11', '5': '9', '6': '11', '7': '13', '8': '10', '9': '11',
                      '10': '10', '11': '9', '12': '11', '13': '10', '14': '8', '15': '9', '16': '9', '17': '10',
                      '18': '9', '19': '9', '20': '10', '21': '8', '22': '10', 'X': '10'}
    version = chr_to_version[chr]
    if chr != 'X' and int(chr) < 10:
        chr = '0' + chr
    if len(ref) == 1 and len(alt) == 1:
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + pos + ref + '>' + alt
    elif len(ref) > 1 and len(alt) == 1:
        deletion = ref[1:]
        if len(deletion) == 1:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+1) + 'del' + deletion
        else:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+1) + '_' + str(int(pos)+len(deletion)) +\
                   'del' + deletion
    elif len(ref) == 1 and len(alt) > 1:
        insertion = alt[1:]
        hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)) + '_' + str(int(pos)+1) + 'ins' + insertion
    elif (len(ref) > 1) and (len(alt) > 1) and (len(ref) > len(alt)):
        deletion = ref[len(alt):]
        if len(deletion) == 1:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+len(alt)) + 'del' + deletion
        else:
            hgvs = 'NC_0000' + chr + '.' + version + ':g.' + str(int(pos)+len(alt)) + '_' +\
                   str(int(pos)+len(alt)+len(deletion)-1) + 'del' + deletion
    else:
        hgvs = '-'
    return hgvs


# get info of VCF from gatk, samtools
def get_info(info, soft):
    vcf_parameters = {'GATK': ['AC', 'AF', 'AN', 'BaseQRankSum', 'ClippingRankSum', 'DP', 'DS', 'END', 'ExcessHet',
                               'FS', 'InbreedingCoeff', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'RAW_MQ',
                               'ReadPosRankSum', 'SOR'],
                      'samtools': ['INDEL', 'IDV', 'IMF', 'DP', 'VDB', 'RPB', 'MQB', 'BQB', 'MQSB', 'SGB',
                                   'MQ0F', 'ICB', 'HOB', 'AC', 'AN', 'DP4', 'MQ']}
    subparas = info.split(';')
    subparas_pre = list(map(lambda tag: tag[:tag.find('=')], subparas))
    paras = []
    for para in vcf_parameters[soft]:
        if para in subparas_pre:
            tag = subparas[subparas_pre.index(para)]
            paras.append(tag[tag.find('=')+1:])
        else:
            paras.append('-')
    return paras


# filter vcf by DP for strelka, samtools, varscan2
def comp_filter_varscan(limists, value):
    result = []
    count = 0
    for para in limists.keys():
        if para == 'DP':
            if value[5] == '-':
                result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[5]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[5]):
                    count += 1
                else:
                    result.append('Low'+para)
    if count == 1:
        return 'PASS'
    else:
        return '|'.join(result)


# filter vcf by DP, QUAL, MQ for strelka, samtools,
def comp_filter_strelka(limists, value):
    result = []
    count = 0
    for para in limists.keys():
        if para == 'DP':
            if value[6] == '-':
                result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[6]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[6]):
                    count += 1
                else:
                    result.append('Low'+para)
        elif para == 'QUAL':
            if value[5] == '-':
                result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[5]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[5]):
                    count += 1
                else:
                    result.append('Low'+para)
        elif para == 'MQ':
            if value[10] == '-':
                result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[10]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[10]):
                    count += 1
                else:
                    result.append('Low'+para)
    if count == 3:
        return 'PASS'
    elif count == 2 and ('QUAL' not in limists.keys()):
        return 'PASS'
    elif count == 1 and ('QUAL' not in limists.keys()) and ('MQ' not in limists.keys()):
        return 'PASS'
    else:
        return '|'.join(result)


# read vcf filter expression
def read_vcf_filter(snp_filter, indel_filter):
    snp = snp_filter.split(' || ')
    snp_limit = {}
    for limit in snp:
        para, char, value = limit.split(" ")
        snp_limit[para] = [char, value]

    indel = indel_filter.split(' || ')
    indel_limit = {}
    for limit in indel:
        para, char, value = limit.split(" ")
        indel_limit[para] = [char, value]
    return snp_limit, indel_limit


# get the values of strelka2 vcf parameters form format and detail
def get_detail_from_strelka_vcf(format, detail):
    detail_dict = {}
    formats = format.split(':')
    details = detail.split(':')
    strelka2 = ['GT', 'GQ', 'GQX', 'DP', 'DPF', 'MIN_DP', 'AD', 'ADF', 'ADR', 'FT', 'DPI', 'PL', 'PS', 'SB']
    for i in range(0, len(strelka2)):
        detail_dict[strelka2[i]] = '-'
        for j in range(0, len(formats)):
            if formats[j] == strelka2[i]:
                detail_dict[strelka2[i]] = details[j]
    return list(detail_dict.values())


# split mutation. such as, ref:CTTT  alt:C,CT
def split_variant_strelka(line):
    chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
    mq = info.split(';')[len(info.split(';'))-1].split('=')[1]
    gt, gq, gqx, dp, dpf, min_dp, ad, adf, adr, ft, dpi, pl, ps, sb = get_detail_from_strelka_vcf(format, detail)
    if dpi == '-':
        var_type = "SNP"
    else:
        var_type = "INDEL"
    return [[chrom, pos, ref, alt, var_type, filter, qual, dp, gt, dpf, min_dp, ad, adf, adr, ft, dpi, mq, sb,  pl, ps, gq, gqx]]


# split mutation. such as, ref:CTTT  alt:C,CT
def split_variant(line, calling):
    if calling == 'GATK':
        chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
        (ac, af, an, baseqranksum, clippingranksum, dp, ds, end, excesshet, fs,
         inbreedingcoeff, mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor) = get_info(info, calling)
        num_mt = len(alt.split(','))
        if len(detail.split(':')) > 5:
            gt, ad, dp = detail.split(':')[0:3]
            pl = detail.split(':')[len(detail.split(':'))-1]
            gq = ':'.join(detail.split(':')[3:len(detail.split(':'))])
        else:
            gt, ad, dp, gq, pl = detail.split(':')
        if len(ref) ==1 and len(alt) ==1:
            var_type = "SNP"
        elif len(ref) ==1 and num_mt >1:
            if len(alt.split(',')[0]) ==1 and len(alt.split(',')[1]) ==1:
                var_type = "SNP"
            else:
                var_type = "INDEL"
        else:
            var_type = "INDEL"
        return [[chrom, pos, ref, alt, var_type, filter, qual, dp, af, ad, ac, an, baseqranksum, clippingranksum, ds,
                 end, excesshet, fs, inbreedingcoeff, mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum,
                 sor, gt, gq, pl]]

    elif calling == 'samtools':
        chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
        indel, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb, sgb, mq0f, icb, hob, ac, an, dp4, mq = get_info(info, calling)
        gt, pl = detail.split(':')
        ref_f, ref_r, alt_f, alt_r = dp4.split(',')
        vf = str(round((float(alt_f)+float(alt_r))/(float(alt_f)+float(alt_r)+float(ref_f)+float(ref_r)), 4))
        if idv == '-':
            var_type = "SNP"
        else:
            var_type = "INDEL"
        return[[chrom, pos, ref, alt, qual, var_type, filter, vf, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb, sgb,
                mq0f, icb, hob, ac, an, dp4, mq, ref_f, ref_r, alt_f, alt_r, gt, pl]]

    elif calling == 'varscan2':
        # get the values of samtools vcf parameters form format and detail
        chrom, pos, ref, alt, poolcall, strandfilt, sref, shet, shom, snc, scalls = line.strip().split('\t')
        cons, cov, read1, read2, freq, pvalue = poolcall.split(':')
        vf = str(float(freq.strip('%'))/100)
        sf = strandfilt.split(':')
        if len(sf) == 6:
            strandfilter, r1f, r1r, r2f, r2r, pval = sf
        else:
            strandfilter, other, r1f, r1r, r2f, r2r, pval = sf
        if int(shet) == 1:
            gt = '0,1'
        elif int(shom) == 1:
            gt = '1,1'
        alt = alt.lstrip('=')
        if alt[0] == '-':
            ref1 = ref + alt[1:]
            alt = ref
            var_type = "INDEL"
        elif alt[0] == '+':
            ref1 = ref
            alt = ref + alt[1:]
            var_type = "INDEL"
        else:
            ref1 = ref
            var_type = "SNP"
        return[[chrom, pos, ref1, alt, var_type, gt, cov, read1, read2, vf, pvalue, strandfilter,
                r1f, r1r, r2f, r2r, pval]]


def read_vcf(variant_vcf, sample_name, calling, snp_filter, indel_filter):
    key = ''
    change1 = ''
    if calling == 'varscan2':
        variant_vcf_snp = os.path.dirname(variant_vcf) + '/' + sample_name + '.raw_varscan2.snp.vcf'
        variant_vcf_indel = os.path.dirname(variant_vcf) + '/' + sample_name + '.raw_varscan2.indel.vcf'
        snp = open(variant_vcf_snp, 'r')
        indel = open(variant_vcf_indel, 'r')

        var = []
        for line in snp.readlines():
            var.append(line.strip())
        for line in indel.readlines():
            var.append(line.strip())

        # output = open(annotated_csv, 'w')
        # output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'DP', 'VF', 'GT', 'DP_read1',
        #                         'DP_read2', 'Fisherexact_var', 'read1F', 'read1R', 'read2F', 'read2R', 'StrandBias',
        #                         'Gene_ID', 'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID', 'Mutation_Description',
        #                         'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
        #                         'Mutation_Strand', 'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score',
        #                         'Mutation_Somatic_Status',  'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF',
        #                         'AMR_AF', 'SAS_AF', 'AFR_AF\n']))
        output = {}
        # print(variant_vcf)
        for line in var:
            if not line.startswith('Chrom'):
                for spl in split_variant(line, calling):
                    (chrom, pos, ref, alt, var_type, gt, dp, read1, read2, vf, pvalue, filter,
                     r1f, r1r, r2f, r2r, pval) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, dp, vf]
                    # upper the ref and alt
                    ref = ref.upper()
                    alt = alt.upper()
                    if len(ref) == len(alt) and len(alt) == 1:
                        change = ref + '>' + alt
                        change1 = base_paired[ref] + '>' + base_paired[alt]
                        filter = comp_filter_varscan(snp_filter, value)
                    elif (len(ref) > len(alt)) and (len(alt) == 1):
                        change = 'del' + ref[1:]
                        filter = comp_filter_varscan(indel_filter, value)
                    elif len(ref) < len(alt) and len(ref) == 1:
                        change = 'ins' + alt[1:]
                        filter = comp_filter_varscan(indel_filter, value)
                    else:
                        change = 'del' + ref + 'ins' + alt
                        filter = comp_filter_varscan(indel_filter, value)
                    key = chrom + ',' + pos + ',' + change
                    key1 = chrom + ',' + pos + ',' + change1
                    key2 = 'chr' + chrom + '_' + pos
                    value = [sample_name, chrom, pos, ref, alt, var_type, filter, dp, vf, gt, read1, read2, pvalue,
                             r1f, r1r, r2f, r2r, pval]
                    yield [key2, key, key1, value]
    elif calling == 'strelka2':
        variant_vcf = os.path.dirname(variant_vcf) + '/results/variants/variants.vcf.gz'
        if not os.path.isfile(os.path.dirname(variant_vcf) + '/variants.vcf'):
            command1 = 'gunzip' + ' -k ' + variant_vcf
            stdout, stderr = stdout_err(command1)
        # output = open(annotated_csv, 'w')
        # output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'QUAL', 'DP', 'GT', 'AD', 'VF', 'MQ',
        #                         'SB', 'Gene_ID', 'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID',
        #                         'Mutation_Description',
        #                         'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
        #                         'Mutation_Strand',
        #                         'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score', 'Mutation_Somatic_Status',
        #                         'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF', 'AMR_AF', 'SAS_AF', 'AFR_AF\n']))
        var = open(os.path.dirname(variant_vcf) + '/variants.vcf', 'r')
        for line in var:
            if not line.startswith('#'):
                for spl in split_variant_strelka(line):
                    (chrom, pos, ref, alt, var_type, filter, qual, dp, gt, dpf, min_dp,
                     ad, adf, adr, ft, dpi, mq, sb, pl, ps, gq, gqx) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, qual, dp, gt, ad, dpf, mq, sb]
                    if len(ref) == len(alt) and len(alt) == 1:
                        change = ref + '>' + alt
                        change1 = base_paired[ref] + '>' + base_paired[alt]
                        if filter == 'PASS':
                            filter = comp_filter_strelka(snp_filter, value)
                    elif (len(ref) > len(alt)) and (len(alt) == 1):
                        change = 'del' + ref[1:]
                        if filter == 'PASS':
                            filter = comp_filter_strelka(indel_filter, value)
                    elif len(ref) < len(alt) and len(ref) == 1:
                        change = 'ins' + alt[1:]
                        if filter == 'PASS':
                            filter = comp_filter_strelka(indel_filter, value)
                    else:
                        change = 'del' + ref + 'ins' + alt
                        if filter == 'PASS':
                            filter = comp_filter_strelka(indel_filter, value)
                    key = chrom + ',' + pos + ',' + change
                    key1 = chrom + ',' + pos + ',' + change1
                    key2 = 'chr' + chrom + '_' + pos
                    value = [sample_name, chrom, pos, ref, alt, var_type, filter, qual, dp, gt, dpf, min_dp,
                             ad, adf, adr, ft, dpi, mq, sb, pl, ps, gq, gqx]
                    yield [key2, key, key1, value]
    elif calling == 'GATK':
        var = open(variant_vcf, 'r')
        # output = open(annotated_csv, 'w')
        # output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'QUAL', 'DP', 'GT1', 'AD', 'VF', 'AC',
        #                         'AN', 'Baseqranksum', 'ClippingRankSum', 'DS', 'END', 'ExcessHet', 'FS',
        #                         'InbreedingCoeff', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'RAW_MQ',
        #                         'ReadPosRankSum', 'SOR', 'GT', 'GQ', 'PL',
        #                         'Gene_ID', 'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID',
        #                         'Mutation_Description',
        #                         'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
        #                         'Mutation_Strand', 'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score',
        #                         'Mutation_Somatic_Status', 'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF', 'AMR_AF',
        #                         'SAS_AF', 'AFR_AF\n']))
        for line in var:
            if not line.startswith('#'):
                for spl in split_variant(line, calling):
                    (chrom, pos, ref, alt, var_type, filter, qual, dp, af, ad, ac, an, baseqranksum, clippingranksum,
                     ds, end, excesshet, fs, inbreedingcoeff,
                     mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor, gt, gq, pl) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, qual, dp, ad, pl, baseqranksum, fs, inbreedingcoeff,
                             mq, mqranksum, qd, readposranksum, sor]
                    # print(value)
                    if len(ref) == len(alt) and len(alt) == 1:
                        change = ref + '>' + alt
                        change1 = base_paired[ref] + '>' + base_paired[alt]
                        filter = comp_filter(snp_filter, value)
                    elif len(ref) > len(alt) and len(alt) == 1:
                        change = 'del' + ref[1:]
                        filter = comp_filter(indel_filter, value)
                    elif len(ref) < len(alt) and len(ref) == 1:
                        change = 'ins' + alt[1:]
                        filter = comp_filter(indel_filter, value)
                    else:
                        change = 'del' + ref + 'ins' + alt
                        filter = comp_filter(indel_filter, value)
                    key = chrom + ',' + pos + ',' + change
                    key1 = chrom + ',' + pos + ',' + change1
                    key2 = 'chr' + chrom + '_' + pos
                    value = [sample_name, chrom, pos, ref, alt, var_type, filter,
                             qual, dp, af, ad, ac, an, baseqranksum,
                             clippingranksum, ds, end, excesshet, fs, inbreedingcoeff,
                             mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor, gt, gq, pl]
                    yield [key2, key, key1, value]
    elif calling == 'samtools':
        variant_vcf = os.path.dirname(variant_vcf) + '/' + sample_name + '.raw_samtools.vcf'
        var = open(variant_vcf, 'r')
        # output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'QUAL', 'DP', 'VF', 'INDEL', 'IDV',
        #                         'IMF', 'VDB', 'RPB', 'MQB', 'BQB', 'MQSB', 'SGB', 'MQ0F', 'ICB', 'HOB',
        #                         'AC', 'AN', 'DP4', 'MQ', 'ref_f', 'ref_r', 'alt_f', 'alt_r', 'GT', 'PL', 'Gene_ID',
        #                         'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID', 'Mutation_Description',
        #                         'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
        #                         'Mutation_Strand',
        #                         'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score', 'Mutation_Somatic_Status',
        #                         'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF', 'AMR_AF', 'SAS_AF', 'AFR_AF\n']))
        # print(variant_vcf)
        for line in var:
            if not line.startswith('#'):
                for spl in split_variant(line, calling):
                    (chrom, pos, ref, alt, qual, var_type, filter, vf, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb,
                     sgb, mq0f, icb, hob, ac, an, dp4, mq, ref_f, ref_r, alt_f, alt_r, gt, pl) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, qual, dp, vf, imf, idv, mq, vdb]
                    # print(value)
                    # upper the ref and alt
                    ref = ref.upper()
                    alt = alt.upper()
                    if len(ref) == len(alt) and len(alt) == 1:
                        change = ref + '>' + alt
                        change1 = base_paired[ref] + '>' + base_paired[alt]
                        filter = comp_filter_strelka(snp_filter, value)
                    elif len(ref) > len(alt) and len(alt) == 1:
                        change = 'del' + ref[1:]
                        filter = comp_filter_strelka(indel_filter, value)
                    elif len(ref) < len(alt) and len(ref) == 1:
                        change = 'ins' + alt[1:]
                        filter = comp_filter_strelka(indel_filter, value)
                    else:
                        change = 'del' + ref + 'ins' + alt
                        filter = comp_filter_strelka(indel_filter, value)

                    key = chrom + ',' + pos + ',' + change
                    key1 = chrom + ',' + pos + ',' + change1
                    key2 = 'chr' + chrom + '_' + pos
                    value = [sample_name, chrom, pos, ref, alt, var_type, filter, qual, dp, vf, idv, imf,
                             vdb, rpb,
                             mqb, bqb, mqsb, sgb, mq0f, icb, hob,
                             ac, an, dp4, mq, ref_f, ref_r, alt_f, alt_r, gt, pl]
                    yield [key2, key, key1, value]


def comp_filter(limists, value):
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


def annotation(variant_vcf, output, sample_name, callings, snp_filter, indel_filter, other_limit):
    varscan = {}
    gatk = {}
    strelka = {}
    samtools ={}
    if "varscan2" in callings:
        calling = "varscan2"
        for keys in zip(read_vcf(variant_vcf, sample_name, calling, other_limit, other_limit)):
            varscan[keys[0][0]] = keys[0][3]
    if "strelka2" in callings:
        calling = "strelka2"
        for keys in zip(read_vcf(variant_vcf, sample_name, calling, other_limit, other_limit)):
            strelka[keys[0][0]] = keys[0][3]
    if "samtools" in callings:
        calling = "samtools"
        for keys in zip(read_vcf(variant_vcf, sample_name, calling, other_limit, other_limit)):
            samtools[keys[0][0]] = keys[0][3]
    if "GATK" in callings:
        calling = "GATK"
        for keys in zip(read_vcf(variant_vcf, sample_name, calling, snp_filter, indel_filter)):
            gatk[keys[0][0]] = keys[0][3]

    muts = np.unique(list(varscan.keys()) + list(gatk.keys()) + list(samtools.keys()) + list(samtools.keys()))
    annotated_csv_add = output.replace(".txt", "_PASS.txt")
    vcfsout = open(output, 'w')
    vcfsout1 = open(annotated_csv_add, 'w')
    vcfsout.write('\t'.join(['Mut', 'ga_CHROM', 'ga_POS', 'ga_REF', 'ga_ALT', 'ga_TYPE',
                              'ga_FILTER', 'ga_QUAL', 'ga_DP', 'ga_AF',
                              'ga_AD', 'ga_AC', 'ga_AN', 'ga_BASEQRANKSUM', 'ga_CLIPPINGRANKSUM', 'ga_DS', 'ga_END',
                              'ga_EXCESSHET', 'ga_FS', 'ga_INBREEDINGCOEFF', 'ga_MLEAC', 'ga_MLEAF', 'ga_MQ',
                              'ga_MQRANKSUM', 'ga_QD', 'ga_RAE_MQ', 'ga_READPOSRANKSUM', 'ga_SOR', 'ga_GT', 'ga_GQ',
                              'ga_PL',
                              'sm_CHROM', 'sm_POS', 'sm_REF', 'sm_ALT', 'sm_TYPE', 'sm_FILTER', 'sm_QUAL', 'sm_DP',
                              'sm_VF', 'sm_IDV', 'sm_IMF', 'sm_VDB', 'sm_RPB', 'sm_MQB', 'sm_BQB',
                              'sm_MQSB', 'sm_SGB', 'sm_MQ0F', 'sm_ICB', 'sm_HOB', 'sm_AC', 'sm_AN', 'sm_DP4', 'sm_MQ',
                              'sm_REF_F', 'sm_REF_R', 'sm_ALT_F', 'sm_ALT_R', 'sm_GT', 'sm_PL',
                              'st_CHROM', 'st_POS', 'st_REF', 'st_ALT', 'st_TYPE', 'st_FILTER', 'st_QUAL', 'st_DP',
                              'st_GT', 'st_DPF', 'st_MIN_DP', 'st_AD', 'st_ADF', 'st_ADR', 'st_FT', 'st_DPI', 'st_MQ',
                              'st_SB', 'st_PL', 'st_PS', 'st_GQ', 'st_GQX',
                              'va_CHROM', 'va_POS', 'va_REF', 'va_ALT', 'va_TYPE', 'va_FILTER', 'va_DP', 'va_VF',
                              'va_GT', 'va_READ1', 'va_READ2', 'va_PVALUE',
                              'va_R1F', 'va_R1R', 'va_R2F', 'va_R2R', 'va_PVAL'
                             ]) + '\n')
    vcfsout1.write('\t'.join(['Mut', 'ga_CHROM', 'ga_POS', 'ga_REF', 'ga_ALT', 'ga_TYPE',
                              'ga_FILTER', 'ga_QUAL', 'ga_DP', 'ga_AF',
                              'ga_AD', 'ga_AC', 'ga_AN', 'ga_BASEQRANKSUM', 'ga_CLIPPINGRANKSUM', 'ga_DS', 'ga_END',
                              'ga_EXCESSHET', 'ga_FS', 'ga_INBREEDINGCOEFF', 'ga_MLEAC', 'ga_MLEAF', 'ga_MQ',
                              'ga_MQRANKSUM', 'ga_QD', 'ga_RAE_MQ', 'ga_READPOSRANKSUM', 'ga_SOR', 'ga_GT', 'ga_GQ',
                              'ga_PL',
                              'sm_CHROM', 'sm_POS', 'sm_REF', 'sm_ALT', 'sm_TYPE', 'sm_FILTER', 'sm_QUAL', 'sm_DP',
                              'sm_VF', 'sm_IDV', 'sm_IMF', 'sm_VDB', 'sm_RPB', 'sm_MQB', 'sm_BQB',
                              'sm_MQSB', 'sm_SGB', 'sm_MQ0F', 'sm_ICB', 'sm_HOB', 'sm_AC', 'sm_AN', 'sm_DP4', 'sm_MQ',
                              'sm_REF_F', 'sm_REF_R', 'sm_ALT_F', 'sm_ALT_R', 'sm_GT', 'sm_PL',
                              'st_CHROM', 'st_POS', 'st_REF', 'st_ALT', 'st_TYPE', 'st_FILTER', 'st_QUAL', 'st_DP',
                              'st_GT', 'st_DPF', 'st_MIN_DP', 'st_AD', 'st_ADF', 'st_ADR', 'st_FT', 'st_DPI', 'st_MQ',
                              'st_SB', 'st_PL', 'st_PS', 'st_GQ', 'st_GQX',
                              'va_CHROM', 'va_POS', 'va_REF', 'va_ALT', 'va_TYPE', 'va_FILTER', 'va_DP', 'va_VF',
                              'va_GT', 'va_READ1', 'va_READ2', 'va_PVALUE',
                              'va_R1F', 'va_R1R', 'va_R2F', 'va_R2R', 'va_PVAL'
                              ]) + '\n')
    for mut in muts:
        if mut in gatk.keys():
            vcf = mut + '\t' + '\t'.join(gatk[mut][1:]) + '\t'
        else:
            vcf = mut + '\t' + '-\t'*29 + '-' + '\t'
        if mut in samtools.keys():
            vcf += '\t'.join(samtools[mut][1:]) + '\t'
        else:
            vcf += '-\t'*29 + '-' + '\t'
        if mut in strelka.keys():
            vcf += '\t'.join(strelka[mut][1:]) + '\t'
        else:
            vcf += '-\t' * 21 + '-' + '\t'
        if mut in varscan.keys():
            vcf += '\t'.join(varscan[mut][1:]) + '\n'
        else:
            vcf += '-\t' * 16 + '-' + '\n'
        vcfsout.write(vcf)
    vcfsout.close()

    for mut in muts:
        if mut in gatk.keys() and 'PASS' in gatk[mut]:
            vcf1 = mut + '\t' + '\t'.join(gatk[mut][1:]) + '\t'
        else:
            vcf1 = mut + '\t' + '-\t' * 29 + '-' + '\t'
        if mut in samtools.keys() and 'PASS' in samtools[mut]:
            vcf1 += '\t'.join(samtools[mut][1:]) + '\t'
        else:
            vcf1 += '-\t' * 29 + '-' + '\t'
        if mut in strelka.keys() and 'PASS' in strelka[mut]:
            vcf1 += '\t'.join(strelka[mut][1:]) + '\t'
        else:
            vcf1 += '-\t' * 21 + '-' + '\t'
        if mut in varscan.keys() and 'PASS' in varscan[mut]:
            vcf1 += '\t'.join(varscan[mut][1:]) + '\n'
        else:
            vcf1 += '-\t' * 16 + '-' + '\n'
        if 'PASS' in vcf1:
            vcfsout1.write(vcf1)
    vcfsout1.close()


#match genename,ENSG and ENST from ensembl.


def main():
    time_start = time.time()
    #(source, sample_name_file) = sys.argv[1:]
    cosmic = '/home/dell/Works/Projects/Datasets/annonDB/COSMIC_variant.csv'
    clinvar = '/home/dell/Works/Projects/Datasets/annonDB/breast_cancer_variant_clinvar.csv'
    g1000 = '/home/dell/Works/Projects/Datasets/annonDB/breast_cancer_variant_1000genomes.csv'
    ref_ens = '/home/dell/Works/Projects/Datasets/annonDB/geneid_cancer_qiagen.csv'
    variant_vcf = '/home/dell/Works/Projects/SampleSheet_180801/SampleSheet_180910/18083102_S23_L001/germline_vc/18083102_S23_L001.raw_variants.vcf'
    sample_name = '18083102_S23_L001'
    output = '/home/dell/Works/Projects/SampleSheet_180801/SampleSheet_180910/18083102_S23_L001/germline_vc/18083102_S23_L001.mutSummary.txt'
    callings = 'GATK,strelka2,samtools,varscan2'.split(',')

    snp_filter ='DP < 20 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
    indel_filter = 'DP < 20 || QD < 2.0 || FS > 200 || ReadPosRankSum < -3.0 || SOR > 10.0'
    other_filter = 'DP < 20'

    # #--
    # dict_cos, dict_clin, dict_g1000 = read_database(cosmic,clinvar,g1000)
    # #--
    snp_limit, indel_limit = read_vcf_filter(snp_filter, indel_filter)
    other_limit1, other_filter2 = read_vcf_filter(other_filter, other_filter)
    annotation(variant_vcf, output, sample_name, callings, snp_limit, indel_limit, other_limit1)
    # sample_name_list=open(sample_name_file,'r')
    # # for sample_name in sample_name_list:
    # #     sample_name = sample_name.strip()
    # #     variant_vcf = source + sample_name + '_L001.raw_variants.vcf'
    # #     annotated_csv = source + sample_name + '.annotated.txt'
    # #     annotated_csv_add = source + sample_name + '.annotated_add.txt'
    # #     non_rs = source + sample_name + '.non_rs.csv'
    # #     non_cos = source + sample_name + '.non_cos.csv'
    # #     stats_file = source + sample_name + '.annotate_stats.txt'
    # #     annotation(dict_cos,dict_clin,dict_g1000,variant_vcf,annotated_csv,stats_file, snp_limit, indel_limit,sample_name)
        #fill_table(annotated_csv,annotated_csv_add,non_rs,non_cos,ref_ens)
    print('The time of used annotation is %s minutes.' % str((time.time() - time_start) / 60))

if __name__ == '__main__':
    main()
