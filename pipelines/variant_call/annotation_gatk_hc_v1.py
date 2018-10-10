from __future__ import barry_as_FLUFL

import os
import sys
import time
import gzip
import pandas as pd
import subprocess
import shlex
sys.path.append("..")
from pipelines.log.log_v1 import store_annotation_logs

__all__ = ['cosmic', 'clinvar', 'g1000',
           'ref_ens', 'vcf', 'sample', 'output', 'snp_filter', 'indel_filter',
           'logger_annotation_process', 'logger_annotation_errors', 'calling']
__version__ = '1.0'
__author__ = 'Wang Xian'


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


# filter vcf by DP, QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum
def comp_filter(limists, value):
    result = []
    count = 0
    for para in limists.keys():
        if para == 'DP':
            if value[6] == '-':
                count += 1
                # result.append('No'+para)
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
        # elif para == 'QUAL':
        #    if value[5] == '-':
        #        count+=1
                # result.append('No'+para)
        #    elif limists[para][0] =='>':
        #        if float(limists[para][1]) >= float(value[5]):
        #            count+=1
        #        else:
        #            result.append('High'+para)
        #    else:
        #        if float(limists[para][1]) <= float(value[5]):
        #            count+=1
        #        else:
        #            result.append('Low'+para)
        elif para == 'QD':
            if value[14] == '-':
                count += 1
                # result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[14]):
                    count += 1
                else:
                    result.append('High'+para)
            elif limists[para][0] == '<':
                if float(limists[para][1]) <= float(value[14]):
                    count += 1
                else:
                    result.append('Low'+para)
        elif para == 'FS':
            if value[10] == '-':
                count += 1
                # result.append('No'+para)
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
        elif para == 'MQ':
            if value[12] == '-':
                count += 1
                # result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[12]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[12]):
                    count += 1
                else:
                    result.append('Low'+para)
        elif para == 'SOR':
            if value[16] == '-':
                count += 1
                # result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[16]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[16]):
                    count += 1
                else:
                    result.append('Low'+para)
        elif para == 'MQRankSum':
            if value[13] == '-':
                count += 1
                # result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[13]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[13]):
                    count += 1
                else:
                    result.append('Low'+para)
        elif para == 'ReadPosRankSum':
            if value[15] == '-':
                count += 1
                # result.append('No'+para)
            elif limists[para][0] == '>':
                if float(limists[para][1]) >= float(value[15]):
                    count += 1
                else:
                    result.append('High'+para)
            else:
                if float(limists[para][1]) <= float(value[15]):
                    count += 1
                else:
                    result.append('Low'+para)
    if count == len(limists.keys()):
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
        if strelka2[i] not in detail_dict.keys():
            for j in range(0, len(formats)):
                if formats[j] == strelka2[i]:
                    detail_dict[strelka2[i]] = details[j]
    return detail_dict


# split mutation. such as, ref:CTTT  alt:C,CT
def split_variant_strelka(line):
    chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
    mq = info.split(';')[len(info.split(';'))-1].split('=')[1]
    detail_dict = get_detail_from_strelka_vcf(format, detail)
    gt = detail_dict['GT']
    ad = detail_dict['AD']
    if "DPI" in format:
        sb = '-'
        dp = detail_dict['DPI']
    elif "DPF" in format:
        sb = detail_dict['SB']
        dp = detail_dict['DP']

    num_mt = len(alt.split(','))
    if num_mt is 1:
        dp0, dp1 = ad.split(',')
        dp0 = int(dp0)
        dp1 = int(dp1)
        if len(gt) < 3:
            af = gt
        else:
            af = gt[0]+','+gt[2]
        if dp0+dp1 == 0:
            # print(detail)
            vf = '-'
        else:
            vf = str(round(dp1/(dp0+dp1), 4))
        return [[chrom, pos, ref, alt, filter, qual, dp, af, ad, vf, mq, sb]]
    elif num_mt is 2:
        alt1, alt2 = alt.split(',')
        # af1, af2 = gt.split('/')
        # print(detail)
        dp0, dp1, dp2 = ad.split(',')
        dp0 = int(dp0)
        dp1 = int(dp1)
        dp2 = int(dp2)
        if dp0 == 0:
            # print(gt)
            af1 = '1,1'
            af2 = '2,2'
        else:
            af1 = '0,' + gt[0]
            af2 = '0,' + gt[1]
        ad1 = str(dp0) + ',' + str(dp1)
        ad2 = str(dp0) + ',' + str(dp2)
        if dp1 == 0:
            vf1 = 0
        else:
            vf1 = str(round(dp1 / (dp0 + dp1), 4))
        if dp2 == 0:
            vf2 = 0
        else:
            vf2 = str(round(dp2 / (dp0 + dp2), 4))
        return [
            [chrom, pos, ref, alt1, filter, qual, dp, af1, ad1, str(vf1), mq, sb],
            [chrom, pos, ref, alt2, filter, qual, dp, af2, ad2, str(vf2), mq, sb]
            ]
    elif num_mt is 3:
        return 0


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
        if num_mt is 1:
            dp0, dp1 = ad.split(',')
            dp0 = int(dp0)
            dp1 = int(dp1)
            # avoid the gt sep by / or | ,that is not good to show in  txt or csv files
            # af=gt.replace('/',',')
            af = gt[0]+','+gt[2]
            vf = str(round(dp1/(dp0+dp1), 4))
            return [[chrom, pos, ref, alt, filter, qual, dp, af, ad, vf, ac, an, baseqranksum, clippingranksum, ds,
                     end, excesshet, fs, inbreedingcoeff, mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum,
                     sor, gt, gq, pl]]
        elif num_mt is 2:
            alt1, alt2 = alt.split(',')
            dp0, dp1, dp2 = ad.split(',')
            dp0 = int(dp0)
            dp1 = int(dp1)
            dp2 = int(dp2)
            if dp0 == 0:
                # print(gt)
                af1 = '1,1'
                af2 = '2,2'
            else:
                af1 = '0,' + gt[0]
                af2 = '0,' + gt[1]
            ad1 = str(dp0) + ',' + str(dp1)
            ad2 = str(dp0) + ',' + str(dp2)
            vf1 = str(round(dp1/(dp0+dp1), 4))
            vf2 = str(round(dp2/(dp0+dp2), 4))
            return [[chrom, pos, ref, alt1, filter, qual, dp, af1, ad1, vf1, ac, an, baseqranksum,
                     clippingranksum,  ds, end, excesshet, fs, inbreedingcoeff,  mleac, mleaf, mq,
                     mqranksum, qd, rae_mq, readposranksum, sor, gt, gq, pl],
                    [chrom, pos, ref, alt2, filter, qual, dp, af2, ad2, vf2, ac, an, baseqranksum,
                     clippingranksum,  ds, end, excesshet, fs, inbreedingcoeff, mleac, mleaf, mq,
                     mqranksum, qd, rae_mq, readposranksum, sor, gt, gq, pl]]
    elif calling == 'samtools':
        chrom, pos, id, ref, alt, qual, filter, info, format, detail = line.strip().split('\t')
        indel, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb, sgb, mq0f, icb, hob, ac, an, dp4, mq = get_info(info, calling)
        gt, pl = detail.split(':')
        ref_f, ref_r, alt_f, alt_r = dp4.split(',')
        vf = str(round((float(alt_f)+float(alt_r))/(float(alt_f)+float(alt_r)+float(ref_f)+float(ref_r)), 4))
        num_mt = len(alt.split(','))
        if num_mt == 1:
            return[[chrom, pos, ref, alt, qual, filter, vf, indel, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb, sgb,
                    mq0f, icb, hob, ac, an, dp4, mq, ref_f, ref_r, alt_f, alt_r, gt, pl]]
        elif num_mt == 2:
            alt1, alt2 = alt.split(',')
            return [[chrom, pos, ref, alt1, qual, filter, vf, indel, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb, sgb,
                     mq0f, icb, hob, ac, an, dp4, mq, ref_f, ref_r, alt_f, alt_r, gt, pl],
                    [chrom, pos, ref, alt2, qual, filter, vf, indel, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb, sgb,
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
        elif alt[0] == '+':
            ref1 = ref
            alt = ref + alt[1:]
        else:
            ref1 = ref
        return[[chrom, pos, ref1, alt, gt, cov, read1, read2, vf, pvalue, strandfilter,
                r1f, r1r, r2f, r2r, pval]]
    elif calling == 'smcounter':
        # get the values of samtools vcf parameters form format and detail
        chrom, pos, ref, alt, type, dp, mt, umt, pi, thr, vmt, vmf, vsm, filter= line.strip().split('\t')
        if filter == 'LSM':
            filter = 'PASS'
        return[[chrom, pos, ref, alt, type, dp, mt, umt, pi, thr, vmt, vmf, vsm, filter]]


def read_vcf(variant_vcf, output, sample_name, calling, snp_filter, indel_filter):
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
        output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'DP', 'VF', 'GT', 'DP_read1',
                                'DP_read2', 'Fisherexact_var', 'read1F', 'read1R', 'read2F', 'read2R', 'StrandBias',
                                'Gene_ID', 'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID', 'Mutation_Description',
                                'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
                                'Mutation_Strand', 'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score',
                                'Mutation_Somatic_Status',  'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF',
                                'AMR_AF', 'SAS_AF', 'AFR_AF\n']))
        # print(variant_vcf)
        for line in var:
            if not line.startswith('Chrom'):
                for spl in split_variant(line, calling):
                    (chrom, pos, ref, alt, gt, dp, read1, read2, vf, pvalue, filter,
                     r1f, r1r, r2f, r2r, pval) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, dp, vf]
                    # print(value)
                    # upper the ref and alt
                    ref = ref.upper()
                    alt = alt.upper()
                    # ref = sub('[+-]', '', ref)
                    # alt = sub('[+-]', '', ref)
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
                    value = [sample_name, chrom, pos, ref, alt, filter, dp, vf, gt, read1, read2, pvalue,
                             r1f, r1r, r2f, r2r, pval]
                    yield [key, key1, value]
    elif calling == 'strelka2':
        variant_vcf = os.path.dirname(variant_vcf) + '/results/variants/variants.vcf.gz'
        if not os.path.isfile(os.path.dirname(variant_vcf) + '/variants.vcf'):
            command1 = 'gunzip' + ' -k ' + variant_vcf
            stdout, stderr = stdout_err(command1)
        # output = open(annotated_csv, 'w')
        output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'QUAL', 'DP', 'GT', 'AD', 'VF', 'MQ',
                                'SB', 'Gene_ID', 'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID',
                                'Mutation_Description',
                                'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
                                'Mutation_Strand',
                                'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score', 'Mutation_Somatic_Status',
                                'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF', 'AMR_AF', 'SAS_AF', 'AFR_AF\n']))
        var = open(os.path.dirname(variant_vcf) + '/variants.vcf', 'r')
        for line in var:
            if not line.startswith('#'):
                for spl in split_variant_strelka(line):
                    chrom, pos, ref, alt, filter, qual, dp, af, ad, vf, mq, sb = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, qual, dp, af, ad, vf, mq, sb]
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
                    value = [sample_name, chrom, pos, ref, alt, filter, qual, dp, af, ad, vf, mq, sb]
                    yield [key, key1, value]
    elif calling == 'GATK':
        var = open(variant_vcf, 'r')
        # output = open(annotated_csv, 'w')
        output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'QUAL', 'DP', 'GT1', 'AD', 'VF', 'AC',
                                'AN', 'Baseqranksum', 'ClippingRankSum', 'DS', 'END', 'ExcessHet', 'FS',
                                'InbreedingCoeff', 'MLEAC', 'MLEAF', 'MQ', 'MQRankSum', 'QD', 'RAW_MQ',
                                'ReadPosRankSum', 'SOR', 'GT', 'GQ', 'PL',
                                'Gene_ID', 'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID',
                                'Mutation_Description',
                                'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
                                'Mutation_Strand', 'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score',
                                'Mutation_Somatic_Status', 'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF', 'AMR_AF',
                                'SAS_AF', 'AFR_AF\n']))
        for line in var:
            if not line.startswith('#'):
                for spl in split_variant(line, calling):
                    (chrom, pos, ref, alt, filter, qual, dp, af, ad, vf, ac, an, baseqranksum, clippingranksum,
                     ds, end, excesshet, fs, inbreedingcoeff,
                     mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor, gt, gq, pl) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, qual, dp, ad, vf, baseqranksum, fs, inbreedingcoeff,
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
                    value = [sample_name, chrom, pos, ref, alt, filter, qual, dp, af, ad, vf, ac, an, baseqranksum,
                             clippingranksum, ds, end, excesshet, fs, inbreedingcoeff,
                             mleac, mleaf, mq, mqranksum, qd, rae_mq, readposranksum, sor, gt, gq, pl]
                    yield [key, key1, value]
    elif calling == 'samtools':
        variant_vcf = os.path.dirname(variant_vcf) + '/' + sample_name + '.raw_samtools.vcf'
        var = open(variant_vcf, 'r')
        output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'QUAL', 'DP', 'VF', 'INDEL', 'IDV',
                                'IMF', 'VDB', 'RPB', 'MQB', 'BQB', 'MQSB', 'SGB', 'MQ0F', 'ICB', 'HOB',
                                'AC', 'AN', 'DP4', 'MQ', 'ref_f', 'ref_r', 'alt_f', 'alt_r', 'GT', 'PL', 'Gene_ID',
                                'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID', 'Mutation_Description',
                                'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
                                'Mutation_Strand',
                                'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score', 'Mutation_Somatic_Status',
                                'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF', 'AMR_AF', 'SAS_AF', 'AFR_AF\n']))
        # print(variant_vcf)
        for line in var:
            if not line.startswith('#'):
                for spl in split_variant(line, calling):
                    (chrom, pos, ref, alt, qual, filter, vf, indel, idv, imf, dp, vdb, rpb, mqb, bqb, mqsb,
                     sgb, mq0f, icb, hob, ac, an, dp4, mq, ref_f, ref_r, alt_f, alt_r, gt, pl) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, qual, dp, vf, indel, idv, mq, vdb]
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
                    value = [sample_name, chrom, pos, ref, alt, filter, qual, dp, vf, indel, idv, imf, vdb, rpb,
                             mqb, bqb, mqsb, sgb, mq0f, icb, hob,
                             ac, an, dp4, mq, ref_f, ref_r, alt_f, alt_r, gt, pl]
                    yield [key, key1, value]
    elif calling == 'smcounter':
        variant_vcf = os.path.dirname(variant_vcf) + '/' + sample_name + '.smCounter.cut.txt'
        var = open(variant_vcf, 'r')
        output.write('\t'.join(['Sample', 'CHR', 'POS', 'REF', 'ALT', 'FILTER', 'TYPE', 'DP', 'MT', 'UMT', 
                                'PI', 'THR', 'VMT', 'VMF', 'VSM','Gene_ID',
                                'RS_ID', 'CLNDN', 'HGVS', 'CLNSIG', 'COSMIC_ID', 'Mutation_Description',
                                'Feature_ID', 'Gene_Name', 'Gene_CDS_Length', 'Mutation_Zygosity', 'LOH',
                                'Mutation_Strand',
                                'HGVS.c', 'HGVS.p', 'FATHMM_Prediction', 'FATHMM_Score', 'Mutation_Somatic_Status',
                                'Gene_Name1', 'RS_ID1', 'EAS_AF', 'EUR_AF', 'AMR_AF', 'SAS_AF', 'AFR_AF\n']))
        # print(variant_vcf)
        for line in var:
            if not line.startswith('CHROM'):
                for spl in split_variant(line, calling):
                    (chrom, pos, ref, alt, type, dp, mt, umt, pi, thr, vmt, vmf, vsm, filter) = spl
                    chrom = chrom[3:]
                    value = [chrom, pos, ref, alt, filter, type, dp, mt, umt, pi, thr, vmt, vmf, vsm]
                    # print(value)
                    if len(ref) == len(alt) and len(alt) == 1:
                        change = ref + '>' + alt
                        change1 = base_paired[ref] + '>' + base_paired[alt]
                        if filter == 'PASS':
                            filter = comp_filter_strelka(snp_filter, value)
                    elif len(ref) > len(alt) and len(alt) == 1:
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
                    value = [sample_name, chrom, pos, ref, alt, filter, type, dp, mt, umt, pi, thr, vmt, vmf, vsm]
                    yield [key, key1, value]


def annotation_v(dict_cos, dict_clin, dict_g1000, variant_vcf, annotated_csv,
                 stats_file, snp_filter, indel_filter, sample_name, logger_annotation_process, calling):
    key_list = []
    num_in_clinvar = 0
    num_in_cosmic = 0
    num_in_g1000 = 0
    num_unmatch = 0
    output = open(annotated_csv, 'w')
    for keys in zip(read_vcf(variant_vcf, output, sample_name, calling, snp_filter, indel_filter)):
        # print(keys)
        key = keys[0][0]
        key1 = keys[0][1]
        value = keys[0][2]
        # print(value)
        unmatch = 0
        # drop duplicate variant
        if key in key_list:
            continue
        if key in dict_clin:
            new = '\t'.join(value) + '\t' + dict_clin[key].replace(',', '\t') + '\t'
            num_in_clinvar += 1
        else:
            new = '\t'.join(value) + '\t' + '-\t'*3 + define_hgvs(value[1], value[2], value[3], value[4]) + '\t-\t'
            unmatch += 1
        if key in dict_cos:
            new += dict_cos[key].replace(',', '\t') + '\t'
            num_in_cosmic += 1
        elif key1 in dict_cos:
            new += dict_cos[key1].replace(',', '\t') + '\t'
            num_in_cosmic += 1
        else:
            new += '-\t'*13
            unmatch += 1
        if key in dict_g1000:
            new += dict_g1000[key].replace(',', '\t') + '\n'
            num_in_g1000 += 1
        else:
            new += '-\t'*6 + '-\n'
        # 3 databases both unmatch
        if unmatch == 3:
            num_unmatch += 1
        else:
            output.write(new)
            key_list.append(key)
    output.close()
    store_annotation_logs(logger_annotation_process,
                          'null', '{0} has {1} variants.'.format(sample_name, key_list))
    store_annotation_logs(logger_annotation_process,
                          'null', '{0} has {1} variants in COSMIC database'.format(sample_name, num_in_cosmic))
    store_annotation_logs(logger_annotation_process,
                          'null', '{0} has {1} variants in Clinvar database'.format(sample_name, num_in_clinvar))
    store_annotation_logs(logger_annotation_process,
                          'null', '{0} has {1}variants in G1000 database'.format(sample_name, num_in_g1000))
    store_annotation_logs(logger_annotation_process, 'null',
                          '{0} has {1} variants unmatch in cosmic and clinvar.'.format(sample_name, num_unmatch))
    stats_out = open(stats_file, 'w')
    stats_out.write('#The sample has %s variants.\n' % len(key_list))
    stats_out.write('#Type\tvariants\n')
    stats_out.write('Variants\t%s\n' % len(key_list))
    stats_out.write('COSMIC\t%s\n' % num_in_cosmic)
    stats_out.write('Clinvar\t%s\n' % num_in_clinvar)
    stats_out.write('G1000\t%s\n' % num_in_g1000)
    stats_out.write('unmatch\t%s\n' % num_unmatch)
    stats_out.close()


# match genename,ENSG and ENST from ensembl.
def fill_table(annotated_csv, annotated_csv_add, ref_ens):
    n2g = {}
    g2n = {}
    for line in open(ref_ens, 'r').readlines():
        name, ensg = line.strip().split(',')
        n2g[name] = ensg
        g2n[ensg] = name
    # df = pd.read_csv(annotated_csv)
    df = pd.read_table(annotated_csv)
    # n,g,t,n1
    subframe = df[['Gene_Name', 'Gene_ID', 'Gene_Name1', 'RS_ID', 'RS_ID1']]
    # for name,id,transcript in subframe.iterrows():
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
    df.insert(6, 'Gene_Name', name)
    df.insert(7, 'Gene_ID', ensg)
    df.insert(8, 'RS_ID', rs)
    df.to_csv(annotated_csv, index=False, sep='\t')
    df[(False ^ df['FILTER'].isin(['PASS']))].to_csv(annotated_csv_add, index=False, sep='\t')


# -annotation main
def annotationmain(cosmic, clinvar, g1000, 
                   ref_ens,
                   vcf, sample, snp_filter, indel_filter,
                   output, logger_annotation_process, logger_annotation_errors, calling):
    callings = calling.split('\t')
    dict_cos, dict_clin, dict_g1000 = read_database(cosmic, clinvar, g1000)
    for call in callings:
        if 'raw_variants_SNP.vcf' in vcf:
            annotated_csv = output + '/' + sample + '.raw_SNP.' + call + '.txt'
            annotated_csv_add = output + '/' + sample + '.raw_SNP.' + call + '_PASS.txt'
            non_rs = output + '/' + sample + '.raw_SNP_non_rs.txt'
            non_cos = output + '/' + sample + '.raw_SNP_non_cos.txt'
            stats_file = output + '/' + sample + '.raw_SNP_' + call + '_stats.txt'
        elif 'raw_variants_indel.vcf' in vcf:
            annotated_csv = output + '/' + sample + '.raw_indel.' + call + '.txt'
            annotated_csv_add = output + '/' + sample + '.raw_indel.' + call + '_PASS.txt'
            non_rs = output + '/' + sample + '.raw_indel_non_rs.txt'
            non_cos = output + '/' + sample + '.raw_indel_non_cos.txt'
            stats_file = output + '/' + sample + '.raw_indel_' + call + '_stats.txt'
        elif 'filter_SNP.vcf' in vcf:
            annotated_csv = output + '/' + sample + '.filter_SNP.' + call + '.txt'
            annotated_csv_add = output + '/' + sample + '.filter_SNP.' + call + '_PASS.txt'
            non_rs = output + '/' + sample + '.filter_SNP_non_rs.txt'
            non_cos = output + '/' + sample + '.filter_SNP_non_cos.txt'
            stats_file = output + '/' + sample + '.filter_SNP_' + call + '_stats.txt'
        elif 'filter_indel.vcf' in vcf:
            annotated_csv = output + '/' + sample + '.filter_indel.' + call + '.txt'
            annotated_csv_add = output + '/' + sample + '.filter_indel.' + call + '_PASS.txt'
            non_rs = output + '/' + sample + '.filter_indel_non_rs.txt'
            non_cos = output + '/' + sample + '.filter_indel_non_cos.txt'
            stats_file = output + '/' + sample + '.filter_indel_' + call + '_stats.txt'
        elif 'raw_variants.vcf' in vcf:
            annotated_csv = output + '/' + sample + '.raw_variants.' + call + '.txt'
            annotated_csv_add = output + '/' + sample + '.raw_variants.' + call + '_PASS.txt'
            non_rs = output + '/' + sample + '.raw_variants_non_rs.txt'
            non_cos = output + '/' + sample + '.raw_variants_non_cos.txt'
            stats_file = output + '/' + sample + '.raw_variants_' + call + '_stats.txt'
        # - read the annotation database
        if not os.path.isfile(vcf) and call == 'GATK':
            store_annotation_logs('null', logger_annotation_errors, vcf + " does not exist!\n")
            # print(vcf + ' does not exist!')
        else:
            # --annotation
            annotation_v(dict_cos, dict_clin, dict_g1000, vcf, annotated_csv, stats_file, snp_filter,
                         indel_filter, sample, logger_annotation_process, call)
            # --add the annotation
            fill_table(annotated_csv, annotated_csv_add, ref_ens)

