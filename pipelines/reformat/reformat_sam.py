from __future__ import barry_as_FLUFL

__all__ = ['alignment_sam' , 'output_sam' , 'logger_reformat_process' , 'logger_reformat_errors']
__version__ = '1.0'
__author__ = 'Maggie Ruimin Sun'

import logging
import os
import sys
import time

hg19_chr_length = {
        'chr1': 249250621,
        'chr2': 243199373,
        'chr3': 198022430,
        'chr4': 191154276,
        'chr5': 180915260,
        'chr6': 171115067,
        'chr7': 159138663,
        'chr8': 146364022,
        'chr9': 141213431,
        'chr10': 135534747,
        'chr11': 135006516,
        'chr12': 133851895,
        'chr13': 115169878,
        'chr14': 107349540,
        'chr15': 102531392,
        'chr16': 90354753,
        'chr17': 81195210,
        'chr18': 78077248,
        'chr19': 59128983,
        'chr20': 63025520,
        'chr21': 48129895,
        'chr22': 51304566,
        'chrX': 155270560,
        'chrY': 59373566,
        'chrM': 16571,
    }

"""
The numbers of consolidated sequences in R1 and R2 are different, because those
consolidated sequences with too many N's are removed. 
"""


def reformat_sam(alignment_sam, output_sam, logger_reformat_process, logger_reformat_errors):
    sam = open(alignment_sam)
    sam.readline()
    sam_out = open(output_sam, 'w')
    for chrom in hg19_chr_length:
        sam_out.write('@SQ\tSN:' + chrom + '\tLN:' + str(hg19_chr_length[chrom]) + '\n')
    for row in sam:
        if row[0:3] == '@SQ':
            continue
        if row[0:3] == '@PG':
            sam_out.write(row)
            continue
        qname, flag, rname, pos, mapq, cigar, rmate, pmate = row.strip().split()[0:8]
        chrom, start = rname.split('_')[0:2]
        pos = str(int(pos) + int(start))
        pmate = str(int(pmate) + int(start))
        if len(qname.split('_')) == 2:
            qname, umi = qname.split('_')
            if flag == '99' or flag == '163':
                qname += ':' + chrom + '-0-' + pos + '-' + umi + ':' + umi
            elif flag == '83' or flag == '147':
                qname += ':' + chrom + '-1-' + pos + '-' + umi + ':' + umi
            else:
                continue  # error flags
        elif len(qname.split('_')) == 1:
            if flag == '99' or flag == '163':
                qname += ':' + chrom + '-0-' + pos
            elif flag == '83' or flag == '147':
                qname += ':' + chrom + '-1-' + pos
            else:
                continue  # error flags
        sam_out.write('\t'.join([qname, flag, chrom, pos, mapq,
                                 cigar, rmate, pmate]) + '\t' + '\t'.join(row.split()[8:]) + '\n')
    sam.close()
    sam_out.close()
