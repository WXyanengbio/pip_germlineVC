from __future__ import barry_as_FLUFL

__all__ = ['samtools_dir', 'alignment_sam', 'min_mapq', 'max_soft_clip', 'out_file', 'stats_file',
           'primers_file', 'primer_stats_file', 'max_dist', 'logger_filter_process', 'logger_filter_errors']
__version__ = '1.0'
__author__ = 'Maggie Ruimin Sun'

import os
import sys
import re
import logging
import time
from difflib import SequenceMatcher
sys.path.append("..")
from pipelines.log.log_v1 import store_filter_logs
time_start = time.time()


def filter_alignment_samtools(samtools_dir, alignment_sam, min_mapq,
                              max_soft_clip, out_file, stats_file,
                              logger_filter_process, logger_filter_errors):
    stats_file_tmp = stats_file + '.tmp'
    command_count = '{0} view {1} | cut -f1 | uniq | wc -l >> {2}'.format(
        samtools_dir, alignment_sam, stats_file_tmp)
    store_filter_logs(logger_filter_process, 'null', 'Samtools counts the total number of read pairs.' + '\n')
    os.system(command_count)

    # Count supplimentary alignments
    command_count1 = '{0} view -f 2048 {1} | cut -f1 | uniq | wc -l >> {2}'.format(
        samtools_dir, alignment_sam, stats_file_tmp)
    store_filter_logs(logger_filter_process, 'null', 'Samtools counts the supplimentary alignments.' + '\n')
    os.system(command_count1)

    out_file_tmp1 = out_file + '_tmp1.sam'
    command_samtools = '{0} view -ShF 2048 {1} > {2}'.format(
        samtools_dir, alignment_sam, out_file_tmp1)
    os.system(command_samtools)

    command_count = '{0} view {1} | cut -f1 | uniq | wc -l >> {2}'.format(
        samtools_dir, out_file_tmp1, stats_file_tmp)
    store_filter_logs(logger_filter_process,'null', 'Samtools counts the number of read pairs without 2048.' + '\n')
    os.system(command_count)

    # Count unmapped reads
    command_count2 = '{0} view -f 8 {1} | cut -f1 | uniq | wc -l >> {2}'.format(
        samtools_dir, out_file_tmp1, stats_file_tmp)
    store_filter_logs(logger_filter_process, 'null', 'Samtools counts the total number of unmapped read pairs' + '\n')
    os.system(command_count2)
    # Discard all unmapped read pairs
    out_file_tmp2 = out_file + '_tmp2.sam'
    command_samtools = '{0} view -ShF 8 {1} > {2}'.format(samtools_dir, out_file_tmp1, out_file_tmp2)
    os.system(command_samtools)
    command_samtools = '{0} view -ShF 4 {1} > {2}'.format(samtools_dir, out_file_tmp2, out_file_tmp1)
    os.system(command_samtools)
    # Discard read pairs mapped to different target sequences
    command_samtools1 = "{0} view -Sh {1} \
    | perl -lane 'print if $F[0] =~ /^@/; print if $F[6] =~ /=/;' > {2}".format(
        samtools_dir, out_file_tmp1, out_file_tmp2)
    store_filter_logs(logger_filter_process, 'null', 'Samtools discards the secondary/supplimentary/unmapped'
                      + ' alignments and those read pairs mapped to different target regions.' + '\n')
    os.system(command_samtools1)
    command_count3 = "{0} view {1} | cut -f1 | uniq | wc -l >> {2}".format(
        samtools_dir, out_file_tmp2, stats_file_tmp)
    os.system(command_count3)
    # Discard alignments with MAPQ < min_mapq
    command_count4 = '{0} view -q {1} {2} | cut -f1 | uniq | wc -l >> {3}'.format(
        samtools_dir, min_mapq, out_file_tmp2, stats_file_tmp)
    store_filter_logs(logger_filter_process, 'null', 'Samtools counts the alignments with higher Maps' + '\n')
    os.system(command_count4)
    command_samtools2 = "{0} view -Shq {1} {2} > {3}".format(samtools_dir, min_mapq,
                                                             out_file_tmp2, out_file)
    os.system(command_samtools2)
    store_filter_logs(logger_filter_process, 'null',
                      'Samtools discards the alignments with lower MapQ and mismatched beginning' + '\n')
    command_count5 = '{0} view {1} | cut -f1 | uniq | wc -l >> {2}'.format(samtools_dir, out_file, stats_file_tmp)
    os.system(command_count5)

    # Output the alignment statistics.
    stats = open(stats_file_tmp)
    (total_count, sec_count, rm_sec_count, unmap_count, same_chr_count,
     hmapq_count, final_count) = [int(x.strip()) for x in stats.readlines()[0:7]]
    stats.close()
    stats_out = open(stats_file, 'w')
    stats_out.write('Total number of read pairs == {0}\n'.format(total_count))
    stats_out.write('Number of secondary/supplimentary alignments == {0}\n'.format(sec_count))
    stats_out.write('Number of unmapped read pairs == {0}\n'.format(unmap_count))
    remained = total_count - unmap_count
    stats_out.write('Threshold MAPQ for filtration == {0}\n'.format(min_mapq))
    stats_out.write('Number of low MAPQ alignments (MAPQ<{0}) == {1}\n'.format(min_mapq, same_chr_count - hmapq_count))
    stats_out.write('Number of alignments passed SAMTools filtration == {0}\n'.format(final_count))
    stats_out.write(
        'Percent of alignments passed SAMTools filtration == {0}(%)\n'.format(100 * final_count / total_count))
    stats_out.close()
    os.system('rm ' + stats_file_tmp)
    os.system('rm ' + out_file_tmp1 + ' ' + out_file_tmp2)


def complement_base(base):
    if base == 'A':
        return 'T'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    elif base == 'T':
        return 'A'
    else:
        return base


def reverse_complement(seq):
    seq_rc = ''
    ls = len(seq)
    for i in range(ls - 1, -1, -1):
        seq_rc += complement_base(seq[i])
    return seq_rc


def compare_seqs(primer, seq, strand):
    matcher = SequenceMatcher(None, primer, seq)
    dist = len(primer)
    i = 0
    j = 0
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == 'equal':
            if i == 0:
                dist -= i1
                if strand == 0:
                    match_start = j1
            dist -= (i2-i1)
            j = j2
            i += 1
    if strand == 1:
        match_start = len(seq) - j
    return dist, match_start


def identify_gs_primers(samtools_dir, alignment_sam, primers_file, max_dist,
                        out_file, stats_file, primer_stats_file,
                        logger_filter_process, logger_filter_errors):
    primers = {}
    primer_pos = {}
    csv = open(primers_file)
    # csv.readline()
    for row in csv:
        chrom, strand, seq, start, stop, gene = row.strip().split(',')[0:6]
        if chrom == 'chr':
            continue
        pos_start = int(start)
        pos_stop = int(stop)
        if chrom not in primer_pos:
            primer_pos[chrom] = {}
            
        if strand == '0':
            for pos in range(pos_start-100, pos_stop+101):
                if pos not in primer_pos[chrom]:
                    primer_pos[chrom][pos] = []
                primer_pos[chrom][pos].append((seq, chrom+'_'+start+'_'+stop))
        else:
            for pos in range(pos_stop-100, pos_start+101):
                if pos not in primer_pos[chrom]:
                    primer_pos[chrom][pos] = []
                primer_pos[chrom][pos].append((reverse_complement(seq), chrom+'_'+stop+'_'+start))

        if strand == '0':
            key_primer = chrom+'_'+start+'_'+stop
            primers[key_primer] = {}
            primers[key_primer]['strand'] = '0'
            primers[key_primer]['seq'] = seq
        else:
            key_primer = chrom+'_'+stop+'_'+start
            primers[key_primer] = {}
            primers[key_primer]['strand'] = '1'
            primers[key_primer]['seq'] = reverse_complement(seq)
        primers[key_primer]['gene'] = gene
        primers[key_primer]['chromosome'] = chrom
        primers[key_primer]['start'] = start
        primers[key_primer]['stop'] = stop
        primers[key_primer]['length'] = len(seq)
        primers[key_primer]['on-target'] = 0
        # primers[key_primer]['mis-target'] = []

    csv.close()
    num_primers = len(primers)
    print('Number of primers of interest == '+str(num_primers))
    
    num_total_alignments = 0
    num_short_fragment = 0
    num_on_target = 0
    num_off_target = 0
    num_off_target_wrong_strand = 0
    num_unproper_pairs = 0
    off_target_reads = {}
    
    soft_clips_start = re.compile('^(\d+)S')
    soft_clips_end = re.compile('(\d+)S$')
    flag_eligible = ['83', '99']
    
    sam = open(alignment_sam)
    ready = open(out_file, 'w')
    off_sam = open(alignment_sam + '.off', 'w')
    proper_pair = True
    sam_rows = sam.readlines()
    while True:
        if sam_rows[0][0] == '@':
            ready.write(sam_rows[0])
            sam_rows.remove(sam_rows[0])
        else:
            break
    
    sam_rows = sorted(sam_rows)
    i = 0
    while i < len(sam_rows)-1:
        row = sam_rows[i]
        mate = sam_rows[i+1]
        qname, flag, rname, pos, mapq, cigar, rmn, pmn, insl, seq = row.strip().split()[0:10]
        if qname in mate:
            proper_pair = True
            i += 2
        else:
            proper_pair = False
            i += 1

        if flag not in flag_eligible:
            if int(flag) > 128:
                flag_mate = mate.split()[1]
                if flag_mate in flag_eligible:
                    qname, flag, rname, pos, mapq, cigar, rmn, pmn, insl, seq = mate.strip().split()[0:10]
                    row, mate = mate, row
            else:    
                proper_pair = False
        
        if not proper_pair:
            num_unproper_pairs += 1
            continue
                
        num_total_alignments += 1
        chrom, start, stop = rname.split('_')
        if chrom not in primer_pos:
            num_off_target += 1
            continue
        pos = int(pos) + int(start) - 1
        pos_rv = pos + len(seq) - 1

        if soft_clips_start.search(cigar):
            pos_rv = pos_rv - int(soft_clips_start.search(cigar).groups()[0])
        if soft_clips_end.search(cigar):
            pos_rv = pos_rv-int(soft_clips_end.search(cigar).groups()[0]) 

        off_target = 0
        off = False
        min_dist = max_dist
        match_start0 = len(seq)
        if flag == '83':  # read 1 is mapped to the reversed strand
            if pos_rv not in primer_pos[chrom]:
                off = True
                num_off_target += 1
                num_off_target_wrong_strand += 1
                continue
            for primer_info in primer_pos[chrom][pos_rv]:
                primer_name = primer_info[1]
                primer_seq = primer_info[0]
                lpr = primers[primer_name]['length']
                if primers[primer_name]['strand'] == '0':
                    off_target += 1
                    continue
                dist, match_start = compare_seqs(primer_seq, seq, 1)
                if dist > max_dist:
                    off_target += 1
                elif dist <= min_dist and match_start < match_start0:
                    min_dist = dist
                    match_start0 = match_start
                    primer_on_target = primer_name

            if off_target == len(primer_pos[chrom][pos_rv]):
                off = True
        elif flag == '99':
            if pos in primer_pos[chrom]:
                for primer_info in primer_pos[chrom][pos]:
                    primer_name = primer_info[1]
                    primer_seq = primer_info[0]
                    lpr = primers[primer_name]['length']
                    if primers[primer_name]['strand'] == '1':
                        off_target += 1
                        continue
                    dist, match_start = compare_seqs(primer_seq, seq, 0)
                    if dist > max_dist:
                        off_target += 1
                    elif dist <= min_dist and match_start < match_start0:
                        min_dist = dist
                        match_start0 = match_start
                        primer_on_target = primer_name
                    
                if off_target == len(primer_pos[chrom][pos]):
                    off = True
            else:
                num_off_target += 1
                num_off_target_wrong_strand += 1
                continue
        if off:
            num_off_target += 1
            off_target_reads[rname+'_'+str(pos)] = seq
            off_sam.write(row)
            off_sam.write(mate)
        else:
            num_on_target += 1
            primers[primer_on_target]['on-target'] += 1
            ready.write(row)
            ready.write(mate)
    sam.close()
    ready.close()
    
    ratio_off = 100*num_off_target / num_total_alignments
    ratio_on = 100*num_on_target / num_total_alignments
    # print('Total number of alignments (PE) == ' + str(num_total_alignments))
    # print('Number of reads mapped in unproper pairs == ' + str(num_unproper_pairs))
    # print('Number of off-target alignments == {0} ({1}%)({2})'.format(num_off_target,
    # ratio_off, num_off_target_wrong_strand))
    # print('Number of on-target alignments == {0} ({1}%)'.format(num_on_target, ratio_on))

    stats_out = open(stats_file, 'a')
    stats_out.write('Number of unproperly-paired alignments == '+str(num_unproper_pairs)+'\n')
    stats_out.write('Number of panel primers == '+str(num_primers)+'\n')
    stats_out.write('Number of off-target alignments == {0} ({1}%)\n'.format(num_off_target, ratio_off))
    stats_out.write('Number of on-target alignments == {0} ({1}%)\n'.format(num_on_target, ratio_on))
    stats_out.close()
    
    st_out = open(primer_stats_file, 'w')
    st_out.write('chrm,F/R,start,stop,sequence,#on-target-reads,%on-target-reads\n')
    for pr in sorted(primers):
        ratio = 100*primers[pr]['on-target'] / float(num_on_target)
        
        st_out.write(','.join([primers[pr]['chromosome'], primers[pr]['strand'], 
                               primers[pr]['start'], primers[pr]['stop'],
                               primers[pr]['seq'], str(primers[pr]['on-target']),
                               str(ratio)])+'\n')
    st_out.close()
