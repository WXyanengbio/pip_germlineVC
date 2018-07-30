from __future__ import barry_as_FLUFL

__all__  =  ['bwa_dir' , 'ref_fa_file' , 'ref_index_name','total_ref_fa_file' ,'exome_target_bed' ,'read1' , 'read2' , 'out_file' , 'num_threads' , 'logger_bwa_process' , 'logger_bwa_errors']
__version__  =  '1.0'
__author__  =  'Wang Xian'


import os
import logging
import time
import sys

def filter_ref_fa_by_bed(ref_fa_file,ref_index_name, exome_target_bed,logger_bwa_process, 
                    logger_bwa_errors):
    if not os.path.isfile(exome_target_bed):
        logger_bwa_errors.error('%s does not exist!', exome_target_bed)
        print("Error: cannot find target bed file!")
        exit()
    if not os.path.isfile(ref_fa_file):
        logger_bwa_errors.error('%s does not exist!', ref_fa_file)
        print("Error: cannot find reference genome file!")
        exit()
    if not os.path.isfile(ref_index_name + '.refSeq.fa'):
        fasta= open(ref_fa_file, 'U')
        fasta_dict= {}
        for line in fasta:
            line= line.strip()
            if line == '':
                 continue
            if line.startswith('>'):
                 seqname= line.lstrip('>')
                 #seqname= re.sub('\..*', '', seqname)
                 fasta_dict[seqname]= ''
            else:
                 fasta_dict[seqname] += line
        fasta.close()

        bed= open(exome_target_bed, 'U')
        fasta_bed_out = open(ref_index_name + '.refSeq.fa','w')
        for line in bed:
            if line.startswith('chr'):
                line= line.strip().split('\t')
                outname= line[0] + '_' + line[1]+ '_' + line[2]
                print('>' + outname)
                fasta_bed_out.write('>' + outname + '\n')
                if fasta_dict.get(outname, 0) == 0:
                    print("{0} is not in {1}".format(outname,FASTA))
                else:
                    fasta_bed_out.write(fasta_dict[outname] + '\n')
            else:
                continue
        bed.close()
    return ref_index_name + '.refSeq.fa'

def align_reads_bwa(bwa_dir, ref_fa_file, ref_index_name, exome_target_bed,
                    read1, read2, 
                    out_file, num_threads, logger_bwa_process, 
                    logger_bwa_errors): 
    if not os.path.isfile(read1):
        logger_bwa_errors.error('%s does not exist!', read1)
        print("Error: cannot find NGS read file!")
        return 1
    ref_fa_file_bed = filter_ref_fa_by_bed(ref_fa_file, ref_index_name, exome_target_bed,logger_bwa_process, 
                                           logger_bwa_errors)
    index_file_extensions = ['.pac', '.amb','.ann','.bwt', '.sa']
    genome_indexed = True
    for extension in index_file_extensions:
        if not os.path.isfile(ref_index_name + extension):
            genome_indexed = False
            break
    if not genome_indexed:
        bwa_index_command = '{0} index -p {1} {2}'.format(
            bwa_dir, ref_index_name, ref_fa_file_bed)
        logger_bwa_process.info(bwa_index_command)
        os.system(bwa_index_command)
        print('BWA genome index files have been built.')
    else:
        print('BWA genome index files exist.')

    bwa_align_command = '{0} mem -t {1} {2} {3} {4} > {5}'.format(
        bwa_dir, num_threads, ref_index_name, read1, read2, out_file)
    print(bwa_align_command)
    os.system(bwa_align_command)
    print('BWA alignment has been completed.')
    logger_bwa_process.info('BWA alignment has been completed.')
    return 0

def align_reads_bwa_based_all(bwa_dir, total_ref_fa_file,
                              read1, read2, 
                              out_file, num_threads, 
                              logger_bwa_process, logger_bwa_errors):
    if not os.path.isfile(read1):
        logger_bwa_errors.error('%s does not exist!', read1)
        print("Error: cannot find NGS read file!")
        return 1
    if not os.path.isfile(total_ref_fa_file):
        logger_bwa_errors.error('%s does not exist!', total_ref_fa_file)
        print("Error: cannot find reference fasta!")
        return 1
    ref_index_name = total_ref_fa_file.rstrip('.fasta')
    index_file_extensions = ['.pac', '.amb','.ann','.bwt', '.sa']
    genome_indexed = True
    for extension in index_file_extensions:
        if not os.path.isfile(ref_index_name + extension):
            genome_indexed = False
            break
    if not genome_indexed:
        bwa_index_command = '{0} index -p {1} {2}'.format(
            bwa_dir, ref_index_name, ref_fa_file_bed)
        logger_bwa_process.info(bwa_index_command)
        os.system(bwa_index_command)
        print('BWA genome index files have been built.')
    else:
        print('BWA genome index files exist.')

    bwa_align_command = '{0} mem -t {1} {2} {3} {4} > {5}'.format(
        bwa_dir, num_threads, ref_index_name, read1, read2, out_file)
    print(bwa_align_command)
    os.system(bwa_align_command)
    print('BWA alignment has been completed.')
    logger_bwa_process.info('BWA alignment has been completed.')
    return 0