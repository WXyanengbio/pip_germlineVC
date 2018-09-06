from __future__ import barry_as_FLUFL

__all__  =  ['log_dir']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
#import re
import sys
import time
#import gzip
#import itertools

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

def store_trim_logs(log_dir):
    formatter_trim_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_trim_errors = logging.Formatter("%(asctime)s;%(levelname)s;%(message)s")
    logger_trim_process = setup_logger('Running Messages of read trimming',
                                       log_dir + '/trim_process.log',
                                       formatter_trim_process)
    logger_trim_errors = setup_logger('Errors & Warnings of read trimming',
                                      log_dir + '/trim_errors.log',
                                      formatter_trim_errors)
    return logger_trim_process, logger_trim_errors

def store_filter_logs(log_dir):
    formatter_filter_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_filter_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_filter_process = setup_logger('Post-alignment Filtration Messages',
                                         log_dir + '/filter_process.log',
                                         formatter_filter_process)
    logger_filter_errors = setup_logger('Errors & Warnings of filtration',
                                        log_dir + '/filter_errors.log',
                                        formatter_filter_errors)
    return logger_filter_process, logger_filter_errors


def store_align_logs(log_dir):
    formatter_bwa_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_bwa_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_bwa_process = setup_logger('BWA Running Messages', 
                                      log_dir + '/bwa_process.log', 
                                      formatter_bwa_process)
    logger_bwa_errors = setup_logger('Errors & Warnings of BWA', 
                                    log_dir + '/bwa_errors.log',
                                    formatter_bwa_errors)
    return logger_bwa_process, logger_bwa_errors


def store_cluster_logs(log_dir):
    formatter_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_umi_process = setup_logger('Running Barcode Clustering Messages', 
                                  log_dir + '/umi_tagging_process.log',
                                  formatter_process)
    logger_umi_errors = setup_logger('Errors & Warnings of Barcode Clustering', 
                                 log_dir + '/umi_tagging_errors.log',
                                 formatter_errors)
    return logger_umi_process, logger_umi_errors

def store_reformat_logs(log_dir):
    formatter_reformat_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_reformat_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_reformat_process = setup_logger('Running Messages of reformating sam file',
                                           log_dir + '/reformat_process.log',
                                           formatter_reformat_process)
    logger_reformat_errors = setup_logger('Errors & Warnings of reformating sam file',
                                          log_dir + '/reformat_errors.log',
                                          formatter_reformat_errors)
    return logger_reformat_process, logger_reformat_errors


def store_germline_vc_logs(log_dir):
    formatter_germline_vc_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_germline_vc_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_germline_vc_process = setup_logger('Germline variant calling Messages',
                                         log_dir + '/germline_vc_process.log',
                                         formatter_germline_vc_process)
    logger_germline_vc_errors = setup_logger('Errors & Warnings of Germline variant calling',
                                        log_dir + '/germline_vc_errors.log',
                                        formatter_germline_vc_errors)
    return logger_germline_vc_process, logger_germline_vc_errors

def store_annotation_logs(log_dir):
    formatter_annotation_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_annotation_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_annotation_dir_process = setup_logger('Annotation variant Messages',
                                         log_dir + '/annotation_process.log',
                                         formatter_annotation_process)
    logger_annotation_dir_errors = setup_logger('Errors & Warnings of Annotation variant',
                                        log_dir + '/annotation_errors.log',
                                        formatter_annotation_errors)
    return logger_annotation_dir_process, logger_annotation_dir_errors


def store_statistics_logs(log_dir):
    formatter_pipeline_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_pipeline_errors = logging.Formatter("%(asctime)s;%(levelname)s; %(message)s")
    logger_statistics_process = setup_logger('Running Messages of statistics',
                                       log_dir + '/statistics_process.log',
                                       formatter_pipeline_process)
    logger_statistics_errors = setup_logger('Errors & Warnings of pipeline',
                                      log_dir + '/statistics_errors.log',
                                      formatter_pipeline_errors)
    return logger_statistics_process, logger_statistics_errors


def store_benchmark_logs(log_dir):
    formatter_pipeline_process = logging.Formatter("%(asctime)s;%(message)s")
    formatter_pipeline_errors = logging.Formatter(
        "%(asctime)s;%(levelname)s; %(message)s")
    logger_benchmark_process = setup_logger('Running Messages of Benchmarking',
                                       log_dir + '/benchmark_process.log',
                                       formatter_pipeline_process)
    logger_benchmark_errors = setup_logger('Errors & Warnings of pipeline',
                                      log_dir + '/benchmark__errors.log',
                                      formatter_pipeline_errors)
    return logger_benchmark_process, logger_benchmark_errors