from __future__ import barry_as_FLUFL

__all__  =  ['logger_process','logger_errors','message']
__version__  =  '1.0'
__author__  =  'Wang Xian'

import logging
import os
#import re
import sys
import time
#import gzip
#import itertools

def get_time_stamp():
    ct = time.time()
    local_time = time.localtime(ct)
    data_head = time.strftime("%Y-%m-%d %H:%M:%S", local_time)
    data_secs = (ct - int(ct)) * 1000
    time_stamp = "%s.%03d" % (data_head, data_secs) + ';'
    return time_stamp

def store_pipeline_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/process.log','a+')
        file.write(get_time_stamp()+message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/error.log','a+')
        file.write(get_time_stamp()+message)
    file.close()
	
def store_trim_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/trim_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/trim_errors.log','a+')
        file.write(message)
    file.close()

def store_filter_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/filter_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/filter_errors.log','a+')
        file.write(message)
    file.close()

def store_align_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/bwa_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/bwa_errors.log','a+')
        file.write(message)
    file.close()

def store_cluster_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/umi_tagging_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/umi_tagging_errors.log','a+')
        file.write(message)
    file.close()
	
def store_reformat_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/reformat_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/reformat_errors.log','a+')
        file.write(message)
    file.close()

def store_germline_vc_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/germline_vc_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/germline_vc_errors.log','a+')
        file.write(message)
    file.close()
	
def store_annotation_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/annotation_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/annotation_errors.log','a+')
        file.write(message)
    file.close()

def store_statistics_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/statistics_process.log','a+')
        file.write(get_time_stamp()+message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/statistics_errors.log','a+')
        file.write(get_time_stamp()+message)
    file.close()

def store_benchmark_logs(logger_process,logger_errors,message):
    if logger_process != 'null':
        file=open(logger_process+ '/benchmark_process.log','a+')
        file.write(message)
    elif logger_errors != 'null':
        file=open(logger_errors+ '/benchmark__errors.log','a+')
        file.write(message)
    file.close()