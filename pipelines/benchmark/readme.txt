模块：hap_benchmark.py
功能：对胚系突变检测方法进行评价，分析其召回率和精确度
输入：benchmark_dir：hap.py软件安装路径
　　　truth_vcf    ：真实突变，参考库文件
　　　filter_vcf   ：胚系突变检测结果，待评价方法分析出的突变结果
　　　confident_region_bed： 真实突变的可靠区域，参考库文件
　　　benchmarking_dir：benchmarking结果路径 
　　　total_ref_fa_file ：与真实突变版本相符的参考基因组
　　　exon_interval：胚系突变检测的目标区域
　　　num_threads：使用hap.py软件线程数
　　　logger_benchmark_process：log文件
　　　logger_benchmark_errors：log文件

编程：Python3.6
      Python模块：os 、 sys

函数：hap_py--实现模块功能
