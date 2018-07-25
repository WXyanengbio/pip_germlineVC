模块：variant_call
功能：模块1对预处理sam结果进行胚系突变检测，模块2对检测出的突变进行注释

编程：Python3.6
      Python模块：os 、 sys 、time 、pandas
软件：GATK4.0.5
      samtools     

模块1：g_variantcall1.py
函数：check_variant_exit
     确定VCF中是否存在突变结果，如果只有文件开头，返回0；有突变结果，返回1

函数：sam_to_bem
     对输入的sam文件进行预处理，生成MarkDuplicates.BQSR.bam，以便GATK能进行突变检测

函数：germline_variant_calling
     对输入的MarkDuplicates.BQSR.bam进行HaplotypeCaller突变检测；基于突变过滤参数进行突变过滤


模块2：annotation_gatk_hc.py
函数：read_database
     读取指定的三个突变注释参考库cosmic,clinvar,g1000

函数：get_info
     读取VCF文件中INFO项的各参数

函数：split_variant
     读取VCF文件，并排除模块g_variantcall1.py中通过参数过滤掉的突变结果

函数：annotation
     基于注释库和VCF文件进行突变注释

函数：fill_table
    对注释结果进行ensembl库匹配

函数：annotationmain
     调用其他函数，根据输入数据，实现突变注释功能


