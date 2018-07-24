# 清理工作环境 clean enviroment object
rm(list=ls()) 

# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)

option_list <- list(
  make_option(c("-p", "--prefix_file"), type="character", 
              help="Input table file to read"),
 # make_option(c("-s", "--suffix_file"), type="character",
 #             help="Input table file to read"),
  make_option(c("-g", "--group_name"), type="character",
              help="set the group names of the sam_bam static files"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The prefix file is: ", unlist(strsplit(opts$prefix_file,",")), sep = " "))
print(paste("The output file prefix is ", opts$output, sep = ""))


# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从文件中读取
if (TRUE){
  dat = c()
  for(i in unlist(strsplit(opts$prefix_file,","))){
  dat1 = read.table(i, sep="\t",stringsAsFactors = F)
  dat <- cbind(dat,dat1[,2])
  }
  rownames(dat)<-dat1[,1]
  colnames(dat)<-unlist(strsplit(opts$group_name,","))
 # suf_dat = read.table(opts$suffix_file, sep="\t")
  #group_name= unlist(strsplit(opts$group_name,","))
}
print(head(dat))

# 4. 统计与绘图
if (TRUE){
  description <- c("总序列条数",
"过滤条数",
"剩余序列条数",
"是否过滤",
"read1序列条数",
"read2序列条数",
"比对上基因组的 read 数目",
"paired reads中两条都比对到参考序列上的reads数目",
"没有比对上基因组的 read 数目",
"正确配对的reads数目",
"配对的reads数目",
"PCR or optical duplicate",
"MQ0 的read数目",
"QC失败的read数目",
"非主要比对",
"总碱基数",
"比对上的碱基数",
"cigar比对上的碱基数",
"trimmed的碱基数",
"重复的碱基数",
"错配的碱基数",
"错配率",
"序列平均长度",
"最大长度",
"序列平均质量",
"插入片段平均长度",
"插入片段标准差",
"两端读取相向的成对reads数目",
"两端读取相反的成对reads数目",
"两端读取为其他方向的成对reads数目",
"paired reads中两条分别比对到两条不同的参考序列的reads数目"
)
dat <- data.frame(description, dat)
}

# 5. 保存图表
if (TRUE){
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  #保存统计结果，有waring正常
  write.table(dat, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
  print(paste("The output table is ", opts$output, ".txt",  sep = ""))
  #----
  # 保存图片至文件，pdf方便AI修改成出版级图片
  if (FALSE){
  ggsave(file=paste(opts$output,".pdf",sep=""), plot=p, width = 5, height = 5)
  ggsave(file=paste(opts$output,".tiff",sep=""), plot=p, width = 5, height = 5)
  print(paste("The output figure is ", opts$output, ".pdf",  sep = ""))
  }
}