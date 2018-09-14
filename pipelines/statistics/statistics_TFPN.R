# 清理工作环境 clean enviroment object
rm(list=ls()) 
options(warn=-1)
# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)

option_list <- list(
  make_option(c("-p", "--prefix_file"), type="character", 
              help="Input table file to read"),
  make_option(c("-r", "--defined_vcf"), type="character",  default="/home/dell/Works/Projects/Datasets/celllines/celllines_vcf.txt",
              help="Input table file to read [default %default]"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))
defined_vcf<-read.table(opts$defined_vcf,header=T)
files<-read.table(opts$prefix_file)
defined_vcf_celllines= colnames(defined_vcf)[-1]
defined_vcf_celllines=gsub('\\.','-',defined_vcf_celllines)
defined_vcf_celllines=gsub('^X','',defined_vcf_celllines)

 if(TRUE){
 vcfs=c()
 defined_vcfs=c()
 names=c()
 for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 cellline = as.character(files[i,2])

 dat = read.table(file,sep="\t",header=T)

# remove the total reads, error reads pair fq, cost time
 raw=nrow(dat)
 dat$type=apply(dat,1,function(x){ifelse(nchar(as.character(x[4]))==1 & nchar(as.character(x[5]))==1,'SNP','Indel')})
 cellline_num = which(cellline==defined_vcf_celllines)
 defined_vcfsub=as.data.frame(defined_vcf[,c(1, cellline_num+1)])
 defined_vcfsub$defined_vcf1=apply(defined_vcfsub,1,function(x){ifelse(length(which(x[1]==dat[,3]))==1,'y','n')})
 raw_cell =length(which('y'==defined_vcfsub[,2]))
 call_cell = length(intersect(which('y'==defined_vcfsub[,2]),which('y'==defined_vcfsub[,3])))
 call_wrong = length(intersect(which('n'==defined_vcfsub[,2]),which('y'==defined_vcfsub[,3])))
 defined_vcf2=c(raw_cell,call_cell,raw_cell-call_cell,call_wrong)
   
 tp = length(which(dat[,'Gene_Name']!='-'))
 fp = raw -tp
 snp = length(which(dat[,'type']=='SNP'))
 indel = raw -snp
 snp_tp = length(intersect(which(dat[,'Gene_Name']!='-'), which(dat[,'type']=='SNP')))
 snp_fp = snp -snp_tp
 snp_pre= ifelse(snp!=0,round(snp_tp/snp,4)*100,NA)
 indel_tp = length(intersect(which(dat[,'Gene_Name']!='-'), which(dat[,'type']=='Indel')))
 indel_fp = indel -indel_tp
indel_pre=ifelse(indel!=0,round(indel_tp/indel,4)*100,NA)
 vcfs<-cbind(vcfs,c(raw,tp,fp,snp,snp_tp,snp_fp,snp_pre,indel,indel_tp,indel_fp,indel_pre))
 defined_vcfs<-cbind(defined_vcfs,defined_vcf2)
 name=as.character(dat[1,1])
 names<-c(names,gsub('_L001','',name))
 }
 vcfs<-as.data.frame(vcfs)
 colnames(vcfs)<-names
colnames(defined_vcfs)<-names
 names1<-c()
  vcfs=cbind(read_set=c("All","TP","FP","SNP","SNP_TP","SNP_FP","SNP_PRE","INDEL","INDEL_TP","INDEL_FP","INDEL_PRE"),vcfs)
  defined_vcfs=cbind(type=c("Targeted in region","Calling","no calling","wrong calling"),defined_vcfs)
 }
 

# 5. 保存图表
if (TRUE){
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  #保存统计结果，有waring正常
  write.table(vcfs, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", 
         na = "NA", dec = ".", row.names = F, col.names = T)
  #保存统计结果，有waring正常
  write.table(t(defined_vcfs), file=paste(opts$output,"defined.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", 
         na = "NA", dec = ".", row.names = T, col.names = F)
}