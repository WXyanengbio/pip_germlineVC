# 清理工作环境 clean enviroment object
rm(list=ls()) 
options(warn=-1)
# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)
library(easyGgplot2)
library(VennDiagram)

option_list <- list(
  make_option(c("-p", "--prefix_file"), type="character", 
              help="Input table file to read"),
  make_option(c("-s", "--suffix_file"), type="character",
              help="Input table file to read"),
  #make_option(c("-g", "--group_name"), type="character",
  #            help="set the group of the two input files,such as :before, after"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
#print(paste("The prefix file is ", opts$prefix_file, sep = ""))
#print(paste("The output file prefix is ", opts$output, sep = ""))


# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从文件中读取
if (TRUE){
  groups=unlist(strsplit(opts$suffix_file,','))
  dat1 = read.table(paste(opts$prefix_file,'.',groups[1],"_PASS.txt",sep=""),sep="\t",header=T)
  print(head(dat1))
  dat2 = read.table(paste(opts$prefix_file,'.',groups[2],"_PASS.txt",sep=""),sep="\t",header=T)

  }
 # suf_dat = read.table(opts$suffix_file, sep="\t")
  #group_name= unlist(strsplit(opts$group_name,","))

# 4. 统计与绘图
if (TRUE){
 dat1$type=apply(dat1,1,function(x){ifelse(nchar(as.character(x[4]))==1 & nchar(as.character(x[5]))==1,'SNP','Indel')})
 snp1= dat1[intersect(which(dat1[,'Gene_Name']!='-'), which(dat1[,'type']=='SNP')),3]
 print(unlist(snp1))
 indel1= dat1[intersect(which(dat1[,'Gene_Name']!='-'), which(dat1[,'type']=='Indel')),3]
 dat2$type=apply(dat2,1,function(x){ifelse(nchar(as.character(x[4]))==1 & nchar(as.character(x[5]))==1,'SNP','Indel')})
 snp2= dat2[intersect(which(dat2[,'Gene_Name']!='-'), which(dat2[,'type']=='SNP')),3]
 print(unlist(snp2))
 indel2= dat2[intersect(which(dat2[,'Gene_Name']!='-'), which(dat2[,'type']=='Indel')),3]	
venn.diagram(list(A=unlist(snp1),B=unlist(snp2)), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename=paste(opts$output, as.character(dat1[1,1]),"_snp_venn.tiff",sep=""))
venn.diagram(list(A=unlist(indel1),B=unlist(indel2)), fill=c("red","green"), alpha=c(0.5,0.5), cex=2, cat.fontface=4, fontfamily=3, filename=paste(opts$output, as.character(dat1[1,1]),"_indel_venn.tiff",sep=""))
}

