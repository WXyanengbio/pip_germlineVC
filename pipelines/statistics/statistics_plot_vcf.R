# 清理工作环境 clean enviroment object
rm(list=ls()) 
options(warn=-1)
# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)
library(easyGgplot2)

option_list <- list(
  make_option(c("-p", "--prefix_file"), type="character", 
              help="Input table file to read")
 # make_option(c("-s", "--suffix_file"), type="character",
 #             help="Input table file to read"),
  #make_option(c("-g", "--group_name"), type="character",
  #            help="set the group of the two input files,such as :before, after"),
  #make_option(c("-o", "--output"), type="character", default="output",
  #            help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
#print(paste("The prefix file is ", opts$prefix_file, sep = ""))
#print(paste("The output file prefix is ", opts$output, sep = ""))


# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从文件中读取
if (TRUE){
  dat1 = read.table(paste("dp20_",opts$prefix_file,".vcf_statis.txt",sep=""))
  dat2 = read.table(paste("dp50_",opts$prefix_file,".vcf_statis.txt",sep=""))
  dat3 = read.table(paste("dp100_",opts$prefix_file,".vcf_statis.txt",sep=""))
 # suf_dat = read.table(opts$suffix_file, sep="\t")
  #group_name= unlist(strsplit(opts$group_name,","))
}

# 4. 统计与绘图
if (TRUE){
  #
  snp = data.frame(dp20=t(dat1[4,-1]),dp50=t(dat2[4,-1]),dp100=t(dat3[4,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=c('DP20','DP50','Dp100')
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  p1<-ggplot(data=snp_melt,aes(x=variable, y=value))+
     geom_line(aes(group=sample,colour=sample),size=1.25)+
      xlab("Depth")+ylab("Variants")+theme_bw()

   ggsave(file=paste(opts$prefix_file,"_snp.tiff",sep=""), plot=p1, width = 6.5, height = 5)

  snp = data.frame(dp20=t(dat1[7,-1]),dp50=t(dat2[7,-1]),dp100=t(dat3[7,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=c('DP20','DP50','Dp100')
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  p1<-ggplot(data=snp_melt,aes(x=variable, y=value))+
     geom_line(aes(group=sample,colour=sample),size=1.25)+
      xlab("Depth")+ylab("Precision")+theme_bw()

   ggsave(file=paste(opts$prefix_file,"_snp_pre.tiff",sep=""), plot=p1, width = 6.5, height = 5)
  
	snp = data.frame(dp20=t(dat1[8,-1]),dp50=t(dat2[8,-1]),dp100=t(dat3[8,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=c('DP20','DP50','Dp100')
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  p1<-ggplot(data=snp_melt,aes(x=variable, y=value))+
     geom_line(aes(group=sample,colour=sample),size=1.25)+
      xlab("Depth")+ylab("Variants")+theme_bw()

   ggsave(file=paste(opts$prefix_file,"_indel.tiff",sep=""), plot=p1, width = 6.5, height = 5)

  snp = data.frame(dp20=t(dat1[11,-1]),dp50=t(dat2[11,-1]),dp100=t(dat3[11,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=c('DP20','DP50','Dp100')
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  p1<-ggplot(data=snp_melt,aes(x=variable, y=value))+
     geom_line(aes(group=sample,colour=sample),size=1.25)+
      xlab("Depth")+ylab("Precision")+theme_bw()

   ggsave(file=paste(opts$prefix_file,"_indel_pre.tiff",sep=""), plot=p1, width = 6.5, height = 5)
    
}

