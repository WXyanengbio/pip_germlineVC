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
              help="Input table file to read"),
  make_option(c("-s", "--suffix_file"), type="character",
              help="Input table file to read")
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
  groups=unlist(strsplit(opts$suffix_file,','))
  dat1 = read.table(paste(opts$prefix_file,'_',groups[1],".vcf_statis.txt",sep=""))
  dat2 = read.table(paste(opts$prefix_file,'_',groups[2],".vcf_statis.txt",sep=""))
  for(i in 1:nrow(dat1)){
  for(j in 2:ncol(dat1)){
   if(is.na(dat1[i,j])){
     dat1[i,j]=1
   }
  }
  }
  for(i in 1:nrow(dat2)){
  for(j in 2:ncol(dat2)){
   if(is.na(dat2[i,j])){
     dat2[i,j]=1
   }
  }
  }
 # suf_dat = read.table(opts$suffix_file, sep="\t")
  #group_name= unlist(strsplit(opts$group_name,","))
}

# 4. 统计与绘图
if (TRUE){
  #
  snp = data.frame(dp20=t(dat1[5,-1]),dp50=t(dat2[5,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=groups
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  print(snp_melt)
  p1<-ggplot(data=snp_melt,aes(sample, weight=value))+
     geom_bar(aes(group=variable,colour=variable,fill=variable),position="dodge")+
      xlab("Sample")+ylab("Variants")+theme_bw()+  theme(axis.text.x=element_text(angle=90,size=14,vjust=0.5))

   ggsave(file=paste(opts$prefix_file,'_',groups[1],'_',groups[2],"_snp.tiff",sep=""), plot=p1, width = 5, height = 4)

  snp = data.frame(dp20=t(dat1[7,-1]),dp50=t(dat2[7,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=groups
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  p1<-ggplot(data=snp_melt,aes(x=variable, y=value))+
     geom_line(aes(group=sample,colour=sample),size=1.25)+
      xlab("Softwares")+ylab("Precision")+theme_bw()

   ggsave(file=paste(opts$prefix_file,'_',groups[1],'_',groups[2],"_snp_pre.tiff",sep=""), plot=p1, width = 5, height = 4)
  
	snp = data.frame(dp20=t(dat1[9,-1]),dp50=t(dat2[9,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=groups
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  print(snp_melt)
  p1<-ggplot(data=snp_melt,aes(sample, weight=value))+
     geom_bar(aes(group=variable,colour=variable,fill=variable),position="dodge")+
      xlab("Sample")+ylab("Variants")+theme_bw()+  theme(axis.text.x=element_text(angle=90,size=14,vjust=0.5))

   ggsave(file=paste(opts$prefix_file,'_',groups[1],'_',groups[2],"_indel.tiff",sep=""), plot=p1, width = 5, height = 4)

  snp = data.frame(dp20=t(dat1[11,-1]),dp50=t(dat2[11,-1]))
  rownames(snp)=gsub('^X','',colnames(dat1)[-1])
  colnames(snp)=groups
  snp$sample=rownames(snp)
  snp_melt=melt(snp)
  p1<-ggplot(data=snp_melt,aes(x=variable, y=value))+
     geom_line(aes(group=sample,colour=sample),size=1.25)+
      xlab("Softwares")+ylab("Precision")+theme_bw()

   ggsave(file=paste(opts$prefix_file,'_',groups[1],'_',groups[2],"_indel_pre.tiff",sep=""), plot=p1, width = 5, height = 4)
    
}

