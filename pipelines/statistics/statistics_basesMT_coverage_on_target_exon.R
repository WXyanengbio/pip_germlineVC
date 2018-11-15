# 清理工作环境 clean enviroment object
rm(list=ls()) 
options(warn=-1)
# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)
suppressMessages(library('data.table'))
option_list <- list(
  make_option(c("-p", "--prefix_file"), type="character", 
              help="Input depth file by posi to read"),
  make_option(c("-s", "--suffix_file"), type="character",
              help="Input target exon to read"),
  make_option(c("-t", "--tiff"), type="character", default = FALSE,
              help="set to plot the depth of based in regions [default %default]"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))
if (TRUE){
  dat = fread(opts$prefix_file)
  region = fread(opts$suffix_file)
}
  dat1 = dat[,c('POS','DP','MT','UMT','UFR')]
  a_s=dat1[,1]
  dat1$y1<- dat1$UMT
  dat1$y2<- dat1$UMT
  dat1$exon <- ''
  for(i in 1:nrow(region)){ 
  b1=intersect(which(a_s>=as.numeric(region[i,2])), which(a_s<=as.numeric(region[i,3])))
  dat1$y1[b1]<-0
  dat1$exon[b1] <- as.character(region[i,1])
  }
  dat_exon1 = dat1[which(dat1$y1==0),]
  for(i in 1:nrow(dat_exon1)){
    if(is.na(dat_exon1[i,4])){
     dat_exon1[i,4]<-0
    }
  }
  dat_exon1$y2 = dat_exon1$UMT-dat_exon1$y1
  meanMT = mean(dat_exon1$UMT)
  p = ggplot(dat_exon1, aes(x = POS)) + 
      geom_line(aes(y = UMT), colour = 'green') +
      geom_area(aes(y = pmin(y2, UMT)), fill = 'red',alpha=0.45)+
      geom_hline(yintercept = meanMT)+
      facet_wrap(~exon,scale="free",ncol=4)+
    xlab("Position of Base") + ylab("Depth of MT")+
    theme_bw()
  ggsave(file=paste(opts$output,".pdf",sep=""), plot=p, width = 30, height = 40,limitsize = FALSE)

  minMDB = min(dat_exon1$UMT)
  maxMDB = max(dat_exon1$UMT)
  pre5MDB = round(length(which(dat_exon1$UMT >= 0.05*meanMT))/nrow(dat_exon1),3)
  pre10MDB = round(length(which(dat_exon1$UMT >= 0.1*meanMT))/nrow(dat_exon1),3)
  pre20MDB = round(length(which(dat_exon1$UMT >= 0.2*meanMT))/nrow(dat_exon1),3)
  pre30MDB = round(length(which(dat_exon1$UMT >= 0.3*meanMT))/nrow(dat_exon1),3)
  exon_statis_mtdp<-cbind(Library_name=c("Mean MT depth per targeted base (mean MDB)",
                                         "Minimun MDB",
                                         "Maximun MDB",
                                         "% of target bases with MT depth >= 5% of MDB",
                                         "% of target bases with MT depth >= 10% of MDB",
                                         "% of target bases with MT depth >= 20% of MDB",
                                         "% of target bases with MT depth >= 30% of MDB"),
                          values=c(round(meanMT,0),minMDB, maxMDB,pre5MDB,pre10MDB,pre20MDB,pre30MDB)
                          )
  write.table(exon_statis_mtdp, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)