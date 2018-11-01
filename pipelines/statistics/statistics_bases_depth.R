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
 # make_option(c("-s", "--suffix_file"), type="character",
 #             help="Input table file to read"),
  #make_option(c("-g", "--group_name"), type="character",
  #            help="set the group of the two input files,such as :before, after"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))



# 从文件中读取
if (TRUE){
  dat = read.table(opts$prefix_file, sep="\t",stringsAsFactors = F)
 # suf_dat = read.table(opts$suffix_file, sep="\t")
  #group_name= unlist(strsplit(opts$group_name,","))
}
#print(pre_dat)
# 弹出窗口选择文件
if (FALSE){
  dat = read.table(file.choose(), header=T, row.names = NULL, sep="\t")
}

# 4. 统计与绘图
if (TRUE){
  #print(dat[,2])
  # remove the zero-X depth bases
  #dat = dat [-1,]
  #dat[1,1] = 0.001
  # get the 
  dat[,3]= cumsum(dat[,2])
  dat[,4]= 1-dat[,3]/sum(dat[,2])
  colnames(dat)<-c("depth","bases","totalbases","ratio")
  #--
  limits<-data.frame(x=c(50,100,200),colour=c("50X","100X","200X"),y=c(0.9,0.9,0.9))
  p<- ggplot(data=dat, aes(x=depth, y= ratio)) +
    #geom_smooth(colour="grey70", se=F,
    #            method="glm",
    #            formula=y~ns(x,8),
    #            family=gaussian(link="log"),
    #            show.legend = FALSE,lwd=0.7) +
    #scale_x_log10(breaks = c(0.001,10,50,100,200,1000,round(max(dat[,1])/1000,0)*1000),label= c(0,10, 50,100,200,1000,round(max(dat[,1])/1000,0)*1000))+ 
    ylim(0,1)+
    xlim(0,500)+
    geom_line(colour="black")+
    #geom_point(size= 0.3,alpha=0.3)+
    #geom_vline(xintercept=10)+
    geom_vline(data=limits,aes(xintercept=x,colour=colour),show.legend = FALSE, size = 0.75)+
    geom_text(data=limits,aes(x=x,y=y,label=colour, colour=colour),show.legend = FALSE, size = 4)+
    #geom_vline(xintercept=100)+
    #geom_vline(xintercept=200)+
    xlab("Depth of Bases") + ylab("Ratio of total Bases") + # Set axis labels
    ggtitle("Cumulative frequency plot of sequencing depths") +
    #scale_colour_hue("Depth",breaks=c("50X","100X","200X"))+
    theme_bw()
    #theme(legend.key.size=unit(1,'cm'))
    #scale_x_log10(breaks=c(0.05, 10, 50, 100, 200, 300, 500, 1000, 2000)) 
    
}

# 5. 保存图表
if (TRUE){
  ggsave(file=paste(opts$output,".pdf",sep=""), plot=p, width = 5, height = 5)
  ggsave(file=paste(opts$output,".tiff",sep=""), plot=p, width = 5, height = 5)
  #print(paste("The output figure is ", opts$output, ".pdf",  sep = ""))
}