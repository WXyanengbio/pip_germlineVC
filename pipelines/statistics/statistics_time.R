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
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))


# 从文件中读取
if (TRUE){
  dat = read.table(opts$prefix_file, sep="\t",stringsAsFactors = F)
}

# 统计与绘图
if (TRUE){
  #print(dat[,2])
  type <- apply(dat,1,function(x){unlist(strsplit(unlist(strsplit(as.character(x),"--"))[3]," "))[1]})
  uniquetype = unique(type)
  cost <- apply(dat,1,function(x){unlist(strsplit(unlist(strsplit(as.character(x),"after "))[2]," "))[1]})
  lines=c()
  for(i in 1:length(uniquetype)){
  lines=c(lines, max(which(uniquetype[i]==type)))
  }
  print(lines)
  dat1=data.frame(type=type[lines],cost=as.numeric(cost[lines]))
  dat1$type = factor(dat1$type,levels=as.character(dat1$type))
  dat1[nrow(dat1),2]=sum(dat1[-c(nrow(dat1)),2])
  dat1$ratio = paste(round(100*dat1$cost/dat1[nrow(dat1),2],2),"%",sep="")
  #colnames(dat)<-c("depth","bases","totalbases","ratio")
  p<- ggplot(data=dat1, aes(x=type, y= cost)) +
    geom_bar(colour="black", fill="#DD8888", width=.8, stat="identity") + 
    geom_text(aes(label=ratio),vjust=0)+
    guides(fill=FALSE) +
    xlab("Modules") +
    ylab("Time cost (min)") + # Set axis labels
    ggtitle("Time cost of the pipeline of germline variant calling") +
    theme_bw()+
    theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
    #scale_x_log10(breaks=c(0.05, 10, 50, 100, 200, 300, 500, 1000, 2000)) 
    
}

# 5. 保存图表
if (TRUE){
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  #保存统计结果，有waring正常
  write.table(dat1, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
  #print(paste("The output table is ", opts$output, ".txt",  sep = ""))
  #----
  # 保存图片至文件，pdf方便AI修改成出版级图片
  ggsave(file=paste(opts$output,".pdf",sep=""), plot=p, width = 6, height = 6)
  ggsave(file=paste(opts$output,".tiff",sep=""), plot=p, width = 6, height = 6)
  #print(paste("The output figure is ", opts$output, ".pdf",  sep = ""))
}