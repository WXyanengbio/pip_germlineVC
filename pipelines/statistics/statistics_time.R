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
  #make_option(c("-g", "--group_name"), type="character",
  #            help="set the group of the two input files,such as :before, after"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The prefix file is ", opts$prefix_file, sep = ""))
print(paste("The output file prefix is ", opts$output, sep = ""))


# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

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
  type <- apply(dat,1,function(x){unlist(strsplit(unlist(strsplit(as.character(x),";"))[2]," "))[1]})
  cost <- apply(dat,1,function(x){unlist(strsplit(unlist(strsplit(as.character(x),"after "))[2]," "))[1]})
  dat1=data.frame(type,cost=as.numeric(cost))
  dat1$type = factor(dat1$type,levels=as.character(dat1$type))
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