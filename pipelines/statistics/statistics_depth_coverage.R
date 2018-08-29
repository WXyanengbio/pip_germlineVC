# 清理工作环境 clean enviroment object
rm(list=ls()) 
options(warn=-1)
# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(easyGgplot2)


option_list <- list(
  make_option(c("-p", "--prefix_file"), type="character", 
              help="Input table file to read"),
  make_option(c("-s", "--suffix_file"), type="character",
              help="Input table file to read"),
  make_option(c("-r", "--region"), type="character",default="null",
             help="the target region [default %default] "),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
#print(paste("The prefix file is ", opts$prefix_file, sep = ""))
#print(paste("The suffix file is ", opts$suffix_file, sep = ""))
#print(paste("The output file prefix is ", opts$output, sep = ""))


# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从文件中读取
if (TRUE){
  pre_dat = read.table(opts$prefix_file, sep="\t")
  suf_dat = read.table(opts$suffix_file, sep="\t")
  region = read.table(opts$region)
  if(substr(region[1,1],1,3) != 'chr'){
    region = region[-1,]
  }
  region$region= apply(region,1,function(x){as.numeric(as.character(x[3]))-as.numeric(as.character(x[2]))})
}

# 4. 统计与绘图
if (TRUE){
  coverage_depth = c()
  for (i in c(1:nrow(pre_dat))){
      exit = match(as.character(pre_dat[i,1]),as.character(suf_dat[,1]))
      if(!is.na(exit)){
        #print(suf_dat[exit,2:3])
        coverage_depth <- rbind(coverage_depth, unlist(suf_dat[exit,2:3]))
      }
      else{
        coverage_depth <- rbind(coverage_depth, c(0,0))
      }
    }
  dat <-cbind(pre_dat,coverage_depth)
  #print(dat)
  colnames(dat)<-c("targe region","targe length","mapped reads","non-mapped reads","coverage position length","sum of depth")
  dat$coverage_ratio<- round(dat[,5]/dat[,2],4)
  dat$mean_depth<-round(dat[,6]/dat[,2],4)
  #--
  #--merge the targe regions for 
  if(nchar(as.character(pre_dat[1,1]))>3){
    print("Merge the targe regions!!!")
    chom<-apply(dat,1,function(x){unlist(strsplit(as.character(x[1]),"_"))[1]})
    unique_chom<-unique(chom)
    merge_dat<-c()
    for(i in unique_chom){
       if(i != "*"){
         row =grep(paste("^",i,"$",sep=""),chom,perl = TRUE)
       }
      else{
        row =grep(paste('\\',i,sep=""),chom,perl = TRUE)
      }
      merge_dat<-rbind(merge_dat,colSums(dat[row,c(2,3,4,5,6)]))
    }
    #print(merge_dat)
    #print(class(merge_dat))
    merge_dat<-as.data.frame(merge_dat)
    rownames(merge_dat)<-c(unique_chom[-length(unique_chom)],"other_regions")
    merge_dat<-merge_dat[-nrow(merge_dat),]
    }
    if(nchar(as.character(pre_dat[1,1]))==4){
    chom<-as.character(region[,1])
    unique_chom<-unique(chom)
    merge_region<-c()
    for(i in unique_chom){
       if(i != "*"){
         row =grep(paste("^",i,"$",sep=""),chom,perl = TRUE)
       }
      else{
        row =grep(paste('\\',i,sep=""),chom,perl = TRUE)
      }
      merge_region<-c(merge_region,sum(region[row,4]))
    }
    merge_region<-as.data.frame(merge_region)
    rownames(merge_region)<- unique_chom
    #print(merge_region)
    merge_dat[,1]<- merge_region[rownames(merge_dat),1]
    #print(merge_dat[,1])
    }
    #print(merge_dat)
    merge_dat = na.omit(merge_dat)
    #print(merge_dat)
    merge_dat$coverage_ratio<- round(merge_dat[,4]/merge_dat[,1],4)
    merge_dat$mean_depth<-round(merge_dat[,5]/merge_dat[,1],4)
    #print(merge_dat)
     
    #--
    merge_dat_sub = data.frame(CHOM = rownames(merge_dat),target=merge_dat[,6],mean_depth=merge_dat[,7])
    #merge_dat_sub = data.frame(CHOM = c("chr1","chr2","chr3","chr4","chr5","chr6",
    #                                    "chr7","chr8","chr9","chr10","chr11","chr12",
    #                                    "chr13","chr14","chr15","chr16","chr17","chr18",
    #                                    "chr19","chr20","chr22","chrX","chrY"),
    #                          target=merge_dat[c("chr1","chr2","chr3","chr4","chr5","chr6",
    #                                    "chr7","chr8","chr9","chr10","chr11","chr12",
    #                                    "chr13","chr14","chr15","chr16","chr17","chr18",
    #                                    "chr19","chr20","chr22","chrX","chrY"),6],
    #                          mean_depth=merge_dat[c("chr1","chr2","chr3","chr4","chr5","chr6",
    #                                    "chr7","chr8","chr9","chr10","chr11","chr12",
    #                                    "chr13","chr14","chr15","chr16","chr17","chr18",
    #                                    "chr19","chr20","chr22","chrX","chrY"),7])
    merge_dat_sub = na.omit(merge_dat_sub)
    #print(head(merge_dat_sub))
    merge_dat_sub1 = merge_dat_sub[which(merge_dat_sub$mean_depth>40),]
    merge_dat_sub1$CHOM = factor(as.character(merge_dat_sub1$CHOM),levels =as.character(merge_dat_sub1$CHOM))
    p <- ggplot(data=merge_dat_sub1, aes(x=CHOM, y= target)) +
    geom_bar(colour="black", fill="#DD8888", width=.8, stat="identity") + 
    geom_text(aes(label=mean_depth),vjust=0)+
    guides(fill=FALSE) +
    xlab("chromosome") +
    ylab("Target ratio") + # Set axis labels
    ggtitle("Target ratio and mean depth of target regions") +
    theme_bw()+
    theme(axis.text.x = element_text(angle=0, hjust=1, vjust=1))
    #---plot the target ratio
    rownames(dat)<-dat[,1]
    dat<-dat[,-1]
  #---plot the depth and the coverage
  #if (nrow(dat)>24){
  #dat<- dat[order(dat[,5],decreasing=T),]
 # dat_sub<-dat[which(dat[,5]>=40),]
  #}
}

# 5. 保存图表
if (TRUE){
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(dat, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
  #print(paste("The output table is ", opts$output, ".txt",  sep = ""))
  #----
  if(nchar(as.character(pre_dat[1,1]))>3){
  write.table("\t", file=paste(opts$output,".mergeCHOM.txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(merge_dat, file=paste(opts$output,".mergeCHOM.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
  #print(paste("The output table is ", opts$output, ".mergeCHOM.txt",  sep = ""))
  
  # 保存图片至文件，pdf方便AI修改成出版级图片
  ggsave(file=paste(opts$output,".mergeCHOM.pdf",sep=""), plot=p, width = 5, height = 5)
  ggsave(file=paste(opts$output,".mergeCHOM.jpg",sep=""), plot=p, width = 5, height = 5)
  #print(paste("The output figure is ", opts$output, ".pdf",  sep = ""))
 }
}