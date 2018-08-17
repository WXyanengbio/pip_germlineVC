# 清理工作环境 clean enviroment object
rm(list=ls()) 

# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)
library(data.table)

option_list <- list(
  make_option(c("-p", "--prefix_file"), type="character", 
              help="Input depth file by posi to read"),
  make_option(c("-g", "--group"), type="character",
              help="the statis model"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The prefix file is ", opts$prefix_file, sep = ""))
print(paste("The output file prefix is ", opts$output, sep = ""))


  #dat = read.csv(opts$prefix_file, sep="\t",header = F)
 files<-read.table(opts$prefix_file)
 print(files)
 if(opts$group=='umi'){
 umis=c()
 names=c()
 for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 dat = fread(file)
 print(head(dat))
 reads= sum(dat[,7])
 sums=nrow(dat)
 mean=reads/sums
 quant=quantile(unlist(dat[,7]), probs = c(0,0.25,0.5,0.75,1))
 umis <- rbind(umis,c(reads,sums,mean,quant))
 name=unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/")))]
 names<-c(names,gsub('_deduplicated_per_umi.tsv','',name))
 }
 umis=t(umis)
 colnames(umis)<-names
 
  umis=cbind(read_set=c("read fragments","MTs","read fragments per MT, mean","read fragments per MT, 0th percentile",
               "read fragments per MT, 25th percentile","read fragments per MT, 50th percentile",
              "read fragments per MT, 75th percentile","read fragments per MT, 100th percentile"),umis)
  write.table(umis, file=paste(opts$output,".umis_statis.txt",sep=""),
              append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
 }

 if(opts$group=='trim'){
 trims=c()
 names=c()
 for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 dat = read.table(file,sep="\t")
 dat = dat[-nrow(dat),]
 trim=c()
 for(j in 1:length(dat)){
 trim<-rbind(trim,unlist(strsplit(as.character(dat[j])," == ")))
 }
 trims<-cbind(trims,trim[,2])
 name=unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/")))]
 names<-c(names,gsub('_L001_basic_stats.txt','',name))
 }
 colnames(trims)<-names
 names1<-c()
  trims=cbind(read_set=trim[,1],trims)
  write.table(trims, file=paste(opts$output,".trims_statis.txt",sep=""),
              append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
 }

  if(opts$group=='filter'){
 filters=c()
 names=c()
for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 dat = read.table(file,sep="\t")
 #dat = dat[-c(4,8,13),]
 filter=c()
 for(j in 1:length(dat)){
 if(j>=10){
   subfilter=unlist(strsplit(as.character(dat[j])," == "))
   number=unlist(strsplit(subfilter[2]," "))[1]
   precent=substr(unlist(strsplit(subfilter[2]," "))[2],2,6)
  filter<-rbind(filter,c(subfilter[1],number))
  filter<-rbind(filter,c(gsub("Precent","Number",subfilter[1]),precent))
 }else if(j==7){
   subfilter=unlist(strsplit(as.character(dat[j])," == "))
   number=substr(subfilter[2],1,5)
  filter<-rbind(filter,c(subfilter[1],number))
 }else{
 filter<-rbind(filter,unlist(strsplit(as.character(dat[j])," == ")))
 }
 }
 filters<-cbind(filters,filter[,2])
 name=unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/")))]
 names<-c(names,gsub('_L001_basic_stats.txt','',name))
 }
 colnames(filters)<-names
 filters=cbind(read_set=filter[,1],filters)
  write.table(filters, file=paste(opts$output,".filters_statis.txt",sep=""),
              append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
 }
