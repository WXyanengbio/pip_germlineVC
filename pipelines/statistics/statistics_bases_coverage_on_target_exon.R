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
  make_option(c("-s", "--suffix_file"), type="character",
              help="Input target exon to read"),
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
  #dat = read.csv(opts$prefix_file, sep="\t",header = F)
  dat = fread(opts$prefix_file)
  suf_dat = read.csv(opts$suffix_file,header=F)
  #group_name= unlist(strsplit(opts$group_name,","))
}
print(head(suf_dat))
print(head(suf_dat[1,2]))
# 4. 统计与绘图
if (TRUE){
  region = data.frame(region=unique(dat[,1]))
  region$start=as.numeric(apply(region,1,function(x){unlist(strsplit(as.character(x[1]),"_"))[2]}))
  region$end=as.numeric(apply(region,1,function(x){unlist(strsplit(as.character(x[1]),"_"))[3]}))

  exon=c()
  for(i in 1:nrow(region)){
     a=intersect(which(region[i,3]>=suf_dat[,3] & suf_dat[,3]> region[i,2]), which(suf_dat[,2] >= region[i,2] & suf_dat[,2] <= region[i,3]))
     b=intersect(which(region[i,3]> suf_dat[,3] & suf_dat[,3]> region[i,2]), which(suf_dat[,2] < region[i,3]))
     c=intersect(which( suf_dat[,3]> region[i,3]), which(suf_dat[,2] > region[i,2] & suf_dat[,2] < region[i,3]))
     #print(a)
     suba=c()
     if(length(a)>0){
     for(j in 1:length(a)){
     #print(a)
     #print(suf_dat[a[j],])
     start=suf_dat[a[j],2]-region[i,2]
     end = suf_dat[a[j],3]-region[i,2]
     target = as.character(region[,1])[i]
     suba<-rbind(suba,c(start,end,target,as.character(suf_dat[a[j],1])))
    }
    exon<-rbind(exon,suba)
    }
    subb=c()
    if(length(a)==0 & length(b)>0){
    #print(b)
    for(j1 in 1:length(b)){
    start=1
    end = suf_dat[b[j1],3]-region[i,2]
    target = as.character(region[,1])[i]
    subb<-rbind(subb,c(start,end,target,as.character(suf_dat[b[j1],1])))
    }
    #print(subb)
    exon<-rbind(exon,subb)
    }
    subc=c()
    if(length(a)==0 & length(c)>0){
    #print(c)
    for(j2 in 1:length(c)){
    #print(suf_dat[a[j],])
    start =suf_dat[c[j2],2]-region[i,2]
    end = region[i,3]-region[i,2]
    target = as.character(region[,1])[i]
    subc<-rbind(subc,c(start,end,target,as.character(suf_dat[c[j2],1])))
    }
    #print(subc)
    exon<-rbind(exon,subc)
    }
  }
  #---
  fun_exon_statis<-function(posi,depth){
  len=length(posi)
  total_depth=sum(depth)
  min_depth=min(depth)
  max_depth=max(depth)
  mean_depth=round(total_depth/len,2)
  posi_min=posi[which(depth==min_depth)]
  if(length(posi_min)==1){
    posi_mins=posi_min
     }else{
    for(i in 1:length(posi_min)){
      if(i<=length(posi_min)-1){
        posi_mins=paste(posi_min[i],posi_min[i+1],sep=",")
      }
  }
  }
  
  #print(posi[which(depth==min_depth)])
  
  x50_ratio=length(which(depth>=50))/len
  x100_ratio=length(which(depth>=100))/len
  x200_ratio=length(which(depth>=200))/len
  x500_ratio=length(which(depth>=500))/len
  x50=length(which(depth>=50))
  x100=length(which(depth>=100))
  x200=length(which(depth>=200))
  x500=length(which(depth>=500))
  #print(c(len,total_depth,mean_depth, min_depth, max_depth,posi_min,x50_ratio,x100_ratio,x200_ratio))
  return(c(len,total_depth, mean_depth, min_depth, max_depth,posi_mins,x50_ratio,
          x100_ratio,x200_ratio,x500_ratio,x50,x100,x200,x500))
  }
  #---
  dat$y1=dat[,3]
  colnames(dat)<-c("chr","posi","depth","y1")
  exon_statis<-c()
  for(i in 1:nrow(exon)){
  #for(i in 1:2){
  #print(which(as.character(exon[i,3])==as.character(dat[,1])))
  a= which(as.character(exon[i,3])==dat[,1])
  a_s=dat[a,2]
  #print(exon[i,1:2])
  b=which(a_s>=as.numeric(exon[i,1]) & a_s<=as.numeric(exon[i,2]))
  #print(dat$posi[a[b]])
  #print(dat$depth[a[b]])
  exon_statis<-rbind(exon_statis,fun_exon_statis(dat$posi[a[b]],dat$depth[a[b]]))
  dat$y1[a[b]]<-0
  }
  #----
  brca_statis_len =sum(as.numeric(exon_statis[,1]))
  brca_statis_total =sum(as.numeric(exon_statis[,2]))
  brca_statis_mean_depth=round(brca_statis_total/brca_statis_len)
  brca_statis_min_depth=min(as.numeric(exon_statis[,4]))
  brca_statis_max_depth=max(as.numeric(exon_statis[,5]))
  brca_statis_posi_mins=NA
  brca_statis_x50=sum(as.numeric(exon_statis[,11]))
  brca_statis_x100=sum(as.numeric(exon_statis[,12]))
  brca_statis_x200=sum(as.numeric(exon_statis[,13]))
  brca_statis_x500=sum(as.numeric(exon_statis[,14]))
  brca_statis_x50_ratio=round(brca_statis_x50/brca_statis_len,4)
  brca_statis_x100_ratio=round(brca_statis_x100/brca_statis_len,4)
  brca_statis_x200_ratio=round(brca_statis_x200/brca_statis_len,4)
  brca_statis_x500_ratio=round(brca_statis_x500/brca_statis_len,4)
  exon_statis<-rbind(c(brca_statis_len,brca_statis_total,
                       brca_statis_mean_depth,brca_statis_min_depth,
                      brca_statis_max_depth,brca_statis_posi_mins,
                      brca_statis_x50_ratio,brca_statis_x50_ratio,
                      brca_statis_x100_ratio,brca_statis_x200_ratio,
                      brca_statis_x50,brca_statis_x100,
                      brca_statis_x200,brca_statis_x500),exon_statis)
  rownames(exon_statis)<-c('BRCA',exon[,4])
  colnames(exon_statis)<-c("exon_lenth","total_depth","mean_depth","min_depth",
                           "max_depth","min_depth_posi_on_region","x50_ratio","x100_ratio","x200_ratio","x500_ratio",
                           "x50","x100","x200","x500")
  
  #----
  dat$y1<-dat[,3]-dat$y1
  data_exon<-data.frame(chr=exon[,3],exon=exon[,4],text_x=rowSums(cbind(as.numeric(exon[,1]),as.numeric(exon[,2])))/2,
                        text_y=as.numeric(exon_statis[2:c(nrow(exon)+1),3]))
  #--
  p = ggplot(dat, aes(x = posi)) + 
      geom_line(aes(y = depth), colour = 'green') +
      #geom_line(aes(y = exon), colour = 'red') +
      geom_area(aes(y = pmin(y1, depth)), fill = 'red',alpha=0.45)+
      geom_text(data=data_exon,aes(x=text_x,y=text_y,label=exon))+
      facet_wrap(~chr,scale="free",ncol=4)+
    xlab("Position of Bases") + ylab("Depth of Bases")+
    theme_bw()
    #theme(legend.key.size=unit(1,'cm'))
    #scale_x_log10(breaks=c(0.05, 10, 50, 100, 200, 300, 500, 1000, 2000)) 
    
}

# 5. 保存图表
if (TRUE){
  dir.create(opts$output)
  ggsave(file=paste(opts$output,"/",opts$output,".pdf",sep=""), plot=p, width = 30, height = 40,limitsize = FALSE)
  #ggsave(file=paste(opts$output,".tiff",sep=""), plot=p, width = 25, height = 50)
  exon_statis=cbind(rownames(exon_statis),exon_statis)
  write.table(exon_statis, file=paste(opts$output,"/",opts$output,".exon_statis.txt",sep=""),
              append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
}
if (FALSE){
for(i in 1:nrow(region)){
dat1= subset(dat,chr ==as.character(region[,1])[i])
data_exon1=subset(data_exon,chr ==as.character(region[,1])[i])
if(nrow(data_exon1)>0){
p = ggplot(dat1, aes(x = posi)) + 
      geom_line(aes(y = depth), colour = 'green') +
      geom_text(data=data_exon1,aes(x=text_x,y=text_y,label=exon))+
      geom_area(aes(y = pmin(depth, y1)), fill = 'red',alpha=0.45)+
    xlab(as.character(region[,1])[i]) + ylab("Depth of Bases")+
    theme_bw()
  ggsave(file=paste(opts$output,"/",as.character(region[,1])[i],".pdf",sep=""), plot=p, width = 8, height = 8,limitsize = FALSE)}
}}