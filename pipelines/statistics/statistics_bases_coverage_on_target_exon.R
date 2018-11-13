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

# 显示输入输出确认是否正确
#print(paste("The prefix file is ", opts$prefix_file, sep = ""))
#print(paste("The output file prefix is ", opts$output, sep = ""))
#dir.create(opts$output)

# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从文件中读取
if (TRUE){
  #dat = read.csv(opts$prefix_file, sep="\t",header = F)
  dat = fread(opts$prefix_file)
  suf_dat = read.csv(opts$suffix_file,header=F)
  #group_name= unlist(strsplit(opts$group_name,","))
}
#print(head(suf_dat))
#print(head(suf_dat[1,2]))
#----
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
  x5_men=length(which(depth>=(0.05*mean_depth)))
  x10_men=length(which(depth>=(0.1*mean_depth)))
  x20_men=length(which(depth>=(0.2*mean_depth)))
  x30_men=length(which(depth>=(0.3*mean_depth)))
  x40_men=length(which(depth>=(0.4*mean_depth)))
  x50_men=length(which(depth>=(0.5*mean_depth)))
  x60_men=length(which(depth>=(0.6*mean_depth)))
  x70_men=length(which(depth>=(0.7*mean_depth)))
  x80_men=length(which(depth>=(0.8*mean_depth)))
  x90_men=length(which(depth>=(0.9*mean_depth)))
  x100_men=length(which(depth>=mean_depth))
  #print(c(len,total_depth,mean_depth, min_depth, max_depth,posi_min,x50_ratio,x100_ratio,x200_ratio))
  return(c(len,total_depth, mean_depth, min_depth, max_depth,posi_mins,x50_ratio,
          x100_ratio,x200_ratio,x500_ratio,x50,x100,x200,x500,
          x5_men,x10_men,x20_men,x30_men,x40_men,x50_men,
          x60_men,x70_men,x80_men,x90_men,x100_men))
  }
#----
fun_exon_statis1<-function(posi,depth,mean_depth){
  len=length(posi)
  total_depth=sum(depth)
  min_depth=min(depth)
  max_depth=max(depth)
  mean_depth1=round(total_depth/len,2)
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
  x5_men=length(which(depth>=(0.05*mean_depth)))
  x10_men=length(which(depth>=(0.1*mean_depth)))
  x20_men=length(which(depth>=(0.2*mean_depth)))
  x30_men=length(which(depth>=(0.3*mean_depth)))
  x40_men=length(which(depth>=(0.4*mean_depth)))
  x50_men=length(which(depth>=(0.5*mean_depth)))
  x60_men=length(which(depth>=(0.6*mean_depth)))
  x70_men=length(which(depth>=(0.7*mean_depth)))
  x80_men=length(which(depth>=(0.8*mean_depth)))
  x90_men=length(which(depth>=(0.9*mean_depth)))
  x100_men=length(which(depth>=mean_depth))
  #print(c(len,total_depth,mean_depth, min_depth, max_depth,posi_min,x50_ratio,x100_ratio,x200_ratio))
  return(c(len,total_depth, mean_depth1, min_depth, max_depth,posi_mins,x50_ratio,
          x100_ratio,x200_ratio,x500_ratio,x50,x100,x200,x500,
          x5_men,x10_men,x20_men,x30_men,x40_men,x50_men,
          x60_men,x70_men,x80_men,x90_men,x100_men))
  }
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
  
  #---
  dat$y1=dat[,3]
  colnames(dat)<-c("chr","posi","depth","y1")
  exon_statis<-c()
  exonlist<-c()
  for(i in 1:nrow(exon)){
  a= which(as.character(exon[i,3])==dat[,1])
  a_s=dat[a,2]
  b=which(a_s>=as.numeric(exon[i,1]) & a_s<=as.numeric(exon[i,2]))
  if(length(dat$depth[a[b]])>0){
  exon_statis<-rbind(exon_statis,fun_exon_statis(dat$posi[a[b]],dat$depth[a[b]]))
   exonlist<-c(exonlist, i)
  }
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
  #---use the total mean
   exon_statis_me<-c()

  for(i in 1:nrow(exon)){
  a= which(as.character(exon[i,3])==dat[,1])
  a_s=dat[a,2]
  b=which(a_s>=as.numeric(exon[i,1]) & a_s<=as.numeric(exon[i,2]))
  if(length(dat$depth[a[b]])>0){
  exon_statis_me<-rbind(exon_statis_me,fun_exon_statis1(dat$posi[a[b]],dat$depth[a[b]],brca_statis_mean_depth))
  }
  dat$y1[a[b]]<-0
  }
  #---
  x5_men=round(sum(as.numeric(exon_statis_me[,15]))/brca_statis_len,4)
  x10_men=round(sum(as.numeric(exon_statis_me[,16]))/brca_statis_len,4)
  x20_men=round(sum(as.numeric(exon_statis_me[,17]))/brca_statis_len,4)
  x30_men=round(sum(as.numeric(exon_statis_me[,18]))/brca_statis_len,4)
  x40_men=round(sum(as.numeric(exon_statis_me[,19]))/brca_statis_len,4)
  x50_men=round(sum(as.numeric(exon_statis_me[,20]))/brca_statis_len,4)
  x60_men=round(sum(as.numeric(exon_statis_me[,21]))/brca_statis_len,4)
  x70_men=round(sum(as.numeric(exon_statis_me[,22]))/brca_statis_len,4)
  x80_men=round(sum(as.numeric(exon_statis_me[,23]))/brca_statis_len,4)
  x90_men=round(sum(as.numeric(exon_statis_me[,24]))/brca_statis_len,4)
  x100_men=round(sum(as.numeric(exon_statis_me[,25]))/brca_statis_len,4)
  exon_statis_me<-rbind(c(brca_statis_len,brca_statis_total,
                       brca_statis_mean_depth,brca_statis_min_depth,
                      brca_statis_max_depth,brca_statis_posi_mins,
                      brca_statis_x50_ratio,brca_statis_x100_ratio,
                      brca_statis_x200_ratio,brca_statis_x500_ratio,
                      brca_statis_x50,brca_statis_x100,
                      brca_statis_x200,brca_statis_x500,
                      x5_men,x10_men,x20_men,x30_men,x40_men,x50_men,
                      x60_men,x70_men,x80_men,x90_men,x100_men),exon_statis_me)
  rownames(exon_statis_me)<-c('BRCA',as.character(exon[exonlist,4]))
  colnames(exon_statis_me)<-c("exon_lenth","total_depth","mean_depth","min_depth",
                           "max_depth","min_depth_posi_on_region","x50_ratio","x100_ratio","x200_ratio","x500_ratio",
                           "x50","x100","x200","x500","x5_men","x10_men","x20_men","x30_men","x40_men","x50_men",
                           "x60_men","x70_men","x80_men","x90_men","x100_men")
  
  #----
  dat$y1<-dat[,3]-dat$y1
  exon <- exon[exonlist, ]
  data_exon<-data.frame(chr=exon[,3],exon=exon[,4],text_x=rowSums(cbind(as.numeric(exon[,1]),as.numeric(exon[,2])))/2,
                        text_y=as.numeric(exon_statis_me[2:c(nrow(exon)+1),3]))
  if(opts$tiff){
  p = ggplot(dat, aes(x = posi)) + 
      geom_line(aes(y = depth), colour = 'green') +
      #geom_line(aes(y = exon), colour = 'red') +
      geom_area(aes(y = pmin(y1, depth)), fill = 'red',alpha=0.45)+
      geom_text(data=data_exon,aes(x=text_x,y=text_y,label=exon))+
      facet_wrap(~chr,scale="free",ncol=4)+
    xlab("Position of Bases") + ylab("Depth of Bases")+
    theme_bw()
  ggsave(file=paste(opts$output,".pdf",sep=""), plot=p, width = 30, height = 40,limitsize = FALSE)
  # ggsave(file=paste(opts$output,"/",opts$output,".pdf",sep=""), plot=p, width = 30, height = 40,limitsize = FALSE)
}
}

  exon_statis_me=cbind(rownames(exon_statis_me),exon_statis_me)
  write.table(exon_statis_me, file=paste(opts$output,".exon_statis.txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
