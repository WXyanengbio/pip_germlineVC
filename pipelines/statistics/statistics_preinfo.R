# 清理工作环境 clean enviroment object
rm(list=ls()) 

# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)
library(data.table)



option_list <- list(
  make_option(c("-q", "--qc"), type="character", 
              help="Input qc"),
  make_option(c("-t", "--trim_qc"), type="character", 
              help="Input trim_qc"),
  make_option(c("-c", "--trim_statis"), type="character", 
              help="Input trim_statis"),
  make_option(c("-f", "--filter_statis"), type="character", 
              help="Input filter_statis"),
  make_option(c("-p", "--primer_statis"), type="character", 
              help="Input primer_statis"),
  make_option(c("-u", "--umi_statis"), type="character", 
              help="Input umi_statis"),
  make_option(c("-a", "--align_base"), type="character", 
              help="Input align_base"),
  make_option(c("-e", "--exon"), type="character", 
              help="Input exon"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

file1s <- read.table(opts$qc)
file2s <- read.table(opts$trim_qc)
 if(TRUE){
 qcs=c()
 names=c()
 for(i in 1:nrow(file1s)){
 file1<-as.character(file1s[i,1])
 file2<-as.character(file2s[i,1])
 dat1 = read.table(file1,header=T,sep='\t')
 dat2 = read.table(file2,,header=T,sep='\t')
 c_q20<-round((as.numeric(gsub('%','',dat2[1,9]))+as.numeric(gsub('%','',dat2[2,9])))/2,2)
 c_q30<-round((as.numeric(gsub('%','',dat2[1,10]))+as.numeric(gsub('%','',dat2[2,10])))/2,2)
 qc<-c(dat1[1,3],dat2[1,3],round(dat2[1,3]*100/dat1[1,3],2),dat2[1,4:7],c_q20,c_q30)
 qcs<-cbind(qcs,qc)
 name=unlist(strsplit(file1,"/"))[length(unlist(strsplit(file1,"/")))]
 names<-c(names,gsub('_L001.QC.statistics.txt','',name))
 }
 qcs<-as.data.frame(qcs)
 colnames(qcs)<-names
 qcs=cbind(read_set=c('raw reads','clean reads','Precent of trim','min length of clean reads',
                      'max length of clean reads','GC content of clean reads','mean base qualit of clean reads',
                      'Q20 of clean reads','Q30 of clean reads'),
           qcs)
  #write.table(trims, file=paste(opts$output,".trims_statis.txt",sep=""),
  #            append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
 }
 files<-read.table(opts$trim_statis)
 if(TRUE){
 trims=c()
 names=c()
 for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 dat = read.table(file,sep="\t")
 dat = dat[-c(1,nrow(dat)),]
 trim=c()
 for(j in 1:length(dat)){
 if(length(unlist(strsplit(as.character(dat[j])," == ")))==2){
 trim<-rbind(trim,unlist(strsplit(as.character(dat[j])," == ")))
 }else{trim<-rbind(trim,unlist(strsplit(as.character(dat[j]),"== ")))}
 }
 
 trims<-cbind(trims,trim[,2])
 name=unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/")))]
 names<-c(names,gsub('_L001_basic_stats.txt','',name))
 }
 trims<-as.data.frame(trims)
 colnames(trims)<-names
 names1<-c()
  trims=cbind(read_set=trim[,1],trims)
  #write.table(trims, file=paste(opts$output,".trims_statis.txt",sep=""),
  #            append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
 }

 files<-read.table(opts$filter_statis)
 if(TRUE){
 filters=c()
 names=c()
for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 dat = read.table(file,sep="\t")
 dat = dat[,1]
 #dat = dat[-c(4,8,13),]
 filter=c()
 for(j in 1:length(dat)){
 if(j>=10){
   subfilter=unlist(strsplit(as.character(dat[j])," == "))
   number=unlist(strsplit(subfilter[2]," "))[1]
   precent=substr(unlist(strsplit(subfilter[2]," "))[2],2,6)
  filter<-rbind(filter,c(subfilter[1],number))
  filter<-rbind(filter,c(gsub("Precent","Number",subfilter[1]),precent))
 }else{
   if(j==7){
   subfilter=unlist(strsplit(as.character(dat[j])," == "))
   number=substr(subfilter[2],1,5)
   filter<-rbind(filter,c(subfilter[1],number))
   }else{filter<-rbind(filter,unlist(strsplit(as.character(dat[j])," == ")))
        }
 }
 }

 filters<-cbind(filters,filter[,2])
 name=unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/")))]
 names<-c(names,gsub('_L001_align_stats.txt','',name))
 }
 filters <-as.data.frame(filters)
 #print(names)
 colnames(filters)<-names
 filters=cbind(read_set=filter[,1],filters)
  #write.table(filters, file=paste(opts$output,".filters_statis.txt",sep=""),
  #            append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
 }

 files<-read.table(opts$umi_statis)
 #print(files)
 if(TRUE){
 umis=c()
 names=c()
 for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 dat = fread(file)
# print(head(dat))
 reads= sum(dat[,7])
 sums=nrow(dat)
 mean=reads/sums
 quant=quantile(unlist(dat[,7]), probs = c(0.25,0.5,0.75,1))
 umis <- cbind(umis,c(reads,sums,mean,quant))
 name=unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/")))]
 names<-c(names,gsub('_L001_deduplicated_per_umi.tsv','',name))
 }
 umis<-as.data.frame(umis)
 colnames(umis)<-names
 
  umis=cbind(read_set=c("read fragments","MTs","read fragments per MT, mean",
               "read fragments per MT, 25th percentile","read fragments per MT, 50th percentile",
              "read fragments per MT, 75th percentile","read fragments per MT, 100th percentile"),umis)
  #write.table(umis, file=paste(opts$output,".umis_statis.txt",sep=""),
  #            append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
 }

 files<-read.table(opts$primer_statis)
  if(TRUE){
 primers=c()
 names=c()
for(i in 1:nrow(files)){
 file<-as.character(files[i,1])
 dat = fread(file)
 num<-nrow(dat)
 meanprimer_mt_depth<-round(umis[2,i+1]/num,2)
 x5<-round(length(which(dat[,6]>meanprimer_mt_depth*0.05))*100/num,2)
 x25<-round(length(which(dat[,6]>meanprimer_mt_depth*0.25))*100/num,2)
 x50<-round(length(which(dat[,6]>meanprimer_mt_depth*0.50))*100/num,2)
 x75<-round(length(which(dat[,6]>meanprimer_mt_depth*0.75))*100/num,2)
 x100<-round(length(which(dat[,6]>meanprimer_mt_depth*1))*100/num,2)
 mean_primer_read_depth <-round(as.numeric(qcs[2,i+1])/num,2)
 rx5<-round(length(which(dat[,6]>mean_primer_read_depth*0.05))*100/num,2)
 rx25<-round(length(which(dat[,6]>mean_primer_read_depth*0.25))*100/num,2)
 rx50<-round(length(which(dat[,6]>mean_primer_read_depth*0.50))*100/num,2)
 rx75<-round(length(which(dat[,6]>mean_primer_read_depth*0.75))*100/num,2)
 rx100<-round(length(which(dat[,6]>mean_primer_read_depth*1))*100/num,2)
 primers<-cbind(primers,c(num,meanprimer_mt_depth,x5,x25,x50,x75,x100,mean_primer_read_depth,rx5,rx25,rx50,rx75,rx100))
 name=unlist(strsplit(file,"/"))[length(unlist(strsplit(file,"/")))]
 names<-c(names,gsub('_L001_primer_stats.csv','',name))
 }
   primers<-as.data.frame(primers)
 colnames(primers)<-names
 primers=cbind(read_set=c("primers","mean primer MT depth","% of primers >= 5% of mean MT depth",
               "% of primers >= 25% of mean MT depth","% of primers >= 50% of mean MT depth",
              "% of primers >= 75% of mean MT depth","% of primers >= 100% of mean MT depth",
                        "mean primer read fragment depth","% of primers >= 5% of mean read fragment depth",
              "% of primers >= 25% of mean read fragment depth","% of primers >= 50% of mean read fragment depth",
                         "% of primers >= 75% of mean read fragment depth","% of primers >= 100% of mean read fragment depth"),primers)
  }

#--------------------------
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

  x50_ratio=length(which(depth>=50))/len
  x100_ratio=length(which(depth>=100))/len
  x200_ratio=length(which(depth>=200))/len
  x500_ratio=length(which(depth>=500))/len
  x50=length(which(depth>=50))
  x100=length(which(depth>=100))
  x200=length(which(depth>=200))
  x500=length(which(depth>=500))
  x5_men=length(which(depth>=(0.05*mean_depth)))
  x25_men=length(which(depth>=(0.25*mean_depth)))
  x50_men=length(which(depth>=(0.5*mean_depth)))
  x75_men=length(which(depth>=(0.75*mean_depth)))
  x100_men=length(which(depth>=mean_depth))
  #print(c(len,total_depth,mean_depth, min_depth, max_depth,posi_min,x50_ratio,x100_ratio,x200_ratio))
  return(c(len,total_depth, mean_depth, min_depth, max_depth,posi_mins,x50_ratio,
          x100_ratio,x200_ratio,x500_ratio,x50,x100,x200,x500,
          x5_men,x25_men,x50_men,
         x75_men,x100_men))
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

  x50_ratio=length(which(depth>=50))/len
  x100_ratio=length(which(depth>=100))/len
  x200_ratio=length(which(depth>=200))/len
  x500_ratio=length(which(depth>=500))/len
  x50=length(which(depth>=50))
  x100=length(which(depth>=100))
  x200=length(which(depth>=200))
  x500=length(which(depth>=500))
  x5_men=length(which(depth>=(0.05*mean_depth)))
  x25_men=length(which(depth>=(0.25*mean_depth)))
  x50_men=length(which(depth>=(0.5*mean_depth)))
  x75_men=length(which(depth>=(0.75*mean_depth)))
  x100_men=length(which(depth>=mean_depth))
  #print(c(len,total_depth,mean_depth, min_depth, max_depth,posi_min,x50_ratio,x100_ratio,x200_ratio))
  return(c(len,total_depth, mean_depth, min_depth, max_depth,posi_mins,x50_ratio,
          x100_ratio,x200_ratio,x500_ratio,x50,x100,x200,x500,
          x5_men,x25_men,x50_men,
         x75_men,x100_men))
  }
#---
 files<-read.table(opts$align_base)
if(TRUE){
 exon_statis_mes=c()
 names=c()
 for(i in 1:nrow(files)){
  file<-as.character(files[i,1])
  dat = fread(file)
  region = data.frame(region=unique(dat[,1]))
  region$start=as.numeric(apply(region,1,function(x){unlist(strsplit(as.character(x[1]),"_"))[2]}))
  region$end=as.numeric(apply(region,1,function(x){unlist(strsplit(as.character(x[1]),"_"))[3]}))
  exon<-c()

  if(opts$exon != 'n'){
    suf_dat = read.csv(opts$exon,header=F)
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
  for(i in 1:nrow(exon)){
  #for(i in 1:2){exon_statis_mes
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
  brca_statis_mean_depth=round(brca_statis_total/brca_statis_len,2)
  brca_statis_min_depth=min(as.numeric(exon_statis[,4]))
  brca_statis_max_depth=max(as.numeric(exon_statis[,5]))
  brca_statis_posi_mins=NA
  brca_statis_x50=sum(as.numeric(exon_statis[,11]))
  brca_statis_x100=sum(as.numeric(exon_statis[,12]))
  brca_statis_x200=sum(as.numeric(exon_statis[,13]))
  brca_statis_x500=sum(as.numeric(exon_statis[,14]))
  brca_statis_x50_ratio=round(brca_statis_x50*100/brca_statis_len,2)
  brca_statis_x100_ratio=round(brca_statis_x100*100/brca_statis_len,2)
  brca_statis_x200_ratio=round(brca_statis_x200*100/brca_statis_len,2)
  brca_statis_x500_ratio=round(brca_statis_x500*100/brca_statis_len,2)
  #---use the total mean
   exon_statis_me<-c()
  for(i in 1:nrow(exon)){
  #for(i in 1:2){
  #print(which(as.character(exon[i,3])==as.character(dat[,1])))
  a= which(as.character(exon[i,3])==dat[,1])
  a_s=dat[a,2]
  #print(exon[i,1:2])
  b=which(a_s>=as.numeric(exon[i,1]) & a_s<=as.numeric(exon[i,2]))
  #print(dat$posi[a[b]])
  #print(dat$depth[a[b]])
  exon_statis_me<-rbind(exon_statis_me,fun_exon_statis1(dat$posi[a[b]],dat$depth[a[b]],brca_statis_mean_depth))
  dat$y1[a[b]]<-0
  }
  #---
  x5_men=round(sum(as.numeric(exon_statis_me[,15]))*100/brca_statis_len,2)
  x25_men=round(sum(as.numeric(exon_statis_me[,16]))*100/brca_statis_len,2)
  x50_men=round(sum(as.numeric(exon_statis_me[,17]))*100/brca_statis_len,2)
  x75_men=round(sum(as.numeric(exon_statis_me[,18]))*100/brca_statis_len,2)
  x100_men=round(sum(as.numeric(exon_statis_me[,19])*100)/brca_statis_len,2)
  exon_statis_me<-rbind(c(brca_statis_len,brca_statis_total,
                       brca_statis_mean_depth,brca_statis_min_depth,
                      brca_statis_max_depth,brca_statis_posi_mins,
                      brca_statis_x50_ratio,brca_statis_x100_ratio,
                      brca_statis_x200_ratio,brca_statis_x500_ratio,
                      brca_statis_x50,brca_statis_x100,
                      brca_statis_x200,brca_statis_x500,
                      x5_men,x25_men,x50_men,
                      x75_men,x100_men),exon_statis_me)

  exon_statis_me1<-exon_statis_me[1,-6]
  exon_statis_mes<-cbind(exon_statis_mes,exon_statis_me1)
  }else{
    len=nrow(dat[,2])
    depth=dat[,3]
  total_depth=sum(depth)
  min_depth=min(depth)
  max_depth=max(depth)
  mean_depth1=round(total_depth/len,2)
  #posi_min=posi[which(depth==min_depth)]
  x50_ratio=round(length(which(depth>=50))*100/len,2)
  x100_ratio=round(length(which(depth>=100))*100/len,2)
  x200_ratio=round(length(which(depth>=200))*100/len,2)
  x500_ratio=round(length(which(depth>=500))*100/len,2)
  x50=length(which(depth>=50))
  x100=length(which(depth>=100))
  x200=length(which(depth>=200))
  x500=length(which(depth>=500))
  x5_men=round(length(which(depth>=(0.05*mean_depth1)))*100/len,2)
  x25_men=round(length(which(depth>=(0.25*mean_depth1)))*100/len,2)
  x50_men=round(length(which(depth>=(0.5*mean_depth1)))*100/len,2)
  x75_men=round(length(which(depth>=(0.75*mean_depth1)))*100/len,2)
  x100_men=round(length(which(depth>=mean_depth1))*100/len,2)
  exon_statis_me<-c(len,total_depth,mean_depth1,min_depth,max_depth,x50_ratio,
          x100_ratio,x200_ratio,x500_ratio,x50,x100,x200,x500,
          x5_men,x25_men,x50_men,
         x75_men,x100_men)
  print(exon_statis_me)
  exon_statis_mes<-cbind(exon_statis_mes, exon_statis_me)
 }
}
  exon_statis_mes<-as.data.frame(exon_statis_mes)
  colnames(exon_statis_mes)<-names
  
  exon_statis_mes<-cbind(read_set=c("target exon/region lenth","total bases depth of target exon/region","mean bases depth of target exon/region","min bases depth of target exon/region",
                           "max bases depth of target exon/region","precent of x50 bases depth of target exon/region","precent of x100 bases depth of target exon/region",
                           "precent of x200 bases depth of target exon/region","precent of x500 bases depth of target exon/region",
                           "Number of x50 bases depth of target exon/region","Number of x100 bases depth of target exon/region",
                           "Number of x200 bases depth of target exon/region","Number of x500 bases depth of target exon/region","precent of 5% mean bases depthof target exon/region",
                           "precent of 25% mean bases depth of target exon/region","precent of 50% mean bases depth of target exon/region",
                           "precent of 75% mean bases depth of target exon/region","precent of 100% mean bases depth of target exon/region"),exon_statis_mes)
print(exon_statis_mes)
}
colnames(qcs)<- colnames(qcs)
colnames(trims)<-colnames(qcs)
colnames(filters)<-colnames(qcs)
colnames(primers)<-colnames(qcs)
colnames(exon_statis_mes)<-colnames(qcs)

data <-rbind(qcs,trims)
data<-rbind(data,filters)
data<-rbind(data,primers)
print(data)
print(exon_statis_mes)
data<-rbind(data,exon_statis_mes)
#print(data)
#print(class(data))
colnames(data)<-colnames(qcs)
data<- as.matrix(data)
write.table(data, file=opts$output, append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)