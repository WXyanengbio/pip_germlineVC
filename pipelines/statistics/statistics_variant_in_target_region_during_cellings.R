# 清理工作环境 clean enviroment object
rm(list=ls()) 

# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)
library(data.table)

option_list <- list(
  make_option(c("-c", "--celllines"), type="character", 
              default="/home/dell/Works/Projects/Datasets/celllines/cell-lines_BT-474,LS-180,HC-15,RL-95-2,MDA-MB-436_mutations.csv",
              help="Input celllines variants data [default %default]"),
  #make_option(c("-r", "--region"), type="character",
  #            help="Input target exon to read"),
  make_option(c("-f", "--vcf_list"), type="character",
              help="the vcf and cellline list to annotation and statistics: vcf file, cellline"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))


if (TRUE){
  #dat = read.csv(opts$prefix_file, sep="\t",header = F)
  dat = fread(opts$celllines)
  #region = fread(opts$region)
  vcflist = read.table(opts$vcf_list, sep="\t",header = T)
  #group_name= unlist(strsplit(opts$group_name,","))
}

#targeted_variants = c()
dat1 = dat[order(dat[,12])]
all_statis=c()
names=c()
unique_raw_muts=c()
vcfsummary=c()
for(j in 1:nrow(vcflist)){
  vcf = fread(as.character(vcflist[,1])[j])
  cellline = as.character(vcflist[,2])[j]
  region = fread(as.character(vcflist[,3])[j])
  #--
  nrows=c()
  for(i in 1:nrow(region)){
   nrows<-c(nrows,intersect(which(dat1[,11]==gsub('chr','',region[i,1])),intersect(which(dat1[,12]>=unlist(region[i,2])),which(dat1[,12]<=unlist(region[i,3])))))
}
#print(nrows)
nrows=unique(nrows)
targeted_variants=dat1[nrows,]

targeted_variants$post=as.numeric(targeted_variants$Start)
#big del
del=setdiff(which(targeted_variants[,12]<targeted_variants[,13]),grep('ins',unlist(targeted_variants[,6])))
#small del
del=c(del,grep('del',unlist(targeted_variants[,6])))
targeted_variants$post[del] = as.numeric(targeted_variants$Start)[del]-1
  #-
  targeted_variants_cellline = targeted_variants[which(targeted_variants[,1]==cellline),]
    # write.table(targeted_variants_cellline, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".targeted.txt",sep=""), append = T, 
    #          quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  #print(targeted_variants_cellline[,12])
  raw=unique(targeted_variants_cellline[,12])

  unique_raw=c()
  for(m in 1:nrow(raw)){
  exist=which(targeted_variants_cellline[,12]==unlist(raw[m]))
  
  #if(length(exist)==3){
  #   if(as.character(targeted_variants_cellline[exist[1],7])!='Unknown' | as.character(targeted_variants_cellline[exist[1],9])=='Verified'){
#unique_raw<-c(unique_raw,exist[1])
#}else if(as.character(targeted_variants_cellline[exist[2],7])!='Unknown' | as.character(targeted_variants_cellline[exist[2],9])=='Verified'){
#		unique_raw<-c(unique_raw,exist[2])
#}else if(as.character(targeted_variants_cellline[exist[3],7])!='Unknown' | as.character(targeted_variants_cellline[exist[3],9])=='Verified'){
	#unique_raw<-c(unique_raw,exist[3])
#}else{
#    unique_raw<-c(unique_raw,exist[1])
#}
#}else if(length(exist)==2){
#     if(as.character(targeted_variants_cellline[exist[1],7])!='Unknown' | as.character(targeted_variants_cellline[exist[1],9])!='Verified'){
#unique_raw<-c(unique_raw,exist[1])
#}else if(as.character(targeted_variants_cellline[exist[2],7])!='Unknown' | as.character(targeted_variants_cellline[exist[2],9])!='Verified'){
	#unique_raw<-c(unique_raw,exist[2])
#}else{
#unique_raw<-c(unique_raw,exist[1])
#}
  #}else if(length(exist)==1){
  #    unique_raw<-c(unique_raw,exist[1])
  #}else{
  #    unique_raw<-c(unique_raw,exist[1])
  #}
  unique_raw<-c(unique_raw, exist[1])
  }
  unique_raw_mut=targeted_variants_cellline[unique_raw,]
  mut=c()
  for(k in 1:nrow(unique_raw_mut)){
  exist=which(unlist(unique_raw_mut[k,14])==vcf[,3])
  if(length(exist)==2){
   mut<-c(mut,exist[1])
  }else if(length(exist)==1){
      mut<-c(mut,exist)
  }else{
        mut<-c(mut,NA)
  }
  }
  all = nrow(unique_raw_mut)
  targeted = length(which(!is.na(mut)))
  no_targeted = all- targeted
  wrong_mu = nrow(vcf)-targeted

  if(targeted==0){
  vcf1=vcf[1,]
  vcf1[1,1:ncol(vcf1)]<-NA
  unique_raw_mut<-cbind(unique_raw_mut,vcf1)
  vcfsummary<-rbind(vcfsummary,vcf)
  }else{
  unique_raw_mut<-cbind(unique_raw_mut,vcf[mut,])
   vcfsummary<-rbind(vcfsummary,vcf[-c(mut[which(!is.na(mut))])])
  }
  wrong_mu_withDB=length(which(vcf[,7]!='-'))-targeted
  vf_dist=quantile(unlist(vcf[which(vcf[,7]!='-'),'VF']), probs = c(0,0.25,0.5,0.75,1))
  vf_dist1=quantile(unlist(vcf[which(vcf[,7]=='-'),'VF']), probs = c(0,0.25,0.5,0.75,1))
  vf_dist2=quantile(unlist(vcf[mut[which(!is.na(mut))],'VF']), probs = c(0,0.25,0.5,0.75,1))
  statis=data.frame(type=c('mut_on_region_in_celllines','targeted','Minimum VF of targeted', '25% VF of targeted',
                           '50% VF of targeted', '75% VF of targeted',
                           'Maximun VF of targeted',
                           'no_targeted','not_region_in_celllines','not_region_in_celling_within_DB',
                           'Minimum VF of not_region_in_celling_within_DB', '25% VF of not_region_in_celling_within_DB',
                           '50% VF of not_region_in_celling_within_DB', '75% VF of not_region_in_celling_within_DB',
                           'Maximun VF of not_region_in_celling_within_DB','Minimum VF of not_region_in_celling_without_DB', '25% VF of not_region_in_celling_without_DB',
                           '50% VF of not_region_in_celling_without_DB', '75% VF of not_region_in_celling_without_DB',
                           'Maximun VF of not_region_in_celling_without_DB'),
                    variant=c(all,targeted,vf_dist2,no_targeted,wrong_mu, wrong_mu_withDB,vf_dist,vf_dist1))
  all_statis<-cbind(all_statis,statis[,2])
  names<-c(names,as.character(vcf[1,1]))
  write.table(statis, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".summary.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  write.table(unique_raw_mut, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".mut_targeted.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  unique_raw_mut1=data.frame(sample=as.character(vcf[1,1]),unique_raw_mut)
  unique_raw_muts=rbind(unique_raw_muts,unique_raw_mut1)
}
all_statis<-cbind(type=c('mut_on_region_in_celllines','targeted','Minimum VF of targeted', '25% VF of targeted',
                           '50% VF of targeted', '75% VF of targeted',
                           'Maximun VF of targeted','no_targeted','not_region_in_celllines','not_region_in_celling_within_DB',
                         'Minimum VF of not_region_in_celling_within_DB', '25% VF of not_region_in_celling_within_DB',
                           '50% VF of not_region_in_celling_within_DB', '75% VF of not_region_in_celling_within_DB',
                           'Maximun VF of not_region_in_celling_within_DB','Minimum VF of not_region_in_celling_without_DB', '25% VF of not_region_in_celling_without_DB',
                           '50% VF of not_region_in_celling_without_DB', '75% VF of not_region_in_celling_without_DB',
                           'Maximun VF of not_region_in_celling_without_DB'),all_statis)
colnames(all_statis)<-c('type',names)
  write.table(t(all_statis), file=paste(opts$output,"/","All.summary.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = F)
  write.table(unique_raw_muts, file=paste(opts$output,"/","All.targeted_summary.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  write.table(vcfsummary, file=paste(opts$output,"/","All.no_celllines_sites_summary.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
