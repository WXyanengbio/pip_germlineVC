# 清理工作环境 clean enviroment object
rm(list=ls()) 
options(warn=-1)
# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(splines)
suppressMessages(library(data.table))
# sudo apt-get install r-cran-rjava
# install.packaegs("XLConnectJars")
# install.packaegs("XLConnect")
suppressMessages(library(XLConnect))

option_list <- list(
  make_option(c("-c", "--celllines"), type="character", 
              default="/home/dell/Works/Projects/Datasets/celllines/cell-lines_BT-474,LS-180,HC-15,RL-95-2,MDA-MB-436_mutations.csv",
              help="Input celllines variants data [default %default]"),
  make_option(c("-m", "--model"), type="character",default="1",
              help="the number of softwares to calling vcf : [default %default]"),
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
if(opts$model=='1'){
  for(j in 1:nrow(vcflist)){
    vcf = fread(as.character(vcflist[,1])[j])
    cellline = as.character(vcflist[,2])[j]
    region = fread(as.character(vcflist[,3])[j])

    nrows=c()
    for(i in 1:nrow(region)){
      nrows<-c(nrows,intersect(which(dat1[,11]==gsub('chr','',region[i,1])),
                               intersect(which(dat1[,12]>=unlist(region[i,2])),
                                         which(dat1[,12]<=unlist(region[i,3])))
                              )
              )
    }

    nrows=unique(nrows)
    targeted_variants=dat1[nrows,]
    targeted_variants$post=as.numeric(targeted_variants$Start)
    #big del
    del=setdiff(which(targeted_variants[,12]<targeted_variants[,13]), 
                grep('ins',unlist(targeted_variants[,6])))
    #small del
    del=c(del,grep('del',unlist(targeted_variants[,6])))
    targeted_variants$post[del] = as.numeric(targeted_variants$Start)[del]-1
    #-
    targeted_variants_cellline = targeted_variants[which(targeted_variants[,1]==cellline),]
    if(nrow(targeted_variants_cellline)>0){
      raw=unique(targeted_variants_cellline[,12])
      unique_raw=c()
      for(m in 1:nrow(raw)){
        exist=which(targeted_variants_cellline[,12]==unlist(raw[m]))
        unique_raw<-c(unique_raw, exist[1])
      }
      unique_raw_mut=targeted_variants_cellline[unique_raw,]
      #print(class(unique_raw_mut))
      mut=c()
      for(k in 1:nrow(unique_raw_mut)){
        if(opts$model=='1'){
          exist=which(unlist(unique_raw_mut[k,14])==unlist(vcf[,3]))
          }
        if(opts$model=='4'){
          posti= apply(vcf,1,function(x){unlist(strsplit(as.character(x),"_"))[2]})
          exist=which(unlist(unique_raw_mut[k,14])==posti)
          }
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
      #---
      all_indel = length(grep('ins',unique_raw_mut[,6])) + length(grep('del',unique_raw_mut[,6]))
      all_snp = all - all_indel
      targeted_vcf = unique_raw_mut[which(!is.na(mut)),]
      target_indel = length(grep('ins',targeted_vcf[,6])) + length(grep('del',targeted_vcf[,6]))
      target_snp = nrow(targeted_vcf) - target_indel
    }else{
      all = 0
      targeted = 0
      no_targeted = 0
      wrong_mu = nrow(vcf)-targeted
      #---
      all_indel = 0
      all_snp = 0
      targeted_vcf = 0
      target_indel = 0
      target_snp = 0
      mut = NA
      # unique_raw_mut = data.frame(matrix(NA, nrow=2,ncol=15))
    }
    if(targeted==0){
      vcf1=vcf[1,]
      vcf1[1,1:ncol(vcf1)]<-NA
      vcf <-cbind(sample_cell=rep(cellline,nrow(vcf)),vcf)
      vcfsummary<-rbind(vcfsummary, vcf)
      if(nrow(targeted_variants_cellline)>0){
        unique_raw_mut<-cbind(unique_raw_mut, vcf1)
      }
    }else{
      unique_raw_mut<-cbind(unique_raw_mut,vcf[mut,])
      vcf <-cbind(sample_cell=rep(cellline,nrow(vcf)), vcf)
      vcfsummary<-rbind(vcfsummary,vcf[-c(mut[which(!is.na(mut))])])
    }
    wrong_mu_withDB=length(which(vcf[,8]!='-'))-targeted
    vf_dist=quantile(unlist(vcf[which(vcf[,8]!='-'),'VF']), probs = c(0,0.25,0.5,0.75,1))
    vf_dist1=quantile(unlist(vcf[which(vcf[,8]=='-'),'VF']), probs = c(0,0.25,0.5,0.75,1))
    vf_dist2=quantile(unlist(vcf[mut[which(!is.na(mut))],'VF']), probs = c(0,0.25,0.5,0.75,1))
    calling_snp = sum(apply(vcf,1,function(x){ifelse(all(nchar(as.character(x[5]))==1,nchar(as.character(x[6]))==1),1,0)}))
    calling_indel=nrow(vcf) - calling_snp
    statis=data.frame(type=c('mut_on_region_in_celllines','targeted', 'snp_on_region_in_celllines', 
                           'indel_on_region_in_celllines', 'targeted snp', 'targeted indel', 
                           'calling snp', 'calling indel',
                           'Minimum VF of targeted', '25% VF of targeted',
                           '50% VF of targeted', '75% VF of targeted',
                           'Maximun VF of targeted',
                           'no_targeted','not_region_in_celllines','not_region_in_celling_within_DB',
                           'Minimum VF of not_region_in_celling_within_DB', '25% VF of not_region_in_celling_within_DB',
                           '50% VF of not_region_in_celling_within_DB', '75% VF of not_region_in_celling_within_DB',
                           'Maximun VF of not_region_in_celling_within_DB','Minimum VF of not_region_in_celling_without_DB', 
                           '25% VF of not_region_in_celling_without_DB',
                           '50% VF of not_region_in_celling_without_DB', '75% VF of not_region_in_celling_without_DB',
                           'Maximun VF of not_region_in_celling_without_DB'),
                    variant=c(all,targeted,all_snp, all_indel,target_snp,target_indel,calling_snp,calling_indel, 
                              vf_dist2,no_targeted,wrong_mu, wrong_mu_withDB,vf_dist,vf_dist1))
    all_statis<-cbind(all_statis,statis[,2])
    names<-c(names,as.character(vcf[1,2]))
  #write.table(statis, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".summary.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  #write.table(unique_raw_mut, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".mut_targeted.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
    if(nrow(targeted_variants_cellline)>0){
      unique_raw_mut1=data.frame(sample=as.character(vcf[1,2]),unique_raw_mut)
      unique_raw_muts=rbind(unique_raw_muts, unique_raw_mut1)
    }
  }

  all_statis<-cbind(type=c('mut_on_region_in_celllines','targeted','snp_on_region_in_celllines', 
                         'indel_on_region_in_celllines', 'targeted snp', 'targeted indel', 
                           'calling snp', 'calling indel', 'Minimum VF of targeted', '25% VF of targeted',
                           '50% VF of targeted', '75% VF of targeted',
                           'Maximun VF of targeted','no_targeted','not_region_in_celllines','not_region_in_celling_within_DB',
                         'Minimum VF of not_region_in_celling_within_DB', '25% VF of not_region_in_celling_within_DB',
                           '50% VF of not_region_in_celling_within_DB', '75% VF of not_region_in_celling_within_DB',
                           'Maximun VF of not_region_in_celling_within_DB','Minimum VF of not_region_in_celling_without_DB',
                         '25% VF of not_region_in_celling_without_DB',
                           '50% VF of not_region_in_celling_without_DB', '75% VF of not_region_in_celling_without_DB',
                           'Maximun VF of not_region_in_celling_without_DB'),all_statis)
  colnames(all_statis)<-c('type',names)
  #write.table(t(all_statis), file=paste(opts$output,"/","All.summary.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = F)
  #write.table(unique_raw_muts, file=paste(opts$output,"/","All.targeted_summary.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  #write.table(vcfsummary, file=paste(opts$output,"/","All.no_celllines_sites_summary.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  #--

  cellline_mut <- paste(unique_raw_muts[,2],unique_raw_muts[,15],sep="_")
  unq_cellline_mut<-unique(cellline_mut)
  muts<-matrix(data=NA,nrow=length(unq_cellline_mut), ncol=length(names))
  mutinfos<-c()
  for(i in 1:length(unq_cellline_mut)){
   celline <-unlist(strsplit(unq_cellline_mut[i],'_'))[1]
   mut <-unlist(strsplit(unq_cellline_mut[i],'_'))[2]
   infos=intersect(which(unique_raw_muts[,2]==celline), which(unique_raw_muts[,15]==mut))
   mutinfos<-rbind(mutinfos, unique_raw_muts[infos[1], c(2,12,15,3,4,6,7,11,42,45,46,55,56)])
   for(j in 1:length(names)){
   subda=subset(unique_raw_muts, sample == names[j])
   num<-which(subda[,15]==mut)
   if(length(num)>0){
    if(!is.na(subda[num,23])){
      muts[i,j] <- paste("DP=",subda[num,23], "  ", "VF=",subda[num,25], sep="")
    }
   }
   }
  }
  colnames(muts)<-names
  muts<-cbind(mutinfos, muts)
  # 读取或创建一个XLSX文件，此步相当于建立一个连接
  filenames = paste(opts$output,'.Celllines.summary.xlsx',sep="")
  xls <- loadWorkbook(filenames,create=TRUE) 
  # 创建工作表
  createSheet(xls,name='summary')
  all_statis1<-data.frame(type=colnames(all_statis), t(all_statis))
  writeWorksheet(xls,all_statis1,'summary',startRow=1,startCol=1,header=FALSE)
  createSheet(xls,name='variants')
  writeWorksheet(xls,muts,'variants',startRow=1,startCol=1,header=TRUE)
  createSheet(xls,name='variants_details')
  writeWorksheet(xls,unique_raw_muts,'variants_details',startRow=1,startCol=1,header=TRUE)
  createSheet(xls,name='no_celllines_mut')
  writeWorksheet(xls,vcfsummary,'no_celllines_mut',startRow=1,startCol=1,header=TRUE)
  # 写入数据
 # 存入硬盘，直到此步方才有文档生成
  saveWorkbook(xls)
  print('Done')
}
  
if(opts$model=='4'){
for(j in 1:nrow(vcflist)){

  vcf = fread(as.character(vcflist[,1])[j])
  cellline = as.character(vcflist[,2])[j]
  region = fread(as.character(vcflist[,3])[j])
  #--
  nrows=c()
  for(i in 1:nrow(region)){
   nrows<-c(nrows,intersect(which(dat1[,11]==gsub('chr','',region[i,1])),
                            intersect(which(dat1[,12]>=unlist(region[i,2])),which(dat1[,12]<=unlist(region[i,3])))))
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
  unique_raw<-c(unique_raw, exist[1])
  }
  unique_raw_mut=targeted_variants_cellline[unique_raw,]

  mut=c()
  for(k in 1:nrow(vcf)){
    posti= unlist(strsplit(as.character(vcf[k,1]),"_"))[2]
    exist=which(unlist(unique_raw_mut[,14])==posti)
    
  if(length(exist)==2){
   mut<-c(mut,'yes')
  }else if(length(exist)==1){
      mut<-c(mut,'yes')
  }else{
        mut<-c(mut,'no')
  }
  }
  vcf=cbind(vcf[,1],mut,vcf[,2:ncol(vcf)])
  name=unlist(strsplit(as.character(vcflist[,1])[j],'/'))[length(unlist(strsplit(as.character(vcflist[,1])[j],'/')))]
  name=gsub('.mutSummary_PASS.txt','',name)
  write.table(vcf, file=paste(opts$output,"/",name,"_",cellline,"_targeted.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
}
}



if(opts$model=='sm'){
  for(j in 1:nrow(vcflist)){
    #print(as.character(vcflist[,1])[j])
    vcf = fread(as.character(vcflist[,1])[j])
    cellline = as.character(vcflist[,2])[j]
    region = fread(as.character(vcflist[,3])[j])
    #--
    #print(vcf)
    nrows=c()
    for(i in 1:nrow(region)){
      nrows<-c(nrows,intersect(which(dat1[,11]==gsub('chr','',region[i,1])),
                            intersect(which(dat1[,12]>=unlist(region[i,2])),which(dat1[,12]<=unlist(region[i,3])))))
    }
    nrows=unique(nrows)

    targeted_variants=dat1[nrows,]

    targeted_variants$post=as.numeric(targeted_variants$Start)

    #print(unlist(targeted_variants[,12]))
    #big del
    del=setdiff(which(targeted_variants[,12]<targeted_variants[,13]),grep('ins',unlist(targeted_variants[,6])))
    #small del
    del=c(del,grep('del',unlist(targeted_variants[,6])))
    targeted_variants$post[del] = as.numeric(targeted_variants$Start)[del]-1
    #-
    targeted_variants_cellline = targeted_variants[which(targeted_variants[,1]==cellline),]
    #print(cellline)
    # write.table(targeted_variants_cellline, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".targeted.txt",sep=""), append = T, 
    #          quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
    #print(targeted_variants_cellline[,12])
    if(nrow(targeted_variants_cellline)>0){
      raw=unique(targeted_variants_cellline[,12])
      #print(raw)
      unique_raw=c()
      for(m in 1:nrow(raw)){
        exist=which(targeted_variants_cellline[,12]==unlist(raw[m]))
        unique_raw<-c(unique_raw, exist[1])
      }

      unique_raw_mut=targeted_variants_cellline[unique_raw,]
      mut=c()
      for(k in 1:nrow(unique_raw_mut)){
        if(opts$model=='sm'){
          #print(unlist(unique_raw_mut[k,14]))
          exist=which(unlist(unique_raw_mut[k,14])==vcf[,3])
        }
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
      #---
      all_indel = length(grep('ins',unique_raw_mut[,6])) + length(grep('del',unique_raw_mut[,6]))
      all_snp = all - all_indel
      targeted_vcf = unique_raw_mut[which(!is.na(mut)),]
      #print(nrow(targeted_vcf))
      target_indel = length(grep('ins',targeted_vcf[,6])) + length(grep('del',targeted_vcf[,6]))
      target_snp = nrow(targeted_vcf) - target_indel
    }else{
      all = 0
      targeted = 0
      no_targeted = 0
      wrong_mu = nrow(vcf)-targeted
      #---
      all_indel = 0
      all_snp = 0
      targeted_vcf = 0
      target_indel = 0
      target_snp = 0
      mut = NA
      # unique_raw_mut = data.frame(matrix(NA, nrow=2,ncol=15))
    }
    if(targeted==0){
      vcf1=vcf[1,]
      vcf1[1,1:ncol(vcf1)]<-NA
      vcf <-cbind(sample_cell=rep(cellline,nrow(vcf)),vcf)
      vcfsummary<-rbind(vcfsummary, vcf)
      if(nrow(targeted_variants_cellline)>0){
        unique_raw_mut<-cbind(unique_raw_mut, vcf1)
      }
    }else{
      unique_raw_mut<-cbind(unique_raw_mut,vcf[mut,])
      vcf <-cbind(sample_cell=rep(cellline,nrow(vcf)),vcf)
      vcfsummary<-rbind(vcfsummary,vcf[-c(mut[which(!is.na(mut))])])
    }
    wrong_mu_withDB=length(which(vcf[,8]!='-'))-targeted
    vf_dist=quantile(unlist(vcf[which(vcf[,8]!='-'),'VMF']), probs = c(0,0.25,0.5,0.75,1))
    #print(vf_dist)
    vf_dist1=quantile(unlist(vcf[which(vcf[,8]=='-'),'VMF']), probs = c(0,0.25,0.5,0.75,1))
    vf_dist2=quantile(unlist(vcf[mut[which(!is.na(mut))],'VMF']), probs = c(0,0.25,0.5,0.75,1))
    calling_snp = sum(apply(vcf,1,function(x){ifelse(all(nchar(as.character(x[5]))==1,nchar(as.character(x[6]))==1),1,0)}))
    calling_indel=nrow(vcf) - calling_snp
    statis=data.frame(type=c('mut_on_region_in_celllines','targeted', 'snp_on_region_in_celllines', 
                           'indel_on_region_in_celllines', 'targeted snp', 'targeted indel', 
                           'calling snp', 'calling indel',
                           'Minimum VF of targeted', '25% VF of targeted',
                           '50% VF of targeted', '75% VF of targeted',
                           'Maximun VF of targeted',
                           'no_targeted','not_region_in_celllines','not_region_in_celling_within_DB',
                           'Minimum VF of not_region_in_celling_within_DB', '25% VF of not_region_in_celling_within_DB',
                           '50% VF of not_region_in_celling_within_DB', '75% VF of not_region_in_celling_within_DB',
                           'Maximun VF of not_region_in_celling_within_DB','Minimum VF of not_region_in_celling_without_DB', 
                           '25% VF of not_region_in_celling_without_DB',
                           '50% VF of not_region_in_celling_without_DB', '75% VF of not_region_in_celling_without_DB',
                           'Maximun VF of not_region_in_celling_without_DB'),
                    variant=c(all,targeted,all_snp, all_indel,target_snp,target_indel,calling_snp,
                              calling_indel, vf_dist2,no_targeted,wrong_mu,wrong_mu_withDB,vf_dist,vf_dist1))
    all_statis<-cbind(all_statis,statis[,2])
    names<-c(names, as.character(vcf[1,2]))

  #write.table(statis, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".summary.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  #write.table(unique_raw_mut, file=paste(opts$output,"/",as.character(vcf[1,1]),"_",cellline,".mut_targeted.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  #print(unique_raw_mut)
    if(nrow(targeted_variants_cellline)>0){
      unique_raw_mut1= data.frame(sample=as.character(vcf[1,2]),unique_raw_mut)
      unique_raw_muts=rbind(unique_raw_muts, unique_raw_mut1)
    }
  }
  all_statis<-cbind(type=c('mut_on_region_in_celllines','targeted','snp_on_region_in_celllines', 'indel_on_region_in_celllines', 
                         'targeted snp', 'targeted indel', 
                           'calling snp', 'calling indel', 'Minimum VF of targeted', '25% VF of targeted',
                           '50% VF of targeted', '75% VF of targeted',
                           'Maximun VF of targeted','no_targeted','not_region_in_celllines','not_region_in_celling_within_DB',
                         'Minimum VF of not_region_in_celling_within_DB', '25% VF of not_region_in_celling_within_DB',
                           '50% VF of not_region_in_celling_within_DB', '75% VF of not_region_in_celling_within_DB',
                           'Maximun VF of not_region_in_celling_within_DB','Minimum VF of not_region_in_celling_without_DB', 
                         '25% VF of not_region_in_celling_without_DB',
                           '50% VF of not_region_in_celling_without_DB', '75% VF of not_region_in_celling_without_DB',
                           'Maximun VF of not_region_in_celling_without_DB'),all_statis)
  colnames(all_statis)<-c('type',names)
  #write.table(t(all_statis), file=paste(opts$output,"/","All.summary.txt",sep=""), append = T, 
  #            quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = F)
  write.table(unique_raw_muts, file=paste(opts$output,".","All.targeted_summary.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
  write.table(vcfsummary, file=paste(opts$output,".","All.no_celllines_sites_summary.txt",sep=""), append = T, 
              quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)

  cellline_mut <- paste(unique_raw_muts[,2],unique_raw_muts[,15],sep="_")
  unq_cellline_mut<-unique(cellline_mut)
  muts<-matrix(data=NA,nrow=length(unq_cellline_mut), ncol=length(names))
  mutinfos<-c()
  for(i in 1:length(unq_cellline_mut)){
   celline <-unlist(strsplit(unq_cellline_mut[i],'_'))[1]
   mut <-unlist(strsplit(unq_cellline_mut[i],'_'))[2]
   infos=intersect(which(unique_raw_muts[,2]==celline), which(unique_raw_muts[,15]==mut))
   mutinfos<-rbind(mutinfos, unique_raw_muts[infos[1], c(2,12,15,3,4,6,7,11,22,36,37,46,47)])
   for(j in 1:length(names)){
   subda=subset(unique_raw_muts, sample == names[j])
   num<-which(subda[,15]==mut)
   if(length(num)>0){
    if(!is.na(subda[num,25])){
      muts[i,j] <- paste("UMT=",subda[num,25], "  ", "VF=",subda[num,29], sep="")
    }
   }
   }
  }
  colnames(muts)<-names
  muts<-cbind(mutinfos, muts)
  print(nrow(vcfsummary))
  # 读取或创建一个XLSX文件，此步相当于建立一个连接
  filenames = paste(opts$output,'.Celllines.summary.xlsx',sep="")
  xls <- loadWorkbook(filenames,create=TRUE) 
  # 创建工作表
  createSheet(xls,name='summary')
  all_statis1<-data.frame(type=colnames(all_statis), t(all_statis))
  writeWorksheet(xls,all_statis1,'summary',startRow=1,startCol=1,header=FALSE)
  createSheet(xls,name='variants')
  writeWorksheet(xls,muts,'variants',startRow=1,startCol=1,header=TRUE)
  createSheet(xls,name='variants_details')
  writeWorksheet(xls,unique_raw_muts,'variants_details',startRow=1,startCol=1,header=TRUE)
  #createSheet(xls,name='no_celllines_mut')
  #writeWorksheet(xls,vcfsummary,'no_celllines_mut',startRow=1,startCol=1,header=TRUE)
  # 写入数据
 # 存入硬盘，直到此步方才有文档生成
  saveWorkbook(xls)
  print('Done')
}
  