# 清理工作环境 clean enviroment object
rm(list=ls()) 

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
  make_option(c("-g", "--group_name"), type="character",
              help="set the group of the two input files,such as :before, after"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="output directory or prefix [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# 显示输入输出确认是否正确
print(paste("The prefix file is ", opts$prefix_file, sep = ""))
print(paste("The suffix file is ", opts$suffix_file, sep = ""))
print(paste("The output file prefix is ", opts$output, sep = ""))


# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从文件中读取
if (TRUE){
  pre_dat = read.table(opts$prefix_file, sep="\t")
  suf_dat = read.table(opts$suffix_file, sep="\t")
  group_name= unlist(strsplit(opts$group_name,","))
}
print(pre_dat)
# 弹出窗口选择文件
if (FALSE){
  dat = read.table(file.choose(), header=T, row.names = NULL, sep="\t")
}

# 4. 统计与绘图
if (TRUE){
  # 将宽表格将换为长表格(变量位于同列方便操作)
  dat= data.frame(row = pre_dat[,1], A = pre_dat[,2], B = suf_dat[,2])
  dat$row = factor(dat$row,levels= as.character(dat$row))
  colnames(dat) = c("row", group_name)
  #rownames(dat) = rownames(pre_dat)
  dat_melt =melt(dat)
  print(dat_melt)
  p <-ggplot2.barplot(data=dat_melt, xName='row', yName="value",
      groupName='variable', groupColors=c('#999999','#E69F00'),
      position=position_dodge(),
      #background and line colors
      backgroundColor="white", color="black", 
      xtitle="", ytitle="Numbers of Variants", 
      #mainTitle="Total bill\n per time of day",
      removePanelGrid=TRUE,
      removePanelBorder=FALSE,
      axisLine=c(0.5, "solid", "black")
      ) 
}

# 5. 保存图表
if (TRUE){
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(dat, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
  print(paste("The output table is ", opts$output, ".txt",  sep = ""))
  # 保存图片至文件，pdf方便AI修改成出版级图片
  ggsave(file=paste(opts$output,".pdf",sep=""), plot=p, width = 6.5, height = 5)
  ggsave(file=paste(opts$output,".jpg",sep=""), plot=p, width = 6.5, height = 5)
  print(paste("The output figure is ", opts$output, ".pdf",  sep = ""))
}