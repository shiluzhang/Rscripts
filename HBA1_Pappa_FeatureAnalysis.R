## 07-29-2019:
## Step1: generate feature for pappa genes:
options(scipen = 999)
library(data.table)
chr=9
gene="Pappa"  # "HBA1"
cell="Hmec"
genFileandBarplot(cell="K562",chr=16,gene="HBA1")
genFileandBarplot(cell="Hmec",chr=9,gene="Pappa")
genFileandBarplot=function(cell,chr,gene)
{
  ## pairfeatures
  path=paste0("/royfs_write/szhang/e_p_project/hic_for_ripple/Result/Prediction5kb/SQRTVC/Window_Merge_DepNormPerCell_Update/",cell,"/upto1000kb/chr",chr,"/FeatureAnalysisGene")
  print(path)
  setwd(path)
  d=read.table("testset_featureanalysis_path_f2_gene_filter.txt",header=T)
  count=colSums(d[,2:ncol(d)])
  
  ### replace E with R1!!!!
  Features=gsub("_E_", "_R1_", names(count))
  Features=gsub("_P_", "_R2_", Features)
  
  data=data.frame(Features,RulesCount=count)
  data=data[order(data$RulesCount),]
  #data$Features=as.character(data$Features)
  id=grep("Distance",data$Features)
  data1=data[-id,]
  
  feat=do.call(rbind,strsplit(as.character(data1$Features),split="[.]"))
  data2=data.frame(feat,data1$RulesCount)
  names(data2)=c("F1","F2","Importance")
  write.table(data2,file=paste0("testset_featureanalysis_path_f2_gene_filter_",gene,"_counts.txt"),col.names = T,row.names = F,quote = F,sep="\t")
  
  data3=data1[(nrow(data1)-19):nrow(data1),]
  data3$Features <- factor(data3$Features, levels = data3$Features)
  library(ggplot2)
  
  pdf(paste0(cell,"_chr",chr,"_",gene,"_signf_pairfeatures_RulesCount_top20.pdf"),height = 4,width = 3)
  p<-ggplot(data3, aes(x=Features, y=RulesCount, fill=Features)) + geom_bar(stat="identity", fill="blue",position=position_dodge(),size=0.1,colour="black")+
    scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=6,angle=90,colour="black"),
          axis.text.y  = element_text(size=8,colour="black"))+ coord_flip()+ guides(fill=FALSE)
  plot(p)
  dev.off()
  
  
  ### individualfeatures
  d=read.table("testset_featureanalysis_f1_gene_signf.txt",header=T)
  count=colSums(d[,2:ncol(d)])

  Features=gsub("_E_", "_R1_", names(count))
  Features=gsub("_P_", "_R2_", Features)
  
  data=data.frame(Features,RulesCount=count)
  data=data[order(data$RulesCount),]
  #data$Features=as.character(data$Features)
  id=grep("Distance",data$Features)
  data1=data[-id,]
  #data1=data1[(nrow(data1)-99):nrow(data1),]
  data2=data1[(nrow(data1)-19):nrow(data1),]
  data2$Features <- factor(data2$Features, levels = data2$Features)
  pdf(paste0(cell,"_chr",chr,"_",gene,"_signf_individualfeatures_RulesCount_top20.pdf"),height = 4,width = 2.2)
  p<-ggplot(data2, aes(x=Features, y=RulesCount, fill=Features)) + geom_bar(stat="identity",fill="blue", position=position_dodge(),size=0.1,colour="black")+
    scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=6,angle=90,colour="black"),
          axis.text.y  = element_text(size=8,colour="black"))+ coord_flip()+ guides(fill=FALSE)
  plot(p)
  dev.off()
}

genDegreeNodeFile(cell="Hmec",chr=9,gene="Pappa")
genDegreeNodeFile=function(cell,chr,gene)
{
  path=paste0("/royfs_write/szhang/e_p_project/hic_for_ripple/Result/Prediction5kb/SQRTVC/Window_Merge_DepNormPerCell_Update/",cell,"/upto1000kb/chr",chr,"/FeatureAnalysisGene/",gene,"_signf_f2")
  print(path)
  setwd(path)
  d=read.table(paste0(cell,"_",gene,"_degree.txt"),header=T)
  d1=d[order(d$Degree,decreasing = T),]
  d1$Label=as.character(d1$Name)
  id=which(d1$Degree<=10)
  kep=nrow(d1)-length(id)
  print(kep)
  if(kep<5)
  {
    cf=d1$Degree[5]
    id=which(d1$Degree<cf)
  }else if(kep>10)
  {
    cf=d1$Degree[10]
    id=which(d1$Degree<cf)
  }
  d1$Label[id]=""
  print(nrow(d1)-length(id))
  write.table(d1,file=paste0(cell,"_",gene,"_degree1.txt"),col.names = T,row.names = F,quote = F,sep="\t")
}