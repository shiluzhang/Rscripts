## 08/18/2019: SupFig 9: Feature ranking based on feature usage counting for all pairs vs top 5% lowest error pairs for all five cell lines. 

genRulesCount=function(cell)
{
  path=paste0("/royfs_write/szhang/e_p_project/hic_for_ripple/Result/Prediction5kb/SQRTVC/Window_Merge_DepNormPerCell_Update/",cell,"/upto1000kb/chr17/FeatureAnalysisAllValue/HardPairs/")
  print(path)
  d=read.table(paste0(path,"/testset_featureanalysis_path_f1_all.txt"),header=T,stringsAsFactors = F)
  RulesCount=colSums(d[,2:ncol(d)])
  ### replace E with R1!!!!
  Features=gsub("[.]", "", names(RulesCount))
  Features=gsub("_E_", "_R1_", Features)
  Features=gsub("_P_", "_R2_", Features)
  
  out1=data.frame(Features,RulesCount)
  write.table(out1,file=paste0(path,"featureanalysis_individual_counting_all.txt"),col.names=T,row.names = F,quote = F,sep="\t")
  
  d=read.table(paste0(path,"testset_featureanalysis_path_f1_top5low.txt"),header=T,stringsAsFactors = F)
  RulesCounttop5=colSums(d[,2:ncol(d)])
  Features=gsub("[.]", "", names(RulesCounttop5))
  Features=gsub("_E_", "_R1_", Features)
  Features=gsub("_P_", "_R2_", Features)
  out2=data.frame(Features,RulesCounttop5)
  write.table(out2,file=paste0(path,"featureanalysis_individual_counting_top5low.txt"),col.names=T,row.names = F,quote = F,sep="\t")
  
  #data=merge(out1,out2,by="Features")
  #write.table(data,file=paste0(path,cell,"_chr17_featureanalysis_all.txt"),col.names=T,row.names = F,quote = F,sep="\t")
  #d=read.table(paste0(cell,"_chr17_featureanalysis_all.txt"),header=T)

  d1=data.frame(out1,"All",cell)
  d2=data.frame(out2,"Top5LowErr",cell)
  names(d1)=c("Features","RulesCount","Pairs","cell")
  names(d2)=c("Features","RulesCount","Pairs","cell")
  d1$RulesCount=d1$RulesCount/sum(d1$RulesCount)
  d2$RulesCount=d2$RulesCount/sum(d2$RulesCount)
  d2=d2[order(d2$RulesCount),]
  data=rbind(d1,d2)
  data$Features <- factor(data$Features, levels = d2$Features)
  return(data)
}

Gm12878=genRulesCount("Gm12878")
K562=genRulesCount("K562")
Huvec=genRulesCount("Huvec")
Hmec=genRulesCount("Hmec")
Nhek=genRulesCount("Nhek")

outpath="/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/NewHiCReg/FeatureAnalysis/RFRulesCounting/"
setwd(outpath)
font=4
pdf(paste0("Allcells_chr17_featureimportance_RulesCount_AllvsTop5_sorted_all.pdf"),height = 11,width = 8.5)
p1<-ggplot(Gm12878, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Gm12878"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p2<-ggplot(K562, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("K562"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p3<-ggplot(Huvec, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Huvec"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p4<-ggplot(Hmec, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Hmec"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text(size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p5<-ggplot(Nhek, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Nhek"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
multiplot(p1, p2, p3,p4,p5,cols=5)
dev.off()

###########################################################################################
## SupFig 10: top 20
###########################################################################################

genDatatop20=function(cell)
{
  path=paste0("/royfs_write/szhang/e_p_project/hic_for_ripple/Result/Prediction5kb/SQRTVC/Window_Merge_DepNormPerCell_Update/",cell,"/upto1000kb/chr17/FeatureAnalysisAllValue/HardPairs/")
  print(path)
  out1=read.table(paste0(path,"featureanalysis_individual_counting_all.txt"),header=T)
  out2=read.table(paste0(path,"featureanalysis_individual_counting_top5low.txt"),header=T)
  d1=data.frame(out1,"All",cell)
  d2=data.frame(out2,"Top5LowErr",cell)
  names(d1)=c("Features","RulesCount","Pairs","cell")
  names(d2)=c("Features","RulesCount","Pairs","cell")
  d1$RulesCount=d1$RulesCount/sum(d1$RulesCount)
  d2$RulesCount=d2$RulesCount/sum(d2$RulesCount)
  id=order(d2$RulesCount)
  d2=d2[id,]
  d1=d1[id,]
  data=rbind(d1[192:211,],d2[192:211,])
  data$Features <- factor(data$Features, levels = d2$Features[192:211])
  
  # id1=order(d1$RulesCount)
  # dd1=d1[id1,]
  # print(length(intersect(dd1[192:211,"Features"],d2[192:211,"Features"])))
  return(data)
}
outpath="/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/NewHiCReg/FeatureAnalysis/RFRulesCounting/"
setwd(outpath)
Gm12878=genDatatop20("Gm12878")
K562=genDatatop20("K562")
Huvec=genDatatop20("Huvec")
Hmec=genDatatop20("Hmec")
Nhek=genDatatop20("Nhek")

font=7
pdf(paste0("Allcells_chr17_featureimportance_RulesCount_AllvsTop5_sorted_top20.pdf"),height = 6,width = 8.3)
p1<-ggplot(Gm12878, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Gm12878"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(angle=90,size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p2<-ggplot(K562, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("K562"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(angle=90,size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p3<-ggplot(Huvec, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Huvec"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(angle=90,size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p4<-ggplot(Hmec, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Hmec"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(angle=90,size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
p5<-ggplot(Nhek, aes(x=Features, y=RulesCount, fill=Pairs)) + geom_bar(stat="identity", position=position_dodge(),size=0.1,colour="black")+
  scale_y_continuous(name="Importance")+theme(axis.title.x = element_blank())+guides(fill=FALSE)+ggtitle(paste("Nhek"))+
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(angle=90,size=5,colour="black"),
        axis.text.y  = element_text(size=font,colour="black"))+ coord_flip()
multiplot(p1, p2, p3,p4,p5,cols=5)
dev.off()


