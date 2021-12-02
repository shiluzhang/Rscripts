## 08/16/2019: Matlab RF importance
### multiplot from ggplot2: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

MatlabRFBarPlot=function(chr,outpath)
{
  library(ggplot2)
  genPlot=function(cell,chr)
  {
    path=paste0("/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/Result/Prediction5kb/SQRTVC/Window_Merge_DepNormPerCell_Update/",cell,"/upto1000kb/chr",chr,"/MatlabRF/")
    setwd(path)
    d=read.table("featureimportance_matlab.txt")
    d$V1=gsub("_E_", "_R1_", d$V1)
    d$V1=gsub("_P_", "_R2_", d$V1)
    d1=d[order(d$V2,decreasing = T),]
    d2=d1[1:20,]
    d2=d2[order(d2$V2),]
    d2$V1 <- factor(d2$V1, levels = d2$V1)
    names(d2)=c("Features","FIMatlab")
    return(d2)
  }
  
  Gm12878=genPlot("Gm12878",chr)
  K562=genPlot("K562",chr)
  Huvec=genPlot("Huvec",chr)
  Hmec=genPlot("Hmec",chr)
  Nhek=genPlot("Nhek",chr)
  
  setwd(outpath)
  
  pdf(paste0("Allcells_chr",chr,"_featureimportance_MatlabRF_top20.pdf"),height = 2.7,width = 8.2)
  p1<-ggplot(Gm12878, aes(x=Features, y=FIMatlab)) + geom_bar(stat="identity", fill="#1E90FF",position=position_dodge(), colour="black")+ggtitle(paste0("Gm12878 chr",chr))+
    scale_y_continuous(name="OOB Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=6,colour="black"),
          axis.text.y  = element_text( size=6,colour="black"),plot.title = element_text(size=8, face="bold"))+ guides(fill=FALSE)+
    coord_flip()
  p2<-ggplot(K562, aes(x=Features, y=FIMatlab)) + geom_bar(stat="identity", fill="#1E90FF",position=position_dodge(), colour="black")+ggtitle(paste0("K562","chr",chr))+
    scale_y_continuous(name="OOB Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=6,colour="black"),
          axis.text.y  = element_text( size=6,colour="black"),plot.title = element_text(size=8, face="bold"))+ guides(fill=FALSE)+
    coord_flip()
  p3<-ggplot(Huvec, aes(x=Features, y=FIMatlab)) + geom_bar(stat="identity",fill="#1E90FF", position=position_dodge(), colour="black")+ggtitle(paste0("Huvec","chr",chr))+
    scale_y_continuous(name="OOB Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=6,colour="black"),
          axis.text.y  = element_text( size=6,colour="black"),plot.title = element_text(size=8, face="bold"))+ guides(fill=FALSE)+
    coord_flip()
  p4<-ggplot(Hmec, aes(x=Features, y=FIMatlab)) + geom_bar(stat="identity",fill="#1E90FF", position=position_dodge(), colour="black")+ggtitle(paste0("Hmec","chr",chr))+
    scale_y_continuous(name="OOB Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=6,colour="black"),
          axis.text.y  = element_text( size=6,colour="black"),plot.title = element_text(size=8, face="bold"))+ guides(fill=FALSE)+
    coord_flip()
  p5<-ggplot(Nhek, aes(x=Features, y=FIMatlab )) + geom_bar(stat="identity",fill="#1E90FF", position=position_dodge(), colour="black")+ggtitle(paste0("Nhek","chr",chr))+
    scale_y_continuous(name="OOB Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text( size=8), axis.title.y = element_blank(),axis.text.x  = element_text(size=6,colour="black"),
          axis.text.y  = element_text( size=6,colour="black"),plot.title = element_text(size=8, face="bold"))+ guides(fill=FALSE)+
    coord_flip()
  multiplot(p1, p2, p3,p4,p5,cols=5)
  dev.off()
}

outpath="/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/NewHiCReg/FeatureAnalysis/MatlabRFFeatImportance"
for(chr in c(14,17,19))
{
  MatlabRFBarPlot(chr,outpath)
}

#############################################################################################################################
## 2. prepare files for heatmap of correlation:
#############################################################################################################################

genDatachr=function(cell,outpath)
{
  for(chr in c(14,17,19))
  {
    path=paste0("/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/Result/Prediction5kb/SQRTVC/Window_Merge_DepNormPerCell_Update/",cell,"/upto1000kb/chr",chr,"/MatlabRF/")
    print(paste0(path,"featureimportance_matlab.txt"))
    setwd(path)
    out=read.table("featureimportance_matlab.txt")
    
    if(chr==14)
    {
      data=out
      names(data)[2]=paste0("chr",chr)
    }
    else{
      data=merge(data,out,by="V1")
      names(data)[ncol(data)]=paste0("chr",chr)
    }
  }
  names(data)[1]="Features"
  setwd(outpath)
  write.table(data,file=paste0(cell,"_featureanalysis_MatlabRFOOB_allpairs_3chroms.txt"),col.names=T,row.names = F,quote = F,sep="\t")
}

outpath="/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/NewHiCReg/FeatureAnalysis/MatlabRFFeatImportance"
for(cell in c("Gm12878","K562","Huvec","Hmec","Nhek"))  #c("Gm12878","K562","Huvec","Hmec","Nhek")
{
  genDatachr(cell,outpath)
}


#############################################################################################################################
# 3. Create a ggheatmap
#############################################################################################################################

outpath="/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/NewHiCReg/FeatureAnalysis/MatlabRFFeatImportance/"
setwd(outpath)
library(reshape2)
library(ggplot2)
cells=c("Gm12878","K562","Huvec","Hmec","Nhek")
p=NULL
for(i in 1:5)
{
  cell=cells[i]
  d=read.table(paste0(cell,"_featureanalysis_MatlabRFOOB_allpairs_3chroms.txt"),header = T)
  corMat=cor(d[,2:4],method = "spearman")
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  cormat=get_upper_tri(corMat)
  cormat=round(cormat,digits=2)
  melted_cormat <- melt(cormat,na.rm=T)
  
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+ggtitle(cell)+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+
    coord_fixed()
  
  p[[i]]<-ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position="none")
}
pdf(paste0("Allcells_featureanalysis_individual_counting_allpairs_3chroms_correlation.pdf"),width = 8.2,height = 2)
multiplot(p[[1]], p[[2]], p[[3]],p[[4]],p[[5]],cols=5)
dev.off()





## Make average across five cells plot:
MakeAverageCrossCellPlot=function(chr)
{
  cells=c("Gm12878","K562","Huvec","Hmec","Nhek")
  for(i in 1:5)
  {
    cell=cells[i]
    path=paste0("/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/Result/Prediction5kb/SQRTVC/Window_Merge_DepNormPerCell_Update/",cell,"/upto1000kb/chr",chr,"/MatlabRF/")
    print(path)
    setwd(path)
    out=read.table("featureimportance_matlab.txt")
    
    if(i==1)
    {
      data=out
      names(data)[2]=cell
    }
    else{
      data=merge(data,out,by="V1")
      names(data)[ncol(data)]=cell
    }
  }
  outpath="/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/NewHiCReg/FeatureAnalysis/MatlabRFFeatImportance"
  setwd(outpath)
  write.table(data,file=paste0("featureanalysis_MatlabRFOOB_allpairs_chr",chr,".txt"),col.names=T,row.names = F,quote = F,sep="\t")

  data=data.frame(Features=data$V1,FIMatlab=rowMeans(data[,2:ncol(data)]))
  d1 = data[order(data$FIMatlab, decreasing = T), ]
  d2 = d1[1:20, ]
  d2 = d2[order(d2$FIMatlab), ]
  d2$Features <- factor(d2$Features, levels = d2$Features)
  pdf(paste0("AverageCrossallcells_chr",chr,"_featureimportance_allhardpairs_Matlab_top20.pdf"),height = 5,width =3.5)
  p<-ggplot(d2, aes(x=Features, y=FIMatlab)) + geom_bar(stat="identity", fill="#1E90FF",position=position_dodge(), colour="black")+ggtitle(paste("allcells","chr",chr))+
    scale_y_continuous(name="OOB Importance")+theme(axis.title.x = element_blank())+
    theme(axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(face="bold",angle=90,size=6,colour="black"),
          axis.text.y  = element_text(face="bold", size=6,colour="black"),plot.title = element_text(size=8, face="bold"))+ guides(fill=FALSE)+
    coord_flip()
  plot(p)
  dev.off()
}

outpath="/mnt/dv/wid/projects3/Roy-enhancer-promoter/HiC-Reg/Gluster/hic_for_ripple/RandomForest/NewHiCReg/FeatureAnalysis/MatlabRFFeatImportance"
setwd(outpath)
for(chr in c(14,17,19))
{
  out=read.table(paste0("featureanalysis_MatlabRFOOB_allpairs_chr",chr,".txt"),header=T)
  out=data.frame(Features=out$V1,FIMatlab=rowMeans(out[,2:ncol(out)]))
  if(chr==14)
  {
    data=out
    names(data)[2]=paste0("chr",chr)
  }
  else{
    data=merge(data,out,by="Features")
    names(data)[ncol(data)]=paste0("chr",chr)
  }
}
write.table(data,file=paste0("featureanalysis_MatlabRFOOB_allpairs_AverageCrossallcells_3chroms.txt"),col.names=T,row.names = F,quote = F,sep="\t")

data=data.frame(Features=data$Features,FIMatlab=rowMeans(data[,2:ncol(data)]))
d1 = data[order(data$FIMatlab, decreasing = T), ]
d2 = d1[1:20, ]
d2 = d2[order(d2$FIMatlab), ]
d2$Features <- factor(d2$Features, levels = d2$Features)
pdf(paste0("AverageCrossallcells_averageacrossChroms_featureimportance_allhardpairs_Matlab_top20.pdf"),height = 5,width =3.5)
p<-ggplot(d2, aes(x=Features, y=FIMatlab)) + geom_bar(stat="identity", fill="#1E90FF",position=position_dodge(), colour="black")+ggtitle("Average allcells 3 chrs")+
  scale_y_continuous(name="OOB Importance")+theme(axis.title.x = element_blank())+
  theme(axis.title.x = element_text(face="bold", size=8), axis.title.y = element_blank(),axis.text.x  = element_text(face="bold",angle=90,size=6,colour="black"),
        axis.text.y  = element_text(face="bold", size=6,colour="black"),plot.title = element_text(size=8, face="bold"))+ guides(fill=FALSE)+
  coord_flip()
plot(p)
dev.off()
