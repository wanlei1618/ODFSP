setwd("/home/shpc_06/ODFSP")
require(switchBox)
library(vcdExtra)
library(caret)
library(forestplot)
library("ktspair")
library(pROC)
library(survcomp)
library(survival)
args(SWAP.Train.KTSP)
require(data.table)
library(reportROC)
library(devtools)
library(ggplot2)
library("easyGgplot2")
ICGC_micro=read.csv("/home/shpc_06/ODFSP/ICGCmicro_inter_pfs1.csv",row.names=1)
ICGC_seq=read.csv("/home/shpc_06/ODFSP/ICGCseq_inter_pfs1.csv",row.names=1)
#View(ICGC)
ICGC=rbind(ICGC_micro,ICGC_seq)
ICGCmicro.grp=read.csv("/home/shpc_06/ODFSP/ICGC/ICGCmicro_pfs1_grp.csv",row.names=1)
ICGCseq.grp=read.csv("/home/shpc_06/ODFSP/ICGC/ICGCmicro_pfs1_grp.csv",row.names=1)
ICGC.grp=rbind(ICGCmicro.grp,ICGCseq.grp)

#####Generating 1000 random models by re-shuffling the labels#####
pred <- list()
sel_pred <- list()
count=0;
b_acc <- vector()
i=1
model <- list()
models_no=1000
count=1
selected_model=list()
set.seed(1987)
sel_b_acc=list()

for(i in 1:models_no){
  #set.seed(1)
  no.low=sample(which(ICGC.grp$pfsrisk==0),80,replace=F)
  index.low=rownames(ICGC.grp)[no.low]
  no.high=sample(which(ICGC.grp$pfsrisk==1),80,replace=F)
  index.high=rownames(ICGC.grp)[no.high]
  
  x1=ICGC[c(index.low,index.high),]
  y_index=c(no.low,no.high) # Selecting the classes of re-sampled samples
  y1=ICGC.grp[y_index,'pfsrisk']

  ### k-TSP models
  zzz=paste('classifier',i,sep="") 
  model[[i]]<- SWAP.KTSP.Train(t(x1), as.factor(y1) )
  selected_model[[count]]=model[[i]]
  selected_model[[count]]$TSPs[,1]=sample(colnames(ICGC),length(selected_model[[count]]$TSPs[,1]))
  selected_model[[count]]$TSPs[,2]=sample(colnames(ICGC),length(selected_model[[count]]$TSPs[,1]))  
  count=count+1;
  print(i)
  
  z=setdiff(1:318,c(no.low,no.high))
  index.test=rownames(ICGC.grp)[z]
  test=ICGC[index.test,]   
  test_grp=ICGC.grp[z,'pfsrisk']
    
  ### Predicting on the test samples
  pred[[i]] <- SWAP.KTSP.Classify(t(test),   selected_model[[i]])   ### Testing on out of bag samples
  print(i+1)
  cc=confusionMatrix(pred[[i]], as.factor(test_grp),  mode = "prec_recall")
  b_acc[i]=as.numeric(cc$byClass)[11]  
}
random_gene_model=selected_model
save(random_gene_model,file="/home/shpc_06/ODFSP/Random_gene_model.RData")
###########################################################

source("validation cohorts formatting.R")
source('normalization.R')
load("Random_gene_model.RData")

predict_ktsp = function(val_mat, val_grp){
  val_pred <- list()
  val_pred_freq<- list()
  auc=vector()
  auc_se=vector()
  
  for(i in 1: length(random_gene_model) ){
    
    print(i)
    
    val_pred[[i]] <- SWAP.KTSP.Classify(t(val_mat), random_gene_model[[i]])  
    
    a=reportROC(val_grp, as.numeric(as.character(val_pred[[i]])),plot = FALSE)
    auc[i]=a$AUC
    auc_se[i]=a$AUC.SE
    
  }
  
  ret_list=list(val_pred,auc,auc_se)
  return(ret_list)  
}

###################################################################################################################################################################
######## Predicting probabliliets using random model for all the validation cohorts
source("/home/shpc_06/ODFSP/normalization.R")
TCGA_list=predict_ktsp(TCGA, TCGA.grp$pfsrisk)
GSE10_list=predict_ktsp(GSE102094,GSE10.grp$pfsrisk)
GSE14_list=predict_ktsp(GSE140082, GSE14.grp$pfsrisk)
GSE32_list=predict_ktsp(GSE32062, GSE32.grp$pfsrisk)
GSE49_list=predict_ktsp(GSE49997,GSE49.grp$pfsrisk)

###################### Calculating meta-estimates for all the 1000 models using all the cohorts

meta_auc=list()
seq_auc=list()
microarray_auc=list()
TCGA_list[[2]]=sapply(TCGA_list[[2]],as.numeric);TCGA_list[[3]]=sapply(TCGA_list[[3]],as.numeric);
GSE10_list[[2]]=sapply(GSE10_list[[2]],as.numeric);GSE10_list[[3]]=sapply(GSE10_list[[3]],as.numeric);
GSE14_list[[2]]=sapply(GSE14_list[[2]],as.numeric);GSE14_list[[3]]=sapply(GSE14_list[[3]],as.numeric);
GSE32_list[[2]]=sapply(GSE32_list[[2]],as.numeric);GSE32_list[[3]]=sapply(GSE32_list[[3]],as.numeric);
GSE49_list[[2]]=sapply(GSE49_list[[3]],as.numeric);GSE49_list[[3]]=sapply(GSE49_list[[3]],as.numeric)

for( i in 1:1000){
  meta_auc[[i]] = combine.est(c(TCGA_list[[2]][i], GSE10_list[[2]][i],GSE14_list[[2]][i], GSE32_list[[2]][i],GSE49_list[[2]][i]),
                              c(TCGA_list[[3]][i], GSE10_list[[3]][i],GSE14_list[[3]][i], GSE32_list[[3]][i],GSE49_list[[3]][i]),na.rm=TRUE)$estimate
  seq_auc[[i]] = combine.est(c(TCGA_list[[2]][i], GSE10_list[[2]][i]), c(TCGA_list[[3]][i], GSE10_list[[3]][i]),na.rm=TRUE)$estimate
  microarray_auc[[i]] = combine.est(c(GSE14_list[[2]][i], GSE32_list[[2]][i],GSE49_list[[2]][i]), c(GSE14_list[[3]][i], GSE32_list[[3]][i],GSE49_list[[3]][i]),na.rm=TRUE)$estimate
  
}

######## Plotting the density plot 
####waiting for more explore
platforms=c( rep("1. Sequencing", 1000),rep("2. Array-based", 1000),rep("3. Overall",1000))
BAC=c(unlist(seq_auc),unlist(microarray_auc), unlist(meta_auc))
dd=data.frame(platforms=platforms, BAC=BAC)

densFindPeak <- function(x){
  td <- density(x)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens],y=td$y[maxDens])
}
seq_peak=densFindPeak(unlist(seq_auc))
micro_peak=densFindPeak(unlist(microarray_auc))
meta_peak=densFindPeak(unlist(meta_auc))
vline.dat =data.frame(platforms=factor(dd$platforms), v1= c(seq_peak$x,micro_peak$x,meta_peak$x))
text.dat=data.frame(platforms=c("1. Sequencing","2. Array-based","3. Overall"),label=c("BAC=0.50","BAC=0.40","BAC=0.50"),
                    x=c(seq_peak$x,micro_peak$x,meta_peak$x),y=c(seq_peak$y,micro_peak$y,meta_peak$y))

pdf("/home/shpc_06/ODFSP/Figure_random_gene_shuffling.pdf")
ggplot2.density(data=dd, xName='BAC', 
                groupName='platforms', legendPosition="top",
                faceting=TRUE, facetingVarNames="platforms",removePanelGrid=TRUE, removePanelBorder=TRUE, showLegend=FALSE, backgroundColor="white",
                fillGroupDensity = TRUE, colorGroupDensityLine = TRUE, 
                xTickLabelFont=c(6, "plain", "black"), yTickLabelFont=c(6, "plain", "black")) + 
  labs(x = "Balanced Accuracy", y="Density")+
  ggtitle("Random gene assignment to k-TSP model") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline( data=vline.dat, aes(xintercept=v1, color= platforms), linetype="dashed",size=0.5) +
  geom_text(data=text.dat, mapping=aes(x=x,y=y,label=label),nudge_y=1,colour='black')

dev.off()



