library(piano)
library(dplyr)
library(tidyverse)
setwd('/home/shpc_06/ODFSP')
gene_selected=vector()
for (i in 1:length(selected_model)){
  gene_selected=append(gene_selected,selected_model[[i]]$TSPs) %>% unique()
  }
gene_selected=gene_selected %>% substr(2,1000000L)


### GO 
genesets= loadGSC("/home/shpc_06/ODFSP/c5.go.bp.v7.5.1.entrez.gmt")
results_hallmark=runGSAhyper(gene_selected,gsc= genesets,  adjMethod="fdr")
results_hallmark$resTab[which(results_hallmark$resTab[,2]<0.05),]
write.table(results$resTab[which(results$resTab[,1]<0.05),],"/home/shpc_06/ODFSP/bp.txt")#pvalue

genesets=loadGSC("/home/shpc_06/ODFSP/c5.go.cc.v7.5.1.entrez.gmt", type="auto")
results=runGSAhyper(gene_selected,gsc= genesets,  adjMethod="fdr")
results$resTab[which(results$resTab[,2]<0.05),]
write.table(results$resTab[which(results$resTab[,1]<0.05),],"/home/shpc_06/ODFSP/cc.txt")

genesets=loadGSC("/home/shpc_06/ODFSP/c5.go.mf.v7.5.1.entrez.gmt", type="auto")
results=runGSAhyper(gene_selected,gsc= genesets,  adjMethod="fdr")
results$resTab[which(results$resTab[,2]<0.1),]
write.table(results$resTab[which(results$resTab[,1]<0.05),],"/home/shpc_06/ODFSP/mf.txt")


##### HALLMARK PATHWAYS
genesets=loadGSC("/home/shpc_06/ODFSP/c1.all.v7.5.1.entrez.gmt", type="auto")
results=runGSAhyper(gene_selected,gsc= genesets,  adjMethod="fdr")
results_p=results$resTab
#View(results_p[which(results_p[,2]<0.05),])#p-adjust
length(results_p[which(results_p[,2]<0.05),])
results$resTab[which(results$resTab[,1]<0.05),]#p
write.table(results$resTab[which(results$resTab[,2]<0.05),],"/home/shpc_06/ODFSP/hallmark.txt")

library("readxl")
library('clusterProfiler')
library("org.Hs.eg.db")
library(ggplot2)
library(tidyr)
library(R.utils)###OC_ktsp_all_0
R.utils::setOption("clusterProfiler.download.method",'auto')
options(clusterProfiler.download.method = "wget")
go=enrichGO(gene_selected, 
            OrgDb = org.Hs.eg.db, 
            ont='ALL',
            pAdjustMethod = 'BH',
            pvalueCutoff = 0.05, 
            qvalueCutoff = 0.2,
            keyType = 'ENTREZID')

go_df<-data.frame(ID=go$ID,
                  Description=go$Description,
                  Count=go$Count,
                  p.adjust=go$p.adjust)
## numbers as data on x axis
go_df$number <- factor(rev(1:nrow(go_df)))
## add GO terms description names
labels <- go_df$Description
names(labels) = rev(1:nrow(go_df))
CPCOLS <- "#eb6bac"
## colors for bar // green, blue, orange
p1 <- ggplot(data=go_df[1:20,], aes(x=number, y=Count,fill=p.adjust)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_gradient(low = "red",high = "blue",space='Lab') + theme_classic() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="black"),legend.position = 'bottom',
        legend.key.height=unit(0.4,'cm'),legend.key.width=unit(1,'cm'),aspect.ratio = 3,
        legend.text=element_text(size=6),legend.justification = c(0.8,1),
        legend.margin=margin(-0.2, 0, 0, 0, "cm"),legend.title=element_text(size=10,vjust=0.8))
plot(p1)
ggsave("go_all.pdf", units="in",width = 12, height = 8,dpi=400)

### KEGG ####
kegg=enrichKEGG(gene_selected, 
                    organism = 'hsa',
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,
                    keyType = 'kegg')

dotplot(kegg,showCategory=10)+theme_classic()+theme(plot.title=element_text(size=14))
ggsave('/home/shpc_06/ODFSP/KEGG.pdf',width=10, height=8,dpi=400)
#goplot_grid=barplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
ggplot(go, aes(Count, Description, fill=ONTOLOGY), split='ONTOLOGY')+ geom_col()+
  facet_grid(ONTOLOGY~., scales="free_y")+theme_classic()+ 
  theme(axis.text=element_text(face = "bold", color="black"),legend.position = 'right',
        legend.key.height=unit(0.4,'cm'),legend.key.width=unit(1,'cm'),aspect.ratio = 3,
        legend.text=element_text(size=8),legend.title=element_text(size=10)) 
ggsave("go_grid.pdf",units="in", width=12, height=8,dpi=400)

### GSEA
TCGA_lr=TCGA[which(TCGA.grp$pfsrisk==1),colnames(TCGA) %>% substr(2,1000000L)%in%gene_selected]
TCGA_hr=TCGA[which(TCGA.grp$pfsrisk==0),colnames(TCGA) %>% substr(2,1000000L)%in%gene_selected]
gsa_list=log2(colMeans(TCGA_hr)/colMeans(TCGA_lr))
names(gsa_list)=gene_selected
gsa_list=sort(gsa_list,decreasing=T)
genesets=read.gmt("/home/shpc_06/ODFSP/msigdb.v7.5.1.entrez.gmt")
GSEA(gsa_list, minGSSize = 10,maxGSSize = 1000,TERM2GENE = genesets)
        
