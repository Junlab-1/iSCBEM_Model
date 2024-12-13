# 12.01.2024
# Try nature method human embryo model. Zhao et al., 2024 (DOI: 10.1038/s41592-024-02493-2)]
# modified from https://zenodo.org/records/12189592

rm(list=ls())
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(data.table)
  library(tibble)
  library(miloR)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(igraph)
  library(batchelor)
  library(uwot)
  library(ggpubr)
})
source("~/Code/LJ.function.R")
setwd("/home/DX6/jwulab/S227384/data/seiya/Rrun/NMmodel_run")
options(Seurat.object.assay.version = "v3")
source("main.function.R")
## loading data
# day0
day0data<-Read10X(data.dir="/path/iSCBEMd0_offical/outs/filtered_feature_bc_matrix/")
day0<-CreateSeuratObject(counts = day0data, min.cells = 3, min.features = 200,project = "iSCBEM")
day0@meta.data$tech<-"10xseq"
day0@meta.data$sampletype<-"Day0"
day0@meta.data$datatype<-"iSCBEM_day0"
day0[["percent.mt"]] <- PercentageFeatureSet(day0, pattern = "^MT-")
plotQC(day0)
day0 <- subset(x = day0, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000)
plotQC(day0,"cut")
# rm doublets
day0<-RMdoublets(day0)

# day3
day3data<-Read10X(data.dir="/path/iSCBEMd3_offical/outs/filtered_feature_bc_matrix/")
day3<-CreateSeuratObject(counts = day3data, min.cells = 3, min.features = 200,project = "iSCBEM")
day3@meta.data$tech<-"10xseq"
day3@meta.data$sampletype<-"Day3"
day3@meta.data$datatype<-"iSCBEM_day3"
day3[["percent.mt"]] <- PercentageFeatureSet(day3, pattern = "^MT-")
plotQC(day3)
day3 <- subset(x = day3, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA < 60000)
plotQC(day3,"cut")
# rm doublets
day3<-RMdoublets(day3)

# day5
day5data<-Read10X(data.dir="/path/iSCBEMd5_offical/outs/filtered_feature_bc_matrix/")
day5<-CreateSeuratObject(counts = day5data, min.cells = 3, min.features = 200,project = "iSCBEM")
day5@meta.data$tech<-"10xseq"
day5@meta.data$sampletype<-"Day5"
day5@meta.data$datatype<-"iSCBEM_day5"
day5[["percent.mt"]] <- PercentageFeatureSet(day5, pattern = "^MT-")
plotQC(day5)
day5 <- subset(x = day5, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 10 & nCount_RNA < 70000)
plotQC(day5,"cut")
day5<-RMdoublets(day5)

# day7
library(Seurat)
library(DoubletFinder)
day7data<-Read10X(data.dir="/path/iSCBEMd7_offical/outs/filtered_feature_bc_matrix/")
day7<-CreateSeuratObject(counts = day7data, min.cells = 3, min.features = 200,project = "iSCBEM")
day7@meta.data$tech<-"10xseq"
day7@meta.data$sampletype<-"Day7"
day7@meta.data$species<-"Human"
day7[["percent.mt"]] <- PercentageFeatureSet(day7, pattern = "^MT-")
plotQC(day7)
day7 <- subset(x = day7, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 10 & nCount_RNA < 90000)
plotQC(day7,"cut")
# rm doublets
day7<-RMdoublets(day7)

# parameters
para.list <- list()
para.list$cellGroup <- TRUE ## whether providing group information for cells
para.list$runMiloR <- TRUE  ## whether run miloR aggregation 
para.list$cor.cutoff <- 0.5 ## correlation cutoff

alldata<-merge(x=day0,y = c(day3,day5,day7))
save(alldata,file = "iSCBEM_all.rdata")
# row is gene, column is cell name , gene name must be "CCR7" type not ENS ID
alldata.counts<-alldata@assays$RNA@counts
alldata.counts.meta <-data.frame(alldata@meta.data)
alldata.counts.meta$cell<-rownames(alldata.counts.meta)
alldata.counts.meta$EML <-"query"
alldata.counts.meta$pj  <-"iSCBEM"
alldata.counts.meta$group<-alldata.counts.meta$sampletype
rownames(alldata.counts.meta)<-NULL

# calculated Representative cells
if (para.list$runMiloR) {
  milo_out <- FunMiloCal(alldata.counts,alldata.counts.meta, temp.cal=TRUE)
}else{
  milo_out <- FunMiloCal(alldata.counts,alldata.counts.meta, temp.cal=FALSE)
}

# loading reference data and generating comparable size factor for query datasets
sf_out <- FunCalSF(milo_out)  

# calculating correlation with reference cells
sf_out$query.sce.cor.out <- FunCalCor(sf_out) 
sf_out$query.sce.cor.out.mean <- sf_out$query.sce.cor.out%>% gather(ref_cell,cor,-query_cell)%>%
  group_by(query_cell) %>% top_n(20,cor)%>% summarise(cor_top_mean=mean(cor))

# generate downsampling samples
sf_out$NWIN <- FunNWIN(sf_out)

# MNN calculation for each datasets
mnn.pairs.list <- list()
for (ref_name in c("SPH2016","D3post","Meistermann_2021","CS7","nBGuo","Yan2013")) {
  print(ref_name)
  mnn.pairs.list[[ref_name]] <-FunCalMNN_each(sf_out,ref_name)
}

# generating the 2D projection and 20D embeddings in latent space
predict_out <-  FunProjCal(mnn.pairs.list,
                           query.sce.ob=sf_out$query.sce.ob,
                           query.sce.cor.out=sf_out$query.sce.cor.out,
                           temp.max=sf_out$NWIN$temp.max,
                           D2_umap_model=ref.umap,
                           Dmulti_umap_model=ref_umap_nDim_model,
                           cor.cutoff=para.list$cor.cutoff)

# including raw meta information
predict_out$HS <- milo_out$HS
predict_out$raw.meta <- milo_out$raw.meta
predict_out$query.sce.cor.out.mean <- sf_out$query.sce.cor.out.mean

# cell identities prediction
predict_out <- FunPredAnno(predict_out,cor.cutoff=para.list$cor.cutoff)
predict_out$full.anno %>% group_by(pred_EML) %>% summarise(nCell=n_distinct(query_cell)) %>% arrange(desc(nCell))

# add celltype information to predict_out$umap
predict_out$umap2<-predict_out$umap
predict_out$umap2$pred_EML<-NA
predict_out$umap2$sub_pred_EML<-NA
predict_out$umap2$pred_EML<-predict_out$anno$pred_EML[match(predict_out$umap2$cell,predict_out$anno$query_cell)]
predict_out$umap2$sub_pred_EML<-predict_out$anno$sub_pred_EML[match(predict_out$umap2$cell,predict_out$anno$query_cell)]

timeline<-c("Day0","Day3","Day5","Day7")
pdf(file = "iSCBM_nmmodel3.pdf",width = 4.8,height = 4)
ggplot()+geom_point(predict_out$umap2 %>% filter(!devTime%in%timeline),
                    mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",size=0.1,alpha=1)+
  geom_point(predict_out$umap2 %>% filter(devTime%in%timeline),mapping=aes(x=UMAP_1,y=UMAP_2,color=devTime),
             size=0.2,alpha=0.8)+theme_bw()
dev.off()

# draw ref time line color
timeline_color<-c("Zygote"="#F1F1F1","E1"="#E5D1D2","E2"="#D9B5B6","E2"="#FFEDA0",
                  "E3"="#FED976", "E4"= "#FEB24C",    "E5"="#FD8D3C",     "E6"="#FC4E2A",
                  "E7"="#98393E",     "E8"="#BD0026",     "E9"="#E31A1C", "E10"="#BD0026",
                  "E12"="#8C2129",    "E14"="#800026",    "CS7"="#7F000D")

pdf(file = "iSCBM_nmmodel_with_timecolor_merge_Ref.pdf",width = 4.8,height = 4)
p<-ggplot()+geom_point(predict_out$umap2 %>% filter(!devTime%in%timeline),
                       mapping=aes(x=UMAP_1,y=UMAP_2,color=devTime),size=0.2,alpha=0.8)+
  theme_bw()+scale_color_manual(values = timeline_color)
print(p)
dev.off()

# show celltype
ctcolor<-c("Zygote"="#999999","2–4 cell"="#F39B7E","8 cell"="#F8766C","AdvMes"="#E9842C",
           "Amnion"="#E71F18","Axial_Mes"="#084334","CTB"="#9CA700","DE"="#A3C8DD",
           "Epiblast"="#00B813","Erythroblasts"="#393A79","EVT"="#53B885",
           "ExE_Mes"="#0085ED","HEP"="#00EBEC","Hypoblast"="#7B95FF","ICM"="#BB81FF",
           "Morula"="#FF7E0E","Mesoderm"="#00C1A7","Prelineage"="#C49C93","PriS"="#F862DF",
           "STB"="#F7B6D2","TE"="#00B5ED","YSE"="#EA618E","Ambiguous"="#666666","low_cor"="#666666")
pdf(file = "iSCBM_nmmodel_predcelltype_noambigugous.pdf",width = 4.8,height = 4)
ggplot()+geom_point(predict_out$umap2 %>% filter(!devTime%in%timeline),
                    mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",size=0.1,alpha=1)+
  geom_point(predict_out$umap2 %>% filter(devTime%in%timeline) %>%filter(!pred_EML%in%c("Ambiguous","low_cor")),
             mapping=aes(x=UMAP_1,y=UMAP_2,color=pred_EML),
             size=0.2,alpha=0.8)+theme_bw()+scale_color_manual(values = ctcolor)
dev.off()
timeline<-c("Day0","Day3","Day5","Day7")
for (i in 1:length(timeline)) {
  pdf(file = paste0("iSCBM_nmmodel_predcelltype_noambigugous",timeline[i],".pdf"),width = 4.8,height = 4)
  p<-ggplot()+geom_point(predict_out$umap2 %>% filter(!devTime%in%timeline),
                         mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",size=0.1,alpha=1)+
    geom_point(predict_out$umap2 %>% filter(devTime%in%timeline[i]) %>%filter(!pred_EML%in%c("Ambiguous","low_cor")),
               mapping=aes(x=UMAP_1,y=UMAP_2,color=pred_EML),
               size=0.2,alpha=0.8)+theme_bw()+scale_color_manual(values = ctcolor)
  print(p)
  dev.off()
}

# Transgenes annotation
day0newdata<- Read10X(data.dir="/path/iSCBEMd0_new10/outs/raw_feature_bc_matrix")
d0ng <- day0newdata[,colnames(day0@assays$RNA)]
colnames(d0ng)<-paste0(colnames(d0ng),"_1")
day3newdata<- Read10X(data.dir="/path/iSCBEMd3_new10/outs/raw_feature_bc_matrix")
d3ng <- day3newdata[,colnames(day3@assays$RNA)]
colnames(d3ng)<-paste0(colnames(d3ng),"_2")
day5newdata<- Read10X(data.dir="/path/iSCBEMd5_new10/outs/raw_feature_bc_matrix")
d5ng <- day5newdata[,colnames(day5@assays$RNA)]
colnames(d5ng)<-paste0(colnames(d5ng),"_3")
day7newdata<- Read10X(data.dir="/path/iSCBEMd7_new10/outs/raw_feature_bc_matrix")
d7ng <- day7newdata[,colnames(day7@assays$RNA)]
colnames(d7ng)<-paste0(colnames(d7ng),"_4")

# Transgenes mean exp
new_out<-milo_out$HS
new_out$sample<-NA
new_out$sample[new_out$Hcell%in%colnames(d0ng)]<-"Day0"
new_out$sample[new_out$Hcell%in%colnames(d3ng)]<-"Day3"
new_out$sample[new_out$Hcell%in%colnames(d5ng)]<-"Day5"
new_out$sample[new_out$Hcell%in%colnames(d7ng)]<-"Day7"
allng<-t(cbind(d0ng,d3ng,d5ng,d7ng))
pjlist<-rownames(predict_out$umap3[predict_out$umap3$pj=="iSCBEM",])
predict_out$umap3$new_iGATA6_mean<-0
predict_out$umap3$new_iYAP_mean<-0
predict_out$umap3$new_RTTA_mean<-0
predict_out$umap3$new_PUROR_mean<-0
for (i in pjlist) {
  icell<-new_out$Scell[new_out$Hcell==i]
  predict_out$umap3[i,17:20]<-colMeans(allng[icell,])
}
# do normalization
meanexp<-t(predict_out$umap3[,17:20])
scaling_factor <- 1
normalized_mean <- sweep(
  meanexp, 
  2, 
  ifelse(colSums(meanexp) == 0, 1, colSums(meanexp)),
  FUN = "/"
) * scaling_factor
log_normalized_mean <- log2(normalized_mean + 1)
predict_out$umap3<-cbind(predict_out$umap3,t(log_normalized_mean))
colnames(predict_out$umap3)[21:24]<-c(paste0(colnames(predict_out$umap3[,17:20]),"_lognormal"))
# draw figures
Featureplot_NMmodel(predict_out$umap3,colnames(predict_out$umap3)[21:24],"newgene_featureplot_mean_lognormal2.pdf")

#cell lineage trace for sample cells
predict_out$umap3<-cbind(predict_out$umap3,t(normalized_mean))
colnames(predict_out$umap3)[25:28]<-c(paste0(colnames(predict_out$umap3[,17:20]),"_normal"))
predict_out$umap3$lineage<-"unlabeled.cell"
predict_out$umap3$lineage[!predict_out$umap3$pj=="iSCBEM"]<-"reference.cell"
# need more filter 
predict_out$umap3$lineage[predict_out$umap3$new_PUROR_mean_normal>=0.6]<-"ESC.cell"
predict_out$umap3$lineage[(predict_out$umap3$new_iGATA6_mean_normal>=0.6|
                             (predict_out$umap3$new_iGATA6_mean_normal>=0.3&predict_out$umap3$new_RTTA_mean_normal>=0.3))&
                            predict_out$umap3$new_PUROR_mean_normal<0.1&predict_out$umap3$new_iYAP_mean_normal<0.1&
                            predict_out$umap3$devTime%in%c("Day0","Day3")
]<-"iGATA6.cell"
predict_out$umap3$lineage[(predict_out$umap3$new_iYAP_mean_normal>=0.6|
                             (predict_out$umap3$new_iYAP_mean_normal>=0.3&predict_out$umap3$new_RTTA_mean_normal>=0.3))&
                            predict_out$umap3$new_PUROR_mean_normal<0.1&predict_out$umap3$new_iGATA6_mean_normal<0.1&
                            predict_out$umap3$devTime%in%c("Day0","Day3")
]<-"iYAP.cell"
predict_out$umap3$lineage[predict_out$umap3$new_RTTA_mean_normal>=0.6&
                            predict_out$umap3$devTime%in%c("Day5","Day7")
]<-"iYAP_iGATA6.cell"

lineage_bar<-cbind(as.matrix(predict_out$umap3$devTime),as.matrix(predict_out$umap3$pred_EML),as.matrix(predict_out$umap3$lineage))
lineage_bar<-lineage_bar[!lineage_bar[,3]%in%c("reference.cell","unlabeled.cell"),]
lineage_bar<-lineage_bar[!lineage_bar[,2]%in%c("Ambiguous","low_cor"),]
for (k in unique(lineage_bar[,3])) {
  kbar<-lineage_bar[lineage_bar[,3]==k,]
  for (i in unique(kbar[,1])) {
    sign<-table(kbar[kbar[,1]==i,2])
    singletype_table<-cbind(rep(i,length(sign)),names(sign),as.matrix(sign))
    rownames(singletype_table)<-NULL
    if (i==unique(kbar[,1])[1]) {
      itable_type<-singletype_table
    }else{
      itable_type<-rbind(itable_type,singletype_table)
    }
  }
  itable_type<-cbind(k,itable_type)
  if (k==unique(lineage_bar[,3])[1]) {
    ktable<-itable_type
  }else{
    ktable<-rbind(ktable,itable_type)
  }
}
lineagetable<-ktable
colnames(lineagetable)<-c("LineageType","SampleType","CellType","CellNumber")
write.table(lineagetable,file = "IndexCell_lineageinfo.txt",col.names = T,row.names = F,sep = "\t")
# use html show this figure

# show new gene percentage
ngpm<-predict_out$umap3[!predict_out$umap3$lineage%in%c("reference.cell","unlabeled.cell"),25:29]
library(ggpubr)

ngpm_long <- pivot_longer(
  ngpm,
  cols = starts_with("new_"),
  names_to = "Feature",
  values_to = "Value"
)
ngpm_long$lineage<-factor(ngpm_long$lineage, levels = c("iGATA6.cell","iYAP.cell","iYAP_iGATA6.cell","ESC.cell"))
ngpm_long$Feature<-factor(ngpm_long$Feature, levels = c("new_iGATA6_mean_normal","new_iYAP_mean_normal",
                                                        "new_RTTA_mean_normal","new_PUROR_mean_normal"))
p1<-ggboxplot(
  data = ngpm_long,
  x = "lineage",
  y = "Value",
  color = "Feature",
  outlier.shape = NA
) +
  theme_bw() +
  labs(
    title = "Boxplot of Features by Lineage",
    x = "Lineage",
    y = "normalized gene counts"
  )+theme(legend.position = "bottom")
ggsave(p1,filename = "IndexCell_lineageboxplot.pdf",height = 4,width = 8)
