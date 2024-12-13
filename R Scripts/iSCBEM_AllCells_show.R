# 12.01.2024
# Integration and Visualization of All Cells in Seurat (v5.0.1)
rm(list=ls())
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(data.table)
  library(patchwork)
  library(ggplot2)
  library(harmony)
  library(pheatmap)
  library(plyr)
  library(AnnotationHub)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(GenomicFeatures)
  library(DESeq2)
  library(ggpubr)
  library(DoubletFinder)
  library(ReactomePA)
  library(tidyverse)
  library(biomaRt)
  library(enrichplot)
})
source("LJ.function.R")
setwd("~/path/NMmodel_run")
options(Seurat.object.assay.version = "v3")
alldata<-readRDS("iSCBEM_all.rdata") # load rds data which is saved from NMH_model_run.R
predict_out<-readRDS("Predict_out.rdata")
# load cell types annotated by  NMH_model_run.R
alldata@meta.data$ct_nmmodel<-NA
rownames(predict_out$full.anno)<-predict_out$full.anno$query_cell
ctmat<-predict_out$full.anno[rownames(alldata@meta.data),c("query_cell","pred_EML")]
if (identical(rownames(alldata@meta.data),rownames(ctmat))) {
  alldata@meta.data$ct_nmmodel<-ctmat$pred_EML
}
ctcolor<-c("Zygote"="#999999","2–4 cell"="#F39B7E","8 cell"="#F8766C","AdvMes"="#E9842C",
           "Amnion"="#E71F18","Axial_Mes"="#084334","CTB"="#9CA700","DE"="#A3C8DD",
           "Epiblast"="#00B813","Erythroblasts"="#393A79","EVT"="#53B885",
           "ExE_Mes"="#0085ED","HEP"="#00EBEC","Hypoblast"="#7B95FF","ICM"="#BB81FF",
           "Morula"="#FF7E0E","Mesoderm"="#00C1A7","Prelineage"="#C49C93","PriS"="#F862DF",
           "STB"="#F7B6D2","TE"="#00B5ED","YSE"="#EA618E","Ambiguous"="#666666","low_cor"="#666666",
           "nb_failed"="#666666")
# remove unknown cells draw again
alldata@meta.data$nmctindex<-!alldata$ct_nmmodel%in%c("Ambiguous","low_cor","nb_failed")
alldata2<-subset(alldata,subset= nmctindex==T)
alldata2.list <- SplitObject(alldata2, split.by = "sampletype")
alldata2new.list<-list()
for (i in 1:length(x = alldata2.list)) {
  alldata2new.list[[i]]<-CreateSeuratObject(counts = alldata2.list[[i]]@assays$RNA@counts,
                                            meta.data = alldata2.list[[i]]@meta.data,project = "iSCBEM")
  alldata2new.list[[i]] <- NormalizeData(object = alldata2new.list[[i]], verbose = FALSE)
  alldata2new.list[[i]] <- FindVariableFeatures(object = alldata2new.list[[i]],
                                                selection.method = "vst",  nfeatures = 3000, verbose = FALSE)
  print(i)
}
alldata2.anchors <- FindIntegrationAnchors(object.list = alldata2new.list, dims = 1:100,k.anchor = 5,k.filter = 20)
alldata2 <- IntegrateData(anchorset = alldata2.anchors, dims=1:100)
alldata2 <- ScaleData(alldata2)
alldata2 <- RunPCA(object = alldata2,npcs=200)
alldata2 <- RunHarmony(alldata2, reduction="pca", group.by.vars = "sampletype", reduction.save="harmony")
# Parameter optimization
for (i in 1:20) {
  pcause=i*10
  alldata2 <- FindNeighbors(object = alldata2, dims = 1:pcause)
  alldata2 <- FindClusters(object = alldata2, resolution = 1)
  alldata2 <- RunUMAP(object = alldata2, reduction="harmony", dims = 1:pcause)
  pdf(file = paste("./umap/alldata2_harmony_cluster", pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = alldata2,reduction = "umap", label = T,pt.size = 0.3))
  dev.off()
  pdf(file = paste("./umap/alldata2_harmony_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
  print(DimPlot(object = alldata2,reduction = "umap", group.by = "sampletype",pt.size = .3))
  dev.off()
  pdf(file = paste("./umap/alldata2_harmony_nmcelltype", pcause,".pdf",sep = ""),width = 8,height = 7)
  print(DimPlot(object = alldata2,reduction = "umap", group.by = "ct_nmmodel",pt.size = .1)+
          scale_color_manual(values = ctcolor))
  dev.off()
}
pcause=90
alldata2 <- FindNeighbors(object = alldata2, dims = 1:pcause)
alldata2 <- FindClusters(object = alldata2, resolution = 1)
alldata2 <- RunUMAP(object = alldata2, reduction="harmony", dims = 1:pcause)
pdf(file = paste("./alldata2_harmony_cluster", pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = alldata2,reduction = "umap", label = T,pt.size = 0.3))
dev.off()
pdf(file = paste("./alldata2_harmony_sample",pcause,".pdf",sep = ""),width = 10,height = 8)
print(DimPlot(object = alldata2,reduction = "umap", group.by = "sampletype",pt.size = .3))
dev.off()
pdf(file = paste("./alldata2_harmony_nmcelltype", pcause,".pdf",sep = ""),width = 8,height = 7)
print(DimPlot(object = alldata2,reduction = "umap", group.by = "ct_nmmodel",shape.by = "sampletype",pt.size = .3)+
        scale_color_manual(values = ctcolor)+
        scale_shape_manual(values = c("Day0"=15,"Day3"=18,"Day5"=16,"Day7"=17)))
dev.off()

ctcolor<-c("Zygote"="#999999","2–4 cell"="#F39B7E","8 cell"="#F8766C","AdvMes"="#E9842C",
           "Amnion"="#E71F18","Axial_Mes"="#084334","CTB"="#9CA700","DE"="#A3C8DD",
           "Epiblast"="#00B813","Erythroblasts"="#393A79","EVT"="#53B885",
           "ExE_Mes"="#0085ED","HEP"="#00EBEC","Hypoblast"="#7B95FF","ICM"="#BB81FF",
           "Morula"="#FF7E0E","Mesoderm"="#00C1A7","Prelineage"="#C49C93","PriS"="#F862DF",
           "STB"="#F7B6D2","TE"="#00B5ED","YSE"="#EA618E","Ambiguous"="#666666","low_cor"="#666666",
           "nb_failed"="#666666")
# 1.2 Only day7 all cell UMAP
D7cells<-names(alldata2$sampletype[alldata2$sampletype=="Day7"])
alldata2@meta.data$D7_show<-alldata2$ct_nmmodel
alldata2@meta.data$D7_show[!alldata2$sampletype=="Day7"]="othercells"

pdf(file = "iSCBEM_alldata2_harmony_nmcelltypedim90_onlyD7_gray.pdf",width = 8,height = 7)
print(DimPlot(object = alldata2, reduction = "umap", group.by = "D7_show",pt.size = .1)+
        scale_color_manual(values = c(ctcolor,"othercells"="gray")))
dev.off()
# Only D7
pdf(file = "iSCBEM_alldata2_harmony_nmcelltypedim90_onlyD7.pdf",width = 8,height = 7)
print(DimPlot(object = alldata2,cells = D7cells, reduction = "umap", group.by = "ct_nmmodel",pt.size = .1)+
        scale_color_manual(values = c(ctcolor,"othercells"="gray"))+
        scale_shape_manual(values =17))
dev.off()

# 2. celltype percent
predict_out$full.anno$sample<-NA
predict_out$full.anno$sample[grepl("_1",predict_out$full.anno$query_cell)]<-"Day0"
predict_out$full.anno$sample[grepl("_2",predict_out$full.anno$query_cell)]<-"Day3"
predict_out$full.anno$sample[grepl("_3",predict_out$full.anno$query_cell)]<-"Day5"
predict_out$full.anno$sample[grepl("_4",predict_out$full.anno$query_cell)]<-"Day7"
merge_bar_all<-cbind(as.matrix(predict_out$full.anno$sample),as.matrix(predict_out$full.anno$pred_EML))

for (i in unique(predict_out$full.anno$sample)) {
  sign<-table(merge_bar_all[merge_bar_all[,1]==i,2])
  singletype_table<-cbind(rep(i,length(sign)),names(sign),as.matrix(sign))
  rownames(singletype_table)<-NULL
  if (i==unique(predict_out$full.anno$sample)[1]) {
    table_type<-singletype_table
  }else{
    table_type<-rbind(table_type,singletype_table)
  }
}
alltable<-table_type

colnames(alltable)<-c("Stage","Celltype","Cellnumber")
table_sample_type<-as.data.frame(alltable)
table_sample_type$Stage<-factor(table_sample_type$Stage)
table_sample_type$Celltype<-as.factor(table_sample_type$Celltype)                                                             ##
table_sample_type$Cellnumber<-as.numeric(as.matrix(table_sample_type$Cellnumber))                                                 ##
library(plyr)
alltable_percent<-ddply(table_sample_type,"Stage",transform,percent_weight=Cellnumber/sum(Cellnumber)*100)
alltable_percent$Celltype2 <- alltable_percent$Celltype
alltable_percent$Celltype2 <- plyr::mapvalues(x=alltable_percent$Celltype2, 
                                              from=c("CTB","TE","STB",
                                                     "Amnion","Epiblast","ICM","Prelineage",
                                                     "AdvMes","ExE_Mes","Mesoderm","PriS",
                                                     "YSE","DE","HEP","Hypoblast",
                                                     "Ambiguous","low_cor","nb_failed"),
                                              to=c("Trophectoderm","Trophectoderm","Trophectoderm",
                                                   rep("Ectoderm",4),
                                                   rep("Mesoderm2",4),
                                                   rep("Endoderm",4),
                                                   rep("Unknown",3)))
alltable_percent$Celltype2<-factor(alltable_percent$Celltype2,levels = c("Ectoderm","Mesoderm2","Endoderm","Trophectoderm","Unknown"))
write.csv(alltable_percent,file = "NMmodel_iSCBEM_CelltypePercent.csv",row.names = F,col.names = T)
alltable_percent_D7<-alltable_percent[alltable_percent$Stage=="Day7",]
p_bar2<-ggplot(alltable_percent_D7,aes(x=Stage,y=percent_weight,fill=Celltype))+geom_bar(stat = "identity",width= 0.6)+
  theme_bw()+ scale_fill_manual(values = ctcolor)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
  facet_grid(rows = vars(Celltype2),scales = "free_y",space = "free")+
  coord_cartesian(clip = 'off')
ggsave("all_celltype_percent_D7.pdf", p_bar2, width=4 ,height=6)

# 4. markergenes heatmap
TopMatrix<-alldata2@assays$RNA@data[,]
celltype_specieslist<-unique(alldata2$ct_nmmodel)
for (i in 1:length(celltype_specieslist)) {
  icell<-TopMatrix[,alldata2$ct_nmmodel==celltype_specieslist[i]]
  if (i==1) {
    MeanMatrix<-rowMeans(icell)
  }else{
    MeanMatrix<-cbind(MeanMatrix,rowMeans(icell))
  }
}
colnames(MeanMatrix)<-celltype_specieslist
mgl<-read.table("human_model_top5_genes_order.txt",header = F)
mgl<-mgl[-c(1:10),] #Remove cell types that contain fewer than 100 cells.
MarkerMeanMatrix<-MeanMatrix[intersect(mgl,rownames(MeanMatrix)),!colnames(MeanMatrix)%in%c("Ambiguous","low_cor","nb_failed")]
MMM_order<-MarkerMeanMatrix[,c("Epiblast","PriS",
                               "Amnion","Mesoderm","AdvMes","ExE_Mes",
                               "Hypoblast","DE","YSE","TE","CTB","STB","HEP")]
annotation_col<-data.frame(Groups=colnames(MMM_order))
rownames(annotation_col)<-colnames(MMM_order)
annotation_col$Groups<-as.factor(annotation_col$Groups)
mycolorsp <- list(Groups = ctcolor)
library(pheatmap)
pdf(file = "Alldata_Markergenemean_heatmap_ordered_rmcelltype100.pdf",width = 8,height = 10)
pheatmap(mat = MMM_order,scale = "row",cluster_cols = F,cluster_rows = F,
         annotation_col = annotation_col, annotation_colors = mycolorsp,
         color = colorRampPalette(c("blue","white","red"))(100),border_color = "white")
dev.off()

# 4.2 changed heatmap
# add new genes to alldata2@meta.data
day0newdata<- Read10X(data.dir="/path/iSCBEMd0_new10/outs/raw_feature_bc_matrix")
colnames(day0newdata)<-paste0(colnames(day0newdata),"_1")
d0ng <- day0newdata[,names(alldata2$sampletype[alldata2$sampletype=="Day0"])]

day3newdata<- Read10X(data.dir="/path/iSCBEMd3_new10/outs/raw_feature_bc_matrix")
colnames(day3newdata)<-paste0(colnames(day3newdata),"_2")
d3ng <- day3newdata[,names(alldata2$sampletype[alldata2$sampletype=="Day3"])]

day5newdata<- Read10X(data.dir="/path/iSCBEMd5_new10/outs/raw_feature_bc_matrix")
colnames(day5newdata)<-paste0(colnames(day5newdata),"_3")
d5ng <- day5newdata[,names(alldata2$sampletype[alldata2$sampletype=="Day5"])]

day7newdata<- Read10X(data.dir="/path/iSCBEMd7_new10/outs/raw_feature_bc_matrix")
colnames(day7newdata)<-paste0(colnames(day7newdata),"_4")
d7ng <- day7newdata[,names(alldata2$sampletype[alldata2$sampletype=="Day7"])]

allnew<-cbind(d0ng,d3ng,d5ng,d7ng)
scaling_factor <- 1
allnew_mean <- sweep(
  allnew, 
  2, 
  ifelse(colSums(allnew) == 0, 1, colSums(allnew)), # 避免除以 0
  FUN = "/"
) * scaling_factor
rownames(allnew_mean)<-c("newiGATA6_norm","newiYAP_norm","newRTTA_norm","newPUROR_norm")
if (identical(colnames(allnew),rownames(alldata2@meta.data))) {
  alldata2@meta.data<-cbind(alldata2@meta.data,t(allnew))
}
if (identical(colnames(allnew_mean),rownames(alldata2@meta.data))) {
  alldata2@meta.data<-cbind(alldata2@meta.data,t(allnew_mean))
}
# PUROR positive heatmap
alldata2@meta.data$ct_NMLT<-"unknown"
indexPUROR<-which(alldata2$newPUROR_norm>=0.6)
alldata2@meta.data$ct_NMLT[indexPUROR]<-paste0(alldata2$ct_nmmodel[indexPUROR],"_PUROR")
# RTTA positive heatmap
index1       <-which(alldata2$newRTTA_norm>=0.6)
index2_iYAP  <-which(alldata2$newRTTA_norm>=0.3&alldata2$newiYAP_norm>=0.3 & alldata2$sampletype%in%c("Day0","Day3"))
index2_iGATA6<-which(alldata2$newRTTA_norm>=0.3&alldata2$newiGATA6_norm>=0.3&alldata2$sampletype%in%c("Day0","Day3"))
indexRTTA<-c(index1,index2_iGATA6,index2_iYAP)
alldata2@meta.data$ct_NMLT[indexRTTA]<-paste0(alldata2$ct_nmmodel[indexRTTA],"_RTTA")

TopMatrix<-alldata2@assays$RNA@data[,]
celltype_specieslist<-unique(alldata2$ct_NMLT)
for (i in 1:length(celltype_specieslist)) {
  icell<-TopMatrix[,alldata2$ct_NMLT==celltype_specieslist[i]]
  if (i==1) {
    MeanMatrix<-rowMeans(icell)
  }else{
    MeanMatrix<-cbind(MeanMatrix,rowMeans(icell))
  }
}
colnames(MeanMatrix)<-celltype_specieslist
mgl<-read.table("human_model_top5_genes_order.txt",header = F)
MarkerMeanMatrix<-MeanMatrix[intersect(mgl[,1],rownames(MeanMatrix)),!colnames(MeanMatrix)=="unknown"]
# 4.2.1 RTTA heatmap
MMM_order_RTTA<-MarkerMeanMatrix[,paste0(c("Prelineage","ICM","Epiblast","PriS",
                                           "Amnion","Mesoderm","AdvMes","ExE_Mes",
                                           "Hypoblast","DE","YSE","TE","CTB","STB","HEP"),"_RTTA")]
annotation_col<-data.frame(Groups=colnames(MMM_order_RTTA))
rownames(annotation_col)<-colnames(MMM_order_RTTA)
annotation_col$Groups<-as.factor(annotation_col$Groups)

ctcolorrtta<-ctcolor
names(ctcolorrtta)<-paste0(names(ctcolor),"_RTTA")
mycolorsp <- list(Groups = ctcolorrtta)
library(pheatmap)
pdf(file = "Alldata_Markergenemean_heatmap_ordered_RTTA.pdf",width = 8,height = 10)
pheatmap(mat = MMM_order_RTTA,scale = "row",cluster_cols = F,cluster_rows = F,
         annotation_col = annotation_col, annotation_colors = mycolorsp,
         color = colorRampPalette(c("blue","white","red"))(100),border_color = "white")
dev.off()
# 4.2.2 PUROR heatmap
MMM_order_PUROR<-MarkerMeanMatrix[,paste0(c("Prelineage","ICM","Epiblast","PriS",
                                            "Amnion","Mesoderm","AdvMes","ExE_Mes",
                                            "Hypoblast","DE","YSE","TE","CTB","STB","HEP"),"_PUROR")]
annotation_col<-data.frame(Groups=colnames(MMM_order_PUROR))
rownames(annotation_col)<-colnames(MMM_order_PUROR)
annotation_col$Groups<-as.factor(annotation_col$Groups)
ctcolorpuror<-ctcolor
names(ctcolorpuror)<-paste0(names(ctcolor),"_PUROR")
mycolorsp <- list(Groups = ctcolorpuror)
library(pheatmap)
pdf(file = "Alldata_Markergenemean_heatmap_ordered_PUROR.pdf",width = 8,height = 10)
pheatmap(mat = MMM_order_PUROR,scale = "row",cluster_cols = F,cluster_rows = F,
         annotation_col = annotation_col, annotation_colors = mycolorsp,
         color = colorRampPalette(c("blue","white","red"))(100),border_color = "white")
dev.off()

# 5 marker genes feature plots.
markerall<-read.table("Public_marker_gene_list.txt",header = T,sep = "\t")
for (i in 1:ncol(markerall)) {
  imarker<-cbind(markerall[,i],rep(colnames(markerall)[i],length(markerall[,i])))
  if (i==1) {
    markers<-imarker
  }else{
    markers<-rbind(markers,imarker)
  }
}
markers<-na.omit(markers)
markers<-as.matrix(markers)
markers<-unique(markers)
markers<-markers[markers[,1]!="",]
library(viridis)
alldata2_fp<-alldata2
DefaultAssay(alldata2_fp) <- "RNA"
markers<-markers[markers[,1]%in%rownames(alldata2_fp),]
for (i in 1:nrow(markers)) {
  pdf(file = paste("GeneUMAP/Alldata_FeaturePlot",markers[i,1],markers[i,2],".pdf",sep = "_"),width =5.5 ,height = 5)
  print(FeaturePlot(alldata2_fp,slot = "data", features = markers[i,1],
                    ncol = 1,alpha=1,order = T)+scale_colour_gradientn(colours =viridis_pal()(100))
  )
  dev.off()
}
