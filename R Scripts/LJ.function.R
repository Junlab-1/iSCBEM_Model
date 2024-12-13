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
  library(Seurat)
  library(DoubletFinder)
})

# QC function
plotQC<-function(seurat_obj,cutinfo="0"){
  # 提取对象名称
  if (!dir.exists("./QC")) dir.create("./QC", recursive = TRUE)
  obj_name <- deparse(substitute(seurat_obj))
  pdf(file = paste0("./QC/",obj_name,"_count_detail",cutinfo,".pdf"),width = 10,height = 8)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  feature_count.data<-cbind(VlnPlot(seurat_obj, features = "nFeature_RNA")$data,
                            VlnPlot(seurat_obj, features = "nCount_RNA")$data)[,-4]
  fc.plot<-ggplot(feature_count.data,aes(x=nCount_RNA,y=nFeature_RNA))+geom_point()+theme_bw()+
    xlab("All counts per cell")+ylab("Gene numbers per cell")
  ggsave(fc.plot,filename = paste0("./QC/",obj_name,"_CountGenepercent",cutinfo,".pdf"),width = 8,height = 8)
}

# rm doublelets function
RMdoublets<-function(seurat_obj){
  datai_doublet <- seurat_obj
  obj_name <- deparse(substitute(seurat_obj))
  datai_doublet <- NormalizeData(object = datai_doublet, verbose = FALSE)
  datai_doublet <- FindVariableFeatures(object = datai_doublet,selection.method = "vst",nfeatures = 3000, verbose = FALSE)
  datai_doublet <- ScaleData(datai_doublet)
  datai_doublet <- RunPCA(object = datai_doublet,npcs=100)
  datai_doublet <- FindNeighbors(object = datai_doublet, dims = 1:50)
  datai_doublet <- FindClusters(object = datai_doublet, resolution = 0.5)
  datai_doublet <- RunUMAP(object = datai_doublet, dims = 1:50)
  sweep.res.list <- paramSweep_v3(datai_doublet, PCs = 1:10, sct = FALSE)
  for(j in 1:length(sweep.res.list)){
    if(length(sweep.res.list[[j]]$pANN[is.nan(sweep.res.list[[j]]$pANN)]) != 0){
      if(j != 1){
        sweep.res.list[[j]] <- sweep.res.list[[j - 1]]
      }else{
        sweep.res.list[[j]] <- sweep.res.list[[j + 1]]
      }
    }
  }
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk_v <- as.numeric(as.character(bcmvn$pK))
  pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
  # set a expect doublet cells, 10 percent is ok
  nExp_poi <- round(0.03*ncol(datai_doublet))
  #doublet table
  datai_doublet <- doubletFinder_v3(datai_doublet, PCs = 1:10, pN = 0.25, 
                                    pK = pk_good, nExp = nExp_poi, 
                                    reuse.pANN = FALSE, sct = FALSE)
  colnames(datai_doublet@meta.data)[ncol(datai_doublet@meta.data)]="Doublet_type"
  # show doublets cells
  pdf(file = paste("./QC/",obj_name,"_UMAP_doublets.pdf",sep = ""),width = 7,height = 6)
  print(DimPlot(datai_doublet,reduction = "umap",group.by = "Doublet_type"))
  dev.off()
  if (identical(colnames(seurat_obj),colnames(datai_doublet))) {
    seurat_obj@meta.data$Doublet_type=datai_doublet@meta.data$Doublet_type
  }else(
    print("Cell Name is not same!")
  )
  seurat_obj<-subset(x = seurat_obj,subset= Doublet_type == "Singlet")
  return(seurat_obj)
}

# 4 in one self feature plot
library(patchwork)
library(viridis)
library(rlang)
Featureplot_NMmodel <- function(mat,genelist,filename = "NMmodel_featureplot.pdf") {
  p<-list()
  for (i in 1:length(genelist)) {
    # arrange by values
    sorted_data <- mat %>% 
      filter(devTime %in% timeline) %>% 
      arrange(!!sym(genelist[i]),) 
    
    p[[i]]<- ggplot()+geom_point(
      mat %>% filter(!devTime%in%timeline),
      mapping=aes(x=UMAP_1,y=UMAP_2),
      color="grey",
      size=0.5,
      alpha=1)+
      geom_point(
        sorted_data,
        mapping=aes(x=UMAP_1,y=UMAP_2,color=!!sym(genelist[i])),
        size=0.5,
        alpha=0.7)+
      theme_bw()+theme(legend.position = "bottom")+
      scale_color_viridis_c(option = "viridis")
  }
  # 保存为 PDF
  combined_plot <- (p[[1]] + p[[2]]) / (p[[3]] + p[[4]])
  pdf(file = filename, width = 10, height = 12)
  print(combined_plot)
  dev.off()
}

