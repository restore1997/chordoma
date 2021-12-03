rm(list=ls())

library(Seurat)
library(ggplot2)

load("data_total.rds")

data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & nCount_RNA >3 & percent.mt < 25)

pancreas.list <- SplitObject(data, split.by = "data")

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


reference.list <- pancreas.list
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
save(pancreas.integrated,file = "pancreas.integrated.rds")