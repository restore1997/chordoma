rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
load("Monocyte.rda")

Monocyte <- subset(Monocyte,subset = seurat_clusters1%in%c(5,10,13,18))

Monocyte <- RunPCA(Monocyte, features = VariableFeatures(object = Monocyte))
DimPlot(Monocyte, reduction = "pca")
DimHeatmap(Monocyte, dims = 1:15, cells = 500, balanced = TRUE)

Monocyte <- JackStraw(Monocyte, num.replicate = 100)
Monocyte <- ScoreJackStraw(Monocyte, dims = 1:20)

JackStrawPlot(Monocyte, dims = 1:15)
ElbowPlot(Monocyte)
Monocyte <- FindNeighbors(Monocyte, dims = 1:10)
Monocyte <- FindClusters(Monocyte, resolution = 0.2)
Monocyte <- RunUMAP(Monocyte, dims = 1:10)
Monocyte <- RunTSNE(Monocyte, dims = 1:10)

DimPlot(Monocyte, reduction = "umap",label = T,pt.size = 1.6)+labs(title = "umap")
DimPlot(Monocyte, reduction = "tsne",label = T,pt.size = 1.6)+labs(title = "tsne")
DimPlot(Monocyte, reduction = "umap",group.by = "seurat_clusters1",label = T,pt.size = 1.6)
DimPlot(Monocyte, reduction = "umap",group.by = "data",pt.size = 1.6)+labs(title = "data")

DimPlot(Monocyte, reduction = "umap",group.by = "sex",pt.size = 1.6)+labs(title = "sex")
DimPlot(Monocyte, reduction = "umap",group.by = "celltype_hpca",label = T,pt.size = 1.6)

Monocyte$location2 <- "NA"
Monocyte$location2[Monocyte$location%in%c("basis_cranii","intracranial")] <- "skull base"
Monocyte$location2[Monocyte$location%in%"sacrococcyx"] <- "sacrum"
Monocyte$location2[Monocyte$location%in%"spinal_canal"] <- "vertebrae"
DimPlot(Monocyte,reduction = "umap",group.by = "location2",pt.size = 1.6)+labs(title = "location")

Monocyte.markers <- FindAllMarkers(Monocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Monocyte.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
save(Monocyte,Monocyte.markers,file = "Monocyte.rda")

top5 <- Monocyte.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Monocyte, features = top5$gene) #+ NoLegend()

top10 <- Monocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Monocyte, features = top10$gene) #+ NoLegend()


features <- c("CD68","CD163","MSR1","C1QA","C1QB","C1QC")
DoHeatmap(Monocyte, features = features) #+ NoLegend()

##
VlnPlot(Monocyte,features = c("CCL3","CXCL2","CXCL3"),pt.size = 0) ##cluster0,1
VlnPlot(Monocyte,features = c("S100A8","S100A9","S100A12"),pt.size = 0) ##cluster2
VlnPlot(Monocyte,features = c("CD274","PDCD1LG2"),pt.size = 0) ##免疫检查点

write.table(top10,file = "monocyte_top10.txt")
features <- unique(top5$gene)
DotPlot(Monocyte, features = features) + RotatedAxis()


load("Monocyte.rda")

VlnPlot(Monocyte, features = c("CD68"), slot = "counts", log = TRUE,pt.size = 0) ##CD11b
VlnPlot(Monocyte, features = c("ITGAM"), slot = "counts", log = TRUE,pt.size = 0) ##CD11b
VlnPlot(Monocyte, features = c("HLA-C"), slot = "counts", log = TRUE,pt.size = 0) ##MHC
VlnPlot(Monocyte, features = c("MRC1"), slot = "counts", log = TRUE,pt.size = 0) ##CD206

##M1
VlnPlot(Monocyte, features = c("CD68"), slot = "counts", log = TRUE,pt.size = 0)
VlnPlot(Monocyte, features = c("CD80"), slot = "counts", log = TRUE,pt.size = 0)
VlnPlot(Monocyte, features = c("CD86"), slot = "counts", log = TRUE,pt.size = 0)
VlnPlot(Monocyte, features = c("SOCS3"), slot = "counts", log = TRUE,pt.size = 0)

##M2
VlnPlot(Monocyte, features = c("MRC1"), slot = "counts", log = TRUE,pt.size = 0) #CD206
VlnPlot(Monocyte, features = c("MSR1"), slot = "counts", log = TRUE,pt.size = 0)  #CD204
VlnPlot(Monocyte, features = c("CD163"), slot = "counts", log = TRUE,pt.size = 0)
VlnPlot(Monocyte, features = c("SOCS3"), slot = "counts", log = TRUE,pt.size = 0)

DoHeatmap(Monocyte,features = c("MRC1","MSR1","CD163","SOCS3",
                                "SIGLEC1","C1QA","C1QB","C1QC"))

VlnPlot(Monocyte, features = c("HLA-DRA"),pt.size = 0)
VlnPlot(Monocyte, features = c("HLA-DRB1"),pt.size = 0)
VlnPlot(Monocyte, features = c("CD86"),pt.size = 0)
VlnPlot(Monocyte, features = c("ITGAM"),pt.size = 0) ##CD11B
VlnPlot(Monocyte, features = c("MHC1"),pt.size = 0) ##CD11B

Monocyte$MM3 <- "NA"
Monocyte$MM3[Monocyte$seurat_clusters%in%"0"] <- "TAM1"
Monocyte$MM3[Monocyte$seurat_clusters%in%"1"] <- "TAM2"
Monocyte$MM3[Monocyte$seurat_clusters%in%"2"] <- "TAM3"
Monocyte$MM3[Monocyte$seurat_clusters%in%"3"] <- "Monocyte"
Monocyte$MM3[Monocyte$seurat_clusters%in%"4"] <- "TAM4"
Monocyte$MM3[Monocyte$seurat_clusters%in%"5"] <- "TAM5"
DimPlot(Monocyte, reduction = "umap",group.by = "MM3",label = T,pt.size = 1.6)+labs(title = "celltype")

library(monocle)
counts.data <- as(as.matrix(Monocyte@assays$RNA@data), 'sparseMatrix')
pheno.data <- Monocyte@meta.data
feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))

pd <- new("AnnotatedDataFrame", data = pheno.data)
fd <- new("AnnotatedDataFrame", data = feature.data)

HSMM <- newCellDataSet(counts.data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
print(head(pData(HSMM)))
table(HSMM$seurat_clusters)

HSMM <- reduceDimension(HSMM,
                        norm_method = 'log',
                        reduction_method = 'tSNE',
                        verbose = T)
HSMM <- clusterCells(HSMM)
#plot_clusters(HSMM)
ordering_genes <- row.names (subset(Monocyte.markers, p_val_adj < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)


plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "MM3")
plot_cell_trajectory(HSMM, color_by = "data")


save(Monocyte,Monocyte.markers,HSMM,file = "Monocyte.rda")

Monocyte$MM2 <- "NA"
Monocyte$MM2[Monocyte$seurat_clusters%in%"0"] <- "TAM1"
Monocyte$MM2[Monocyte$seurat_clusters%in%"1"] <- "TAM2"
Monocyte$MM2[Monocyte$seurat_clusters%in%"2"] <- "Monocyte1"
Monocyte$MM2[Monocyte$seurat_clusters%in%"3"] <- "TAM3"
Monocyte$MM2[Monocyte$seurat_clusters%in%"4"] <- "Monocyte2"
Monocyte$MM2[Monocyte$seurat_clusters%in%"5"] <- "TAM4"
Monocyte$MM2[Monocyte$seurat_clusters%in%"6"] <- "TAM5"
Monocyte$MM2[Monocyte$seurat_clusters%in%"7"] <- "M1"

DimPlot(Monocyte, reduction = "umap",group.by = "MM3",pt.size = 1.6)
DimPlot(Monocyte, reduction = "umap",pt.size = 1.6,label = T)

FeaturePlot(Monocyte,features = "C1QA")
VlnPlot(Monocyte, features = c("CD68"), slot = "counts", log = TRUE,pt.size = 0)+labs(title = "Macro-CD68") ##Macrophage
VlnPlot(Monocyte, features = c("CD163"), slot = "counts", log = TRUE,pt.size = 0)+labs(title = "TAM-CD163") ##M2
VlnPlot(Monocyte, features = c("C1QA"), slot = "counts", log = TRUE,pt.size = 0)+labs(title = "Macro-C1QA") ##Macrophage