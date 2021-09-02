rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
load("Malignant.rda")
Malignant$seurat_clusters1 <- Malignant$seurat_clusters

Malignant <- RunPCA(Malignant, features = VariableFeatures(object = Malignant))
DimPlot(Malignant, reduction = "pca")
DimHeatmap(Malignant, dims = 1:15, cells = 500, balanced = TRUE)

Malignant <- JackStraw(Malignant, num.replicate = 100)
Malignant <- ScoreJackStraw(Malignant, dims = 1:20)

JackStrawPlot(Malignant, dims = 1:15)
ElbowPlot(Malignant)
Malignant <- FindNeighbors(Malignant, dims = 1:10)
Malignant <- FindClusters(Malignant, resolution = 0.2)
Malignant <- RunUMAP(Malignant, dims = 1:10)
Malignant <- RunTSNE(Malignant, dims = 1:10)

DimPlot(Malignant, reduction = "umap",label = T,pt.size = 1.6)
DimPlot(Malignant, reduction = "umap",group.by = "seurat_clusters1",label = T,pt.size = 2)
DimPlot(Malignant, reduction = "umap",group.by = "data",pt.size = 1.6)+labs(title = "data")
Malignant$location2 <- "NA"
Malignant$location2[Malignant$location%in%c("basis_cranii","intracranial")] <- "skull base"
Malignant$location2[Malignant$location%in%"sacrococcyx"] <- "sacrum"
Malignant$location2[Malignant$location%in%"spinal_canal"] <- "vertebrae"
DimPlot(Malignant, reduction = "umap",group.by = "location2",pt.size = 1.6)+labs(title = "location")

Malignant.markers <- FindAllMarkers(Malignant, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Malignant.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
save(Malignant,Malignant.markers,file = "Malignant.rda")

top5 <- Malignant.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Malignant, features = top5$gene) #+ NoLegend()

top10 <- Malignant.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Malignant, features = top10$gene) #+ NoLegend()


VlnPlot(Malignant, features = c("THY1"), slot = "counts", log = TRUE,pt.size = 0) ##C0
VlnPlot(Malignant, features = c("SPP1"), slot = "counts", log = TRUE,pt.size = 0) ##C1
VlnPlot(Malignant, features = c("IGLC3"), slot = "counts", log = TRUE) ##C2
VlnPlot(Malignant, features = c("ID4"), slot = "counts", log = TRUE,pt.size = 0) ##C3
VlnPlot(Malignant, features = c("MAD2L1"), slot = "counts", log = TRUE,pt.size = 0) ##C4
VlnPlot(Malignant, features = c("HLA-DRA"), slot = "counts", log = TRUE,pt.size = 0) ##C5

VlnPlot(Malignant, features = c("CDKN2A"), slot = "counts", log = TRUE,pt.size = 0) 
VlnPlot(Malignant, features = c("SMARCB1"), slot = "counts", log = TRUE,pt.size = 0) ##C5

VlnPlot(Malignant, features = c("EGFR"), slot = "counts", log = TRUE,pt.size = 0) 
VlnPlot(Malignant, features = c("ERBB2"), slot = "counts", log = TRUE,pt.size = 0) ##HER2
VlnPlot(Malignant, features = c("PDGFRB"), slot = "counts", log = TRUE,pt.size = 0) ##PDGFR
VlnPlot(Malignant, features = c("KDR"), slot = "counts", log = TRUE,pt.size = 0) ##VEGFR
VlnPlot(Malignant, features = c("EZH2"), slot = "counts", log = TRUE,pt.size = 0) ##VEGFR

VlnPlot(Malignant, features = c("CD274"), slot = "counts", log = TRUE,pt.size = 0) ##PD-L1
VlnPlot(Malignant, features = c("PDCD1"), slot = "counts", log = TRUE,pt.size = 0) ##PD-1
VlnPlot(Malignant, features = c("CTLA4"), slot = "counts", log = TRUE,pt.size = 0) ##CTLA4
VlnPlot(Malignant, features = c("HAVCR2"), slot = "counts", log = TRUE,pt.size = 0) ##TIMI3
VlnPlot(Malignant, features = c("TIGIT"), slot = "counts", log = TRUE,pt.size = 0) ##TIGIT



library(monocle)
counts.data <- as(as.matrix(Malignant@assays$RNA@data), 'sparseMatrix')
pheno.data <- Malignant@meta.data
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
ordering_genes <- row.names (subset(Malignant.markers, p_val_adj < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)

load("HSMM.rda")

plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "data")
plot_cell_trajectory(HSMM, color_by = "CC")

plot_cell_trajectory(HSMM, color_by = "MM2")

DimPlot(Malignant, reduction = "umap",label = T,group.by = "CC")


save(Malignant,Malignant.markers,HSMM,total.markers2,file = "Malignant.rda")

Malignant$MM <- Malignant$CC
Malignant$MM[Malignant$MM%in%"CC1"] <- "M1"
Malignant$MM[Malignant$MM%in%"CC2"] <- "M2"
Malignant$MM[Malignant$MM%in%"CC3"] <- "M3"
Malignant$MM[Malignant$MM%in%"CC4"] <- "M4"
Malignant$MM[Malignant$MM%in%"CC5"] <- "M5"
Malignant$MM[Malignant$MM%in%"CC6"] <- "M6"
Malignant$MM[Malignant$MM%in%"CC7"] <- "M7"
Malignant$MM[Malignant$MM%in%"CC8"] <- "M8"

save(Malignant,Malignant.markers,HSMM,total.markers2,file = "Malignant.rda")

load("Malignant.rda")
##MM分群
Malignant$MM2 <- Malignant$MM
Malignant$MM <- "Mali"
DimPlot(Malignant, reduction = "umap",group.by = "MM2")

VlnPlot(Malignant,features = "EGFR",pt.size = 0)
VlnPlot(Malignant,features = "ERBB2",pt.size = 0) ##HER2
VlnPlot(Malignant,features = "PDGFRB",pt.size = 0) ##PDGFR
VlnPlot(Malignant,features = "VEGFA",pt.size = 0) ##VEGF
VlnPlot(Malignant,features = "VEGFA",pt.size = 0) ##VEGF
VlnPlot(Malignant,features = "KDR",pt.size = 0) ##VEGFR2
VlnPlot(Malignant,features = "KIT",pt.size = 0) ##KIT

VlnPlot(Malignant,features = "HLA-DRA",pt.size = 0,slot = "counts")
VlnPlot(Malignant,features = "HLA-DRB1",pt.size = 0,slot = "counts")
VlnPlot(Malignant,features = "CD274",pt.size = 0,slot = "counts")
VlnPlot(Malignant,features = "PDCD1LG2",pt.size = 0,slot = "counts")

