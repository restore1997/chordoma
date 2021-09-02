rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
load("Fibroblast.rda")
Fibroblast$seurat_clusters1 <- Fibroblast$seurat_clusters

Fibroblast <- RunPCA(Fibroblast, features = VariableFeatures(object = Fibroblast))
DimPlot(Fibroblast, reduction = "pca")
DimHeatmap(Fibroblast, dims = 1:15, cells = 500, balanced = TRUE)

Fibroblast <- JackStraw(Fibroblast, num.replicate = 100)
Fibroblast <- ScoreJackStraw(Fibroblast, dims = 1:20)

JackStrawPlot(Fibroblast, dims = 1:15)
ElbowPlot(Fibroblast)
Fibroblast <- FindNeighbors(Fibroblast, dims = 1:10)
Fibroblast <- FindClusters(Fibroblast, resolution = 0.2)
Fibroblast <- RunUMAP(Fibroblast, dims = 1:10)
Fibroblast <- RunTSNE(Fibroblast, dims = 1:10)

DimPlot(Fibroblast, reduction = "umap",label = T,pt.size = 1.6)
DimPlot(Fibroblast, reduction = "umap",group.by = "seurat_clusters1",label = T,pt.size = 1.6)
DimPlot(Fibroblast, reduction = "umap",group.by = "data",pt.size = 1.6)+labs(title = "data")

Fibroblast$location2 <- "NA"
Fibroblast$location2[Fibroblast$location%in%c("basis_cranii","intracranial")] <- "skull base"
Fibroblast$location2[Fibroblast$location%in%"sacrococcyx"] <- "sacrum"
Fibroblast$location2[Fibroblast$location%in%"spinal_canal"] <- "vertebrae"
DimPlot(Fibroblast, reduction = "umap",group.by = "location2",pt.size = 1.6)+labs(title = "location")

DimPlot(Fibroblast, reduction = "umap",group.by = "sex",pt.size = 1.6)+labs(title = "sex")
DimPlot(Fibroblast, reduction = "umap",group.by = "celltype_hpca",label = T,pt.size = 1.6)

Fibroblast.markers <- FindAllMarkers(Fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top5 <- Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Fibroblast, features = top5$gene) #+ NoLegend()

top10 <- Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Fibroblast, features = top10$gene) #+ NoLegend()

write.csv(Fibroblast.markers,file = "Fibroblast_markers.csv")

load("Fibroblast.rda")

VlnPlot(Fibroblast, features = c("HSPA6"), slot = "counts", log = TRUE,pt.size = 0) ##C0
VlnPlot(Fibroblast, features = c("SPP1"), slot = "counts", log = TRUE,pt.size = 0) ##C1
VlnPlot(Fibroblast, features = c("IGLC3"), slot = "counts", log = TRUE) ##C2
VlnPlot(Fibroblast, features = c("ID4"), slot = "counts", log = TRUE,pt.size = 0) ##C3
VlnPlot(Fibroblast, features = c("MAD2L1"), slot = "counts", log = TRUE,pt.size = 0) ##C4
VlnPlot(Fibroblast, features = c("HLA-DRA"), slot = "counts", log = TRUE,pt.size = 0) ##C5

VlnPlot(Fibroblast, features = c("ACTA2"), slot = "counts", log = TRUE,pt.size = 0) 
VlnPlot(Fibroblast, features = c("COL1A2"), slot = "counts", log = TRUE,pt.size = 0) 
VlnPlot(Fibroblast, features = c("PDGFRB"), slot = "counts", log = TRUE,pt.size = 0) 

VlnPlot(Fibroblast, features = c("TGFB1"), slot = "counts", log = TRUE,pt.size = 0) 


#####细胞群注释
library(SingleR)

load("../hpca.rda")
refdata <- ref.data

testdata <- GetAssayData(Fibroblast, slot="data")
clusters <- Fibroblast@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

Fibroblast@meta.data$celltype_hpca = "NA"
for(i in 1:nrow(celltype)){
  Fibroblast@meta.data[which(Fibroblast@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_hpca'] <- celltype$celltype[i]}
DimPlot(Fibroblast, group.by="celltype_hpca", repel=T, label=T, reduction='umap')


###
library(monocle)
counts.data <- as(as.matrix(Fibroblast@assays$RNA@data), 'sparseMatrix')
pheno.data <- Fibroblast@meta.data
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
ordering_genes <- row.names (subset(Fibroblast.markers, p_val_adj < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)


plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "data")
plot_cell_trajectory(HSMM, color_by = "location2")


save(Fibroblast,Fibroblast.markers,HSMM,file = "Fibroblast.rda")

####Featureplot
library(ggplot2)
FeaturePlot(Fibroblast,features = "COL1A2",pt.size = 1.2)+labs(title = "pan-COL1A2") ##7*5
FeaturePlot(Fibroblast,features = "ACTA2",pt.size = 1.2)+labs(title = "pan-ACTA2")
FeaturePlot(Fibroblast,features = "PDGFRB",pt.size = 1.2)+labs(title = "pan-PDGFRB")

FeaturePlot(Fibroblast,features = "COL1A1",pt.size = 1.2,min.cutoff = 0)+labs(title = "C0-COL1A1")
FeaturePlot(Fibroblast,features = "SFRP2",pt.size = 1.2,min.cutoff = 0)+labs(title = "C0-SFRP2")
FeaturePlot(Fibroblast,features = "THBS2",pt.size = 1.2,min.cutoff = 0)+labs(title = "C0-THBS2")

FeaturePlot(Fibroblast,features = "CHRDL1",pt.size = 1.2,min.cutoff = 0)+labs(title = "C1-CHRDL1")
FeaturePlot(Fibroblast,features = "GADD45G",pt.size = 1.2,min.cutoff = 0)+labs(title = "C1-GADD45G")
FeaturePlot(Fibroblast,features = "CCL3",pt.size = 1.2,min.cutoff = 0)+labs(title = "C1-CCL3")

FeaturePlot(Fibroblast,features = "FNDC1",pt.size = 1.2,min.cutoff = 0)+labs(title = "C2-FNDC1")
FeaturePlot(Fibroblast,features = "WWTR1",pt.size = 1.2,min.cutoff = 0)+labs(title = "C2-WWTR1")
FeaturePlot(Fibroblast,features = "CRABP2",pt.size = 1.2,min.cutoff = 0)+labs(title = "C2-CRABP2")

FeaturePlot(Fibroblast,features = "ISG15",pt.size = 1.2,min.cutoff = 0)+labs(title = "C3-ISG15")
FeaturePlot(Fibroblast,features = "TTR",pt.size = 1.2,min.cutoff = 0)+labs(title = "C3-TTR")
FeaturePlot(Fibroblast,features = "CST1",pt.size = 1.2,min.cutoff = 0)+labs(title = "C3-CST1")
FeaturePlot(Fibroblast,features = "SFRP4",pt.size = 1.2,min.cutoff = 0)+labs(title = "C3-SFRP4")
