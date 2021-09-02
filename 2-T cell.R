rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)
load("T_NK.rda")
T_NK$seurat_clusters1 <- T_NK$seurat_clusters

T_NK <- RunPCA(T_NK, features = VariableFeatures(object = T_NK))
DimPlot(T_NK, reduction = "pca")
DimHeatmap(T_NK, dims = 1:15, cells = 500, balanced = TRUE)

T_NK <- JackStraw(T_NK, num.replicate = 100)
T_NK <- ScoreJackStraw(T_NK, dims = 1:20)

JackStrawPlot(T_NK, dims = 1:15)
ElbowPlot(T_NK)
T_NK <- FindNeighbors(T_NK, dims = 1:10)
T_NK <- FindClusters(T_NK, resolution = 0.2)
T_NK <- RunUMAP(T_NK, dims = 1:10)
T_NK <- RunTSNE(T_NK, dims = 1:10)

DimPlot(T_NK, reduction = "umap",label = T,pt.size = 1.6)
DimPlot(T_NK, reduction = "umap",group.by = "seurat_clusters1",label = T,pt.size = 1.6)
DimPlot(T_NK, reduction = "umap",group.by = "data",pt.size = 1.6)+labs(title = "data")
DimPlot(T_NK, reduction = "umap",group.by = "sex",pt.size = 1.6)+labs(title = "sex")
DimPlot(T_NK, reduction = "umap",group.by = "location",pt.size = 1.6)+labs(title = "sex")

T_NK$location2 <- "NA"
T_NK$location2[T_NK$location%in%c("basis_cranii","intracranial")] <- "skull base"
T_NK$location2[T_NK$location%in%"sacrococcyx"] <- "sacrum"
T_NK$location2[T_NK$location%in%"spinal_canal"] <- "vertebrae"
DimPlot(T_NK,reduction = "umap",group.by = "location2",pt.size = 1.6)+labs(title = "location")

T_NK.markers <- FindAllMarkers(T_NK, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T_NK.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
save(T_NK,T_NK.markers,file = "T_NK.rda")

top5 <- T_NK.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(T_NK, features = top5$gene) #+ NoLegend()

top10 <- T_NK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(T_NK, features = top10$gene) #+ NoLegend()

library(clusterProfiler)
library(org.Hs.eg.db)

genelist <- T_NK.markers$gene
genelist <- unique(genelist)
EntrezID <- bitr(genelist,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(EntrezID)

ego <- enrichGO(gene = EntrezID$ENTREZID,
                OrgDb='org.Hs.eg.db'
                ,ont ='ALL'
                ,pAdjustMethod ='BH'
                ,pvalueCutoff =0.05
                ,qvalueCutoff =0.2
                ,keyType ='ENTREZID'
)
table(ego@result$ONTOLOGY);dim(ego)
barplot(ego,showCategory = 10,drop=T,title = "GO")
dotplot(ego)

gene <- EntrezID$ENTREZID
kk <- enrichKEGG(gene= gene,
                 organism     = 'hsa',
                 #universe     = gene_all,
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.2)
barplot(kk,showCategory = 10,drop=T)
dotplot(kk)

load("T_NK.rda")

VlnPlot(T_NK, features = c("CTLA4"), slot = "counts", log = TRUE,pt.size = 0) ##Treg
VlnPlot(T_NK, features = c("FOXP3"), slot = "counts", log = TRUE,pt.size = 0) ##Treg
VlnPlot(T_NK, features = c("CD8A"), slot = "counts", log = TRUE,pt.size = 0) ##CD8+T
VlnPlot(T_NK, features = c("CD3E"), slot = "counts", log = TRUE,pt.size = 0) ##T

##C0
VlnPlot(T_NK, features = c("IL7R"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("GPR183"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("LTB"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("S100A4"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("CD40LG"), slot = "counts", log = TRUE,pt.size = 0) ##T

##C1
VlnPlot(T_NK, features = c("CD8A"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("CCL4"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("NKG7"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("GZMK"), slot = "counts", log = TRUE,pt.size = 0) ##cytotoxic marker
VlnPlot(T_NK, features = c("PRF1"), slot = "counts", log = TRUE,pt.size = 0) ##cytotoxic marker
VlnPlot(T_NK, features = c("IFNG"), slot = "counts", log = TRUE,pt.size = 0) ##cytotoxic marker

##C2
VlnPlot(T_NK, features = c("FGFBP2"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("GNLY"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("KLRD1"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("TYROBP"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("NKG7"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("KLRF1"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("GZMK"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("GZMH"), slot = "counts", log = TRUE,pt.size = 0) ##T

##C3
VlnPlot(T_NK, features = c("TNFRSF18"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("SCX"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("LST1"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("SCN1B"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("IL4I1"), slot = "counts", log = TRUE,pt.size = 0) ##T

##C4
VlnPlot(T_NK, features = c("CTLA4"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("BATF"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("TIGIT"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("FOXP3"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("LAIR2"), slot = "counts", log = TRUE,pt.size = 0) ##T

##标记细胞群
VlnPlot(T_NK, features = c("CD4"), slot = "counts", log = TRUE,pt.size = 0) ##CD4 T
VlnPlot(T_NK, features = c("FOXP3"), slot = "counts", log = TRUE,pt.size = 0) ##Treg
VlnPlot(T_NK, features = c("CD8A"), slot = "counts", log = TRUE) ##CD8 T
VlnPlot(T_NK, features = c("KLRF1"), slot = "counts", log = TRUE,pt.size = 0) ##NK
VlnPlot(T_NK, features = c("CD3E"), slot = "counts", log = TRUE,pt.size = 0) ##T
VlnPlot(T_NK, features = c("FCGR3A"), slot = "counts", log = TRUE,pt.size = 0) ##NK
VlnPlot(T_NK, features = c("CD160"), slot = "counts", log = TRUE) ##
VlnPlot(T_NK, features = c("CD160"), slot = "counts", log = TRUE) ##

##C3
VlnPlot(T_NK, features = c("CD3E"), slot = "counts", log = TRUE,pt.size = 0) ##

VlnPlot(T_NK, features = c("SCX"), slot = "counts", log = TRUE,pt.size = 0) ##
VlnPlot(T_NK, features = c("VEGFA"), slot = "counts", log = TRUE,pt.size = 0) ##
VlnPlot(T_NK, features = c("SCX"), slot = "counts", log = TRUE,pt.size = 0) ##


T_NK$MM2 <- "NA"
T_NK$MM2[T_NK$seurat_clusters%in%"0"] <- "CD4T"
T_NK$MM2[T_NK$seurat_clusters%in%"1"] <- "CD8T"
T_NK$MM2[T_NK$seurat_clusters%in%"2"] <- "NK"
T_NK$MM2[T_NK$seurat_clusters%in%"3"] <- "T0"
T_NK$MM2[T_NK$seurat_clusters%in%"4"] <- "Treg"

DimPlot(T_NK, reduction = "umap",group.by = "MM2",pt.size = 1.6)

library(SingleR)
load("../hpca.rda")
refdata <- ref.data

##C3
testdata3 <- subset(T_NK,subset=seurat_clusters%in%"3")
testdata3 <- GetAssayData(testdata3, slot="data")
cellpred3 <- SingleR(test = testdata3, ref = refdata, labels = refdata$label.fine,
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
sort(table(cellpred3$first.labels))
c3 <- T_NK.markers[T_NK.markers$cluster%in%"3",]
C3 <- c3$gene
write.table(C3,file = "c3_marker.txt")

##pseudotime
library(monocle)
counts.data <- as(as.matrix(T_NK@assays$RNA@data), 'sparseMatrix')
pheno.data <- T_NK@meta.data
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
ordering_genes <- row.names (subset(T_NK.markers, p_val_adj < 0.01))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
HSMM <- orderCells(HSMM)


plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "seurat_clusters")
plot_cell_trajectory(HSMM, color_by = "data")
plot_cell_trajectory(HSMM, color_by = "MM2")

VlnPlot(T_NK, features = c("EGFR"), slot = "counts", log = TRUE,pt.size = 0) 
VlnPlot(T_NK, features = c("ERBB2"), slot = "counts", log = TRUE,pt.size = 0) ##HER2
VlnPlot(T_NK, features = c("PDGFRB"), slot = "counts", log = TRUE,pt.size = 0) ##PDGFR
VlnPlot(T_NK, features = c("KDR"), slot = "counts", log = TRUE,pt.size = 0) ##VEGFR
VlnPlot(T_NK, features = c("EZH2"), slot = "counts", log = TRUE,pt.size = 0) ##VEGFR

VlnPlot(T_NK, features = c("CD274"), slot = "counts", log = TRUE,pt.size = 0) ##PD-L1
VlnPlot(T_NK, features = c("PDCD1"), slot = "counts", log = TRUE,pt.size = 0) ##PD-1
VlnPlot(T_NK, features = c("CTLA4"), slot = "counts", log = TRUE,pt.size = 0) ##CTLA4
VlnPlot(T_NK, features = c("HAVCR2"), slot = "counts", log = TRUE,pt.size = 0) ##TIMI3
VlnPlot(T_NK, features = c("TIGIT"), slot = "counts", log = TRUE,pt.size = 0) ##TIGIT

features <- unique(top5$gene)
DotPlot(T_NK, features = features) + RotatedAxis()


FeaturePlot(T_NK, features ="CD8A",slot = "counts",pt.size = 1.6) ##CD8
FeaturePlot(T_NK, features ="KLRF1",slot = "counts",pt.size = 1.6) ##NK
FeaturePlot(T_NK, features ="FOXP3",slot = "counts",pt.size = 1.6) ##Treg
FeaturePlot(T_NK, features ="CD4",slot = "counts",pt.size = 1.6)  ##CD4
FeaturePlot(T_NK, features ="LAG3",slot = "counts",pt.size = 1.6) ##exhausted T cell
FeaturePlot(T_NK, features ="TOP2A",slot = "counts",pt.size = 1.6) ##Proliferating T cells



###C1
C1 <- Malignant.markers[Malignant.markers$cluster%in%"1",]
library(clusterProfiler)
library(org.Hs.eg.db)

genelist <- C1$gene
genelist <- unique(genelist)
EntrezID <- bitr(genelist,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(EntrezID)

###GO
ego <- enrichGO(gene = EntrezID$ENTREZID,
                OrgDb='org.Hs.eg.db'
                ,ont ='ALL'
                ,pAdjustMethod ='BH'
                ,pvalueCutoff =0.05
                ,qvalueCutoff =0.2
                ,keyType ='ENTREZID'
)
table(ego@result$ONTOLOGY);dim(ego)
barplot(ego,showCategory = 10,drop=T,title = "GO")
dotplot(ego)

gene <- EntrezID$ENTREZID
kk <- enrichKEGG(gene= gene,
                 organism     = 'hsa',
                 #universe     = gene_all,
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.2)
barplot(kk,showCategory = 10,drop=T)
dotplot(kk)
save(ego,kk,file = "C1_GO.rda")