rm(list = ls())

library(Seurat)
load("T_NK.rda")
table(T_NK$MM2)

load("../Malignant/Malignant.rda")
table(Malignant$MM2)
Malignant$MM2 <- "Malignant"

total <- merge(T_NK,Malignant)
table(total$MM2)

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)

total

data.input  <- total@assays$RNA@data
identity = data.frame(group =total$MM2, row.names = names(total$MM2)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels

cellchat <- createCellChat(data = data.input)
cellchat
summary(cellchat)

str(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human

colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

unique(CellChatDB$interaction$annotation)

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat) 
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat@netP$pathways

head(cellchat@LR$LRsig)

save(cellchat,file = "cellchat2.rda")

load("cellchat.rda")
cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4) # a numeric vector
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
# check the order of cell identity to set suitable vertex.receiver
#cellchat@LR$LRsig$pathway_name
#cellchat@LR$LRsig$antagonist
pathways.show <- "TGFb"
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)   # 原函数

netVisual_aggregate(cellchat, signaling = "TGFb", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "GDF", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "MSTN", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "GDNF", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

##"TGFb"       "BMP"        "GDF"        "MSTN"       "GDNF"       "ACTIVIN"    "EGF"  
##"NRG"        "FGF"        "PDGF"       "VEGF"       "IGF"        "CCL"        "CXCL"      
#"MIF"        "CX3C"       "IL2"        "IL4"        "IL6"        "LIFR"       "OSM"       
#"IL1"        "IL17"       "CSF"        "IL16"       "EPO"        "GH"         "CSF3"      
#"IFN-II"     "TNF"        "LT"         "LIGHT"      "FASLG"      "VEGI"       "TRAIL"     
#"NGF"        "RANKL"      "TWEAK"      "CD30"       "CD137"      "OX40"       "BAFF"      
#"CD40"       "SPP1"       "VISFATIN"   "ANGPTL"     "MK"         "PTN"        "COMPLEMENT"
#"NMU"        "PARs"       "NPR2"       "KIT"        "NT"         "FLT3"       "HGF"       
#"SEMA3"      "GAS"        "GRN"        "PROS"       "PSAP"       "BTLA"       "BAG" 

netVisual_aggregate(cellchat, signaling = "PDGF", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "VEGF", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "CCL", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

netVisual_aggregate(cellchat, signaling = "MIF", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "IL4", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "IL6", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

netVisual_aggregate(cellchat, signaling = "IL1", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "EPO", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "GH", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

#"IFN-II"     "TNF"        "LT"         "LIGHT"      "FASLG"      "VEGI"       "TRAIL"     
#"NGF"        "RANKL"      "TWEAK"      "CD30"       "CD137"      "OX40"       "BAFF"      

netVisual_aggregate(cellchat, signaling = "IFN-II", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "TNF", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "LT", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "LIGHT", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "FASLG", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "TRAIL", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

netVisual_aggregate(cellchat, signaling = "TWEAK", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "CD30", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "CD137", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "OX40", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "BAFF", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

#"CD40"       "SPP1"       "VISFATIN"   "ANGPTL"     "MK"         "PTN"        "COMPLEMENT"
#"NMU"        "PARs"       "NPR2"       "KIT"        "NT"         "FLT3"       "HGF"       
#"SEMA3"      "GAS"        "GRN"        "PROS"       "PSAP"       "BTLA"       "BAG" 

netVisual_aggregate(cellchat, signaling = "CD40", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "PTN", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "COMPLEMENT", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

netVisual_aggregate(cellchat, signaling = "PARs", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)
netVisual_aggregate(cellchat, signaling = "NPR2", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

netVisual_aggregate(cellchat, signaling = "BAG", layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

netAnalysis_contribution(cellchat, signaling = pathways.show)
cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(cellchat, signaling = pathways.show, width = 12, height = 2.5, font.size = 10)


nPatterns = 5 
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# Visualize the communication pattern using river plot
netAnalysis_river(cellchat, pattern = "outgoing")


# Visualize the communication pattern using dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

