
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)

setwd("E:\\分析过程\\COVID-19\\data")
pbmc <- readRDS(file = "pbmc_asc_RNA.rds")

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")

#对数据进行标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
## #PCA降维之前的标准预处理步骤
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
pbmc=ScaleData(pbmc)

pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA分析

pcSelect=15  #根据前面PCA的图像，选择p-value小于0.05的15个PC


pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)               
pbmc <- FindClusters(object = pbmc, resolution = 0.3)                  #对细胞分组,优化标准模块化

pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)                      #TSNE聚类

# 绘图
TSNEPlot(object = pbmc, pt.size = 1, label = TRUE)    #TSNE可视化

sample <- c("Control","Control","Control","COVID-19","COVID-19","COVID-19")
names(sample) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, sample)
pbmc[["sample"]] <-pbmc@active.ident


TSNEPlot(object = pbmc, pt.size = 0.5, label = TRUE)

celltype <- c("EPCAM+ cells","T cells","Myeloid cells","NK cells","Fibroblasts cells","Fibroblasts cells","Plasma & B cells","Endothelial cells","Myeloid cells","EPCAM+ cells","EPCAM+ cells","Mast cells","Endothelial cells")
names(celltype) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, celltype)
pbmc[["celltype"]] <-pbmc@active.ident


TSNEPlot(object = pbmc, pt.size = 0.1, label = TRUE)


pdf(file="5.TSNE_Cell_2.pdf",width=8,height=6)
TSNEPlot(object = pbmc, pt.size = 0.001, label = TRUE,cols =c("#F5BAC5","#3778AC","#4EA649","#8E4C99","#EE7C1C","#E6DF84","#A15427","#A2CDE6"))
dev.off()

pdf(file="5.TSNE_Cell_3.pdf",width=8,height=6)
TSNEPlot(object = pbmc, pt.size = 0.001,cols =c("#F16D7A","#F65066","#B0706E","#C86F67","#F1B8E3","#F1CCB9","#B9F1CC","#E7DAC9"))
dev.off()


##配色
TSNEPlot(object = pbmc, pt.size = 0.001, label = TRUE,cols =c("#F5BAC5","#B9F1CC","#4EA649","#F1B8E3","#F1CCB9","#E6DF84","#E7DAC9","#A2CDE6"))
pdf(file="5.TSNE_Cell_colour.pdf",width=8,height=6)
TSNEPlot(object = pbmc, pt.size = 0.001, label = TRUE,cols =c("#F5BAC5","#B9F1CC","#4EA649","#F1B8E3","#F1CCB9","#E6DF84","#E7DAC9","#A2CDE6"))
dev.off()


subcelltype <- c("EPCAM+ cells_1","T cells","Myeloid cells_1","NK cells","Fibroblasts cells_1","Fibroblasts cells_2","Plasma & B cells","Endothelial cells_1","Myeloid cells_2","EPCAM+ cells_2","EPCAM+ cells_3","Mast cells","Endothelial cells_2")
names(subcelltype) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, subcelltype)
pbmc[["subcelltype"]] <-pbmc@active.ident


TSNEPlot(object = pbmc, pt.size = 0.1, label = TRUE)

pdf(file="5.TSNE_subcelltyp.pdf",width=8,height=6)
TSNEPlot(object = pbmc, pt.size = 0.1, label = TRUE)
dev.off()


pdf(file="5.TSNE_patient_colour.pdf",width=8,height=6)
TSNEPlot(object = pbmc, pt.size = 0.001,cols =c("#489994","#DF5759"))
dev.off()

####亚型颜色设置

TSNEPlot(object = pbmc, pt.size = 0.1, label = TRUE,cols =c("#4E7AA7","#9ECCE7","#F28D26","#FCBD7A","#8DCA7D","#B49830","#489994","#E15759","#F69A99","#78706E","#D37194","#F7BFD2","#B079A0"))

pdf(file="5.TSNE_subcelltyp_colour.pdf",width=8,height=6)
TSNEPlot(object = pbmc, pt.size = 0.1, label = TRUE,cols =c("#4E7AA7","#9ECCE7","#F28D26","#FCBD7A","#8DCA7D","#B49830","#489994","#E15759","#F69A99","#78706E","#D37194","#F7BFD2","#B079A0"))
dev.off()


saveRDS(pbmc,file="pbmc_2asc_RNA.rds")

