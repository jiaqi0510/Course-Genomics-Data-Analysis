library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

# 读取数据
AIS <- readRDS("E:/course/大三秋/基因组学数据分析/final/AIS.Rdata")
MIA <- readRDS("E:/course/大三秋/基因组学数据分析/final/MIA.Rdata")
IAC <- readRDS("E:/course/大三秋/基因组学数据分析/final/IAC.Rdata")
AIS[["RNA"]]<-JoinLayers(AIS[["RNA"]])
MIA[["RNA"]]<-JoinLayers(MIA[["RNA"]])
IAC[["RNA"]]<-JoinLayers(IAC[["RNA"]])

# 识别高变基因
ais1 <- FindVariableFeatures(AIS, selection.method = "vst", nfeatures = 3000)
mia1 <- FindVariableFeatures(MIA, selection.method = "vst", nfeatures = 3000)
iac1 <- FindVariableFeatures(IAC, selection.method = "vst", nfeatures = 3000)
# 查看高变基因
head(VariableFeatures(ais1))

# 进行 PCA 降维
ais1 <- ScaleData(ais1, features = VariableFeatures(ais1))
ais1<- RunPCA(ais1,features = VariableFeatures(ais1),npcs=50)

mia1 <- ScaleData(mia1, features = VariableFeatures(mia1))
mia1<- RunPCA(mia1,features = VariableFeatures(mia1),npcs=50)

iac1 <- ScaleData(iac1, features = VariableFeatures(iac1))
iac1<- RunPCA(iac1,features = VariableFeatures(iac1),npcs=50)
# 可视化前两个主成分
DimPlot(ais1, reduction = "pca")

# 使用前30个PC进行Louvain聚类
ais1<- FindNeighbors(ais1, dims = 1:30)
ais1<- FindClusters(ais1, resolution = 0.4)

mia1<- FindNeighbors(mia1, dims = 1:30)
mia1<- FindClusters(mia1, resolution = 0.4)

iac1<- FindNeighbors(iac1, dims = 1:30)
iac1<- FindClusters(iac1, resolution = 0.4)

# 进行UMAP降维
ais1 <- RunUMAP(ais1, dims = 1:50,min.dist=0.01,spread=1)
mia1 <- RunUMAP(mia1, dims = 1:50,min.dist=0.01,spread=1)
iac1 <- RunUMAP(iac1, dims = 1:50,min.dist=0.01,spread=1)

# 可视化簇
DimPlot(ais1, reduction = "umap", label = TRUE,pt.size = 1)
DimPlot(mia1, reduction = "umap", label = TRUE,pt.size = 1)
DimPlot(iac1, reduction = "umap", label = TRUE,pt.size = 1)

#识别每个聚类的差异基因
all.markers <- FindAllMarkers(ais1, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
significant.markers <- all.markers [all.markers $p_val_adj < 0.05 & all.markers $avg_log2FC > 0.25, ]
#writecsv保存到本地
write.csv(significant.markers, file = "E:/course/大三秋/基因组学数据分析/final/significant.markers.csv")

markers <- c('SCGB1A1','GZMK','CD3D','GNLY','NKG7','CD79A','JCHAIN','IGHG1','LYZ','TPSB2','AIF1','HLA-DRA','HLA-DRB1','CPA3','CLDN5','COL1A1','DCN')
#B cells
markers1 <- c('JCHAIN','IGHM','IGHG1', 'IGHA1','IGHG4', 'IGHA2','IGHG3','IGHG2','IGKC','IGHD','MT1X','GZMB','GNLY','SFTPC')
#endothelial cells
markers2 <- c('CLDN5','RAMP2','CAV1','VWF','PECAM1','CDH5','EMCN','CD34','CD31','THBD') 
#fibroblasts        
markers3 <- c('DCN','LUM','COL1A2','COL3A1','COL1A1','MMP11','FN1','COL6A2','COL6A3','ACTA2','COL5A2','COL6A1','COL4A2') 
#Myeloid    
markers4 <- c('TPSB2','CCL3L1','FABP4','FCN1','CPA3','COL1A1','DCN','LYZ','CTSD','CCL18','CD68','CD14','TNF','CD63','CD163','FCGR2A') 
#T
markers5 <- c('LTB','CXCL13','GZMK','CCL4L2','GZMB', 'FOXP3', 'IL2RA')
#NK
markers7 <- c('GNLY', 'NKG7', 'KLRD1', 'KLRC1', 'TIGIT', 'LTB')
#epithelial
markers6 <- c('CAPS', 'SCGB1A1', 'WFDC2', 'KRT8', 'KRT19', 'SCGB3A1', 'SCGB3A2', 'KRT18', 'EPCAM', 'MUC1', 'SFTPC', 'SFTPA1', 'SFTPA2', 'SFTPB', 'PGC', 'AGER', 'CAV1', 'CYP4B1')
#plasma
plasma_cell_markers <- c('IGHG1', 'MZB1', 'SDC1', 'CD79A', 'CD27', 'CD38', 'CD81')
#MALT B
markers8 <- c('GZMB','IGHD')

DotPlot(mia_ann,features = plasma_cell_markers)+coord_flip()+theme(axis.text.x = element_text(size = 9))



#UMAP图
markers<- c('EPCAM','NKG7','LYZ','CD79A','CLDN5','DCN')
FeaturePlot(ais1,features = markers)
FeaturePlot(mia_ann,features = markers)
FeaturePlot(iac_ann,features = markers)

ais1 <- subset(ais1, idents = c('16','18','19','20','21'), invert = TRUE)#去掉低质量细胞群
mia1 <- subset(mia1, idents = c('12','14','16','17','18','19'), invert = TRUE)
iac1 <- subset(iac1, idents = c('11','13','16'), invert = TRUE)

ais.cluster.ids <- c("0"="T&NK cells",  
                     "1"="Epithelial cells",  
                     "2"="Mast cells", 
                     "3"="T&NK cells",
                     "4"="Myeloid cells", 
                     "5"="B cells", 
                     "6"="T&NK cells", 
                     '7'="Epithelial cells",
                     '8'="Myeloid cells",
                     '9'="Myeloid cells",
                     "10"="Fibroblast",
                     "11"="Plasma cells",
                     "12"="Plasma cells",
                     "13"="Endothelial cells",
                     "14"="Myeloid cells",
                     "15"="Fibroblast",
                     '16'="others",
                     '17'="Myeloid cells",
                     "18"="others",
                     '19'="others",
                     '20'="others",
                     '21'="others",
                     "22"="MALT B"
)

mia.cluster.ids <- c("0"="T&NK cells",  
                     "1"="T&NK cells",  
                     "2"="B cells", 
                     "3"="Myeloid cells",
                     "4"="Myeloid cells", 
                     "5"="Myeloid cells", 
                     "6"="Endothelial cells", 
                     '7'="Fibroblast",
                     '8'="Endothelial cells",
                     '9'="Myeloid cells",
                     "10"="Epithelial cells",
                     "11"="Myeloid cells",
                     "12"="others",
                     "13"="Epithelial cells",
                     "14"="others",
                     "15"="Epithelial cells",
                     '16'="others",
                     '17'="others",
                     "18"="others",
                     '19'="others"
                     
)

iac.cluster.ids <- c("0"="T&NK cells",  
                     "1"="T&NK cells",  
                     "2"="B cells", 
                     "3"="T&NK cells",
                     "4"="Epithelial cells", 
                     "5"="Myeloid cells", 
                     "6"="Myeloid cells", 
                     '7'="Myeloid cells",
                     '8'="Epithelial cells",
                     '9'="T&NK cells",
                     "10"="Epithelial cells",
                     "11"="others",
                     "12"="Plasma cells",
                     "13"="others",
                     "14"="Myeloid cells",
                     "15"="Fibroblast",
                     '16'="others",
                     '17'="Endothelial cells",
                     "18"="Epithelial cells",
                     '19'="Myeloid cells",
                     '20'='MALT B'
                     
)

#给ais1加celltype注释，获得注释后结果ais_ann
ais_ann <- RenameIdents(ais1, ais.cluster.ids)                        
ais_ann$celltype <- ais_ann@active.ident

mia_ann <- RenameIdents(mia1, mia.cluster.ids)                        
mia_ann$celltype <- mia_ann@active.ident

iac_ann <- RenameIdents(iac1, iac.cluster.ids)                        
iac_ann$celltype <- iac_ann@active.ident

#保存注释结果
saveRDS(ais_ann, file = "E:/course/大三秋/基因组学数据分析/final/ais_ann.rds")
saveRDS(mia_ann, file = "E:/course/大三秋/基因组学数据分析/final/mia_ann.rds")
saveRDS(iac_ann, file = "E:/course/大三秋/基因组学数据分析/final/iac_ann.rds")

#可视化
DimPlot(ais_ann, group.by = "celltype",label = 'TRUE')
DimPlot(mia_ann, group.by = "celltype",label = 'TRUE')
DimPlot(iac_ann, group.by = "celltype",label = 'TRUE')

#保存细胞注释图片
save_path="E:/course/大三秋/基因组学数据分析/final/"
ggsave(filename = paste0(save_path, "celltype", ".", "png"), plot = p, width = 10, height = 8)

#分别计算AIS与MIA和IAC的各类细胞比例
ais_ann <- readRDS("E:/course/大三秋/基因组学数据分析/final/ais_ann.rds")
mia_ann <- readRDS("E:/course/大三秋/基因组学数据分析/final/mia_ann.rds")
iac_ann <- readRDS("E:/course/大三秋/基因组学数据分析/final/iac_ann.rds")
# 获取总细胞数
total_cells_count_ais <- ncol(ais_ann)
total_cells_count_ais
total_cells_count_mia <- ncol(mia_ann)
total_cells_count_mia
total_cells_count_iac <- ncol(iac_ann)
total_cells_count_iac

# 计算B细胞的比例
b_cells_count_ais <- sum(ais_ann$celltype == "B cells")
b_cells_ratio_ais <- b_cells_count_ais / total_cells_count_ais

b_cells_count_mia <- sum(mia_ann$celltype == "B cells")
b_cells_ratio_mia <- b_cells_count_mia / total_cells_count_mia

b_cells_count_iac <- sum(iac_ann$celltype == "B cells")
b_cells_ratio_iac <- b_cells_count_iac / total_cells_count_iac

# 计算T&NK细胞的比例
TNK_cells_count_ais <- sum(ais_ann$celltype == "T&NK cells")
TNK_cells_ratio_ais <- TNK_cells_count_ais / total_cells_count_ais

TNK_cells_count_mia <- sum(mia_ann$celltype == "T&NK cells")
TNK_cells_ratio_mia <- TNK_cells_count_mia / total_cells_count_mia

TNK_cells_count_iac <- sum(iac_ann$celltype == "T&NK cells")
TNK_cells_ratio_iac <- TNK_cells_count_iac / total_cells_count_iac

# 计算Myeloid细胞的比例
Myeloid_cells_count_ais <- sum(ais_ann$celltype == "Myeloid cells")
Myeloid_cells_ratio_ais <- Myeloid_cells_count_ais / total_cells_count_ais

Myeloid_cells_count_mia <- sum(mia_ann$celltype == "Myeloid cells")
Myeloid_cells_ratio_mia <- Myeloid_cells_count_mia / total_cells_count_mia

Myeloid_cells_count_iac <- sum(iac_ann$celltype == "Myeloid cells")
Myeloid_cells_ratio_iac <- Myeloid_cells_count_iac / total_cells_count_iac

# 计算Fibroblast细胞的比例
Fibroblast_cells_count_ais <- sum(ais_ann$celltype == "Fibroblast")
Fibroblast_cells_ratio_ais <- Fibroblast_cells_count_ais / total_cells_count_ais

Fibroblast_cells_count_mia <- sum(mia_ann$celltype == "Fibroblast")
Fibroblast_cells_ratio_mia <- Fibroblast_cells_count_mia / total_cells_count_mia

Fibroblast_cells_count_iac <- sum(iac_ann$celltype == "Fibroblast")
Fibroblast_cells_ratio_iac <- Fibroblast_cells_count_iac / total_cells_count_iac

# 计算Epithelial细胞的比例
Epithelial_cells_count_ais <- sum(ais_ann$celltype == "Epithelial cells")
Epithelial_cells_ratio_ais <- Epithelial_cells_count_ais / total_cells_count_ais

Epithelial_cells_count_mia <- sum(mia_ann$celltype == "Epithelial cells")
Epithelial_cells_ratio_mia <- Epithelial_cells_count_mia / total_cells_count_mia

Epithelial_cells_count_iac <- sum(iac_ann$celltype == "Epithelial cells")
Epithelial_cells_ratio_iac <- Epithelial_cells_count_iac / total_cells_count_iac

# 计算Endothelial细胞的比例
Endothelial_cells_count_ais <- sum(ais_ann$celltype == "Endothelial cells")
Endothelial_cells_ratio_ais <- Endothelial_cells_count_ais / total_cells_count_ais

Endothelial_cells_count_mia <- sum(mia_ann$celltype == "Endothelial cells")
Endothelial_cells_ratio_mia <- Endothelial_cells_count_mia / total_cells_count_mia

Endothelial_cells_count_iac <- sum(iac_ann$celltype == "Endothelial cells")
Endothelial_cells_ratio_iac <- Endothelial_cells_count_iac / total_cells_count_iac

cell_proportions_ais <- data.frame(celltype = c("B cells","T&NK cells","Myeloid cells","Fibroblast","Epithelial cells","Endothelial cells"), 
                               proportion = c(b_cells_ratio_ais,TNK_cells_ratio_ais,Myeloid_cells_ratio_ais,Fibroblast_cells_ratio_ais,Epithelial_cells_ratio_ais,Endothelial_cells_ratio_ais))
cell_proportions_mia <- data.frame(celltype = c("B cells","T&NK cells","Myeloid cells","Fibroblast","Epithelial cells","Endothelial cells"), 
                                   proportion = c(b_cells_ratio_mia,TNK_cells_ratio_mia,Myeloid_cells_ratio_mia,Fibroblast_cells_ratio_mia,Epithelial_cells_ratio_mia,Endothelial_cells_ratio_mia))
cell_proportions_iac <- data.frame(celltype = c("B cells","T&NK cells","Myeloid cells","Fibroblast","Epithelial cells","Endothelial cells"), 
                                  proportion = c(b_cells_ratio_iac,TNK_cells_ratio_iac,Myeloid_cells_ratio_iac,Fibroblast_cells_ratio_iac,Epithelial_cells_ratio_iac,Endothelial_cells_ratio_iac))

#合并数据
cell_proportions <- rbind(
  cell_proportions_mia %>% mutate(sample = 'MIA'),
  cell_proportions_ais %>% mutate(sample = 'AIS')
)

cell_proportions <- rbind(
  cell_proportions_mia %>% mutate(sample = 'MIA'),
  cell_proportions_iac %>% mutate(sample = 'IAC')
)

cell_proportions <- rbind(
  cell_proportions_ais %>% mutate(sample = 'AIS'),
  cell_proportions_mia %>% mutate(sample = 'MIA'),
  cell_proportions_iac %>% mutate(sample = 'IAC')
)
cell_proportions$sample <- factor(cell_proportions$sample, levels = c("AIS", "MIA", "IAC"))
cell_proportions$sample <- factor(cell_proportions$sample, levels = c( "IAC", "MIA"))
# 使用ggplot2绘制条形图
p <- ggplot(cell_proportions, aes(x = celltype, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") + # 使用dodge使条形分开
  labs(title = "Proportion of Cells in AIS and MIA") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # 旋转X轴标签以便更好地显示

p <- ggplot(cell_proportions, aes(x = sample, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") + # 使用dodge使条形分开
  labs(title = "Proportion of Cells in IAC and MIA",x = "Sample", y = "Proportion",fill = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_text( hjust = 0.5),,plot.title = element_text(hjust = 0.5))+
  facet_wrap(~celltype, scales = "free")

p <- ggplot(cell_proportions, aes(x = sample, y = proportion, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge") + # 使用dodge使条形分开
  labs(title = "Proportion of Cells in AIS, MIA and IAC",x = "Sample", y = "Proportion",fill = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_text( hjust = 0.5),plot.title = element_text(hjust = 0.5))+
  facet_wrap(~celltype, scales = "free")
# 打印图表
print(p)

#对epithelial做细胞亚型注释，用于inferCNV识别癌细胞
#AIS
#提取epithelial
epithelial_subset <- ais_ann[,ais_ann@meta.data$celltype %in% c('Epithelial cells')]
epithelial_subset<- NormalizeData(epithelial_subset, normalization.method = "LogNormalize", scale.factor = 1e4)
epithelial_subset <- FindVariableFeatures(epithelial_subset, selection.method = "vst", nfeatures = 2000)
epithelial_subset <- ScaleData(epithelial_subset, features = VariableFeatures(epithelial_subset))
epithelial_subset <- RunPCA(epithelial_subset,features = VariableFeatures(epithelial_subset),npcs=50)
DimPlot(epithelial_subset, reduction = "umap", label = TRUE,pt.size = 1)

# 使用前30个PC进行Louvain聚类
epithelial_subset<- FindNeighbors(epithelial_subset, dims = 1:30)
epithelial_subset<- FindClusters(epithelial_subset, resolution = 0.4)

# 进行UMAP降维
epithelial_subset <- RunUMAP(epithelial_subset, dims = 1:30,min.dist=0.001,spread=1)

# 可视化簇
DimPlot(epithelial_subset, reduction = "umap", label = TRUE,pt.size = 1)+ 
  labs(title = "Epithelial clusters")+ 
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(hjust = 0.5))

#去掉低质量细胞
epithelial_subset <- subset(epithelial_subset, idents = c('0','1','2','3','4','5','6','7','8'), invert = FALSE)

#识别每个聚类的差异基因
all.markers <- FindAllMarkers(epithelial_subset, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
markers__epithelial <- all.markers [all.markers $p_val_adj < 0.05 & all.markers $avg_log2FC > 0.25, ]
#保存到本地
write.csv(markers__epithelial, file = "E:/course/大三秋/基因组学数据分析/final/markers_epithelial.csv")

markers <- c('CAV1','CLDN18','UBE2C','SCGB1A1','TM4SF1','CAPS','TMSF1')
VlnPlot(epithelial_subset,features = markers, pt.size = 0)+theme(axis.text.x = element_text(size = 9))

#UMAP图
FeaturePlot(epithelial_subset,features = markers)
epithelial_subset <- subset(epithelial_subset, idents = c('0','1','3','4'), invert = TRUE)#去掉低质量细胞群
epithelial_subset.cluster.ids <- c('2'='AT2',
                     '5'='Clara cells',
                     '6'='Clara-like cancer',
                     '7'='Ciliated cells',
                     '8'='Clara-like cancer'
                     )

#加celltype注释
epithelial_subset <- RenameIdents(epithelial_subset, epithelial_subset.cluster.ids)                        
epithelial_subset$celltype <- epithelial_subset@active.ident

#可视化
DimPlot(epithelial_subset, group.by = "celltype",label = 'TRUE')+ 
  labs(title = "Epithelial cell subtype")+ 
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(hjust = 0.5))

#保存注释结果
saveRDS(epithelial_subset, file = "E:/course/大三秋/基因组学数据分析/final/epithelial_ann.rds")



#MIA
epithelial_subset <- mia_ann[, mia_ann@meta.data$celltype %in% c('Epithelial cells')]
epithelial_subset=CreateSeuratObject(counts = epithelial_subset@assays$RNA@counts,
                       meta.data = epithelial_subset@meta.data) 

epithelial_subset <- subset(mia_ann, idents = "Epithelial cells")
epithelial_subset<- NormalizeData(epithelial_subset, normalization.method = "LogNormalize", scale.factor = 1e4)
epithelial_subset<- ScaleData(epithelial_subset)
epithelial_subset <- FindVariableFeatures(epithelial_subset, selection.method = "vst", nfeatures = 2000)
epithelial_subset <- ScaleData(epithelial_subset, features = VariableFeatures(epithelial_subset))
epithelial_subset <- RunPCA(epithelial_subset,features = VariableFeatures(epithelial_subset),npcs=50)
DimPlot(epithelial_subset, reduction = "umap", label = TRUE,pt.size = 1)
# 使用前30个PC进行Louvain聚类
epithelial_subset<- FindNeighbors(epithelial_subset, dims = 1:30)
epithelial_subset<- FindClusters(epithelial_subset, resolution = 0.4)

# 进行UMAP降维
epithelial_subset <- RunUMAP(epithelial_subset, dims = 1:30,min.dist=0.001,spread=1)

# 可视化簇
DimPlot(epithelial_subset, reduction = "umap", label = TRUE,pt.size = 1)+ 
  labs(title = "MIA Epithelial clusters")+ 
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(hjust = 0.5))

#去掉低质量细胞
epithelial_subset <- subset(epithelial_subset, idents = c('0','1','2','3','4','5','6','7','8'), invert = FALSE)

#识别每个聚类的差异基因
all.markers <- FindAllMarkers(epithelial_subset, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
markers__epithelial <- all.markers [all.markers $p_val_adj < 0.05 & all.markers $avg_log2FC > 0.25, ]
#保存到本地
write.csv(markers__epithelial, file = "E:/course/大三秋/基因组学数据分析/final/mia_markers_epithelial.csv")

markers <- c('CAV1','CLDN18','CRABP2','SCGB1A1','TM4SF1','CAPS','UBE2C')
DotPlot(epithelial_subset,features = markers, pt.size = 0)+theme(axis.text.x = element_text(size = 9))

#UMAP图
FeaturePlot(epithelial_subset,features = markers)
epithelial_subset <- subset(epithelial_subset, idents = c('5','9'), invert = TRUE)#去掉低质量细胞群
epithelial_subset.cluster.ids <- c('0'='Cancer cells',
                                   '1'='AT1',
                                   '2'='Cancer cells',
                                   '3'='AT2',
                                   '4'='Epi4',
                                   '6'='Epi6',
                                   '7'='Epi7',
                                   '8'='Epi8'
)

#加celltype注释
epithelial_subset <- RenameIdents(epithelial_subset, epithelial_subset.cluster.ids)                        
epithelial_subset$celltype <- epithelial_subset@active.ident

#可视化
DimPlot(epithelial_subset, group.by = "celltype",label = 'TRUE')+ 
  labs(title = "Epithelial cell subtype")+ 
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(hjust = 0.5))

#保存注释结果
saveRDS(epithelial_subset, file = "E:/course/大三秋/基因组学数据分析/final/mia_epithelial_ann.rds")

endothelial_subset <- subset(mia_ann, idents = "Endothelial cells")

#IAC
iac_ann <- readRDS("E:/course/大三秋/基因组学数据分析/final/iac_ann.rds")
epithelial_subset <- subset(iac_ann, idents = "Epithelial cells")
epithelial_subset <- NormalizeData(epithelial_cells, normalization.method = "LogNormalize", scale.factor = 1e4)
epithelial_subset <- FindVariableFeatures(epithelial_subset, selection.method = "vst", nfeatures = 2000)
epithelial_subset <- ScaleData(epithelial_subset, features = VariableFeatures(epithelial_subset))
epithelial_subset <- RunPCA(epithelial_subset,features = VariableFeatures(epithelial_subset),npcs=50)
DimPlot(epithelial_subset, reduction = "umap", label = TRUE,pt.size = 1)
# 使用前30个PC进行Louvain聚类
epithelial_subset<- FindNeighbors(epithelial_subset, dims = 1:30)
epithelial_subset<- FindClusters(epithelial_subset, resolution = 0.4)

# 进行UMAP降维
epithelial_subset <- RunUMAP(epithelial_subset, dims = 1:30,min.dist=0.001,spread=1)

# 可视化簇
DimPlot(epithelial_subset, reduction = "umap", label = TRUE,pt.size = 1)+ 
  labs(title = "IAC Epithelial clusters")+ 
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(hjust = 0.5))

#去掉低质量细胞
epithelial_subset <- subset(epithelial_subset, idents = c('0','1','2','3','4','5','6','7','8'), invert = FALSE)

#识别每个聚类的差异基因
all.markers <- FindAllMarkers(epithelial_subset, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
markers__epithelial <- all.markers [all.markers $p_val_adj < 0.05 & all.markers $avg_log2FC > 0.25, ]
#保存到本地
write.csv(markers__epithelial, file = "E:/course/大三秋/基因组学数据分析/final/iac_markers_epithelial.csv")

markers <- c('CAV1','CLDN18','CRABP2','SCGB1A1','TM4SF1','CAPS','TMSF1')
VlnPlot(epithelial_subset,features = markers, pt.size = 0)+theme(axis.text.x = element_text(size = 9))

#UMAP图
FeaturePlot(epithelial_subset,features = markers)
epithelial_subset <- subset(epithelial_subset, idents = c('0','1','3','4'), invert = TRUE)#去掉低质量细胞群
epithelial_subset.cluster.ids <- c('0'='Cancer cells 2',
                                   '1'='Cancer cells 2',
                                   '2'='Cancer cells 1',
                                   '3'='Ciliated cells',
                                   '4'='AT1',
                                   '5'='Epi5',
                                   '6'='Epi6',
                                   '7'='Ciliated cells',
                                   '8'='Cancer cells 1',
                                   '9'='AT2',
                                   '10'='Cancer cells 1',
                                   '11'='Cancer cells 1'
)

#加celltype注释
epithelial_subset <- RenameIdents(epithelial_subset, epithelial_subset.cluster.ids)                        
epithelial_subset$celltype <- epithelial_subset@active.ident

#可视化
DimPlot(epithelial_subset, group.by = "celltype",label = 'TRUE')+ 
  labs(title = "Epithelial cell subtype")+ 
  theme(axis.text.x = element_text(size = 9),
        plot.title = element_text(hjust = 0.5))

#保存注释结果
saveRDS(epithelial_subset, file = "E:/course/大三秋/基因组学数据分析/final/iac_epithelial_ann.rds")

endothelial_subset <- subset(iac_ann, idents = "Endothelial cells")
