library(Seurat)
library(monocle)
library(dplyr)
library(future)

#读取数据
AIS<-readRDS("E:/course/大三秋/基因组学数据分析/final/epithelial_ann.rds")
markers <- FindAllMarkers(AIS, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
significant.markers <- markers [markers $p_val_adj < 0.05 & markers $avg_log2FC > 0.25, ]
MIA<-readRDS("E:/course/大三秋/基因组学数据分析/final/mia_ann.rds")


#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(AIS@assays$RNA@layers$counts), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = AIS@meta.data)
fData <- data.frame(gene_short_name = row.names(AIS@assays$RNA@features), row.names = row.names(AIS@assays$RNA@features))
fd <- new('AnnotatedDataFrame', data = fData)

data <- as(as.matrix(MIA@assays$RNA@layers$counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = MIA@meta.data)
fData <- data.frame(gene_short_name = row.names(MIA@assays$RNA@features), row.names = row.names(MIA@assays$RNA@features))
fd <- new('AnnotatedDataFrame', data = fData)

#构建S4对象，CellDataSet
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

plan("multisession", workers = 4)
#估计sizefactor和dispersion
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

#过滤低质量细胞
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
head(pData(HSMM))

significantGenes <- read.csv("E:/course/大三秋/基因组学数据分析/final/AIS_annotation/significant.markers.csv")

# 使用DEGs排序细胞
HSMM <- setOrderingFilter(HSMM, significant.markers)

# 降维和排序细胞
HSMM <- reduceDimension(HSMM, method = "DDRTree",max_components = 2)
HSMM <- orderCells(HSMM)

# 绘制UMAP或tSNE图来可视化细胞排序结果
plot_cell_trajectory(HSMM, color_by = "celltype")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "state")
