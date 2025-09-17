library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)

MIA<-readRDS("E:/course/大三秋/基因组学数据分析/final/mia_epithelial_ann.rds")
IAC<-readRDS("E:/course/大三秋/基因组学数据分析/final/iac_epithelial_ann.rds")

MIA[["RNA"]]<-JoinLayers(MIA[["RNA"]])
IAC[["RNA"]]<-JoinLayers(IAC[["RNA"]])

MIA <- SetIdent(MIA, value = "MIA")
IAC <- SetIdent(IAC, value = "IAC")

cancer_mia <- MIA[, MIA@meta.data$celltype %in% c('Cancer cells')]
cancer_mia=CreateSeuratObject(counts = cancer_mia@assays$RNA@layers$counts,
                                     meta.data = cancer_mia@meta.data)
cancer_mia[["RNA"]]<-JoinLayers(cancer_mia[["RNA"]])
cancer_mia <- SetIdent(cancer_mia, value = "MIA")
cancer_iac <- IAC[, IAC@meta.data$celltype %in% c('Cancer cells 1','Cancer cells 2')]
cancer_iac=CreateSeuratObject(counts = cancer_iac@assays$RNA@layers$counts,
                              meta.data = cancer_iac@meta.data)
cancer_iac[["RNA"]]<-JoinLayers(cancer_iac[["RNA"]])
cancer_iac <- SetIdent(cancer_iac, value = "IAC")

combined_seurat <- merge(cancer_mia, y = cancer_iac, add.cell.ids = c("MIA", "IAC"))
combined_seurat[["RNA"]]<-JoinLayers(combined_seurat[["RNA"]])
combined_seurat<- NormalizeData(combined_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
deg_results <- FindMarkers(combined_seurat, ident.1 = "MIA", ident.2 = "IAC", test.use = "wilcox",, logfc.threshold = 0)

# 查看结果
head(deg_results)

# 保存结果
write.csv(deg_results, file = "E:/course/大三秋/基因组学数据分析/final/deg_results.csv")

#Volcano
deg_results <- as.data.frame(deg_results)
fold_change_threshold <- 1
EnhancedVolcano(deg_results,
                x = "avg_log2FC",
                y = "p_val_adj",
                title = 'Volcano Plot of Differential Expression',
                pCutoff = 0.1,     # 设置p值的阈值
                FCcutoff = fold_change_threshold, # 设置差异倍数阈值
                lab = rownames(deg_results),  #设置标签名称为基因名
                pointSize = 0.3,   #设置散点大小
                labSize = 4   #设置标签大小
)   
