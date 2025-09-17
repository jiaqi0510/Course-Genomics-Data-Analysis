#inferCNV 环境配置
Sys.setenv(JAGS_HOME="C:/Users/lenovo/AppData/Local/Programs/JAGS/JAGS-4.3.1")
library(rjags)
install.packages("rjags")
BiocManager::install("infercnv")
library(infercnv)
library(ggplot2)
library(Seurat)
library(data.table)
library(devtools)

#读取数据
AIS<-readRDS("E:/course/大三秋/基因组学数据分析/final/epithelial_subset.rds")
AIS <- merge(epithelial_subset, y = endothelial_subset, add.cell.ids = c("epithelial", "endothelial"), project = "MergedAIS")
AIS[["RNA"]]<-JoinLayers(AIS[["RNA"]])

MIA<-readRDS("E:/course/大三秋/基因组学数据分析/final/epithelial_subset.rds")
MIA <- merge(epithelial_subset, y = endothelial_subset, add.cell.ids = c("epithelial", "endothelial"), project = "MergedMIA")
MIA[["RNA"]]<-JoinLayers(MIA[["RNA"]])

IAC<-readRDS("E:/course/大三秋/基因组学数据分析/final/epithelial_subset.rds")
IAC <- merge(epithelial_subset, y = endothelial_subset, add.cell.ids = c("epithelial", "endothelial"), project = "MergedIAC")
IAC[["RNA"]]<-JoinLayers(IAC[["RNA"]])
#制作3个文件用于后续读取
dat <- GetAssayData(IAC,
                    layer = 'counts',assay = 'RNA')
groupinfo <- data.frame(v1 = colnames(dat),
                        v2 = IAC@meta.data$celltype)

library(AnnoProbe)
geneInfor <- annoGene(rownames(dat),"SYMBOL","human")
geneInfor <- geneInfor[!geneInfor$chr %in% c("chrM", "chrX", "chrY"), ]
#提取chr后后面的数字并转化为num，从而按这个num排序
#sub函数用于把chr替换为空
geneInfor$chr_num <- as.numeric(sub("chr", "", geneInfor$chr))
colnames(geneInfor)

geneInfor <- geneInfor[with(geneInfor,order(chr_num,start)),c(1,4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

#保留行名在geneInfor第一列中存在的行。
dat <- dat[rownames(dat) %in% geneInfor[,1],]
#match(x, table)：match函数返回x中每个元素在table中的位置索引。
#获得位置后对dat进行重新排序使其跟geneInfor中的顺序一致
dat <- dat[match(geneInfor[,1],rownames(dat)),]

#统一把”-“改成”_“,因为后续运行inferCNV的时候需要读取保存文档，若不修改则会出现错误。
colnames(dat) <- gsub("-", "_", colnames(dat))

#groupinfo数据处理
groupinfo$v1 <- gsub("-", "_", groupinfo$v1)
# 然后将第一列设置为行名
rownames(groupinfo) <- groupinfo[, 1]
groupinfo <- groupinfo[, -1]
groupinfo <- as.data.frame(groupinfo)
rownames(groupinfo) <-colnames(dat)


#geneInfor数据处理
head(geneInfor)
rownames(geneInfor) <- geneInfor[, 1]
rownames(geneInfor) <- gsub("-", "_", rownames(geneInfor))
geneInfor <- geneInfor[, -1]


#ref_group_names 这里的细胞是正常对照，然后跟其他的细胞比较
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = dat,
                                     annotations_file = groupinfo,
                                     delim = "\t",
                                     gene_order_file = geneInfor,
                                     ref_group_names = c("Endothelial cells")
)

infercnv_obj2 <- infercnv::run(infercnv_obj,
                               cutoff =  0.1, #smart-seq选择1,10X选择0.1
                               out_dir = "E:/course/大三秋/基因组学数据分析/final/iac_infercnv_output", # dir is auto
                               # 是否根据细胞注释文件的分组
                               # 对肿瘤细胞进行分组
                               # 影响read.dendrogram, 如果有多个细胞类型，且设置为TRUE，
                               # 后续的read.dendrogram无法执行
                               cluster_by_groups =  T, 
                               hclust_method = "ward.D2",# ward.D2 方法进行层次聚类
                               analysis_mode = "subclusters", # 默认是samples，推荐是subclusters
                               denoise = TRUE, # 去噪音
                               HMM = F,  ##特别耗时间,是否要去背景噪音
                               plot_steps = F, #不在每个步骤后生成图形。
                               leiden_resolution = "auto", #可以手动调参数
                               num_threads = 4 #多线程工作，加快速度
)
