#导入相关包
library(Seurat)
library(circlize)
library(nichenetr)
library(tidyverse)

#读取数据
AIS<-readRDS("E:/course/大三秋/基因组学数据分析/final/ais_ann.rds")
MIA<-readRDS("E:/course/大三秋/基因组学数据分析/final/mia_ann.rds")
IAC<-readRDS("E:/course/大三秋/基因组学数据分析/final/iac_ann.rds")
lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS("F:/ligand_target_matrix_nsga2r_final.rds")
weighted_networks<- readRDS("F:/weighted_networks_nsga2r_final.rds")

#查看数据
ligand_target_matrix[1:5,1:5]
head(lr_network)
head(weighted_networks$lr_sig) 
head(weighted_networks$gr)
AIS@meta.data %>% head()
AIS@meta.data$celltype %>% table()

#设置分组类型
Idents(AIS) <- 'celltype'
Idents(MIA) <- 'celltype'
Idents(IAC) <- 'celltype'

#receiver
receiver = c('T&NK cells')
expressed_genes_receiver = get_expressed_genes(receiver, MIA, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#sender
sender_celltypes = c('Myeloid cells', 'B cells', 'Endothelial cells','Epithelial cells','Fibroblast')
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, MIA, 0.10)
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, IAC, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()


markers_sig <- read.csv('E:/course/大三秋/基因组学数据分析/final/mia.significant.markers.csv', check.names = F)
#markers_sig <-markers_sig[which(markers_sig$cluster %in% ''),]

geneset_oi<-markers_sig$gene 
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands) 
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

#ligand activity DotPlot 可视化配体活性分析中排名靠前的配体在不同细胞中的表达情况
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(MIA, features = best_upstream_ligands %>% rev(), cols = 'RdYlBu') + RotatedAxis()+coord_flip()
DotPlot(IAC, features = best_upstream_ligands %>% rev(), cols = 'RdYlBu') + RotatedAxis()+coord_flip()

#ligand_pearson heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames('Pearson')
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot('Prioritized ligands','Ligand activity', color = 'darkorange',legend_position = 'top', x_axis_position = 'top', legend_title = 'Pearson correlation coefficient\ntarget gene prediction ability') + theme(legend.text = element_text(size = 9))
p_ligand_pearson

#ligand activity heatmap
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot('Prioritized ligands','Predicted target genes', color = 'purple',legend_position = 'top', x_axis_position = 'top',legend_title = 'Regulatory potential') + theme(axis.text.x = element_text(face = 'italic')) + scale_fill_gradient2(low = 'whitesmoke', high = 'purple', breaks = c(0,0.01,0.02,0.03,0.04),labels = c('0','0.01','0.02','0.03','0.04'))
p_ligand_target_network

#ligand receptor heatmap


