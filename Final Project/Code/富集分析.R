#富集分析——GO
#使用DAVID网站对差异表达基因significant.markers.csv做富集分析，并将结果保存为excel
#加载包
library (dplyr)
library (ggplot2)  
library(tidyverse)
library(openxlsx)

#由细胞注释annotation文件中保存差异基因csv

#识别每个聚类的差异基因
all.markers <- FindAllMarkers(ais1, only.pos = T, min.pct = 0.5, logfc.threshold = 0.25)
significant.markers <- all.markers [all.markers $p_val_adj < 0.05 & all.markers $avg_log2FC > 0.25, ]
#writecsv保存到本地
write.csv(significant.markers, file = "E:/course/大三秋/基因组学数据分析/final/significant.markers.csv")
#然后使用DAVID网站上传significant.markers文件，得到富集通路，分为BP CC MF三大类，在网站中保存结果GOTERM_BP.xlsx,GOTERM_CC.xlsx,GOTERM_MF.xlsx


#导入富集通路数据
BP = read.xlsx('E:/course/大二春/生信实验/自主实验/Marker&Annotation&GO/PD_GO/GOTERM_BP_PD.xlsx',sep=',')
CC = read.xlsx('E:/course/大二春/生信实验/自主实验/Marker&Annotation&GO/PD_GO/GOTERM_CC_PD.xlsx',sep=',')
MF = read.xlsx('E:/course/大二春/生信实验/自主实验/Marker&Annotation&GO/PD_GO/GOTERM_MF_PD.xlsx',sep=',')
head(BP)

#DAVID获得的数据term列是GO:0005515~protein binding形式的，故需要将其拆分
BP = separate(BP,Term, sep="~",into=c("ID","Description"))
CC = separate(CC,Term, sep="~",into=c("ID","Description"))
MF = separate(MF,Term, sep="~",into=c("ID","Description"))

#提取各组数据需要展示的数量
display_number = c(15, 10, 5)  
go_result_BP = as.data.frame(BP)[1:display_number[1], ]
go_result_CC = as.data.frame(CC)[1:display_number[2], ]
go_result_MF = as.data.frame(MF)[1:display_number[3], ]

#将提取的各组数据进行整合
go_enrich = data.frame(
  ID=c(go_result_BP$ID, go_result_CC$ID, go_result_MF$ID),  #指定ego_result_BP、ego_result_CC、ego_result_MFID为ID                        
  Description=c(go_result_BP$Description,go_result_CC$Description,go_result_MF$Description),
  GeneNumber=c(go_result_BP$Count, go_result_CC$Count, go_result_MF$Count), #指定ego_result_BP、ego_result_CC、ego_result_MF的Count为GeneNumber
  type=factor(c(rep("Biological Process", display_number[1]), #设置biological process、cellular component、molecular function 的展示顺序
                rep("Cellular Component", display_number[2]),
                rep("Molecular Function", display_number[3])),
              levels=c("Biological Process", "Cellular Component","Molecular Function" )))

#设置GO term名字
for(i in 1:nrow(go_enrich)){
  description_splite=strsplit(go_enrich$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ")  
  go_enrich$Description[i]=description_collapse
  go_enrich$Description=gsub(pattern = "NA","",go_enrich$Description)  
}

#转成因子，防止重新排列
#go_enrich$type_order=factor(rev(as.integer(rownames(go_enrich))),labels=rev(go_enrich$Description))
go_enrich$type_order = factor(go_enrich$Description,levels=go_enrich$Description,ordered = T)
head(go_enrich)

#绘制柱状图
ggplot(go_enrich,
       aes(x=type_order,y=GeneNumber, fill=type)) +  
  geom_bar(stat="identity", width=0.8) +  
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + 
  xlab("GO term") + 
  ylab("Gene_Number") +  
  labs(title = "GO Terms Enrich")+ 
  theme_bw() +
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))

#分组气泡图
colnames(go_result_BP)[5] <- "GeneRatio"
ggplot(go_result_BP,  
       aes(x=GeneRatio,y=Description,color= PValue)) +  
  geom_point(aes(size=Count))  +
  ylim(rev(go_result_BP$Description)) +  
  labs(size="Genecounts",x="GeneRatio",y="GOTerm",title="GO_BP",color=expression(-log[10](PValue))) + 
  scale_color_gradient(low="red",high ="blue") +   
  theme(axis.text=element_text(size=10,color="black"),  
        axis.title = element_text(size=16),title = element_text(size=13))

colnames(go_result_CC)[5] <- "GeneRatio"
ggplot(go_result_CC,  
       aes(x=GeneRatio,y=Description,color= PValue)) +  
  geom_point(aes(size=Count))  +
  ylim(rev(go_result_CC$Description)) +  
  labs(size="Genecounts",x="GeneRatio",y="GOTerm",title="GO_CC",color=expression(-log[10](PValue))) + 
  scale_color_gradient(low="red",high ="blue") +   
  theme(axis.text=element_text(size=10,color="black"),  
        axis.title = element_text(size=16),title = element_text(size=13))

colnames(go_result_MF)[5] <- "GeneRatio"
ggplot(go_result_MF,  
       aes(x=GeneRatio,y=Description,color= PValue)) +  
  geom_point(aes(size=Count))  +
  ylim(rev(go_result_MF$Description)) +  
  labs(size="Genecounts",x="GeneRatio",y="GOTerm",title="GO_MF",color=expression(-log[10](PValue))) + 
  scale_color_gradient(low="red",high ="blue") +   
  theme(axis.text=element_text(size=10,color="black"),  
        axis.title = element_text(size=16),title = element_text(size=13))

