library(tidyverse)
library(ggtern)
library(hrbrthemes)
library(ggtext)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(scales)
library(tidyr)
#install.packages("ggrepel")
#install.packages("matrixStats")
library(matrixStats)
library(ggrepel)
require(pals)
library(rdist)
library(rpart)
library(ggthemes)
library(ComplexHeatmap)
library(ggdendro)
#install.packages('ggdendro')
#导入检查过的文???
setwd("D:/ZLH/文章图表")
C0<-read.csv("C0_abd_triads_R1.csv",header=T)
cluster0<-cbind(C0[,4],C0[,5],C0[,6])
  head(cluster0)
  sc0<-rowSums(cluster0)
  max0 <- rowMaxs(cluster0)
  head(max0)
  aai<-cluster0[,1]
  bbi<-cluster0[,2]
  ddi<-cluster0[,3]
  clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
  head(clu0)
  row.names(clu0)<-C0[,1]
  colnames(clu0)<-c("AA","BB","DD", "size")
  clu0 <- clu0[rowSums(clu0) > 0.5,]
  head(clu0)
  write.csv(clu0,"C0.cbindout.csv")

#计算欧氏距离
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
#expectation_distance<-rdist(clu0,centers)
#?rdist

tridat<-clu0[,1:3]
rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]

head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))

#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"C0_tridat_clusterR1.csv")
tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
#tridats<-read.csv("C0_tridat_cluster.csv",header=T, row.names=1)
#head(tridats)
#tridat$group <- factor(tridat$group, levels = names(centers))
#绘图

head(tridat01)
#label=ifelse(tridat[,1] > 0.5 | tridat[,2] > 0.5 | tridat[,3] > 0.5, row.names(tridat),"")
#head(label)
#画图
#triadpos <- centers
#triadcol <- setNames(brewer.dark2(8)[-7], rownames(triadpos))
#triadcol['Central'] <- 'lightgrey'
#显示两边的线
#选择颜色
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
require(pals)
cubebasis <- pal.compress(brewer.accent(8))
show_col(cubebasis)
#coll<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#CCCCCC","#FF7F00","#F781BF","#A65628")
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#plot_data1=tridat01[which(tridat[,1]>0.6),]
#plot_data2=tridat01[which(tridat[,2]>0.6),]
#plot_data3=tridat01[which(tridat[,3]>0.6),]
tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
plot_data<-tridat01[which(tridat01[,4]=="Central"),]
plot_data1<-tridat01[which(tridat01[,4]!="Central"),]
#每个group都着色
tridat01<-as.data.frame(tridat01)
ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  # legend
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))


#balance vs. snRNA.gene
platte<-c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
"#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
"#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
"#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#cluster0单独画
AA<-read.csv("AA.clusterR1.csv",header=T,row.names=1)
BB<-read.csv("BB.clusterR1.csv",header=T,row.names=1)
DD<-read.csv("DD.clusterR1.csv",header=T,row.names=1)
clustmarker<-tridat01[which(tridat01[,4]=="Central"),]
head(clustmarker)
clu0<-read.csv(paste0("cluster0.cbindout.csv"), header=T, row.names=1)
row.names(clu0)<-row.names(AA)
aacluster<-merge(clustmarker, clu0, by="row.names",all.x=FALSE)
#提取cluster0的数据
cluster0<-cbind(aacluster$AA.y,aacluster$BB.y,aacluster$DD.y,aacluster$size.y)
colnames(cluster0)<-c("AA","BB","DD","size")
row.names(cluster0)<-aacluster[,1]
head(cluster0)
dim(cluster0)
#计算平衡情况
cluster0<-cluster0[,1:3]
rownames(cluster0) <- rownames(cluster0)
cluster0 <- data.frame(cluster0)
cluster0$group <- apply(cluster0, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
cluster0$size <- aacluster$size.y
head(cluster0)
cluster001<-matrix(unlist(cluster0[rownames(cluster0),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=cluster0$group, size=cluster0$size) %>% 
  magrittr::set_rownames(rownames(cluster0))
cluster001<-na.omit(cluster001)
#cluster001 <- apply(cluster001,2,as.character)
head(cluster001)
rowname<-row.names(cluster001)
head(rowname)
cluster001 <- apply(cluster001,2,as.character)
rownames(cluster001) <- rowname
write.csv(cluster001,"cluster0_C0balance0.csv")

#绘图
head(cluster001)
a0<-read.csv("cluster0_C0balance0.csv",header=T, row.names=1)
plot_data<-a0[which(a0[,4]=="Central"),]
plot_data1<-a0[which(a0[,4]!="Central"),]
p0<-ggtern(as.data.frame(a0), mapping = aes(AA, BB, DD, size)) +
  geom_mask() +
  ##geom_point(data=plot_data, aes(size=size), color="lightgrey", show.legend=F)  +
  geom_point(data=plot_data, color="lightgrey", size=0.1, show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  ##geom_point(data=plot_data1, aes(size=size), color="#DC143C", show.legend=T)  +
  geom_point(data=plot_data1, color="#DC143C", size=0.1, show.legend=T)  +
  scale_colour_manual(values=platte) +
  ##scale_size(0.1) +
  # legend
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
  axis.ticks=element_blank()) +
  ggtitle(paste0("Cluster 0")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0)) +
  theme_ggtern(base_size = 8, base_family = "")

#画堆积图统计数量
#准备数据
data<-read.csv("cluster_C0balance0.csv",header=T, row.names=1)
tree1<-rpart(group ~.,
             data=data,
             method="class",
             parms=list(split="gini"));#自行指定参数建立决策树【consume_class：输出变量名称】
group.number<-as.data.frame(table(data$group))
group.number$type<-replicate(nrow(group.number), "cluster0")
x1<-group.number[5,]
y1<-group.number[1:4,]
z1<-group.number[6:7,]
xyz1<-rbind.data.frame(x1,y1,z1)
colnames(xyz1)<-c("group","Number","type")
rownames(xyz1)<-c("1","2","3","4","5","6","7")
xyz1$group<-factor(xyz$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
head(xyz1)
#作图
coll<-c("#A6CEE3","#F781BF","#A65628","#4DAF4A","#984EA3","#FF7F00","#E41A1C","#377EB8")
d0<-ggplot(xyz1,aes(type,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="stack",show.legend=F, width=0.1) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  theme() +
  xlab(NULL) + ylab("Number (×1000)")+
  guides(fill=guide_legend(title=NULL))
#组图  
#plot<-ggtern::grid.arrange(p0, d0, nrow=1, ncol=2)


#循环来做其它几个cluster不对称图
#准备文件
clustmarker<-tridat01[which(tridat01[,4]=="Central"),]
head(clustmarker)
#准备数据
for(i in 1:21)
{
  clu0<-read.csv(paste0("cluster",i,".cbindout.csv"), header=T, row.names=1)
  row.names(clu0)<-row.names(AA)
  aacluster<-merge(clustmarker, clu0, by="row.names",all.x=FALSE)
  #提取cluster0的数据
  cluster0<-cbind(aacluster$AA.y,aacluster$BB.y,aacluster$DD.y,aacluster$size.y)
  colnames(cluster0)<-c("AA","BB","DD","size")
  row.names(cluster0)<-aacluster[,1]
  head(cluster0)
  dim(cluster0)
  #计算平衡情况
  cluster0<-cluster0[,1:3]
  rownames(cluster0) <- rownames(cluster0)
  cluster0 <- data.frame(cluster0)
  cluster0$group <- apply(cluster0, 1, function(x){
    rownames(centers)[which.min(cdist(t(x), centers))]
  })
  cluster0$size <- aacluster$size.y
  head(cluster0)
  cluster001<-matrix(unlist(cluster0[rownames(cluster0),c('AA','BB','DD')]), ncol=3)  %>%
    as.data.frame() %>% 
    magrittr::set_colnames(c('AA','BB','DD')) %>% 
    mutate(group=cluster0$group, size=cluster0$size) %>% 
    magrittr::set_rownames(rownames(cluster0))
  cluster001<-na.omit(cluster001)
  #cluster001 <- apply(cluster001,2,as.character)
  head(cluster001)
  rowname<-row.names(cluster001)
  cluster001 <- apply(cluster001,2,as.character)
  rownames(cluster001) <- rowname
  write.csv(cluster001,paste0("../cluster_C0balance",i,".csv"))
}

#画不对称图
"#DC143C",
colo<-c("#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
        "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
        "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
        "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
colo[1]
plot_list=list()
for(i in 1:21)
{
  a<-read.csv(paste0("cluster_C0balance",i,".csv"), header=T, row.names=1)
  #  label=ifelse(a[,7] > 0.5 | a[,8] > 0.5 | a[,9] > 0.5, row.names(a),"")
  plot_data<-a[which(a[,4]=="Central"),]
  plot_data1<-a[which(a[,4]!="Central"),]
  p<- ggtern(as.data.frame(a), mapping = aes(AA, BB, DD, size)) +
    geom_mask() +
    ##geom_point(data=plot_data, aes(size=size), color="lightgrey", show.legend=F)  +
    geom_point(data=plot_data, size=0.1, color="lightgrey", show.legend=F)  +
    #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
    ##geom_point(data=plot_data1, aes(size=size), color=colo[i], show.legend=T)  +
    geom_point(data=plot_data1, size=0.1, color=colo[i], show.legend=T)  +
    scale_colour_manual(values=platte) +
    #scale_size(range=c(0, 0.1)) +
    # legend
    guides(size="none") +
    theme_bw() +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank()) +
    ggtitle(paste0("Cluster ",i)) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0)) +
    theme_ggtern(base_size = 8, base_family = "")
  #  print(p)
  #  assign(paste0("plot_list",i), p)
  plot_list[[i]] <- p
  #ggsave(paste0(i,".pdf"))
}

#循环来做其它几个cluster柱形图
#画堆积图统计数量
#准备数据
dlot_list=list()
for(i in 1:21)
{data0<-read.csv(paste0("cluster_C0balance",i,".csv"),header=T, row.names=1)
tree<-rpart(group ~.,
             data=data0,
             method="class",
             parms=list(split="gini"));#自行指定参数建立决策树【consume_class：输出变量名称】
group.number0<-as.data.frame(table(data0$group))
group.number0$type<-replicate(nrow(group.number0), paste0("cluster",i))
x1<-group.number0[5,]
y1<-group.number0[1:4,]
z1<-group.number0[6:7,]
xyz1<-rbind.data.frame(x1,y1,z1)
colnames(xyz1)<-c("group","Number","type")
rownames(xyz1)<-c("1","2","3","4","5","6","7")
xyz1$group<-factor(xyz1$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
head(xyz1)
str(xyz1)
#作图
#coll<-c("#A6CEE3","#F781BF","#A65628","#4DAF4A","#984EA3","#FF7F00","#E41A1C","#377EB8")
d<-  ggplot(xyz1,aes(type,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="stack",show.legend=F, width=0.1) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  theme() +
  xlab(NULL) + ylab("Number (×1000)")+
  guides(fill=guide_legend(title=NULL))
dlot_list[[i]] <- d
}

#画图注
i=21
data0<-read.csv(paste0("cluster_C0balance",i,".csv"),header=T, row.names=1)
tree<-rpart(group ~.,
            data=data0,
            method="class",
            parms=list(split="gini"));#自行指定参数建立决策树【consume_class：输出变量名称】
group.number0<-as.data.frame(table(data0$group))
group.number0$type<-replicate(nrow(group.number0), paste0("cluster",i))
x1<-group.number0[5,]
y1<-group.number0[1:4,]
z1<-group.number0[6:7,]
xyz1<-rbind.data.frame(x1,y1,z1)
colnames(xyz1)<-c("group","Number","type")
rownames(xyz1)<-c("1","2","3","4","5","6","7")
xyz1$group<-factor(xyz1$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
head(xyz1)
str(xyz1)
d21<-ggplot(xyz1,aes(type,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="stack",show.legend=T, width=0.1) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  theme() +
  xlab(NULL) + ylab("Number (×1000)")+
  guides(fill=guide_legend(title=NULL))

#拼图
ggtern::grid.arrange(p0, d0,plot_list[[1]], dlot_list[[1]], plot_list[[2]], 
                           dlot_list[[2]],plot_list[[3]], dlot_list[[3]],
                           plot_list[[4]], dlot_list[[4]],
                     nrow=2, ncol=6)  %>%  ggsave("plot0.pdf",.,width=210,height=297, units="mm")
ggtern::grid.arrange(plot_list[[5]],dlot_list[[5]],
                           plot_list[[6]], dlot_list[[6]],plot_list[[7]], dlot_list[[7]],
                           plot_list[[8]], dlot_list[[8]],plot_list[[9]], dlot_list[[9]],
                           plot_list[[10]], dlot_list[[10]],plot_list[[11]], dlot_list[[11]],
                     nrow=3, ncol=6)  %>%  ggsave("plot1.pdf",.,width=210,height=297, units="mm")
ggtern::grid.arrange(plot_list[[12]], dlot_list[[12]],plot_list[[13]], dlot_list[[13]],
                           plot_list[[14]], dlot_list[[14]],plot_list[[15]], dlot_list[[15]],
                           plot_list[[16]], dlot_list[[16]],plot_list[[17]], dlot_list[[17]],
                     nrow=3, ncol=6)  %>%  ggsave("plot2.pdf",.,width=210,height=297, units="mm")
ggtern::grid.arrange(plot_list[[18]], dlot_list[[18]],plot_list[[19]], dlot_list[[19]],
                           plot_list[[20]], dlot_list[[20]],plot_list[[21]], dlot_list[[21]],d21,
                           #plot_list[[22]],
                     nrow=3, ncol=6)  %>%  ggsave("plot3.pdf",.,width=210,height=297, units="mm")
#ggsave(file="x.pdf",width=210,height=297, units="mm")


#marker基因与C0平衡基因
#读入marker基因列表
marker<- read.csv("D:/ZLH/snRNAseq/Cell_ranger/3ajbkR1/ajbk/新建文件夹/1dim49-2dim49-0.8-all.marker.csv",header=T, row.names=1)
head(marker)
#提取cluster的marker基因的表达量
total<-read.csv("ABD111_v1.1.csv",header=T)
head(total)
for(i in 0:21)
{clust0<-marker[which(marker[,6]==i),]
x0<-read.csv(paste0("cluster",i,".cbindout.csv"), header=T, row.names=1)
row.names(x0)<-row.names(AA)
test<-merge(clust0, x0, by="row.names",all.x=FALSE)
test<-cbind(test[,c(1,7,9,10,11,12)])
colnames(test)<-c("A","cluster","AA","BB","DD","size")
head(test)
row.names(x0)<-row.names(BB)
test0<-merge(clust0, x0, by="row.names",all.x=FALSE)
row.names(total)<-total[,2]
row.names(test0)<-test0[,1]
test0<-distinct(test0[,-1])
head(test0)
test0<-merge(total,test0, by="row.names",all.x=FALSE)
head(test0)
test0<-cbind(test0[,c(2,10,12,13,14,15)])
head(test0)
row.names(x0)<-row.names(DD)
test00<-merge(clust0, x0, by="row.names",all.x=FALSE)
row.names(total)<-total[,3]
row.names(test00)<-test00[,1]
test00<-distinct(test00[,-1])
head(test00)
test00<-merge(total,test00, by="row.names",all.x=FALSE)
test00<-cbind(test00[,c(2,10,12,13,14,15)])
head(test00)
cluste0<-rbind(test, test0, test00)
head(cluste0)
cluste0<- filter(cluste0,!duplicated(cluste0[,1]))
row.names(cluste0)<-cluste0[,1]
head(cluste0)
dim(cluste0)
#与C0平衡基因求交集
clustmarker<-tridat01[which(tridat01[,4]=="Central"),]
head(clustmarker)
xm<-merge(clustmarker, cluste0, by="row.names",all.x=FALSE)
#提取cluster0的数据
xn<-cbind(xm$AA.y,xm$BB.y,xm$DD.y,xm$size.y)
colnames(xn)<-c("AA","BB","DD","size")
row.names(xn)<-xm[,1]
head(xn)
dim(xn)
#计算平衡情况
size<-xn[,4]
xn<-xn[,1:3]
rownames(xn) <- rownames(xn)
xn <- data.frame(xn)
xn$group <- apply(xn, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
xn$size <- size
head(xn)
xn001<-matrix(unlist(xn[rownames(xn),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=xn$group, size=xn$size) %>% 
  magrittr::set_rownames(rownames(xn))
xn001<-na.omit(xn001)
#xn001 <- apply(xn001,2,as.character)
head(xn001)
rowname<-row.names(xn001)
head(rowname)
xn001 <- apply(xn001,2,as.character)
rownames(xn001) <- rowname
write.csv(xn001,paste0("marker_cluster",i,"_C0balance0.csv"))
}
#绘制图片
plot_marker=list()
for(i in 0:21)
{mm<-read.csv(paste0("marker_cluster",i,"_C0balance0.csv"),header=T, row.names=1)
mm$gene<-row.names(mm)
head(mm)
plot_data<-mm[which(mm[,4]=="Central"),]
plot_data1<-mm[which(mm[,4]!="Central"),]
#cluster0作图
pm0<-ggtern(as.data.frame(mm), mapping = aes(AA, BB, DD, size, gene)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  geom_point(data=plot_data, aes(size=size), color="#A6CEE3", alpha=0.8, show.legend=F) +
  geom_point(data=plot_data1, aes(size=size, color=group), show.legend=F)  +
  scale_colour_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
  #scale_size(range=c(0, 6)) +
  scale_size(range=c(0, 4)) +
  # legend
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("cluster",i)) +
  theme_ggtern(base_size = 8, base_family = "")
plot_marker[[i+1]] <- pm0
}
#组图，留出画堆积比例分布图的位置
ggtern::grid.arrange( plot_marker[[1]],  plot_marker[[2]], 
                            plot_marker[[3]],
                            plot_marker[[4]],  plot_marker[[5]],
                            plot_marker[[6]],  plot_marker[[7]],
                            plot_marker[[8]],  plot_marker[[9]],
                            plot_marker[[10]],  plot_marker[[11]],
                            plot_marker[[12]],  plot_marker[[13]],
                            plot_marker[[14]],  plot_marker[[15]],
                            plot_marker[[16]],  plot_marker[[17]],
                            plot_marker[[18]],  plot_marker[[19]],
                            plot_marker[[20]],  plot_marker[[21]],
                            plot_marker[[22]],
                           nrow=8, ncol=3)  %>%  ggsave("cluster marker gene without barplot.pdf",.,width=210,height=297, units="mm")


#将所有marker基因合并在一个数据框里
#循环
newdata=c()
for(i in 0:21)
{ a<-read.csv(paste0("marker_cluster",i,"_C0balance0.csv"), header=T, row.names=1)
  a$type<-replicate(nrow(a), paste0("cluster",i))
  newdata<-rbind(newdata,a)
}

#cluster0marker基因表达情况及其在其它类群中的表达情况
zxt=list()
for(i in 0:21)
{
  
  i=5
a<-read.csv(paste0("marker_cluster",i,"_C0balance0.csv"), header=T, row.names=1)
a$type<-replicate(nrow(a), paste0("cluster",i))
#cluster0不平衡的在其它类群中的表达情况
au<-a[which(a[,4]!="Central"),]
if(nrow(au)!=0)
{
#其它cluster循环做
neww=c()
for(j in 0:21)
{
if(j!=i){
b<-read.csv(paste0("cluster",j,".cbindout.csv"),header=T,row.names=1)
b$type<-replicate(nrow(b), paste0("cluster",j))
row.names(b)<-row.names(AA)
b0<-merge(au, b, by="row.names",all.x=FALSE)
bb0<-cbind(b0$AA.y,b0$BB.y,b0$DD.y,b0$size.y, b0$type.y)
colnames(bb0)<-c("AA","BB","DD","size","type")
row.names(bb0)<-b0[,1]
#计算平衡情况
bb0<-bb0[,1:3]
bb0<-data.frame(bb0)
bb0 <- as.data.frame(lapply(bb0,as.numeric))
bb0$group <- apply(bb0, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
bb0$size <- b0$size.y
bb0$type<-replicate(nrow(bb0), paste0("cluster",j))
bb0$gene<-b0[,1]
head(bb0)
neww<-rbind(neww,bb0)
}
  }
  #作图
  ##在cluster0中特异不平衡表达的基因
  head(neww)
  dim(neww)
  #统计数量
  #去掉不表达的行
  neww<-na.omit(neww)
  neww$group<-unlist(neww$group)
  neww_num<-table(neww$gene, neww$group)
  #加上基因在当前cluster中的类型
  aa0<-data.frame(row.names(au),au$group)
  row.names(aa0)<-row.names(au)
  colnames(aa0)<-c("gene","type")
  head(aa0)
  #通过合并的方法加上这些基因在cluster0中的类型
  neww_num<-as.data.frame(neww_num)
  colnames(neww_num)<-c("gene","group","Freq")
  aaa0<-merge(neww_num, aa0, by="gene",all.x=FALSE)
  aaa0$type<-factor(aaa0$type,levels=c("A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
  #绘制堆积柱形图
  #数据转换
  #neww_num<-as.data.frame(neww_num)
  #排序
  xx0<-aaa0[order(aaa0[,4], aaa0[,1]),]
  xx0$group<-factor(xx0$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
  xx0$type<-factor(xx0$type,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
  #neww_num$Var2<-factor(neww_num$Var2,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
  pp0<-ggplot(xx0,mapping = aes(x=gene,y=Freq,fill=group))+
    geom_bar(stat='identity',position='fill',show.legend = F) +
    labs(y = 'Frequnency (100%)') +
    scale_fill_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggtitle(paste0("other cluster",i)) +
    theme(panel.background = element_blank())#去除背景
  my<-replicate(nrow(xx0),"1")
  len_gene<-nrow(au)
  pq0<-ggplot(xx0,mapping=aes(x=gene,y=my, fill=type)) +
    geom_bar(stat='identity',position='fill',show.legend = F) +
    scale_fill_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggtitle(paste0("cluster",i," unbalance marker",len_gene)) +
    theme(panel.background = element_blank())#去除背景
  zxt[[i+1]]<-grid.arrange(pp0,pq0,nrow=2) 
}
}

#组图
grid.arrange(zxt[[1]],zxt[[2]], 
             zxt[[3]],
             zxt[[5]],zxt[[6]],
             zxt[[7]],zxt[[8]],
             zxt[[10]],
             zxt[[11]],zxt[[12]],
             zxt[[13]],zxt[[14]],
             zxt[[15]],zxt[[16]],
             zxt[[17]],zxt[[18]],
             zxt[[19]],zxt[[20]],
             zxt[[21]],zxt[[22]],
             nrow=8,ncol=3)     %>%  ggsave("cluster_marker_expression_in_other_cluster.pdf",.,width=210,height=297, units="mm")

#画总图
#pa<-newdata[which(newdata[,4]=="Central"),]
#pa0<-newdata[which(newdata[,4]!="Central"),]
#ggtern(as.data.frame(newdata), mapping = aes(AA, BB, DD, size, group,type)) +
#  geom_mask() +
#  geom_point(data=as.data.frame(pa), aes(size=size), color="#A6CEE3", show.legend=F)  +
#  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
#  geom_point(data=as.data.frame(pa0), aes(size=size, color=type), show.legend=F)  +
#  scale_colour_manual(values=platte) +
# # scale_size(range=c(0, 6)) +
#  # legend
#  guides(size="none") +
#  theme_bw() +
#  theme(axis.text=element_blank(),
#        axis.ticks=element_blank()) +
#  ggtitle(paste0("cluster expression")) +
#  theme(plot.title = element_text(hjust = 0.5, vjust = -2))
#cluster0平衡表达的基因在其它cluster中的表达情况作图
#总图
#plotb<-neww[which(neww[,4]=="Central"),]
#plotb0<-neww[which(neww[,4]!="Central"),]
#pn0<-ggtern(as.data.frame(neww), mapping = aes(AA, BB, DD, size, group,type)) +
#  geom_mask() +
#  geom_point(data=as.data.frame(plotb), aes(size=size), color="#CCCCCC", show.legend=F)  +
#  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
#  geom_point(data=as.data.frame(plotb0), aes(size=size, color=type), show.legend=F)  +
#  scale_colour_manual(values=platte[c(2:24)]) +
#  # scale_size(range=c(0, 6)) +
#  # legend
#  guides(size="none") +
#  theme_bw() +
#  theme(axis.text=element_blank(),
#        axis.ticks=element_blank()) +
#  ggtitle(paste0("cluster expression")) +
#  theme(plot.title = element_text(hjust = 0.5, vjust = -2))



#只画marker基因，不考虑C0的情况
#读入marker基因列表
marker<- read.csv("D:/ZLH/snRNAseq/Cell_ranger/3ajbkR1/ajbk/新建文件夹/1dim49-2dim49-0.8-all.marker.csv",header=T, row.names=1)
head(marker)
#提取cluster的marker基因的表达量
total<-read.csv("ABD111_v1.1.csv",header=T)
head(total)
for(i in 0:21)
{
i=3  
clust0<-marker[which(marker[,6]==i),]
x0<-read.csv(paste0("cluster",i,".cbindout.csv"), header=T, row.names=1)
row.names(x0)<-row.names(AA)
test<-merge(clust0, x0, by="row.names",all.x=FALSE)
test<-cbind(test[,c(1,7,9,10,11,12)])
colnames(test)<-c("A","cluster","AA","BB","DD","size")
head(test)
row.names(x0)<-row.names(BB)
test0<-merge(clust0, x0, by="row.names",all.x=FALSE)
row.names(total)<-total[,2]
row.names(test0)<-test0[,1]
test0<-distinct(test0[,-1])
head(test0)
test0<-merge(total,test0, by="row.names",all.x=FALSE)
head(test0)
test0<-cbind(test0[,c(2,10,12,13,14,15)])
head(test0)
row.names(x0)<-row.names(DD)
test00<-merge(clust0, x0, by="row.names",all.x=FALSE)
row.names(total)<-total[,3]
row.names(test00)<-test00[,1]
test00<-distinct(test00[,-1])
head(test00)
test00<-merge(total,test00, by="row.names",all.x=FALSE)
test00<-cbind(test00[,c(2,10,12,13,14,15)])
head(test00)
cluste0<-rbind(test, test0, test00)
head(cluste0)
cluste0<- filter(cluste0,!duplicated(cluste0[,1]))
row.names(cluste0)<-cluste0[,1]
head(cluste0)
dim(cluste0)
#计算平衡情况
size<-cluste0[,6]
xn<-cluste0[,3:5]
rownames(xn) <- rownames(cluste0)
xn <- data.frame(xn)
xn$group <- apply(xn, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
xn$size <- size
head(xn)
xn001<-matrix(unlist(xn[rownames(xn),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=xn$group, size=xn$size) %>% 
  magrittr::set_rownames(rownames(xn))
xn001<-na.omit(xn001)
#xn001 <- apply(xn001,2,as.character)
head(xn001)
rowname<-row.names(xn001)
head(rowname)
xn001 <- apply(xn001,2,as.character)
rownames(xn001) <- rowname
write.csv(xn001,paste0("marker_cluster",i,".csv"))
}
#绘制图片
plot_marker=list()
for(i in 0:21)
{mm<-read.csv(paste0("marker_cluster",i,".csv"),header=T, row.names=1)
mm$gene<-row.names(mm)
head(mm)
plot_data<-mm[which(mm[,4]=="Central"),]
plot_data1<-mm[which(mm[,4]!="Central"),]
#cluster0作图
pm0<-ggtern(as.data.frame(mm), mapping = aes(AA, BB, DD, size, gene)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  geom_point(data=plot_data, aes(size=size), color="#A6CEE3", alpha=0.8, show.legend=F) +
  geom_point(data=plot_data1, aes(size=size, color=group), show.legend=F)  +
  scale_colour_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
  #scale_size(range=c(0, 6)) +
  scale_size(range=c(0, 4)) +
  # legend
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("cluster",i)) +
  theme_ggtern(base_size = 8, base_family = "")
plot_marker[[i+1]] <- pm0
}
#组图，留出画堆积比例分布图的位置
ggtern::grid.arrange( plot_marker[[1]],  plot_marker[[2]], 
                      plot_marker[[3]],
                      plot_marker[[4]],  plot_marker[[5]],
                      plot_marker[[6]],  plot_marker[[7]],
                      plot_marker[[8]],  plot_marker[[9]],
                      plot_marker[[10]],  plot_marker[[11]],
                      plot_marker[[12]],  plot_marker[[13]],
                      plot_marker[[14]],  plot_marker[[15]],
                      plot_marker[[16]],  plot_marker[[17]],
                      plot_marker[[18]],  plot_marker[[19]],
                      plot_marker[[20]],  plot_marker[[21]],
                      plot_marker[[22]],
                      nrow=8, ncol=3)  %>%  ggsave("cluster marker gene without C0 and barplot.pdf",.,width=210,height=297, units="mm")
#堆积柱形图统计数量
dlot_list=list()
neww=c()
for(i in 0:21)
{
data0<-read.csv(paste0("marker_cluster",i,".csv"),header=T, row.names=1)
tree<-rpart(group ~.,
            data=data0,
            method="class",
            parms=list(split="gini"));#自行指定参数建立决策树【consume_class：输出变量名称】
group.number0<-as.data.frame(table(data0$group))
group.number0$type<-replicate(nrow(group.number0), paste0("cluster",i))
x1<-group.number0[5,]
y1<-group.number0[1:4,]
z1<-group.number0[6:7,]
xyz1<-rbind.data.frame(x1,y1,z1)
colnames(xyz1)<-c("group","Number","type")
rownames(xyz1)<-c("1","2","3","4","5","6","7")
xyz1$group<-factor(xyz1$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
head(xyz1)
neww<-rbind(neww,xyz1)
}
write.csv(neww,"cluster_marker_number.csv")
#作图
coll<-c("#A6CEE3","#F781BF","#A65628","#4DAF4A","#984EA3","#FF7F00","#E41A1C","#377EB8")
x<-as.data.frame(gsub('[cluster]', '', neww$type))
x<-x[,1]
neww$type1<-x
neww$type1<-factor(neww$type1,levels=c("0","1","2","3","4","5","6",
                                       "7","8","9","10","11","12","13",
                                       "14","15","16","17","18","19","20","21"))

d<-  ggplot(neww,aes(type1,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="stack",show.legend=T, width=0.1) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  theme() +
  xlab(NULL) + ylab("Gene Number (×1000)")+
  guides(fill=guide_legend(title=NULL))
#作百分比堆积图
#write.csv(neww,"cluster_marker_number.csv")
#excel中操作
neww<-read.csv("cluster_marker_number.csv",header=T, row.names=1)
head(neww)
neww$group<-factor(neww$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
ggplot(neww,aes(type1,rate*100,fill=as.factor(group))) +
  geom_bar(stat="identity",position="stack",show.legend=T, width=0.1) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  theme() +
  xlab(NULL) + ylab("Rate of gene Number (%)")+
  guides(fill=guide_legend(title=NULL))

##作每cluster不平衡表达基因在其它cluster中的情况
#cluster0marker基因表达情况及其在其它类群中的表达情况，不考虑C0
zxt=list()
for(i in 0:21)
{
  a<-read.csv(paste0("marker_cluster",i,".csv"), header=T, row.names=1)
  a$type<-replicate(nrow(a), paste0("cluster",i))
  #cluster0不平衡的在其它类群中的表达情况
  au<-a[which(a[,4]!="Central"),]
  if(nrow(au)!=0)
  {
    #其它cluster循环做
    neww=c()
    for(j in 0:21)
    {
        if(j!=i){
        b<-read.csv(paste0("cluster",j,".cbindout.csv"),header=T,row.names=1)
        b$type<-replicate(nrow(b), paste0("cluster",j))
        row.names(b)<-row.names(AA)
        b0<-merge(au, b, by="row.names",all.x=FALSE)
        bb0<-cbind(b0$AA.y,b0$BB.y,b0$DD.y,b0$size.y, b0$type.y)
        colnames(bb0)<-c("AA","BB","DD","size","type")
        row.names(bb0)<-b0[,1]
        #计算平衡情况
        bb0<-bb0[,1:3]
        if(ncol(as.data.frame(bb0))==3){
        bb0<-data.frame(bb0)}
        bb0 <- as.data.frame(lapply(bb0,as.numeric))
        bb0$group <- apply(bb0, 1, function(x){
          rownames(centers)[which.min(cdist(t(x), centers))]
        })
        bb0$size <- b0$size.y
        bb0$type<-replicate(nrow(bb0), paste0("cluster",j))
        bb0$gene<-b0[,1]
        head(bb0)
        neww<-rbind(neww,bb0)
      }
    }
    #作图
    ##在cluster0中特异不平衡表达的基因
    head(neww)
    dim(neww)
    #统计数量
    #去掉不表达的行
    neww<-na.omit(neww)
    neww$group<-unlist(neww$group)
    neww_num<-table(neww$gene, neww$group)
    #加上基因在当前cluster中的类型
    aa0<-data.frame(row.names(au),au$group)
    row.names(aa0)<-row.names(au)
    colnames(aa0)<-c("gene","type")
    head(aa0)
    #通过合并的方法加上这些基因在cluster0中的类型
    neww_num<-as.data.frame(neww_num)
    colnames(neww_num)<-c("gene","group","Freq")
    aaa0<-merge(neww_num, aa0, by="gene",all.x=FALSE)
    aaa0$type<-factor(aaa0$type,levels=c("A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
    #绘制堆积柱形图
    #数据转换
    #neww_num<-as.data.frame(neww_num)
    #排序
    xx0<-aaa0[order(aaa0[,4], aaa0[,1]),]
    xx0$group<-factor(xx0$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
    xx0$type<-factor(xx0$type,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
    write.csv(xx0,paste0(i,"_unique.csv"))
  }
}
    
    #neww_num$Var2<-factor(neww_num$Var2,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
    pp0<-ggplot(xx0,mapping = aes(x=gene,y=Freq,fill=group))+
      geom_bar(stat='identity',position='fill',show.legend = F) +
      labs(y = 'Frequnency (100%)') +
      scale_fill_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      ggtitle(paste0("other cluster",i)) +
      theme(panel.background = element_blank())#去除背景
    my<-replicate(nrow(xx0),"1")
    len_gene<-nrow(au)
    pq0<-ggplot(xx0,mapping=aes(x=gene,y=my, fill=type)) +
      geom_bar(stat='identity',position='fill',show.legend = F) +
      scale_fill_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      ggtitle(paste0("cluster",i," unbalance marker",len_gene)) +
      theme(panel.background = element_blank())#去除背景
    zxt[[i+1]]<-grid.arrange(pp0,pq0,nrow=2) 
  }
}

#组图
grid.arrange(zxt[[1]],zxt[[2]], 
             zxt[[3]],zxt[[4]],
             zxt[[5]],zxt[[6]],
             zxt[[7]],zxt[[8]],
             zxt[[9]],zxt[[10]],
             zxt[[11]],zxt[[12]],
             zxt[[13]],zxt[[14]],
             zxt[[15]],zxt[[16]],
             zxt[[17]],zxt[[18]],
             zxt[[19]],zxt[[20]],
             zxt[[21]],zxt[[22]],
             nrow=8,ncol=3)     %>%  ggsave("without_C0_cluster_marker_expression_in_other_cluster.pdf",.,width=210,height=297, units="mm")

#附表数据
nw <- data.frame()
for (i in 0:21) {
  x <- read.csv(paste0(i,"_unique.csv"))
  x$cluster <- rep(i,nrow(x))
  nw <- rbind(x,nw)
}
write.csv(nw,"附表7-3.csv")

#挑C0不平衡的基因在不同clsuter中的不平衡表达情况，不看marker
c0<-read.csv("cluster_C0balance0.csv",header=T,row.names=1)
c1<-read.csv("cluster_C0balance1.csv",header=T,row.names=1)
c2<-read.csv("cluster_C0balance2.csv",header=T,row.names=1)
c3<-read.csv("cluster_C0balance3.csv",header=T,row.names=1)
c4<-read.csv("cluster_C0balance4.csv",header=T,row.names=1)
c5<-read.csv("cluster_C0balance5.csv",header=T,row.names=1)
c6<-read.csv("cluster_C0balance6.csv",header=T,row.names=1)
c7<-read.csv("cluster_C0balance7.csv",header=T,row.names=1)
c8<-read.csv("cluster_C0balance8.csv",header=T,row.names=1)
c9<-read.csv("cluster_C0balance9.csv",header=T,row.names=1)
c10<-read.csv("cluster_C0balance10.csv",header=T,row.names=1)
c11<-read.csv("cluster_C0balance11.csv",header=T,row.names=1)
c12<-read.csv("cluster_C0balance12.csv",header=T,row.names=1)
c13<-read.csv("cluster_C0balance13.csv",header=T,row.names=1)
c14<-read.csv("cluster_C0balance14.csv",header=T,row.names=1)
c15<-read.csv("cluster_C0balance15.csv",header=T,row.names=1)
c16<-read.csv("cluster_C0balance16.csv",header=T,row.names=1)
c17<-read.csv("cluster_C0balance17.csv",header=T,row.names=1)
c18<-read.csv("cluster_C0balance18.csv",header=T,row.names=1)
c19<-read.csv("cluster_C0balance19.csv",header=T,row.names=1)
c20<-read.csv("cluster_C0balance20.csv",header=T,row.names=1)
c21<-read.csv("cluster_C0balance21.csv",header=T,row.names=1)
c0$type <- rep("c0",nrow(c0))
c1$type <- rep("c1",nrow(c1))
c2$type <- rep("c2",nrow(c2))
c3$type <- rep("c3",nrow(c3))
c4$type <- rep("c4",nrow(c4))
c5$type <- rep("c5",nrow(c5))
c6$type <- rep("c6",nrow(c6))
c7$type <- rep("c7",nrow(c7))
c8$type <- rep("c8",nrow(c8))
c9$type <- rep("c9",nrow(c9))
c10$type <- rep("c10",nrow(c10))
c11$type <- rep("c11",nrow(c11))
c12$type <- rep("c12",nrow(c12))
c13$type <- rep("c13",nrow(c13))
c14$type <- rep("c14",nrow(c14))
c15$type <- rep("c15",nrow(c15))
c16$type <- rep("c16",nrow(c16))
c17$type <- rep("c17",nrow(c17))
c18$type <- rep("c18",nrow(c18))
c19$type <- rep("c19",nrow(c19))
c20$type <- rep("c20",nrow(c20))
c21$type <- rep("c21",nrow(c21))
ct <- rbind(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21)
str(ct)
write.csv(ct,"total_cluster_C0balance.csv")

#挑cluster 不平衡的基因在不同clsuter中的不平衡表达情况，不看marker
c0<-read.csv("marker_cluster0.csv",header=T,row.names=1)
c1<-read.csv("marker_cluster1.csv",header=T,row.names=1)
c2<-read.csv("marker_cluster2.csv",header=T,row.names=1)
c3<-read.csv("marker_cluster3.csv",header=T,row.names=1)
c4<-read.csv("marker_cluster4.csv",header=T,row.names=1)
c5<-read.csv("marker_cluster5.csv",header=T,row.names=1)
c6<-read.csv("marker_cluster6.csv",header=T,row.names=1)
c7<-read.csv("marker_cluster7.csv",header=T,row.names=1)
c8<-read.csv("marker_cluster8.csv",header=T,row.names=1)
c9<-read.csv("marker_cluster9.csv",header=T,row.names=1)
c10<-read.csv("marker_cluster10.csv",header=T,row.names=1)
c11<-read.csv("marker_cluster11.csv",header=T,row.names=1)
c12<-read.csv("marker_cluster12.csv",header=T,row.names=1)
c13<-read.csv("marker_cluster13.csv",header=T,row.names=1)
c14<-read.csv("marker_cluster14.csv",header=T,row.names=1)
c15<-read.csv("marker_cluster15.csv",header=T,row.names=1)
c16<-read.csv("marker_cluster16.csv",header=T,row.names=1)
c17<-read.csv("marker_cluster17.csv",header=T,row.names=1)
c18<-read.csv("marker_cluster18.csv",header=T,row.names=1)
c19<-read.csv("marker_cluster19.csv",header=T,row.names=1)
c20<-read.csv("marker_cluster20.csv",header=T,row.names=1)
c21<-read.csv("marker_cluster21.csv",header=T,row.names=1)
c0$type <- rep("c0",nrow(c0))
c1$type <- rep("c1",nrow(c1))
c2$type <- rep("c2",nrow(c2))
c3$type <- rep("c3",nrow(c3))
c4$type <- rep("c4",nrow(c4))
c5$type <- rep("c5",nrow(c5))
c6$type <- rep("c6",nrow(c6))
c7$type <- rep("c7",nrow(c7))
c8$type <- rep("c8",nrow(c8))
c9$type <- rep("c9",nrow(c9))
c10$type <- rep("c10",nrow(c10))
c11$type <- rep("c11",nrow(c11))
c12$type <- rep("c12",nrow(c12))
c13$type <- rep("c13",nrow(c13))
c14$type <- rep("c14",nrow(c14))
c15$type <- rep("c15",nrow(c15))
c16$type <- rep("c16",nrow(c16))
c17$type <- rep("c17",nrow(c17))
c18$type <- rep("c18",nrow(c18))
c19$type <- rep("c19",nrow(c19))
c20$type <- rep("c20",nrow(c20))
c21$type <- rep("c21",nrow(c21))
mt <- rbind(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21)
str(mt)
write.csv(mt,"total_cluster_marker.csv")

coll=c("#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
       "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
       "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
       "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

#筛选原位杂交探针验证亚基因组表达不平衡——可视化
setwd("D:/ZLH/文章图表/验证")
bluk1 <- read.csv("marker.csv")
head(bluk1)
tridat01<-as.data.frame(bluk1[which(bluk1$cluster == 1),c(1,3,4,5,6,7)])
ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, gene_ID)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=gene_ID), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  # legend
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))


#TraesCS4D02G079500在cluster3和cluster6中的不对称性
#TraesCS4A02G234500	TraesCS4B02G080700	TraesCS4D02G079500
#cluster3
plot_marker=list()
for(i in 0:20)
{
c3 <- read.csv(paste0("cluster",i,".cbindout.csv"),row.names=1,header=T)
DD<-read.csv("DD.clusterR1.csv")
head(c3)
head(DD)
row.names(c3) <- DD[,1]
d<-c3["TraesCS4D02G079500",]
d
#分组
tridat<-d[,1:3]
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- d[,4]
head(tridat)
#cluster0作图
plot_data<-tridat[which(tridat[,4]=="Central"),]
plot_data1<-tridat[which(tridat[,4]!="Central"),]
#cluster0作图
pm0<-ggtern(as.data.frame(tridat), mapping = aes(AA, BB, DD, size, gene)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  geom_point(data=plot_data, aes(size=size), color="#A6CEE3", alpha=0.8, show.legend=F) +
  geom_point(data=plot_data1, aes(size=size, color=group), show.legend=F)  +
  scale_colour_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
  #scale_size(range=c(0, 6)) +
  scale_size(range=c(0, 4)) +
  # legend
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("cluster",i)) +
  theme_ggtern(base_size = 8, base_family = "")
plot_marker[[i+1]] <- pm0

}
i=21
c3 <- read.csv(paste0("cluster",i,".cbindout.csv"),row.names=1,header=T)
DD<-read.csv("DD.clusterR1.csv")
head(c3)
head(DD)
row.names(c3) <- DD[,1]
d<-c3["TraesCS4D02G079500",]
d
#分组
tridat<-d[,1:3]
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- d[,4]
head(tridat)
#cluster0作图
plot_data<-tridat[which(tridat[,4]=="Central"),]
plot_data1<-tridat[which(tridat[,4]!="Central"),]
pm0<-ggtern(as.data.frame(tridat), mapping = aes(AA, BB, DD, size, gene)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  geom_point(data=plot_data, aes(size=size), color="#A6CEE3", alpha=0.8, show.legend=T) +
  geom_point(data=plot_data1, aes(size=size, color=group), show.legend=F)  +
  scale_colour_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
  #scale_size(range=c(0, 6)) +
  scale_size(range=c(0, 4)) +
  # legend
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("cluster",i)) +
  theme_ggtern(base_size = 8, base_family = "")
plot_marker[[i+1]] <- pm0
#组图
grid.arrange(plot_marker[[1]],plot_marker[[2]], 
             plot_marker[[3]],plot_marker[[4]],
             plot_marker[[5]],plot_marker[[6]],
             plot_marker[[7]],plot_marker[[8]],
             plot_marker[[9]],plot_marker[[10]],
             plot_marker[[11]],plot_marker[[12]],
             plot_marker[[13]],plot_marker[[14]],
             plot_marker[[15]],plot_marker[[16]],
             plot_marker[[17]],plot_marker[[18]],
             plot_marker[[19]],plot_marker[[20]],
             plot_marker[[21]],plot_marker[[22]],
             nrow=8,ncol=3) %>%  ggsave("TraesCS4D02G079500.pdf",.,width=210,height=297, units="mm")


###marker基因不对称表达情况画环状柱形堆叠图，应该是所有基因，不只是marker基因
#准备数据
library( "tidyr")
library(tidyverse)
library(hrbrthemes)
library(kableExtra)
options(knitr.table.format = "html")
library(viridis)
coll <- c("#A6CEE3","#F781BF","#A65628","#4DAF4A","#984EA3","#FF7F00","#E41A1C")
dc=data.frame()
for(i in 0:21)
{
  data0<-read.csv(paste0("total_cluster",i,"0.csv"),header=T, row.names=1)
  #data0<-read.csv(paste0("cluster_C0balance",i,".csv"),header=T, row.names=1)
  data0<-na.omit(data0)
  tree<-rpart(group ~.,
              data=data0,
              method="class",
              parms=list(split="gini"));#自行指定参数建立决策树【consume_class：输出变量名称】
  group.number0<-as.data.frame(table(data0$group))
  group.number0$type<-replicate(nrow(group.number0), paste0("cluster",i))
  x1<-group.number0[5,]
  y1<-group.number0[1:4,]
  z1<-group.number0[6:7,]
  xyz1<-rbind.data.frame(x1,y1,z1)
  colnames(xyz1)<-c("group","Number","type")
  rownames(xyz1)<-c("1","2","3","4","5","6","7")
  xyz1$group<-factor(xyz1$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
  head(xyz1)
  str(xyz1)
  dc <- rbind(dc,xyz1)
}
head(dc)
dc$type <- factor(dc$type ,levels=c("cluster0","cluster1","cluster2","cluster3","cluster4","cluster5",
                                    "cluster6","cluster7","cluster8","cluster9","cluster10","cluster11",
                                    "cluster12","cluster13","cluster14","cluster15","cluster16","cluster17",
                                    "cluster18","cluster19","cluster20","cluster21"), ordered=TRUE)
p <- ggplot(dc,aes(type,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="stack",show.legend=T, width=0.8) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  theme() +
  xlab(NULL) + ylab("Number (×1000)")+
  guides(fill=guide_legend(title=NULL)) +
  ylim(-10, 15)
pdf("total_triad_环状柱形图.pdf")
p + theme(legend.position = "none",
          #axis.text = element_blank,
          #axis.title = element_blank,
          #panel.grid = element_blank,
          plot.margin = unit(rep(- 1, 4), "cm")) +
  coord_polar(start = 0) +
  annotate( "text", x = rep(max(dc$type), 4), y = c(0, 4, 8, 12), label = c( "0", "4", "8", "12") ,
            color= "grey", size= 2, angle= 0, fontface= "bold", hjust= 1)
dev.off()
#百分比堆积图
p <- ggplot(dc,aes(type,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="fill",show.legend=T, width=0.8) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  xlab(NULL) + ylab("Rate of gene Number (%)")+
  theme(panel.grid=element_blank()) +
  guides(fill=guide_legend(title=NULL))
  #coord_cartesian(ylim = c(-2, 1))
pdf("total_triad_环状柱形图R1.pdf")
p +  ylim(-0.5,20) +
  coord_polar(start = 0) +
  annotate( "text", x = rep(max(dc$type), 3), y = c(0, 0.5, 1), label = c( "0","0.5", "1") ,
            color= "grey", size= 2, angle= 0, fontface= "bold", hjust= 1)
dev.off()

write.csv(dc,"total_traid_ratio")


p + coord_polar() +
