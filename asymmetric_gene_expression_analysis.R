library(tidyverse)
library(ggtern)
library(hrbrthemes)
library(ggtext)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(scales)
library(tidyr)
library(matrixStats)
library(ggrepel)
require(pals)
library(rdist)
library(rpart)
library(ggthemes)
library(ComplexHeatmap)
library(ggdendro)
###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
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
head(tridat01)
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
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")
tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
plot_data<-tridat01[which(tridat01[,4]=="Central"),]
plot_data1<-tridat01[which(tridat01[,4]!="Central"),]
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
###2. Group the bluk cells balance expressed triads in the 21 cell clusters
#generate the triads
cluster<-read.csv("pbmc.pesedo.bulk.csv",header=T,row.names=1)
homo<-read.csv("ABD111_v1.1.csv",header=T)
head(cluster)
head(homo)
a<-cbind(homo[,1],homo[,1])
head(a)
rownames(a)<-a[,1]
head(a)
aa<-merge(cluster, a, by="row.names",all.x=FALSE)
head(aa)
write.csv(aa,"AA.cluster.csv")
b<-cbind(homo[,2],homo[,2])
head(b)
rownames(b)<-b[,2]
head(b)
bb<-merge(cluster, b, by="row.names",all.x=FALSE)
head(bb)
write.csv(bb,"BB.cluster.csv")
c<-cbind(homo[,3],homo[,3])
head(c)
rownames(c)<-c[,2]
head(c)
dd<-merge(cluster, c, by="row.names",all.x=FALSE)
head(dd)
write.csv(dd,"DD.cluster.csv")
#use excel to check these datasets and reload
AA<-read.csv("AA.clusterR1.csv",header=T,row.names=1)
BB<-read.csv("BB.clusterR1.csv",header=T,row.names=1)
DD<-read.csv("DD.clusterR1.csv",header=T,row.names=1)
head(AA)
head(BB)
head(DD)
for(i in 0:21)
{
  cluster0<-cbind(AA[,i+2],BB[,i+2],DD[,i+2])
  head(cluster0)
  sc0<-rowSums(cluster0)
  max0 <- rowMaxs(cluster0)
  head(max0)
  aai<-cluster0[,1]
  bbi<-cluster0[,2]
  ddi<-cluster0[,3]
  clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
  head(clu0)
  row.names(clu0)<-AA[,1]
  colnames(clu0)<-c("AA","BB","DD", "size")
  head(clu0)
  write.csv(clu0,paste0("cluster",i,".cbindout.csv"))
}
#group the bluk cells balance expressed triads in the 21 cell clusters
clustmarker<-tridat01[which(tridat01[,4]=="Central"),]
head(clustmarker)
for(i in 0:21)
{
  clu0<-read.csv(paste0("cluster",i,".cbindout.csv"), header=T, row.names=1)
  row.names(clu0)<-row.names(AA)
  aacluster<-merge(clustmarker, clu0, by="row.names",all.x=FALSE)
  cluster0<-cbind(aacluster$AA.y,aacluster$BB.y,aacluster$DD.y,aacluster$size.y)
  colnames(cluster0)<-c("AA","BB","DD","size")
  row.names(cluster0)<-aacluster[,1]
  head(cluster0)
  dim(cluster0)
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
  head(cluster001)
  rowname<-row.names(cluster001)
  cluster001 <- apply(cluster001,2,as.character)
  rownames(cluster001) <- rowname
  write.csv(cluster001,paste0("cluster_C0balance",i,".csv"))
 }
#visualize
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
  data0<-read.csv(paste0("cluster_C0balance",i,"0.csv"),header=T, row.names=1)
  data0<-na.omit(data0)
  tree<-rpart(group ~.,
              data=data0,
              method="class",
              parms=list(split="gini"))
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
pdf("total_triad_circos.pdf")
p + theme(legend.position = "none",
          #axis.text = element_blank,
          #axis.title = element_blank,
          #panel.grid = element_blank,
          plot.margin = unit(rep(- 1, 4), "cm")) +
  coord_polar(start = 0) +
  ylim(-0.5,20) +
  coord_polar(start = 0) +
  annotate( "text", x = rep(max(dc$type), 3), y = c(0, 0.5, 1), label = c( "0","0.5", "1") ,
            color= "grey", size= 2, angle= 0, fontface= "bold", hjust= 1)
dev.off()
pdf("total_triad_barplot.pdf")
p <- ggplot(dc,aes(type,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="fill",show.legend=T, width=0.8) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  xlab(NULL) + ylab("Rate of gene Number (%)")+
  theme(panel.grid=element_blank()) +
  guides(fill=guide_legend(title=NULL))
  #coord_cartesian(ylim = c(-2, 1))
dev.off()

###2.draw the ternary diagram for the expression of triads in each cell cluster
AA<-read.csv("AA.clusterR1.csv",header=T,row.names=1)
BB<-read.csv("BB.clusterR1.csv",header=T,row.names=1)
DD<-read.csv("DD.clusterR1.csv",header=T,row.names=1)
dlot_list=list()
neww=c()
for(i in 0:21)
{
clu0<-read.csv(paste0("cluster",i,".cbindout.csv"), header=T, row.names=1)
row.names(clu0)<-row.names(AA)
cluster0 <- clu0
head(cluster0)
dim(cluster0)
#group the triad
cluster0<-cluster0[,1:3]
rownames(cluster0) <- rownames(cluster0)
cluster0 <- data.frame(cluster0)
cluster0$group <- apply(cluster0, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
cluster0$size <- cluster0$size
head(cluster0)
cluster001<-matrix(unlist(cluster0[rownames(cluster0),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=cluster0$group, size=cluster0$size) %>% 
  magrittr::set_rownames(rownames(cluster0))
rowname<-rownames(cluster001)
cluster001 <- apply(cluster001,2,as.character)
write.csv(cluster001,paste0("total_cluster",i,"0.csv"))
cluster001<-read.csv(paste0("total_cluster",i,"0.csv"),header=T,row.names=1)
cluster001<-na.omit(cluster001)
tree<-rpart(group ~.,
            data=cluster001,
            method="class",
            parms=list(split="gini"))
group.number0<-as.data.frame(table(cluster001$group))
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
write.csv(neww,"total_homologygene_cluster_marker_numberR1.csv")
#visualize
neww <- read.csv("total_homologygene_cluster_marker_numberR1.csv",header=T)
head(neww)
coll<-c("#A6CEE3","#F781BF","#A65628","#4DAF4A","#984EA3","#FF7F00","#E41A1C","#377EB8")
x<-as.data.frame(gsub('[cluster]', '', neww$type))
x<-x[,1]
neww$type1<-x
neww$type1<-factor(neww$type1,levels=c("0","1","2","3","4","5","6",
                                       "7","8","9","10","11","12","13",
                                       "14","15","16","17","18","19","20","21"), ordered=TRUE)
neww$group <-factor(neww$group ,levels=c("Central","A.dominant","A.suppressed",
                                         "B.dominant","B.suppressed","D.dominant",
                                         "D.suppressed"), ordered=TRUE)
d1<-  ggplot(neww,aes(type1,Number/1000,fill=as.factor(group))) +
  geom_bar(stat="identity",position="stack",show.legend=T, width=0.8) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  theme() +
  xlab(NULL) + ylab("Gene Number (×1000)")+
  theme(panel.grid=element_blank()) +
  guides(fill=guide_legend(title=NULL))
neww$group<-factor(neww$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
d2<-ggplot(neww,aes(type1,Number,fill=as.factor(group))) +
  geom_bar(stat="identity",position="fill",show.legend=T, width=0.8) +
  theme_bw() +
  scale_fill_manual(values=coll) +
  #scale_fill_solarized()
  xlab(NULL) + ylab("Rate of gene Number (%)")+
  theme(panel.grid=element_blank()) +
  guides(fill=guide_legend(title=NULL))
#merge the figures
grid.arrange(d1,d2,
             nrow=2,ncol=1)     %>%  ggsave("Number_all_homology_in_all_clusterR1.pdf",.,width=210,height=400, units="mm")

###3. Expression of triad homologs that are unbalanced in one cell population in other cell populations
total<-read.csv("ABD111_v1.1.csv",header=T)
head(total)
for(i in 0:21)
{
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
head(xn001)
rowname<-row.names(xn001)
head(rowname)
xn001 <- apply(xn001,2,as.character)
rownames(xn001) <- rowname
write.csv(xn001,paste0("marker_cluster",i,".csv"))
}
zxt=list()
for(i in 0:21)
{
  a<-read.csv(paste0("marker_cluster",i,".csv"), header=T, row.names=1)
  a$type<-replicate(nrow(a), paste0("cluster",i))
  au<-a[which(a[,4]!="Central"),]
  if(nrow(au)!=0)
  {
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
    #visualize
    neww<-na.omit(neww)
    neww$group<-unlist(neww$group)
    neww_num<-table(neww$gene, neww$group)
    aa0<-data.frame(row.names(au),au$group)
    row.names(aa0)<-row.names(au)
    colnames(aa0)<-c("gene","type")
    head(aa0)
    neww_num<-as.data.frame(neww_num)
    colnames(neww_num)<-c("gene","group","Freq")
    aaa0<-merge(neww_num, aa0, by="gene",all.x=FALSE)
    aaa0$type<-factor(aaa0$type,levels=c("A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
    xx0<-aaa0[order(aaa0[,4], aaa0[,1]),]
    xx0$group<-factor(xx0$group,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
    xx0$type<-factor(xx0$type,levels=c("Central","A.dominant","A.suppressed","B.dominant","B.suppressed","D.dominant","D.suppressed"))
    write.csv(xx0,paste0(i,"_unique.csv"))
  }
}
    pp0<-ggplot(xx0,mapping = aes(x=gene,y=Freq,fill=group))+
        geom_bar(stat='identity',position='fill',show.legend = F) +
        labs(y = 'Frequnency (100%)') +
        scale_fill_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
        theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      ggtitle(paste0("other cluster",i)) +
      theme(panel.background = element_blank())
    my<-replicate(nrow(xx0),"1")
    len_gene<-nrow(au)
    pq0<-ggplot(xx0,mapping=aes(x=gene,y=my, fill=type)) +
      geom_bar(stat='identity',position='fill',show.legend = F) +
      scale_fill_manual(values=c("Central"="#A6CEE3","A.dominant"="#F781BF","A.suppressed"="#A65628","B.dominant"="#4DAF4A","B.suppressed"="#984EA3","D.dominant"="#FF7F00","D.suppressed"="#E41A1C")) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      ggtitle(paste0("cluster",i," unbalance marker",len_gene)) +
      theme(panel.background = element_blank())
    zxt[[i+1]]<-grid.arrange(pp0,pq0,nrow=2) 
  }
}
#merge barplot
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
             nrow=8,ncol=3)     %>%  ggsave("cluster_marker_expression_in_other_cluster.pdf",.,width=210,height=297, units="mm")
#visualize using ggttern
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
#merge figures
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
                      nrow=8, ncol=3)  %>%  ggsave("cluster_marker_expression_in_other_cluster_ggttern.pdf",.,width=210,height=297, units="mm")
