# unsupervised learning
library(gplots)
library(RColorBrewer)
library(factoextra)
library(ppclust)
library(ggsci)
library(e1071)
library(corrplot)

fcluster<-function(data, method, rows, whetherlabel, ncenter){
  if (is.na(rows)){print("no rows available")}else{
  ds<-data[rows,]
  if (method=="Hierarchical Clustering"){
  my_palette<-colorRampPalette(c("blue", "white", "red"))(n = 100)
  m<-as.matrix(t(apply(ds,1,scale)))
  colnames(m)<-colnames(ds)
  heatmap.2(m,trace="none",col=my_palette, margins=c(7,5),
               Rowv=FALSE, labRow = whetherlabel,
               density.info = "none", cexRow=0.75, cexCol=0.75,)
  }
  if (method=="k-means"){
    res.km<-kmeans(t(na.omit(ds)),centers=ncenter)
    g<-fviz_cluster(res.km, data=t(na.omit(ds)),
                 palette="npg",
                 geom = c("point", "text"),
                 repel=TRUE,
                 ellipse.type = "convex", 
                 ggtheme = theme_bw()
    )
    return(g)
  }
  if (method=="Fuzzy Clustering"){
    res.fcm <- cmeans(t(na.omit(ds)), centers=ncenter)
    corrplot(t(res.fcm$membership), is.corr = FALSE, 
             tl.cex=0.5,cl.pos="b",cl.length = 5,cl.ratio=0.5,cl.cex=0.5)
  }
  }
}
