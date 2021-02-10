# unsupervised learning
library(gplots)
library(RColorBrewer)
library(factoextra)
library(ppclust)
library(ggsci)
library(e1071)
library(plot.matrix)

fcluster<-function(data, method, rows, whetherlabel, ncenter){
  if (is.null(rows)){return()}else{
  ds<-data[rows,]
  ds<-na.omit(ds)
  if (nrow(ds)<2){return()}
  if (method=="Hierarchical Clustering"){
  my_palette<-colorRampPalette(c("blue", "white", "red"))(n = 100)
  m<-as.matrix(t(apply(ds,1,scale)))
  colnames(m)<-colnames(ds)
  if (whetherlabel==TRUE){
  rowlabs<-rownames(ds)}else{
    rowlabs<-""
  }
  heatmap.2(m,trace="none",col=my_palette, margins=c(7,7),
               Rowv=FALSE, labRow = rowlabs,
               density.info = "none", cexRow=0.75, cexCol=0.75,)
  g<-recordPlot()
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
  }
  if (method=="Fuzzy Clustering"){
    res.fcm <- cmeans(t(na.omit(ds)), centers=ncenter)
    par(mar=c(5, 5, 2, 5))
    plot(t(res.fcm$membership), main="", col= colorRampPalette(brewer.pal(8, "Blues"))(11), cex=0.5, axis.col=list(side=1, las=2),
         breaks=seq(0,1,0.1), xlab="", ylab="Fuzzy Clustering membership")
    g<-recordPlot()
  }
  return(g)
  }
}
