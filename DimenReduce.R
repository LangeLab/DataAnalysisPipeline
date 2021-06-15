# Dimensionality Reduction methods function
library(factoextra)
library(Rtsne)
library(ggplot2)
library(ggsci)
library(M3C)

dimen.reduce<-function(data, DoE, method, colorfactor, tSNEper=NULL, rows){
  if (is.null(rows)){return()}
    data<-data[rows,]
    if (nrow(na.omit(data))<2){return()}
  if (method=="PCA"){
    m<-log2(t(na.omit(data)))
    res.pca <- prcomp(m,scale. = TRUE)
    g<-fviz_pca_ind(res.pca,
                 geom=c("point","text"), palette="npg",
                 habillage= DoE[,colorfactor], # color by groups
                 addEllipses = TRUE, # Concentration ellipses
                 repel = FALSE
    )
  }
  if (method=="t-SNE"){
    m<-log2(t(na.omit(data)))
    if (tSNEper==0){tSNEper<-round((nrow(m)-1)/3)-1}
    tsne <- Rtsne(m, dims = 2, perplexity=tSNEper, 
    theta=0, verbose=TRUE, max_iter = 1500, normalize=FALSE)
    df<-data.frame(tsne$Y)
    colnames(df)<-c("V1","V2")
    df$group<-as.factor(DoE[,colorfactor])
    g<-ggplot(df, aes(x=V1, y=V2, color=group))+
      geom_point()+
      theme_classic()+scale_color_npg()
  }
  if (method=="UMAP"){
    m<-log2(na.omit(data))
    g<-umap(m,labels=as.factor(DoE[,colorfactor]))
  }
  return(g)
}
