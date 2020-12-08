# Dimensionality Reduction methods function
library(factoextra)
library(Rtsne)
library(ggplot2)
library(ggsci)
library(M3C)

dimen.reduce<-function(data, DoE, method, colorfactor){
  if (method=="PCA"){
    m<-log10(t(na.omit(data)))
    res.pca <- prcomp(m,scale. = TRUE)
    g<-fviz_pca_ind(res.pca,
                 geom="point", palette="npg",
                 habillage= DoE[,colorfactor], # color by groups
                 addEllipses = TRUE, # Concentration ellipses
                 repel = FALSE
    )
  }
  if (method=="t-SNE"){
    m<-log10(t(na.omit(data)))
    tsne <- Rtsne(m, dims = 2, perplexity=round((nrow(m)-1)/3)-1, 
    theta=0, verbose=TRUE, max_iter = 1000, normalize=FALSE)
    df<-data.frame(tsne$Y)
    colnames(df)<-c("V1","V2")
    df$group<-as.factor(DoE[,colorfactor])
    g<-ggplot(df, aes(x=V1, y=V2, color=group))+
      geom_point()+
      theme_classic()+scale_color_npg()
  }
  if (method=="UMAP"){
    m<-log10(na.omit(data))
    g<-umap(m,labels=as.factor(DoE[,colorfactor]))
  }
  return(g)
}