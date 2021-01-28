# pre-processing functions
library(MSnbase)
library(limma)
library(reshape2)
impute_data = function(df, width=0.2 ,downshift=3.5) {
  # df = data frame containing filtered 
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  i<-1
  df <- log10(df)
  for (i in 1:ncol(df)){
    if (any(is.na(df[,i]))){
      temp <- df[!is.na(df[,i]),i]
      temp.sd <- width * sd(temp)
      temp.mean = mean(temp) - downshift * sd(temp)
      n.missing = sum(is.na(df[,i]))
      df[is.na(df[,i]),i] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
    }
  }
  10**df
}

preprocessing<-function(dat,DoE,filterlevel, flag.normalize, flag.impute, impute.method){
#filter
dat<-data.frame(dat[apply(dat,1,function(x) sum(!is.na(x))>filterlevel/100*length(x)),])
#normalization
if(flag.normalize==TRUE){
dat<-normalizeBetweenArrays(dat, method="scale")}
#impute functions
na.ind<-NA
if (flag.impute==TRUE){
na.ind<-which(is.na(dat))
if (impute.method=="Down-shifted Normal samples"){
dat<-impute_data(dat)}else{
  phenos<-AnnotatedDataFrame(DoE)
  features<-AnnotatedDataFrame(data.frame(name=rownames(dat)))
  rownames(features)<-rownames(dat)
  object<-MSnSet(as.matrix(log10(dat)), phenoData=phenos,featureData=features)
  msn.impute<- MSnbase::impute(object, method=impute.method)
  msn.impute<-t(as(msn.impute,"data.frame"))
  dat<-10**msn.impute}
}
list(data=data.frame(dat),na.index=na.ind)
}


# violin plot for filtering, normalization and imputation
plotviolin <- function(data){
  data$Protein<-rownames(data)
  df<-reshape2::melt(data, id.vars="Protein")
  names(df)[names(df) == "variable"] <- "Sample"
  names(df)[names(df) == "value"] <- "Intensity"
  g<-ggplot(na.omit(df), aes(x = Sample, y = log2(Intensity))) + 
    geom_violin()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(g)
}

