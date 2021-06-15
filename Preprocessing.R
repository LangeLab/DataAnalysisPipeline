# pre-processing functions
library(MSnbase)
library(limma)
library(reshape2)
impute_data = function(df, condition_group) {
  downshift=3.5
  # df = data frame containing filtered 
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  i<-1
  df <- log10(df)
  for (i in 1:ncol(df)){
    if (any(is.na(df[,i]))){
      ind<-which(condition_group==condition_group[i])
      temp <- df[!is.na(df[,i]),ind]
      temp.sd <-sd(as.vector(as.matrix(data.frame(temp))),na.rm = TRUE)
      temp.mean = mean(as.vector(as.matrix(data.frame(temp))), na.rm=TRUE) - downshift *temp.sd
      n.missing = sum(is.na(df[,i]))
      df[is.na(df[,i]),i] <- rnorm(n.missing, mean = temp.mean, sd = temp.sd)
    }
  }
  10**df
}

normalization_by_condition<-function(proteindata=NULL,protein_anno=NULL,maindata,DoE,FoI){
  i<-1
  if (is.null(proteindata)){
    for (i in 1:nrow(maindata)){
      temp.m<-aggregate(t(maindata[i,]),by=list(DoE[,FoI]),mean,na.action=na.pass, na.rm=TRUE)
      maindata[i,]<-maindata[i,]/unlist(lapply(DoE[,FoI],function(x) sum((x==temp.m[,1])*temp.m[,2])))
    }
  }else{
    for (i in 1:nrow(proteindata)){
      temp.m<-aggregate(t(proteindata[i,]),by=list(DoE[,FoI]),mean,na.action=na.pass, na.rm=TRUE)
      proteindata[i,]<-proteindata[i,]/unlist(lapply(DoE[,FoI],function(x) sum((x==temp.m[,1])*temp.m[,2])))
    }
    i<-1
    for (i in 1:nrow(maindata)){
      temp.m<-aggregate(t(maindata[i,]),by=list(DoE[,FoI]),mean,na.action=na.pass, na.rm=TRUE)
      ind.pro<-which(rownames(proteindata)==protein_anno[i])
      if (length(ind.pro)==0){maindata[i,]<-NA}else{
        maindata[i,]<-maindata[i,]/unlist(lapply(DoE[,FoI],function(x) sum((x==temp.m[,1])*temp.m[,2])))/proteindata[ind.pro,]
      }
      }
    
  }
  ind.na<-which(apply(maindata,1,function(x) all(is.na(x))))
  return(maindata[-ind.na,])
}

preprocessing<-function(type,dat.ls,DoE,filterlevel,normalize.method,impute.method, FoI_norm, FoI_impute,col_protein_anno=NULL){
#filter
dat<-dat.ls[[type]]$data
dat<-data.frame(dat[apply(dat,1,function(x) sum(!is.na(x))>filterlevel/100*length(x)),])
#normalization
if(normalize.method=="median randomization"){
dat<-normalizeBetweenArrays(dat, method="scale")}
if(normalize.method=="normalization over samples under same condition"){
  dat<-normalization_by_condition(proteindata=NULL,protein_anno=NULL,dat,DoE,FoI_norm)
}
if(normalize.method=="normalization over same protein and samples under same condition"){
  dat<-normalization_by_condition(dat.ls[["protein data"]]$data,dat.ls[[type]]$other_annotation[,col_protein_anno],dat,DoE,FoI_norm)
}
#impute functions
na.ind<-NA
if (impute.method!="No imputation"){
na.ind<-which(is.na(dat))
if (impute.method=="Down-shifted Normal samples"){
dat<-impute_data(dat, DoE[,FoI_impute])}else{
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

#aggregatereplica<-function(dat,DoE,replicacol,naind){
#  dt<-aggregate(t(dat),list(DoE[,replicacol]),mean, na.rm=TRUE)
#  df<-data.frame(t(dt[,-1]))
#  colnames(df)<-dt[,1]
#  m<-matrix(0, ncol=ncol(dat),nrow=nrow(dat))
#  m[naind]<-1
#  naind<-which((t(aggregate(t(m),list(DoE[,replicacol]),sum)[,-1]))>0)
#  return(list(data=df,na.index=naind))
#}

# violin plot for filtering, normalization and imputation
plotviolin <- function(data, custom_title){
  data$Protein<-rownames(data)
  df<-reshape2::melt(data, id.vars="Protein")
  names(df)[names(df) == "variable"] <- "Sample"
  names(df)[names(df) == "value"] <- "Intensity"
  g<-ggplot(na.omit(df), aes(x = Sample, y = log2(Intensity))) + 
    geom_violin()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ggtitle(custom_title)
  return(g)
}

