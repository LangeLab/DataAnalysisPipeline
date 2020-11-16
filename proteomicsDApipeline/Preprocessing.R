# pre-processing functions
library(MSnbase)
library(limma)
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

preprocessing<-function(filterlevel, flag.normalize, flag.impute, impute.method){
dat<-fixed_data
#filter
dat<-data.frame(dat[apply(dat,1,function(x) sum(!is.na(x))>filterlevel/100*length(x)),])
#normalization
if(flag.normalize==TRUE){
dat<-normalizeBetweenArrays(dat, method="scale")}
#impute functions
if (flag.impute==TRUE){
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
dat
}
