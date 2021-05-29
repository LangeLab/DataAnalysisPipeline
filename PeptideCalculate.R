# Calculate intensities from peptides
calIntPep<-function(pepdata, sumform, DoE, col.ACC){
  col.sample<-Doe[,1]
  # grouping by ACC
  pepdata[,col.ACC]<-unlist(lapply(pepdata[,col.ACC], function(x) unlist(strsplit(as.character(x), split=";"))[1]))
  pepdata[pepdata=="NaN"]<-NA
  pepdata[,col.sample]<-sapply(pepdata[,col.sample], as.numeric)
  pepdata<-pepdata[!apply(pepdata[,col.sample],1,function(x) all(is.na(x))),]
  ACC<-unique(pepdata[,col.ACC])
  prodata<-data.frame(matrix(NA,nrow=length(ACC),ncol=length(col.sample)+1, 
                             dimnames=list(NULL, c("ProteinACC", colnames(pepdata)[col.sample]))))
  prodata$ProteinACC<-ACC
  i<-1
  for (i in 1:nrow(prodata)){
  ind<-which(pepdata[,col.ACC]==prodata$ProteinACC[i])
  # Sum of all peptide intensities
  if (sumform=="sum"){
  prodata[i,-1]<-colSums(pepdata[ind,col.sample], na.rm=TRUE)}
  # Weighted sum (most intense peptide get highest weight)
  if (sumform=="weighted"){
    weights<-rep(NA,length(ind))
    
  }
  # Sum of top 3 most intense peptides
  if (sumform=="top3"){
    temp.d<-pepdata[ind,col.sample]
    temp.d[is.na(temp.d)]<-0
    temp.d<-unlist(apply(temp.d, 2, function(x) tail(sort(x),3)))
    if (!is.null(nrow(temp.d))){
    prodata[i,-1]<-colSums(temp.d)}else{
    prodata[i,-1]<-temp.d
    }
  }
  }
  prodata[prodata==0]<-NA
  return(prodata)
}