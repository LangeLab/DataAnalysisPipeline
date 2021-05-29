library(openxlsx)

setdoe<-function(meta_data){
  # set up DoE matrix
  DoE<-data.frame(meta_data[!is.na(meta_data[,1]),])
  return(DoE)
}

inputraw<-function(raw_data, DoE, col_id){
if(is.null(raw_data)){return()}
raw_data<-raw_data[!is.na(raw_data[,col_id]),]

rownames(raw_data)<-raw_data[,col_id]
raw_data[,col_id]<-NULL

# fix NaNs
raw_data[raw_data=="NaN"]<-NA
raw_data[raw_data=="Filtered"]<-NA

# fix characters to numeric
fixed_data<-sapply(raw_data[,DoE[,1]], as.numeric)
rownames(fixed_data)<-rownames(raw_data)

# filter out rows that are completely missing
i<-1
ind <- numeric()
for (i in 1:nrow(fixed_data)){ 
  if (all(is.na(fixed_data[i,]))){
    ind <- c(ind, i)
  }
}
if(length(ind)!=0){
  fixed_data <- fixed_data[-ind, ]
}

ind.tf<-(!(colnames(raw_data)%in%DoE[,1]))
return(list(data=data.frame(fixed_data), other_annotation=cbind(rownames(fixed_data),data.frame(raw_data[,ind.tf]))))
}



averagereplica<-function(fixed_data, replica_list){
  if(is.null(fixed_data)){return()}
  ls<-replica_list[match(replica_list[,1],colnames(fixed_data)),2]
  aa<-t(aggregate(t(fixed_data), by=list(ls),mean,na.action=na.pass, na.rm=TRUE))
  df<-data.frame(aa[-1,])
  colnames(df)<-aa[1,]
  df<-sapply(df, as.numeric)
  rownames(df)<-rownames(fixed_data)
  return(df)
}
                   