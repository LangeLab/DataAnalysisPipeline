library(openxlsx)

inputraw<-function(raw_data){
raw_data<-raw_data[!is.na(raw_data[,1]),]
rownames(raw_data)<-raw_data[,1]
raw_data[,1]<-NULL

# fix NaNs
raw_data[raw_data=="NaN"]<-NA
raw_data[raw_data=="Filtered"]<-NA

# fix characters to numeric
fixed_data<-sapply(raw_data, as.numeric)
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
return(data.frame(fixed_data))
}

setdoe<-function(meta_data,fixed_data){
# set up DoE matrix
DoE<-data.frame(meta_data[na.omit(match(meta_data[,1],colnames(fixed_data))),])
return(DoE)
}

averagereplica<-function(fixed_data, replica_list){
  ls<-replica_list[match(replica_list[,1],colnames(fixed_data)),2]
  aa<-t(aggregate(t(fixed_data), by=list(ls),mean,na.action=na.pass, na.rm=TRUE))
  df<-data.frame(aa[-1,])
  colnames(df)<-aa[1,]
  df<-sapply(df, as.numeric)
  rownames(df)<-rownames(fixed_data)
  return(df)
}
                   