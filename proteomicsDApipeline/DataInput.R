library(openxlsx)
raw_data <- read.xlsx("rawdata_example.xlsx",
                      colNames=TRUE, rowNames=FALSE)
raw_data<-raw_data[!is.na(raw_data[,1]),]
rownames(raw_data)<-raw_data[,1]
raw_data[,1]<-NULL
meta_data <- read.xlsx("metadata_example.xlsx",
                       colNames=TRUE)
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

# set up DoE matrix
DoE<-meta_data[match(meta_data[,1],colnames(fixed_data)),]

                   