# create the Samples csv
df<-data.frame(IDs=c(1:nrow(DoE)), shortname=DoE[,1],
               longname=NA)
#write.csv(df,file="Samples.csv",row.names=FALSE)

# create the Group csv
i<-2
n<-0
shortname<-character()
for (i in 2:ncol(DoE)){
n<-length(levels(as.factor(DoE[,i])))+n
shortname<-c(shortname, as.character(levels(as.factor(DoE[,i]))))
}
df1<-data.frame(ID=c(1:n), 
               shortname=shortname,
               longname=NA, description=NA)
#write.csv(df,file="Group.csv",row.names=FALSE)

# create the Group_has_Samples csv
i<-1
df2<-data.frame(Group_ID=NULL, Samples_IDs=NULL)
for (i in 1:nrow(meta_data)){
  s.id<-df$IDs[which(df$shortname==meta_data[i,1])]
  groups<-unname(meta_data[which(df$shortname==meta_data[i,1]),c(2:4)])
  g.id<-df1$ID[which(df1$shortname %in% groups)]  
  temp.df<-data.frame(Group_ID=g.id, Samples_IDs=s.id)
  df2<-rbind(df2,temp.df)
  
}
#write.csv(df2,file="SQL/Group_has_Samples.csv",row.names=FALSE)

# create the replica csv
df3<-data.frame(SampleID=unlist(lapply(df$IDs, function(x){rep(x, length(grep(df$shortname[x], colnames(raw_data))))})))
df3$ID<-c(1:nrow(df))
#write.csv(df,file="SQL/Replica.csv",row.names=FALSE)

# create protein instance csv
df4<-data.frame(ProteinID=NULL, Replica_ID=NULL)
i<-1
for (i in 1:nrow(df)){
  df4<-rbind(df4, data.frame(ProteinID=rep(c(1:nrow(raw_data))),
                             Replica_ID=df3$ID[i]))
}
df4$ID<-c(1:nrow(df4))
#write.csv(df4,file="SQL/ProteinInstance.csv",row.names=FALSE)

# create protein quant csv
df5<-data.frame(ProteinInstanceID=df4$ID,name=NA,intensity=NA)
i<-1
for (i in 1:nrow(df4)){
  ind.col<-as.character(df$shortname[df3$SampleID[df$ID[df4$Replica_ID[i]]]])
  ind.row<-df4$ProteinID[i]
  df5$name[i]<-as.character(rownames(raw_data)[ind.row])
  df5$intensity[i]<-raw_data[ind.row,ind.col]
}
df5$ID<-c(1:nrow(df5))
#write.csv(df5,file="SQL/ProteinQuant.csv",row.names=FALSE)



library(DBI)
library(RSQLite)
library(dplyr)
con<-dbConnect(RSQLite::SQLite(), "sql/aging.db" )
dbListTables(con)

#set up tables
dbSendStatement(con, '
CREATE TABLE Samples(
	IDs INT PRIMARY KEY,
   	shortname VARCHAR(255),
    longname VARCHAR(255)
)')
dbSendStatement(con, '
CREATE TABLE Groupdesign(
	ID INT PRIMARY KEY,
   	shortname VARCHAR(255),
    longname VARCHAR(255),
    descriptions VARCHAR (255)
)')
dbSendStatement(con, '
CREATE TABLE GroupHasSamples(
	Group_ID INT,
    Samples_IDs INT,
	PRIMARY KEY (Group_ID, Samples_IDs),
    FOREIGN KEY (Group_ID) 
      REFERENCES Groupdesign (ID),
	FOREIGN KEY (Samples_IDs) 
      REFERENCES Samples (IDs)
)')
dbSendStatement(con, '
CREATE TABLE Replica(
	ID INT PRIMARY KEY,
   	SampleID INT,
    FOREIGN KEY (SampleID) 
      REFERENCES Samples (IDs)
)')
dbSendStatement(con, '
CREATE TABLE ProteinInstance(
	ID INT PRIMARY KEY,
   	ProteinID INT,
    Replica_ID INT,
     FOREIGN KEY (Replica_ID) 
      REFERENCES Replica (ID)
)')
dbSendStatement(con, '
CREATE TABLE ProteinQuant(
	ID INT PRIMARY KEY,
   	ProteinInstanceID INT,
    name VARCHAR(255),
    intensity double,
     FOREIGN KEY (ProteinInstanceID) 
      REFERENCES ProteinInstance (ID)
)')
#write tables
dbWriteTable(conn = con, name = "Samples", value = df, append=TRUE)  
dbWriteTable(conn = con, name = "Group", value = df1, append=TRUE)  
dbWriteTable(conn = con, name = "Group_has_Samples", value = df2, append=TRUE)  
dbWriteTable(conn = con, name = "Replica", value = df3, append=TRUE)
dbWriteTable(conn = con, name = "ProteinInstance", value = df4, append=TRUE)
dbWriteTable(conn = con, name = "ProteinQuant", value = df5, append=TRUE)

res<-dbSendQuery(con, 'SELECT * FROM Samples')
dbFetch(res)
dbDisconnect(con)
