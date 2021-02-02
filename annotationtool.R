# Annotation tools
library(UniProt.ws)
library(stringr)


# Retrieve protein sequence from uniport based on protein ACC

annotationtool<-function(identifier,strippedseq,numextend){
up <- UniProt.ws(taxId=9606)
columns <- c("ENTRY-NAME","GENES","PROTEIN-NAMES","SEQUENCE")
kt <- "UNIPROTKB"
res <- UniProt.ws::select(up, identifier, columns, kt)

# Match entry sequence to protein sequence
i<-1
df<-data.frame(identifier=identifier, strippedseq=strippedseq,
               Nterm_start=NA,Cterm_end=NA, Nterm_non_prime_sequence=NA, Nterm_prime_sequence=NA,
               Cterm_non_prime_sequence=NA, Cterm_prime_sequence=NA,
               Nterm_seq_window=NA, Cterm_seq_window=NA,
               Preceding_AA=NA, Nterm_AA=NA,
               Cterm_AA=NA, Following_AA=NA)
df$strippedseq<-as.character(df$strippedseq)
for (i in 1:nrow(df)){
pos.v<-str_locate(res$SEQUENCE[i],df$strippedseq[i])
df$Nterm_start[i]<-pos.v[1]
df$Cterm_end[i]<-pos.v[2]
if(any(is.na(pos.v))){next}else{
df$Nterm_non_prime_sequence[i]<-str_sub(res$SEQUENCE[i],max(0,pos.v[1]-10),pos.v[1]-1)
df$Nterm_prime_sequence[i]<-str_sub(res$SEQUENCE[i],pos.v[1]+1,pos.v[1]+10)
df$Cterm_non_prime_sequence[i]<-str_sub(res$SEQUENCE[i],max(0,pos.v[2]-10),pos.v[2]-1)
df$Cterm_prime_sequence[i]<-str_sub(res$SEQUENCE[i],pos.v[2]+1,pos.v[2]+10)
df$Nterm_seq_window[i]<-str_sub(res$SEQUENCE[i], max(0,pos.v[1]-numextend), pos.v[1]+numextend)
df$Cterm_seq_window[i]<-str_sub(res$SEQUENCE[i], max(0,pos.v[2]-numextend), pos.v[2]+numextend)
df$Preceding_AA[i]<-str_sub(res$SEQUENCE[i],max(0,pos.v[1]-1),max(0,pos.v[1]-1))
if(pos.v[1]==1){df$Preceding_AA[i]="O"}
df$Nterm_AA[i]<-str_sub(res$SEQUENCE[i],pos.v[1],pos.v[1])
df$Cterm_AA[i]<-str_sub(res$SEQUENCE[i],pos.v[2],pos.v[2])
df$Following_AA[i]<-str_sub(res$SEQUENCE[i],pos.v[2]+1, pos.v[2]+1)
if (pos.v[2]==nchar(res$SEQUENCE[i])){df$Following_AA[i]="O"}
}
}
return(df)
}
