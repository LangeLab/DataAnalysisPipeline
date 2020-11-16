#tests
library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggsci)
library(reshape2)

eb.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  log2FC <- fit.eb$coefficients[, 2]
  p.mod <- fit.eb$p.value[, 2]
  results.eb <- data.frame(log2FC, p.mod)
  return(results.eb)
}

testingDEmethods<-function(data,methodin,FoI,DoE){
  ind<-which(DoE[,1] %in% colnames(data))
  DoE<-DoE[ind,]
  if (methodin=="limma"){
    design <- model.matrix(~DoE[,FoI])
    res.eb <- eb.fit(log2(data), design)
    pvalues<-res.eb$p.mod
    names<-rownames(data)
    log2FC<-res.eb$log2FC
  }
  if (methodin=="t-test"){
    pvalues<-numeric(nrow(data))
    log2FC<-numeric(nrow(data))
    i<-1
    for (i in 1:nrow(data)){
      ds<-melt(t(log2(data[i,])))[,-1]
      ds[,3]<-DoE[match(ds[,1],DoE[,1]),FoI]
      colnames(ds)<-c("sample","intensity","FoI")
      ts<-t.test(formula=intensity~FoI, data=ds)
      log2FC[i]<-diff(ts[["estimate"]])
      pvalues[i]<-p.adjust(ts$p.value,method="BH")
    }
    names<-rownames(data)
  }
  
  df<-data.frame(log2FC=log2FC,pvalue=pvalues, name=names)
  df<- df %>% mutate(significance=(pvalue<0.05)&(abs(log2FC)>1))
  df
}

testingDE<-function(data,methodin,FoI,flagblock,blockfactor){
  data<-na.omit(data)
  DoE[,FoI]<-as.factor(DoE[,FoI])
  if (flagblock==FALSE){
  df<-testingDEmethods(data,methodin,FoI,DoE)
  ggplot(df, aes(x=log2FC, y=-log10(pvalue),color=significance,
                 fill=significance, alpha=significance, label=name))+
  geom_point()+
  geom_vline(xintercept=1, linetype="dashed", color="darkgrey")+
  geom_vline(xintercept=-1, linetype="dashed", color="darkgrey")+
  geom_hline(yintercept=1.30103, linetype="dashed", color="darkgrey")+
  scale_alpha_manual(values=c(0.2, 1.0)) + 
  scale_color_npg()+
  geom_text_repel(data=filter(df, significance == TRUE & log2FC < 0),
                    force=0.5, nudge_x = 0.5, direction="y",
                    hjust=2,segment.size=0.1,size=3)+
  geom_text_repel(data=filter(df, significance == TRUE & log2FC > 0),
                   force=0.5, nudge_x = 2, direction="y",
                   hjust=0,segment.size=0.1,size=3)+
  labs(x=paste0("log2 fold change of ", FoI,": ",levels(DoE[,FoI])[2]," over ", levels(DoE[,FoI])[1]))+
  theme_classic()
  }else{
   if (blockfactor==FoI){print("Error: blocking factor cannot be the same as testing factor.")}else{
     grouping<-levels(as.factor(DoE[,blockfactor]))
     samples1<-DoE[which(DoE[,blockfactor]==grouping[1]),1]
     samples2<-DoE[which(DoE[,blockfactor]==grouping[2]),1]
     df1<-testingDEmethods(data[,samples1],methodin,FoI,DoE)
     df2<-testingDEmethods(data[,samples2],methodin,FoI,DoE)
     df<-data.frame(log2FC1=df1$log2FC, log2FC2=df2$log2FC,
                    pvalue1=df1$pvalue, pvalue2=df2$pvalue, name=df1$name)
     df <- df %>% mutate(significance=case_when(
       (pvalue1 < 0.05) & (abs(log2FC1)>1) &(!((pvalue2 < 0.05) & (abs(log2FC2)>1))) ~ "significance in block 1", 
       (pvalue2 < 0.05) & (abs(log2FC2)>1) &(!((pvalue1 < 0.05) & (abs(log2FC1)>1))) ~ "significance in block 2",
       (pvalue1 < 0.05) & (abs(log2FC1)>1) & (pvalue2 < 0.05) & (abs(log2FC2)>1) ~ "significance in both blocks"))
     df$significance[is.na(df$significance)]<-"no significance"
     df$significance <-factor(df$significance,levels=(c("significance in block 1", "significance in block 2", "significance in both blocks", "no significance")))
    
     ggplot(df, aes(x=log2FC1,y=log2FC2,color=significance, alpha=significance, label=name))+
       geom_point()+
       geom_vline(xintercept=1, linetype="dashed", color="darkgrey")+
       geom_vline(xintercept=-1, linetype="dashed", color="darkgrey")+
       geom_hline(yintercept=1, linetype="dashed", color="darkgrey")+
       geom_hline(yintercept=-1, linetype="dashed", color="darkgrey")+
       scale_alpha_manual(values=c(1.0, 1.0, 1.0, 0.2)) + 
       scale_color_npg()+
       geom_text_repel(data=filter(df, significance == "significance in both blocks" & log2FC1 < 0),
                       force=0.5, nudge_x = 0.5, direction="y",
                       hjust=2,segment.size=0.1,size=3)+
       geom_text_repel(data=filter(df, significance == "significance in both blocks" & log2FC1 > 0),
                       force=0.5, nudge_x = 2, direction="y",
                       hjust=0,segment.size=0.1,size=3)+
       labs(y="log2 fold change on block 1",x="log2 fold change on block 2")+theme_classic() 
     
    }
    
  }
}
