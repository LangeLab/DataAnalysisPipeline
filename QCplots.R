library(ggplot2)
library(reshape2)
library(dplyr)
library(ggsci)
library(ggpubr)
library(UpSetR)
library(ggcorrplot)
library(patchwork)
cvplots<-function(fixed_data, DoE){
  cvs <- apply(na.omit(fixed_data), 1, function(x) (sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)) * 100)
  g1<-ggplot(data.frame(name = "CV", CV = cvs), aes(x = CV, y = name)) + 
    geom_violin() + 
    geom_vline(xintercept = median(cvs,na.rm=TRUE), color = "red", linetype = "dashed") + 
    geom_text(data = data.frame(x = median(cvs),y = 0), aes(x, y), 
              label = round(median(cvs), digits = 1), vjust = -0.8, hjust = -0.2,color = "red", size = 3.5) + 
    labs(y = "", x = "%CV") + coord_flip() + theme_classic()

  temp.df <- data.frame(protein = "Groups", CV = cvs, number = 1)
  temp.df <- temp.df %>% mutate(range = case_when(CV < 10 ~ "<10%", 
                                                  (CV > 10) & (CV < 20) ~ "10%~20%", 
                                                  (CV > 20) & (CV < 50) ~ "20%~50%", 
                                                  CV > 50 ~ ">50%"))
  temp.df$range <- ordered(temp.df$range, levels = c(">50%", "20%~50%", "10%~20%", "<10%"))
  g2<-ggplot(temp.df, aes(x = protein, y = number, fill = range)) +
    geom_bar(position = "stack", stat = "identity", width = 0.5) + 
    scale_fill_npg() +
    labs(y = "Numbers of available rows",x = "", fill = "%CV") + theme_classic()
  g1+g2
}

distIndProtein<-function(fixed_data, DoE, group){
  num.proteins<-apply(fixed_data, 2, function(x) sum(!is.na(x)))
  temp.df<-data.frame(number.row=num.proteins, sample=colnames(fixed_data),condition=DoE[match(colnames(fixed_data),DoE[,1]),group])
  ggbarplot(temp.df, x="sample",y= "number.row",
            fill = "condition",               # change fill color by condition
            color = "white",            # Set bar border colors to white
            palette = "jco",            # jco journal color palett. see ?ggpar
            sort.val = "asc",          # Sort the value in dscending order
            sort.by.groups = TRUE,     # Don't sort inside each group
            x.text.angle = 90,           # Rotate vertically x axis texts
            ggtheme = theme_pubclean()
  )+
    font("x.text", size = 8, vjust = 0.5)+
    ggtitle("Number of available rows per sample, colored by condition")
}

upsetplot<-function(fixed_data, DoE, group){
  flag.df <- data.frame(fixed_data)
  flag.df <-data.frame(1*(!is.na(flag.df)))
  protein=rownames(flag.df)
  group.flag.df <- list()
  i<-1
  for (i in levels(as.factor(DoE[,group]))){
    group.flag.df[[i]]<-protein[rowSums(flag.df[,DoE[which(DoE[,group]==i),1]])>1]
  }
  upset(fromList(group.flag.df), order.by="freq",decreasing=T,cutoff=0)
}

datacompleteness<-function(fixed_data, DoE){
  percent.samples<-apply(fixed_data, 1, function(x) sum(!is.na(x))/ncol(fixed_data))
  percent.samples<-percent.samples[order(percent.samples,decreasing=TRUE)]
  temp.df<-data.frame(protein=1:nrow(fixed_data),datacompleteness=percent.samples)
  ggplot(temp.df,aes(protein,datacompleteness))+geom_point()+labs(y="Data Completeness",x="Unique feature")+
  geom_vline(xintercept=max(which(temp.df$datacompleteness>=0.99)), linetype="dashed", color = "red")+
  annotate("text", max(which(temp.df$datacompleteness>=0.99)), 1.1, vjust = -0.5, label = "99%", color="red")+
  geom_vline(xintercept=max(which(temp.df$datacompleteness>=0.9)), linetype="dashed", color = "red")+
  annotate("text", max(which(temp.df$datacompleteness>=0.9)), 1.1, vjust = -0.5, label = "90%", color="red")+
  geom_vline(xintercept=max(which(temp.df$datacompleteness>=0.5)), linetype="dashed", color = "red")+
  annotate("text", max(which(temp.df$datacompleteness>=0.5)), 1.1, vjust = -0.5, label = "50%", color="red")+
  theme_classic()
}

corplot<-function(fixed_data){
  fixed_data<-log2(fixed_data)
  M.pearson<-cor(na.omit(fixed_data), method="pearson")
  M.spearman<-cor(na.omit(fixed_data), method="spearman")
  g1<-ggcorrplot(M.pearson,tl.cex=5,type = "lower",outline.col = "white",hc.order=TRUE)+ggtitle("Pearson correlation for all samples")
  g2<-ggcorrplot(M.spearman,tl.cex=5,type = "lower",outline.col = "white",hc.order=TRUE)+ggtitle("Spearman correlation for all samples")
  g1+g2
}
