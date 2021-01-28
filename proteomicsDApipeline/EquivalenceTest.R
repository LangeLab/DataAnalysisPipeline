# Equivalence test tool
library(TOSTER)
library(ggplot2)
library(ggrepel)

#d=="By row, test each protein respectively"
eqtest.row<-function(data, DoE, FoI, FoI1, FoI2, lowbound, upbound){
i<-1
data<-log10(data)
ind1<-DoE[which(DoE[,FoI]==FoI1),1]
ind2<-DoE[which(DoE[,FoI]==FoI2),1]
eq.pvalues<-matrix(NA, nrow=nrow(data),ncol=3)
for (i in 1:nrow(data)){
  m1<-mean(as.numeric(data[i,ind1]),na.rm=TRUE)
  sd1<-sd(as.numeric(data[i,ind1]),na.rm=TRUE)
  m2<-mean(as.numeric(data[i,ind2]),na.rm=TRUE)
  sd2<-sd(as.numeric(data[i,ind2]),na.rm=TRUE)
  if(any(is.na(c(m1,sd1,m2,sd2)))){next}
  n1<-length(na.omit(as.numeric(data[i,ind1])))
  n2<-length(na.omit(as.numeric(data[i,ind2])))
  res<-TOSTtwo(m1,m2,sd1,sd2,n1,n2,
          low_eqbound_d=lowbound, high_eqbound_d=upbound,
          alpha = 0.05, plot = FALSE, var.equal=FALSE)
  eq.pvalues[i,]<-c(res[["TOST_p1"]],res[["TOST_p2"]],m1-m2)
}
colnames(eq.pvalues)<-c("lowerp","upperp","log2FC")

eq.pvalues<-data.frame(eq.pvalues)
eq.pvalues$equivalence<-(apply(eq.pvalues[,1:2],1,max)<0.05)
eq.pvalues$maxp<-(apply(eq.pvalues[,1:2],1,max))
eq.pvalues$name<-rownames(data)

g<-ggplot(eq.pvalues, aes(x=log2FC, y=-log10(maxp),
                          color=equivalence,
                          alpha=equivalence, 
                          label=name))+
  geom_point()+scale_alpha_manual(values=c(0.2, 1.0 ,0))
  return(g)
}

#d=="By column, test all proteins at the same time"
eqtest.all<-function(data, DoE, FoI, FoI1, FoI2, lowbound, upbound){
  data<-log10(data)
  ind1<-DoE[which(DoE[,FoI]==FoI1),1]
  ind2<-DoE[which(DoE[,FoI]==FoI2),1]
  m1<-mean(as.matrix(data[,ind1]),na.rm=TRUE)
  sd1<-sd(as.matrix(data[,ind1]),na.rm=TRUE)
  m2<-mean(as.matrix(data[,ind2]),na.rm=TRUE)
  sd2<-sd(as.matrix(data[,ind2]),na.rm=TRUE)

  n1<-sum(!is.na(as.matrix(data[,ind1])))
  n2<-sum(!is.na(as.matrix(data[,ind2])))
  res<-TOSTtwo(m1,m2,sd1,sd2,n1,n2,
               low_eqbound_d=lowbound, high_eqbound_d=upbound,
               alpha = 0.05,plot=FALSE, var.equal=FALSE)
  
  ptost<-max(res$TOST_p1,res$TOST_p2)
  TOSToutcome<-ifelse(ptost<res$alpha,"significant","non-significant")
  
  title.c<-paste("Equivalence bounds ",round(res$low_eqbound,digits=3),
                 " and ",round(res$high_eqbound,digits=3),
                 "\nMean difference = ",round(res$diff,digits=3),
                 " \n TOST: ", 100*(1-res$alpha*2),"% CI [",round(res$LL_CI_TOST,digits=3),";",round(res$UL_CI_TOST,digits=3),"] ",
                 TOSToutcome," \n NHST: ", 100*(1-res$alpha),"% CI [",round(res$LL_CI_TTEST,digits=3),";",round(res$UL_CI_TTEST,digits=3),"] ")
  
  g<-ggplot() + 
  scale_x_continuous(name="Mean difference", limits=c(min(res$LL_CI_TOST,res$low_eqbound)-max(res$UL_CI_TOST-res$LL_CI_TOST, res$high_eqbound-res$low_eqbound)/10, max(res$UL_CI_TOST,res$high_eqbound)
                                        +max(res$UL_CI_TOST-res$LL_CI_TOST,res$high_eqbound-res$low_eqbound)/10))+ 
  scale_y_continuous(name="", limits=c(0,1)) +
  geom_vline(mapping=aes(xintercept=res$low_eqbound),linetype="twodash")+ 
  geom_vline(mapping=aes(xintercept=res$high_eqbound),linetype="twodash")+
  geom_vline(mapping=aes(xintercept=0),linetype="twodash",color="grey")+
  geom_point(mapping=aes(x=res$diff, y=0.5),size=3, shape=15)+
  geom_segment(mapping=aes(x=res$LL_CI_TOST,y=0.5,xend=res$UL_CI_TOST,yend=0.5))+
  geom_segment(mapping=aes(x=res$LL_CI_TTEST,y=0.5,xend=res$UL_CI_TTEST,yend=0.5))+
  ggtitle(title.c)+
  theme_classic()+ 
  theme(plot.title = element_text(face="bold",size=12))
  
  return(g)
}
