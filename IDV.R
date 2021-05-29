library(ggplot2)
library(reshape2)
library(corrplot)
library(dplyr)
library(patchwork)
IDV_plot<-function(data){
  data$name<-rownames(data)
  df<-melt(data,id.vars="name",variable.name="sample",value.name="intensity")
  g<-ggplot(df,aes(y=name,x=log2(intensity)))+geom_boxplot()+geom_point(aes(color=factor(sample)))
  return(g)
}

corrplot_customize<-function(data, corr_sign, p_threshold, order, ncluster, colorscheme){
  if (colorscheme=="red-white-blue"){colorscheme.val<-colorRampPalette(c("blue", "white", "red"))(100)}
  if (colorscheme=="heat"){colorscheme.val<-heat.colors(100)[c(100:1)]}
  if (colorscheme=="cm"){colorscheme.val<-cm.colors(100)}
  t_data<-t(data)
  M<-cor(t_data,use="pairwise.complete.obs")
  for (i in rownames(M)){
    if (sum(!is.na(M[,i]))<2){
      M<-M[-which(rownames(M)==i),-which(rownames(M)==i)]}
  }
  res1 <- cor.mtest(M, conf.level = .95)
  corrplot(M, type = "lower",p.mat = res1$p , insig = "label_sig",
                sig.level = p_threshold, order = "hclust",addrect = ncluster,
                pch.cex = .9, pch.col = "white",addgrid.col = NA, tl.col = "black", tl.srt = 45,tl.cex = 0.5,col=colorscheme.val)
}

combined_lolipop_plot<-function(selected_protein, PTM.data.ls,ind,proteinACC, col_position, modificationType){
  sub.PTM.data<-PTM.data.ls[["data"]][ind,]
  anno.PTM<-PTM.data.ls[["other_annotation"]][ind,c(proteinACC, col_position, modificationType)]
  df<-melt(cbind(sub.PTM.data, anno.PTM))
  g1<-ggplot(df, aes(x=position, y=log10(value))) +
      geom_segment(aes(x=position, xend=position, y=0, yend=log10(value))) +
      geom_point(size=4, alpha=0.6)+facet_wrap(~modificationType,strip.position = "left")+theme_bw()
  
  prot_data <- get_features(selected_protein)
  prot_data <- feature_to_dataframe(prot_data)
  p <- draw_canvas(prot_data)
  p <- draw_chains(p, prot_data)
  p <- draw_domains(p, prot_data)
  p <- draw_repeat(p, prot_data)
  p <- draw_motif(p, prot_data)
  p <- draw_phospho(p, prot_data, size = 8)
  
  # background and y-axis
  p <- p + theme_bw(base_size = 20) + # white backgnd & change text size
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank()) +
    theme(panel.border = element_blank())
  
  # add titles
  rel_subtitle <- paste0("circles = phosphorylation sites\n",
                         "RHD = Rel Homology Domain\nsource:Uniprot")
  
  p <- p + labs(title = paste0("Schematic of protein", selected_proteins[1]),
                subtitle = rel_subtitle)
  g2<-p
  return(g1+g2+plot_layout(nrow = 2))
}

library(dplyr)
library(reshape2)
library(circlize)
library(RColorBrewer)
circosplot.fun<-function(SIoutput.ls,dat.ls,proteinACC,ds.included){
  colorscheme<-brewer.pal(length(ds.included), "Set3")
  df<-data.frame()
  for(i in c("protein data",ds.included)){
    df<-rbind(df, cbind(SIoutput.ls[[i]][["alldf"]],type=i))
  }
  
  df<- df %>% mutate(change=case_when(
    log2FC<(-0.5) ~ "down", 
    log2FC>(-0.5)&(log2FC<0.5) ~ "no change",
    log2FC>0.5 ~ "up"))
  df<- df %>% mutate(sectionID=paste0(type,"_",change))
  
  circos.par("track.height" = 0.1, cell.padding=c(0,0,0,0))
  circos.initialize(factors = df$sectionID, x = df$log2FC)
  circos.track(track.index=1, sectors = df$sectionID, y = df$pvalue, panel.fun = function(x, y) {
    circos.axis(h='top', labels=TRUE, major.tick=TRUE, labels.cex=1, 
                labels.font=15, direction='outside', minor.ticks=10, lwd=1)
  })
  circos.track(track.index=2, factors = df$sectionID, y = df$pvalue, panel.fun = function(x, y) {})
  circos.track(track.index=3, factors = df$sectionID, y = df$pvalue, panel.fun = function(x, y) {}) 
  circos.trackPoints(df$sectionID,  df$log2FC, df$pvalue, pch = 16, cex = 0.5, col='#03040550')
  
  connection_df<-data.frame()
  for(i in ds.included){
    anno<-dat.ls[[i]][["other_annotation"]]
    anno<-data.frame(cbind(anno[,1],anno[,proteinACC[i]]))
    colnames(anno)<-c("name","proteinACC")
    anno<-anno[anno$proteinACC %in% df$name,]
    connection_df<-rbind(connection_df, cbind(anno,type=i,color=colorscheme[which(ds.included==i)]))
  }
  j<-1
  connection_df$index1<-NA
  connection_df$value1<-NA
  connection_df$index2<-NA
  connection_df$value2<-NA
  for(j in 1:nrow(connection_df)){
    connection_df$index1[j]<-df$sectionID[which(df$name==connection_df$name[j])]
    connection_df$index2[j]<-df$sectionID[which(df$name==connection_df$proteinACC[j])]
    connection_df$value1[j]<-df$log2FC[which(df$name==connection_df$name[j])]
    connection_df$value2[j]<-df$log2FC[which(df$name==connection_df$proteinACC[j])]
  }
  
  for (row in 1:nrow(connection_df)){
    circos.link(as.character(connection_df[row, 'index1']), as.double(connection_df[row, 'value1']), 
                as.character(connection_df[row, 'index2']), as.double(connection_df[row, 'value2']), 
                w=1, h2=0.95, col=as.character(connection_df[row, 'color']))
  }
  sectors<-unique(df$sectionID)
  highlight.sector(sectors[grep("down",sectors)], track.index=2, col="#00AFBB", text='Down', cex=0.7, text.col='black', niceFacing=TRUE)
  highlight.sector(sectors[grep("no change",sectors)], track.index=2, col="#94B0B3", text='', cex=0.7, text.col='black', niceFacing=TRUE)
  highlight.sector(sectors[grep("up",sectors)], track.index=2, col="#FC4E07", text='Up', cex=0.7, text.col='black', niceFacing=TRUE)
  
  colorscheme2<-brewer.pal(12, "Set3")
  
  lapply(c("protein data",ds.included),function(x){highlight.sector(sectors[grep(x,sectors)], 
                                                                    track.index=1,col=sample(colorscheme2,1), text=x, cex=0.8, text.col='white', niceFacing=TRUE)})
  cp<-recordPlot()
  return(cp)
}
