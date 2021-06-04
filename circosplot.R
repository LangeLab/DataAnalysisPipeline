library(tidyverse)
library(reshape2)
library(circlize)
ball_LFC_data <- read_csv('Connection_data/ball_LFC_data.csv', col_types = cols())
connections_data <- read_csv('Connection_data/ball_connections_data_for_full.csv', col_types = cols())
circos.par("track.height" = 0.1, cell.padding=c(0,0,0,0))
circos.initialize(factors = ball_LFC_data$SectionIdx, x = ball_LFC_data$log2FC)

circos.track(track.index=1, factors = ball_LFC_data$SectionIdx, y = ball_LFC_data$y, panel.fun = function(x, y) {
  circos.axis(h='top', labels=TRUE, major.tick=TRUE, labels.cex=1, 
              labels.font=15, direction='outside', minor.ticks=8, lwd=1)
})
circos.track(track.index=2, factors = ball_LFC_data$SectionIdx, y = ball_LFC_data$y, panel.fun = function(x, y) {})

circos.track(track.index=3, factors = ball_LFC_data$SectionIdx, y = ball_LFC_data$y, panel.fun = function(x, y) {}) 

circos.trackPoints(ball_LFC_data$SectionIdx,  ball_LFC_data$log2FC, ball_LFC_data$y, pch = 16, cex = 0.5, col='#03040550')

for (row in 1:nrow(connections_data)){
  
  circos.link(as.character(connections_data[row, 'index1']), as.double(connections_data[row, 'value1']), 
              as.character(connections_data[row, 'index2']), as.double(connections_data[row, 'value2']), 
              w=1, h2=0.95, col=as.character(connections_data[row, 'color']))
}

highlight.sector(c('Sec 01', 'Sec 07', 'Sec 13'), track.index=2, col="#00AFBB", text='Down', cex=0.7, text.col='black', niceFacing=TRUE)
highlight.sector(c('Sec 02', 'Sec 08', 'Sec 14'), track.index=2, col="#94B0B3", text='', cex=0.7, text.col='black', niceFacing=TRUE)
highlight.sector(c('Sec 03', 'Sec 09', 'Sec 15'), track.index=2, col="#FC4E07", text='Up', cex=0.7, text.col='black', niceFacing=TRUE)

highlight.sector(c('Sec 01', 'Sec 02', 'Sec 03'), track.index=1, col='#2274A5', text='Total protein abundance', cex=0.8, text.col='white', niceFacing=TRUE)
highlight.sector(c('Sec 07', 'Sec 08', 'Sec 09'), track.index=1, col='#32936F', text='N Termini', cex=0.8, text.col='white', niceFacing=TRUE)
highlight.sector(c('Sec 13', 'Sec 14', 'Sec 15'), track.index=1, col='#FC6471', text='Phosphopeptides', cex=0.8, text.col='white', niceFacing=TRUE)
