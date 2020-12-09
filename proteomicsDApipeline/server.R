library(shiny)
library(shinydashboard)
source("DataInput.R")
source("QCPlots.R")
source("Preprocessing.R")
source("testingDE.R")
source("DimenReduce.R")
source("UnsupervisedLearning.R")

server <- function(input, output) {
  fixed_data<-reactive({
    req(input$raw_data)
    if(grepl(".xlsx",input$raw_data$name)){
      df <- read.xlsx(input$raw_data$datapath)}
    if (grepl(".csv",input$raw_data$name)){
      df <- read.csv(input$raw_data$datapath, header=TRUE)
    }
    inputraw(df)
  })
  DoE<-reactive({
    req(input$meta_data)
    if(grepl(".xlsx",input$meta_data$name)){
      df <- read.xlsx(input$meta_data$datapath)}
    if (grepl(".csv",input$meta_data$name)){
      df <- read.csv(input$meta_data$datapath, header=TRUE)
    }
    setdoe(df)
  })
  output$group<-renderUI({
    DoE<-DoE()
    selectizeInput(inputId = "group",
                   label="Select group column for grouped graphs",
                   choices = colnames(DoE)[-1])
  })
  output$blockfactor<-renderUI({
    DoE<-DoE()
    selectInput(inputId = "blockfactor", 
                label="Blocking on...(supports 2-level factors only)", 
                choices = c(colnames(DoE)[-1]),
                selected=NULL)
  })
  output$FoI<-renderUI({
    DoE<-DoE()
    selectInput(inputId = "FoI", 
                label="test the factor of...", 
                choices = c(colnames(DoE)[apply(DoE, 2, function(x) length(unique(x)))==2]),
                selected = NULL)
  })
  output$colorfactor<-renderUI({
    DoE<-DoE()
    selectInput(inputId = "colorfactor", 
                label="color by...", 
                choices = c(colnames(DoE)[-1]),
                selected = NULL)
  })
  
  output$graph <- renderPlot({
    if (input$options=="%CV plots for QC samples"){
      g<-cvplots(fixed_data(),DoE(),1)}
    if (input$options=="%CV plots for QC samples stacked bar"){
      g<-cvplots(fixed_data(),DoE(),2)}
    if (input$options=="Distribution of identified proteins"){
      g<-distIndProtein(fixed_data(),DoE(),input$group)}
    if (input$options=="Overlap of identified proteins"){
      g<-upsetplot(fixed_data(),DoE(),input$group)}
    if (input$options=="Data Completeness"){
      g<-datacompleteness(fixed_data(),DoE())}
    g
  })
  
  react_fixed_data<-reactive({preprocessing(fixed_data(), DoE(), input$filterlevel, input$normalization,
                                            input$whetherimpute, input$imputation)})
  output$view <- renderTable({
    sub.head<-head(react_fixed_data(), n=5)
    sub.head[,1:5]
  })
  
  listoutput<-reactive({testingDE(react_fixed_data(), DoE(),input$DEtest, input$FoI, input$whetherblock,input$blockfactor)})
  
  output$volcanoplot <- renderPlot({
    listoutput()[["graph"]]
  })
  
  output$dimenreduction <- renderPlot({
    g<-dimen.reduce(react_fixed_data(), DoE(), input$DRmethod, input$colorfactor)
    g
  })
  
  output$cluster <- renderPlot({
    if (input$rows=="all"){includedrows=c(1:nrow(react_fixed_data()))}else{
      includedrows=listoutput()[["DEdf"]]$name}
    if (input$Cmethod=="Hierarchical Clustering"){
      g<-fcluster(react_fixed_data(), input$Cmethod, includedrows, input$whetherlabel, 0)
    }else{
      g<-fcluster(react_fixed_data(), input$Cmethod, includedrows, FALSE, input$clusternum)
    }
    g
  })
}
