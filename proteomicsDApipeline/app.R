library(shiny)
library(shinydashboard)
source("DataInput.R")
source("QCPlots.R")
source("Preprocessing.R")
source("testingDE.R")
source("DimenReduce.R")
source("UnsupervisedLearning.R")

# Define UI ----
ui <- dashboardPage(
    dashboardHeader(title="Proteomics Data Analysis Pipeline"),
    dashboardSidebar(
        h4("File input"),
        fileInput("raw_data", "please upload .csv/.xlsx file of raw data with quantitative intensities of proteins:",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv", ".xlsx")),
        fileInput("meta_data", "please upload .csv/.xlsx file of meta data with only sample names matched with raw data, and corresponding covariates of interest:",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv", ".xlsx")),
        h4("1. Graphs at QC stage"),
        #QC Plots
        selectizeInput("options", 
        label="Quality Control graph options", 
        choices = c( "%CV plots for QC samples",
                     "%CV plots for QC samples stacked bar",
                      "Distribution of identified proteins",
                      "Overlap of identified proteins",
                      "Data Completeness")),
        conditionalPanel('(input.options == "Distribution of identified proteins")|(input.options == "Overlap of identified proteins")',
                         uiOutput("group")),
        # Normalization, imputation, etc.
        h4("2. Processing"),
        numericInput("filterlevel", 
                     h5("Row filter level (%)"), 
                     value = 0),
        checkboxInput(inputId = "normalization", 
                      "Median Normalization", value = FALSE),
        checkboxInput(inputId = "whetherimpute", 
                      "Imputation", value = FALSE),
        conditionalPanel(condition='input.whetherimpute==1',
        selectInput(inputId = "imputation", 
                    label="Imputation Method", 
                    choices = c( "Down-shifted Normal samples",
                                 "MinProb",
                                 "knn",
                                 "min"),
                    selected = NULL)),
        # tests
        h4("3. Inference on Differential Expression"),
        checkboxInput(inputId = "whetherblock", 
                      "Blocking", value = FALSE),
        conditionalPanel('input.whetherblock==1',
                         uiOutput("blockfactor")),
        selectInput(inputId = "DEtest", 
                    label="Statistical testing for DE", 
                    choices = c("limma",
                                 "t-test"),
                    selected = NULL),
        uiOutput("FoI"),
        
        # Dimensionality Reduction
        h4("4. Dimensionality reduction"),
        selectInput(inputId = "DRmethod", 
                    label="use the method of...",
                    choices=c("PCA","t-SNE", "UMAP")),
        uiOutput("colorfactor"),
        
        # clustering
        h4("5. Clustering"),
        selectInput(inputId = "Cmethod", 
                    label="use the method of...",
                    choices=c("Hierarchical Clustering",
                              "Fuzzy Clustering")),
        selectInput(inputId = "rows", 
                    label="use the rows of...",
                    choices=c("significant","all")),
        conditionalPanel('input.Cmethod == "Hierarchical Clustering"',
            checkboxInput("whetherlabel", 
                         h5("include row labels"), 
                         value = FALSE)),
        conditionalPanel('input.Cmethod != "Hierarchical Clustering"',
            numericInput("clusternum", 
                         h5("the number of centers: "), 
                         value = 2))
        ),
        dashboardBody(
            fluidRow(tabBox(title = "Input & Pre-processing",
                       tabPanel("Graphs at QC stage",
                                plotOutput("graph")),
                       tabPanel("Head", 
                                h5("The first 5 rows and 5 columns of data currently"),
                                tableOutput("view")))
            ),
            fluidRow(tabBox(title="DE test and others",
                tabPanel("Volcano plot",
                           h5("Volcano plot of DE testing results"),
                           plotOutput("volcanoplot")),
                tabPanel("Dimensionality Reduction",
                         h5("Output of dimensionality reduction methods"),
                         plotOutput("dimenreduction")),
                tabPanel("Clustering",
                         h5("Output of clustering"),
                         plotOutput("cluster"))
            )
            )
        )
)
# Define server logic ----
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

# Run the app ----
shinyApp(ui = ui, server = server)