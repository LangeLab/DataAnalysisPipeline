library(shinydashboard)
library(shiny)
library(shinycssloaders)
library(stringr)
library(DT)
source("DataInput.R")
source("QCplots.R")
source("Preprocessing.R")
source("testingDE.R")
source("DimenReduce.R")
source("UnsupervisedLearning.R")
source("annotationtool.R")
source("EquivalenceTest.R")
source("PeptideCalculate.R")
source("batcheffect.R")
source("IDV.R")
options(shiny.maxRequestSize=20*1024^2)
# Define UI ----
ui <- dashboardPage(
    dashboardHeader(title="Proteomics Data Analysis Pipeline", titleWidth = 350),
    dashboardSidebar(width = 250,sidebarMenu(
        menuItem("Input and annotation", tabName = "IA"),
        menuItem("QC" , tabname = "QC", icon = icon("bacon"),
            menuSubItem("Standard based QC", tabName = "StandardQC"),
            menuSubItem("Sample based QC", tabName = "SampleQC")),
        menuItem("Pre-processing", tabName = "PP"),
        menuItem("Statistical Inference", tabName = "StatInf"),
        menuItem("Dimensionality reduction", tabName = "DR"),
        menuItem("Clustering", tabName = "Clustering"),
        menuItem("IV",tabname="IV",icon=icon("bars"),
        menuSubItem("Individual Protein Visualization 1", tabName = "IDV1"),
        menuSubItem("Individual Protein Visualization 2", tabName = "IDV2"),
        menuSubItem("Individual Protein Visualization 3", tabName = "IDV3"),
        menuSubItem("Individual Protein Visualization 4", tabName = "IDV4")),
        menuItem("Generate a report", tabName="GenerateReport"))),
        dashboardBody(
            tabItems(
                # Input and annotation ----
                tabItem(tabName="IA",
                        h4("File input"),
                        fluidRow(
                        box(h4("Meta file input"),
                                fileInput("meta_data", "please upload .csv/.xlsx file of meta data with only sample names matched with raw data, and corresponding covariates of interest:",
                                          multiple = FALSE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv", ".xlsx")),
                                checkboxInput("whether_replica", "Is there a column indicating replica in the meta data file?", value=FALSE),
                                conditionalPanel('input.whether_replica==1', 
                                                 uiOutput("replicacol")),
                            
                        checkboxInput("whether_protein_file", "Upload Protein file", value=FALSE),
                        conditionalPanel('input.whether_protein_file==1',
                                         fileInput("protein_data", "please upload .csv/.xlsx file of raw data with quantitative intensities of proteins:",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx")),
                                         numericInput("protein_col_id","Which column contains unique identifiers of this dataset?",1)),
                        checkboxInput("whether_termini_file", "Upload Termini file", value=FALSE),
                        conditionalPanel('input.whether_termini_file==1',
                                         fileInput("termini_data", "please upload .csv/.xlsx file of raw termini data with quantitative intensities of terminis:",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx")),
                                         numericInput("termini_col_id","Which column contains unique identifiers of this dataset?",1)),
                        checkboxInput("whether_peptide_file","Upload Peptide file", value=FALSE),
                        conditionalPanel('input.whether_peptide_file==1',
                                         fileInput("peptide_data", "please upload .csv/.xlsx file of raw peptide data with quantitative intensities of peptides:",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx")),
                                         numericInput("peptide_col_id","Which column contains unique identifiers of this dataset?",1)),
                        checkboxInput("whether_PTM_file","Upload PTM file", value=FALSE),
                        conditionalPanel('input.whether_PTM_file==1',
                                         fileInput("PTM_data", "please upload .csv/.xlsx file of raw PTM file",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx")),
                                         numericInput("PTM_col_id","Which column contains unique identifiers of this dataset?",1))),
                        
                       
                        box(h4("Annotation tools"),
                            conditionalPanel('input.whether_termini_file==1 || input.whether_PTM_file==1',
                            uiOutput("select_dataset_for_anno"),
                            uiOutput("anno_idf"),
                            uiOutput("anno_seq"),
                            numericInput("lseqwindow", 
                                         h5("Length of sequence window:"), 
                                         value = 1),
                            actionButton("annoButton", "Get annotation"),
                            div(dataTableOutput("annotationpreview")%>% withSpinner(color="#0dc5c1"), style = "font-size:60%"),
                            downloadButton("downloadAnnotation", "Download full annotation"))),
                        
                        box(h4("Protein intensity calculator based on peptide"),
                            conditionalPanel('input.whether_peptide_file==1',
                                             numericInput("peptide_col_ACC",
                                                          "Please input the column of protein accession:",
                                                          value=1),
                            selectInput("sumform", "Choose method of sum:",choices=c("sum","weighted","top3")),
                            actionButton("pepcalButton", "Calculate"),
                            div(dataTableOutput("peptidepreview")%>% withSpinner(color="#0dc5c1"), style = "font-size:60%"),
                            checkboxInput("whether_replace_protein", "Use the dataset calculated from peptides to replace the protein data?", 
                            value=FALSE)))
                        )),
                # Standard QC ----
                tabItem(tabName = "StandardQC",
                        box(fileInput("standardQC_data", "please upload .csv/.xlsx file of raw data with quantitative intensities of standard samples:",
                            multiple = FALSE,
                            accept = c("text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv", ".xlsx")),
                            fileInput("standardQC_meta", "please upload .csv/.xlsx file of meta data for standard samples:",
                                      multiple = FALSE,
                                      accept = c("text/csv",
                                                 "text/comma-separated-values,text/plain",
                                                 ".csv", ".xlsx"))),
                        fluidRow(
                        box(uiOutput("batch_col"),
                            uiOutput("order_col"),
                            uiOutput("factor_cols"),
                            selectizeInput("correctmethod",label="Batch effect correction method",
                                           choices=c("MedianCentering","MeanCentering")),
                            actionButton("standar_batcheffectdButton", "Confirm"),
                            downloadButton("DownloadCorrectedQCData","Download corrected QC data")),
                        box(plotOutput("standard_batcheffectplots")%>% withSpinner(color="#0dc5c1")),
                        box(plotOutput("standard_batcheffect_clustering")%>% withSpinner(color="#0dc5c1")),
                        box(plotOutput("standard_batcheffect_heatmap")%>% withSpinner(color="#0dc5c1"))
                        )),
                # Sample QC ----
                tabItem(tabName = "SampleQC",
                        fluidRow(
                        box(h4("Show the QC graphs of"),
                            uiOutput("selectdata_QC")),
                        tabBox(title = "QC plots",
                               tabPanel("%CV plots for QC samples",
                                        plotOutput("graph1")),
                               tabPanel("Distribution of identified proteins",
                                        uiOutput("group"),
                                        plotOutput("graph3")),
                               tabPanel("Overlap of identified proteins",
                                        uiOutput("group2"),
                                        plotOutput("graph4")),
                               tabPanel("Data Completeness",
                                        plotOutput("graph5")),
                               tabPanel("Correlation Plot",
                                         plotOutput("graph6"))),
                        box(checkboxInput("whethercorrectbatch","Correct the batch effect with the same setting in Standard QC?"),
                            checkboxInput("whether_average_replica",
                                          "Average the intensities across replica?",
                                          value=FALSE))
                            #conditionalPanel("input.whethercorrectbatch==1",
                             #                plotOutput("sample_batch_plot"))
                        )),
                # Pre-processing ----
                tabItem(tabName = "PP",
                        h4("2. Pre-processing"),
                        # Normalization, imputation, etc.
                        fluidRow(
                        box(uiOutput("preprocessing_panel"))
                        )),
                        #plotOutput("viewviolin")%>% withSpinner(color="#0dc5c1"),
                        #downloadButton("DownloadProcessedData", "Download")),
                # Statistical Inference ----
                tabItem(tabName="StatInf",
                        # Equivalence test
                        #h4("3.1 Inference on Equivalence"),
                        uiOutput("selectdata_SI"),
                        # box(
                        #     selectInput(inputId = "dd", 
                        #             label="Specify the way of testing...", 
                        #             choices = c("By column, test all proteins at the same time",
                        #                         "By row, test each protein respectively")),
                        #     numericInput("lowerbound", "Lower bound:", -0.5),
                        #     numericInput("upperbound", "Upper bound:", 0.5),
                        #     uiOutput("eqFoI"),
                        #     uiOutput("eqFoIlevels"),
                        #     conditionalPanel('input.dd=="By row, test each protein respectively"',
                        #     plotOutput("eqtestres1")%>% withSpinner(color="#0dc5c1")),
                        #     conditionalPanel('input.dd=="By column, test all proteins at the same time"',
                        #      plotOutput("eqtestres2"))
                        # ),
                        # DE tests
                        #h4("3.2 Inference on Differential Expression"),
                        fluidRow(
                        box(width=12,
                            checkboxInput(inputId = "whetherstandardDE", 
                                          "Standardization?", 
                                          value=FALSE),
                            #conditionalPanel('input.whetherstandardDE==1',
                            #                 selectInput(inputId = "standardmethodDE", 
                            #                             label="Standardization by", 
                            #                             choices = c("z-score",
                            #                                         "log2FC over median"))),
                            checkboxInput(inputId = "whetherblock", 
                                      "Blocking", value = FALSE),
                            conditionalPanel('input.whetherblock==1',
                                         uiOutput("blockfactor")),
                            checkboxInput(inputId = "whetherweighting", 
                                      "Weighting?", 
                                      value=FALSE),
                            conditionalPanel('input.whetherweighting==1',
                                         numericInput("NAweights", "Weight",value=10e-5)),
                            selectInput(inputId = "DEtest", 
                                    label="Statistical testing for DE", 
                                    choices = c("limma",
                                                "t-test"),
                                    selected = NULL, width="33%"),
                            uiOutput("FoI"),
                            uiOutput("FoIlevels"),
                            h5("Volcano plot of DE testing results"),
                            plotOutput("volcanoplot")%>% withSpinner(color="#0dc5c1"),
                            downloadButton("DownloadDEresults","Download Differential Expression test results"))
                        )),
                # Dimensionality Reduction ----
                tabItem(tabName = "DR",
                        h4("4. Dimensionality reduction"),
                        uiOutput("selectdata_DR"),
                        fluidRow(
                        box(
                        selectInput(inputId = "DRrows", 
                                    label="use the rows of...",
                                    choices=c("all","significant")),
                        uiOutput("colorfactor")),
                        box(
                        selectInput(inputId = "DRmethod", 
                                    label="use the method of...",
                                    choices=c("PCA","t-SNE", "UMAP")),
                        conditionalPanel('input.DRmethod=="t-SNE"',
                                         numericInput("tSNEper","Specify the perplexity for t-SNE", 0)))),
                        box(width=12,plotOutput("dimenreduction"))),
                # clustering ----
                tabItem(tabName="Clustering",
                        h4("5. Clustering"),
                        uiOutput("selectdata_cluster"),
                        fluidRow(
                        box(
                        selectInput(inputId = "Cmethod", 
                                    label="use the method of...",
                                    choices=c("Hierarchical Clustering",
                                              "k-means",
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
                                                      value = 4)))),
                        box(width=12,plotOutput("cluster"))),
                # Individual Protein Visualization ----
                tabItem(tabName="IDV1",
                        uiOutput("selectdata_IDV"),
                        selectInput("proteingroup","Select the group of...",
                                    choices=c("significant in differential expression",
                                              "specify below...")),
                        uiOutput("proteingroupID"),
                        fluidRow(
                        box(width=12,
                            column(width=6,
                                DT::dataTableOutput("IDV_table"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
                            column(width=6,
                                plotOutput("IDV_boxplot"))))
                        ),
                tabItem(tabName="IDV2",
                        selectInput("corr_sign",label="Please select the correlation sign you want to include",
                                                    choices=c("+","-","both")),
                        numericInput("p_threshold", "Please specify p-value threshold",value=0.05),
                        selectInput("corrorder","Please specify the order of rows in correlation plot",
                                                    choices=c("hclust","none")),
                        numericInput("ncluster","Please specify the number of clusters",value=3),
                        selectInput("colorscheme","Please specify the color scheme",
                                                    choices=c("red-white-blue","heat","cm")),
                        fluidRow(
                        box(width=12,
                            column(width=6,DT::dataTableOutput("IDV_corrtable"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
                            column(width=6,plotOutput("IDV_corrplot"))))
                ),
            tabItem(tabName="IDV3",
                    uiOutput("proteinACC"),
                    uiOutput("col_position"),
                    uiOutput("modificationType"),
                    fluidRow(
                    uiOutput("IDVseperate_panel"))
            ),
            tabItem(tabName="IDV4",
                    uiOutput("selectdata_circos"),
                    textInput("groupproACC", "Specify column names of protein accession (separate by comma)", value = ""),
                    fluidRow(box(plotOutput("circosplot")))
            ),
                # generate report ----
                tabItem(tabName="GenerateReport",
                        textInput("title","Customize a title of report...",value="Analysis"),
                        textInput("description", "Customize a description of the report...", value=""),
                        downloadButton("report", "Generate report"))
            )
    )
)
# Define server logic ----
server <- function(input, output) {
    report_vals <- reactiveValues()
    # Data input ----
    protein_data<-reactive({
        req(input$protein_data)
        if(grepl(".xlsx",input$protein_data$name)){
            df <- read.xlsx(input$protein_data$datapath)}
        if (grepl(".csv",input$protein_data$name)){
            df <- data.frame(read.csv(input$protein_data$datapath, header=TRUE, check.names=FALSE))
        }
        df
    })
    
    termini_data<-reactive({
        req(input$termini_data)
        if(grepl(".xlsx",input$termini_data$name)){
            df <- read.xlsx(input$termini_data$datapath)}
        if (grepl(".csv",input$termini_data$name)){
            df <- read.csv(input$termini_data$datapath, header=TRUE, check.names=FALSE)
        }
        data.frame(df)
    })
    
    PTM_data<-reactive({
        req(input$PTM_data)
        if(grepl(".xlsx",input$termini_data$name)){
            df <- read.xlsx(input$termini_data$datapath)}
        if (grepl(".csv",input$termini_data$name)){
            df <- read.csv(input$termini_data$datapath, header=TRUE, check.names=FALSE)
        }
        data.frame(df)
    })
    
    peptide_data<-reactive({
        req(input$peptide_data)
        if(grepl(".xlsx",input$peptide_data$name)){
            df <- read.xlsx(input$peptide_data$datapath)}
        if (grepl(".csv",input$peptide_data$name)){
            df <- read.csv(input$peptide_data$datapath, header=TRUE, check.names=FALSE)
        }
        data.frame(df)
    })
    
    # data integration ----
    DoE0<-reactive({
        req(input$meta_data)
        if(grepl(".xlsx",input$meta_data$name)){
            df <- read.xlsx(input$meta_data$datapath)}
        if (grepl(".csv",input$meta_data$name)){
            df <- read.csv(input$meta_data$datapath, header=TRUE, check.names=FALSE)
        }
        dt<-setdoe(data.frame(df))
        dt
    })
    output$replicacol<-renderUI({
        selectizeInput(inputId = "replicacol",
                       label="Select the column indicates replica",
                       choices = colnames(DoE0())[-1])
    })
    replica_list<-reactive({
        validate(need(input$replicacol, ""))
        df<-cbind(DoE0()[,1],DoE0()[,input$replicacol])
    })
    DoE<-reactive({
        dt<-DoE0()
        if(input$whether_average_replica==TRUE){
            dt[,1]<-dt[,input$replicacol]
            dt[,input$replicacol]<-NULL
            dt<-unique(dt)
        }
        dt
    })
   
    data_collection0<-reactive({
        ls<-list()
        if(input$whether_protein_file==TRUE){ls[["protein data"]]<-inputraw(protein_data(),DoE0(),input$protein_col_id)}
        if(input$whether_termini_file==TRUE){ls[["termini data"]]<-inputraw(termini_data(),DoE0(),input$termini_col_id)}
        if(input$whether_peptide_file==TRUE){ls[["peptide data"]]<-inputraw(peptide_data(),DoE0(),input$peptide_col_id)}
        if(input$whether_PTM_file==TRUE){ls[["PTM data"]]<-inputraw(PTM_data(),DoE0(),input$PTM_col_id)}
        if(input$whether_replace_protein==TRUE){ls[["protein data"]]$data<-peptide_based_fixed_data()}
        ls
    })
    
    data_collection<-reactive({
        ls<-data_collection0()
        if(input$whethercorrectbatch==TRUE){
            ls[["protein data"]]$data<-QCeffectPlot(data_collection0()[["protein data"]]$data,DoE0(),input$batch_col,
                                                         input$order_col, input$factor_cols, input$correctmethod)[["batchcorrectedmatrix"]]
            ls[["termini data"]]$data<-QCeffectPlot(data_collection0()[["termini data"]]$data,DoE0(),input$batch_col,
                                                         input$order_col, input$factor_cols, input$correctmethod)[["batchcorrectedmatrix"]]
            ls[["peptide data"]]$data<-QCeffectPlot(data_collection0()[["peptide data"]]$data,DoE0(),input$batch_col,
                                                         input$order_col, input$factor_cols, input$correctmethod)[["batchcorrectedmatrix"]]
            ls[["PTM data"]]$data<-QCeffectPlot(data_collection0()[["PTM data"]]$data,DoE0(),input$batch_col,
                                                     input$order_col, input$factor_cols, input$correctmethod)[["batchcorrectedmatrix"]]
        }
        
        if(input$whether_average_replica==TRUE){
            ls[["protein data"]]$data<-averagereplica(ls[["protein data"]]$data, replica_list())
            ls[["termini data"]]$data<-averagereplica(ls[["termini data"]]$data, replica_list())
            ls[["peptide data"]]$data<-averagereplica(ls[["peptide data"]]$data, replica_list())
            ls[["PTM data"]]$data<-averagereplica(ls[["PTM data"]]$data, replica_list())
        }
        ls
    })
    
    
    # annotation tool ----
    output$select_dataset_for_anno<-renderUI({
        cc<-c("termini data", "PTM data")[c(input$whether_termini_file, input$whether_PTM_file)]
        selectInput("select_dataset_for_anno",
                    label="Choose the dataset for annotation",
                    choices=cc,
                    selected=NULL)
    })
    
    dataforanno<-reactive({
    validate(need(!is.null(input$select_dataset_for_anno),""))
    if (input$select_dataset_for_anno=="termini data"){df<-termini_data()} 
    if (input$select_dataset_for_anno=="PTM data"){df<-PTM_data()}
    df
    })
    
    output$anno_idf<-renderUI({
        validate(need(!is.null(dataforanno()),""))
        coln<-colnames(dataforanno())
        selectizeInput(inputId = "anno_idf",
                       label="Select the column of identifier for data",
                       choices = coln)
    })
    output$anno_seq<-renderUI({
        validate(need(!is.null(dataforanno()),""))
        coln<-colnames(dataforanno())
        selectizeInput(inputId = "anno_seq",
                       label="Select the column of stripped sequence for data",
                       choices = coln)
    })

    annotation_res<-eventReactive(input$annoButton,{
        annotationtool(dataforanno()[,input$anno_idf], dataforanno()[,input$anno_seq], input$lseqwindow)}
    )
    
    output$downloadAnnotation <- downloadHandler(
        filename = function() {
            paste("AnnotationResults", ".csv", sep = "")
        },
        content = function(file) {
            write.csv(annotation_res(), file, row.names = FALSE)
        }
    )
    
    output$annotationpreview<-renderDataTable({
        annotation_res()[c(1:20),c(1:4,9:10)]
    })
    
    # peptide calculator ----
    peptide_based_fixed_data<-eventReactive(
        input$pepcalButton,{
        calIntPep(peptide_data(),input$sumform, DoE0(),input$peptide_col_ACC)
    })
    
    output$peptidepreview<-renderDataTable({
        peptide_based_fixed_data()[c(1:20),c(1:6)]
    })
    

    # QC plots ----
    standardQC_meta<-reactive({
        req(input$standardQC_meta)
        if(grepl(".xlsx",input$standardQC_meta$name)){
            df <- read.xlsx(input$standardQC_meta$datapath)}
        if (grepl(".csv",input$standardQC_meta$name)){
            df <- read.csv(input$standardQC_meta$datapath, header=TRUE, check.names=FALSE)
        }
        setdoe(data.frame(df))
    })
    
    standardQC_data<-reactive({
        req(input$standardQC_data)
        if(grepl(".xlsx",input$standardQC_data$name)){
            df <- read.xlsx(input$standardQC_data$datapath)}
        if (grepl(".csv",input$standardQC_data$name)){
            df <- read.csv(input$standardQC_data$datapath, header=TRUE, check.names=FALSE)
        }
        inputraw(data.frame(df),standardQC_meta(),1)
    })
    
    output$batch_col<-renderUI({
        validate(need(standardQC_data(),""),need(standardQC_meta(),""))
        selectInput(inputId = "batch_col", 
                    label="select the column indicates batch of samples:",
                    choices=colnames(standardQC_meta())[-1])
    })
    output$order_col<-renderUI({
        validate(need(standardQC_data(),""),need(standardQC_meta(),""))
        selectInput(inputId = "order_col", 
                    label="select the column indicates running order of samples:",
                    choices=colnames(standardQC_meta())[-1])
    })
    output$factor_cols<-renderUI({
        validate(need(standardQC_data(),""),need(standardQC_meta(),""))
        selectInput(inputId = "factor_cols", 
                    label="select the factors you would like to include in annotation", 
                    choices = colnames(standardQC_meta())[-1],
                    selected = NULL, multiple=TRUE)
    })
    
    standardQCoutput<-eventReactive(
        input$standar_batcheffectdButton,{
        QCeffectPlot(standardQC_data(),standardQC_meta(),input$batch_col,
                     input$order_col, input$factor_cols, input$correctmethod)
        })
    
    output$standard_batcheffectplots <- renderPlot({
        grid.arrange(grobs=standardQCoutput()[["plot"]]) 
    })
    output$standard_batcheffect_clustering <- renderPlot({
        standardQCoutput()[["clustering"]]  
    })
    output$standard_batcheffect_heatmap <- renderPlot({
        standardQCoutput()[["heatmap"]]  
    })
    output$DownloadCorrectedQCData <- downloadHandler(
        filename = function() {
            paste("CorrectedQCData", ".csv", sep = "")
        },
        content = function(file) {
            write.csv(standardQCoutput()[["batchcorrectedmatrix"]], file, row.names = TRUE)
        }
    )
    output$selectdata_QC<-renderUI({
        dataset<-c("protein data", "termini data", "peptide data", "PTM data")[
            c(input$whether_protein_file,input$whether_termini_file,input$whether_peptide_file, input$whether_PTM_file)]
        selectInput(inputId = "selectdata_QC", 
                    label="",
                    choices=dataset)
    })
    
    output$group<-renderUI({
        DoE<-DoE()
        selectizeInput(inputId = "group",
                       label="Select group column for grouped graphs",
                       choices = colnames(DoE)[-1])
    })
    output$group2<-renderUI({
        DoE<-DoE()
        selectizeInput(inputId = "group2",
                       label="Select group column for grouped graphs",
                       choices = colnames(DoE)[-1])
    })
    
    output$graph1 <- renderPlot({
        validate(need(input$selectdata_QC,""))
        report_vals$cvviolin<-cvplots(data_collection()[[input$selectdata_QC]]$data,DoE())
        report_vals$cvviolin})
    output$graph3 <- renderPlot({
        validate(need(input$selectdata_QC,""))
        report_vals$distprotein<-distIndProtein(data_collection()[[input$selectdata_QC]]$data,DoE(),input$group)
        report_vals$distprotein})
    output$graph4 <- renderPlot({
        validate(need(input$selectdata_QC,""))
        report_vals$upset<-upsetplot(data_collection()[[input$selectdata_QC]]$data,DoE(),input$group2)
        report_vals$upset})
    output$graph5 <- renderPlot({
        validate(need(input$selectdata_QC,""))
        report_vals$datacompleteness<-datacompleteness(data_collection()[[input$selectdata_QC]]$data,DoE())
        report_vals$datacompleteness})
    output$graph6 <- renderPlot({
        validate(need(input$selectdata_QC,""))
        report_vals$corplot<-corplot(data_collection()[[input$selectdata_QC]]$data)
        report_vals$corplot})
    # pre process ----
    react_data_collection<-reactiveValues()
    react_na_index<-reactiveValues()
    output$preprocessing_panel<-renderUI({
        available_sets<-names(data_collection())
        cc<-c("No normalization",
              "median randomization",
              "normalization over same protein and samples under same condition",
              "normalization over samples under same condition")
        lapply(available_sets,function(x){
            output[[paste0("filter_level_for_",x)]]<-
                renderUI(
                    numericInput(inputId = paste0("filter_level_for_",x),
                                label=paste0("select the row filter level (%) for ", x),
                                value=0)
                        )
        })
        lapply(available_sets,function(x){
            output[[paste0("normalization_method_for_",x)]]<-
                renderUI(
                 selectInput(inputId = paste0("normalization_method_for_",x),
                 label=paste0("select the normalization method for ", x),
                 choices=cc[c(TRUE,TRUE,(x %in% c("termini data","PTM data")), TRUE)])
                 )
        })
        lapply(available_sets, function(x){
            output[[paste0("norm_condition_for_",x)]]<-renderUI({
                validate(need(input[[paste0("normalization_method_for_",x)]],""))
                if (input[[paste0("normalization_method_for_",x)]] %in% c("normalization over same protein and samples under same condition",
                                                                  "normalization over samples under same condition")){
                    cc<-colnames(DoE())[-1]
                }else{cc<-c(NULL)}
                selectInput(inputId = paste0("norm_condition_for_",x),
                            label="select the condition",
                             choices=cc)
                })

        })
        lapply(available_sets,function(x){
            output[[paste0("protein_anno_for_",x)]]<-
                renderUI({
                    validate(need(input[[paste0("normalization_method_for_",x)]]=="normalization over same protein and samples under same condition",""))
                    selectInput(inputId = paste0("protein_anno_for_",x),
                                 label=paste0("select the column cotains protein groups for ", x),
                                 choices=colnames(data_collection()[[x]]$other_annotation))
                })
        })
        lapply(available_sets,function(x){
            output[[paste0("imputation_method_for_",x)]]<-
                renderUI(
                    selectInput(inputId = paste0("imputation_method_for_",x), 
                                label="Imputation Method", 
                                choices = c( "No imputation",
                                             "Down-shifted Normal samples",
                                             "MinProb",
                                             "knn",
                                             "min"),
                                selected = NULL)
                    )
        })
        lapply(available_sets, function(x){
            output[[paste0("impute_condition_for_",x)]]<-renderUI({
                validate(need(input[[paste0("imputation_method_for_",x)]],""))
                if (input[[paste0("imputation_method_for_",x)]]== "Down-shifted Normal samples"){
                    cc<-colnames(DoE())[-1]
                }else{cc<-c(NULL)}
                selectInput(inputId = paste0("impute_condition_for_",x),
                            label="select the condition",
                            choices=cc)
            })
            
        })
        
        lapply(available_sets,function(x){
            output[[paste0("viewviolin_for_",x)]]<-
                renderPlot({
                    validate(need(data_collection(),""))
                    ctitle<-paste0(c("filtered for ",input[[paste0("filter_level_for_",x)]],
                                             "% completeness, normalized by ", input[[paste0("normalization_method_for_",x)]],"\n", 
                                             "imputation by ", input[[paste0("imputation_method_for_",x)]]), collapse="")
                    aa<-preprocessing(x,data_collection(),DoE(),
                                      input[[paste0("filter_level_for_",x)]], 
                                      input[[paste0("normalization_method_for_",x)]],
                                      input[[paste0("imputation_method_for_",x)]],
                                      input[[paste0("norm_condition_for_",x)]],
                                      input[[paste0("impute_condition_for_",x)]],
                                      input[[paste0("protein_anno_for_",x)]])
                    react_data_collection[[x]]<-aa[["data"]]
                    react_na_index[[x]]<-aa[["na.index"]]
                    plotviolin(aa[["data"]],ctitle)
                })
        })
        lapply(available_sets, function(x){
            output[[paste0("DownloadProcessedData_for",x)]]<- downloadHandler(
                filename = function() {
                    paste("ProcessedData_", x, ".csv", sep = "")
                },
                content = function(file) {
                    write.csv(react_data_collection[[x]], file, row.names = TRUE)
                }
            )
        })
        
        myTabs = lapply(available_sets, function(x){
            tabPanel(title=x, 
                     uiOutput(paste0("filter_level_for_",x)),
                     uiOutput(paste0("normalization_method_for_",x)),
                     uiOutput(paste0("norm_condition_for_",x)),
                     uiOutput(paste0("protein_anno_for_",x)),
                     uiOutput(paste0("imputation_method_for_",x)),
                     uiOutput(paste0("impute_condition_for_",x)),
                     plotOutput(paste0("viewviolin_for_",x)),
                     downloadButton(paste0("DownloadProcessedData_for",x),"Download processed data")
            )
        })
        do.call(tabsetPanel, myTabs)
    })

    # statistical inference ----
    output$selectdata_SI<-renderUI({
        dataset<-c("protein data", "termini data", "peptide data", "PTM data")[
            c(input$whether_protein_file,input$whether_termini_file,input$whether_peptide_file, input$whether_PTM_file)]
        selectInput(inputId = "selectdata_SI", 
                    label="",
                    choices=dataset)
    })
    # equivalence test FoI and FoI levels
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
                    choices = colnames(DoE)[-1],
                    selected = NULL)
    })
    output$colorfactor<-renderUI({
        DoE<-DoE()
        selectInput(inputId = "colorfactor", 
                    label="color by...", 
                    choices = c(colnames(DoE)[-1]),
                    selected = NULL)
    })
    
    # output$eqFoI<-renderUI({
    #     DoE<-DoE()
    #     selectInput(inputId = "eqFoI", 
    #                 label="test the factor of...", 
    #                 choices = c(colnames(DoE)[-1]),
    #                 selected = NULL)
    # })
    # output$eqFoIlevels<-renderUI({
    #     DoE<-DoE()
    #     selectInput(inputId = "eqFoIlevels", 
    #                 label="test the level of...", 
    #                 choices = unique(DoE[,input$eqFoI]),
    #                 selected = NULL, multiple=TRUE)
    # })
    # output$eqtestres1<- renderPlot({
    #     validate(need(input$eqFoIlevels[2], "Please select two levels to conduct this test."))
    #     report_vals$eqtest1<-eqtest.row(react_fixed_data(), DoE(), input$eqFoI, 
    #                         input$eqFoIlevels[1], input$eqFoIlevels[2], 
    #                         input$lowerbound, input$upperbound)
    #     report_vals$eqtest1
    # })
    # output$eqtestres2<-renderPlot({
    #     validate(need(input$eqFoIlevels[2], "Please select two levels to conduct this test."))
    #     report_vals$eqtest2<-eqtest.all(react_fixed_data(), DoE(), input$eqFoI, 
    #                         input$eqFoIlevels[1], input$eqFoIlevels[2], 
    #                         input$lowerbound, input$upperbound)
    #     report_vals$eqtest2
    # })
    
    output$FoIlevels<-renderUI({
        DoE<-DoE()
        selectInput(inputId = "FoIlevels", 
                    label="test the level of...", 
                    choices = unique(DoE[,input$FoI]),
                    selected = NULL, multiple=TRUE)
    })
    
    listoutput<-reactiveValues()
 
    output$volcanoplot <- renderPlot({
        validate(need(input$FoIlevels[2], "Please select two levels to conduct this test."))
        listoutput[[input$selectdata_SI]]<-testingDE(react_data_collection[[input$selectdata_SI]], DoE(),input$DEtest, 
                                             input$FoI, input$FoIlevels, input$whetherblock,input$blockfactor, 
                                             input$whetherweighting, input$NAweights, react_na_index[[input$selectdata_SI]], input$whetherstandardDE)
        report_vals$volcanoplot<-listoutput[[input$selectdata_SI]][["graph"]]
        report_vals$volcanoplot
    })
    
    output$DownloadDEresults <- downloadHandler(
        filename = function() {
            paste("StatisticalTestResultsfor_", input$selectdata_SI, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(listoutput[[input$selectdata_SI]][["alldf"]], file, row.names = TRUE)
        }
    )

    # dimensionality reduction ----
    output$selectdata_DR<-renderUI({
        dataset<-c("protein data", "termini data", "peptide data", "PTM data")[
            c(input$whether_protein_file,input$whether_termini_file,input$whether_peptide_file, input$whether_PTM_file)]
        selectInput(inputId = "selectdata_DR", 
                    label="",
                    choices=dataset)
    })
    
    output$dimenreduction <- renderPlot({
        validate(need(input$selectdata_DR,""))
        if (input$DRrows=="all"){includedrows=c(1:nrow(react_data_collection[[input$selectdata_DR]]))}else{
            includedrows=listoutput[[input$selectdata_DR]][["DEdf"]]$name}
        report_vals$dmr<-dimen.reduce(react_data_collection[[input$selectdata_DR]], DoE(), input$DRmethod, input$colorfactor, input$tSNEper, includedrows)
        report_vals$dmr
    })
    
    # clustering ----
    output$selectdata_cluster<-renderUI({
        dataset<-c("protein data", "termini data", "peptide data", "PTM data")[
            c(input$whether_protein_file,input$whether_termini_file,input$whether_peptide_file, input$whether_PTM_file)]
        selectInput(inputId = "selectdata_cluster", 
                    label="",
                    choices=dataset)
    })
    
    output$cluster <- renderPlot({
        validate(need(input$selectdata_cluster,""))
        if (input$rows=="all"){includedrows=c(1:nrow(react_data_collection[[input$selectdata_cluster]]))}else{
            includedrows=listoutput[[input$selectdata_cluster]][["DEdf"]]$name}
        if (input$Cmethod=="Hierarchical Clustering"){
            report_vals$clustering<-fcluster(react_data_collection[[input$selectdata_cluster]], input$Cmethod, includedrows, input$whetherlabel, 0)
        }else{
            report_vals$clustering<-fcluster(react_data_collection[[input$selectdata_cluster]], input$Cmethod, includedrows, FALSE, input$clusternum)
        }
        report_vals$clustering
    })
    
    # Individual Protein Visualization ----
    output$selectdata_IDV<-renderUI({
        dataset<-c("protein data", "termini data", "peptide data", "PTM data")[
            c(input$whether_protein_file,input$whether_termini_file,input$whether_peptide_file, input$whether_PTM_file)]
        selectInput(inputId = "selectdata_IDV", 
                    label="Select dataset:",
                    choices=dataset)
    })

    output$proteingroupID<-renderUI({
        validate(need(input$proteingroup=="specify below...",""))
        textInput("proteingroupID", "Which to visualize? (separate by comma)", value = "")
    })
    
    rows_include<-reactive({
        validate(need(input$selectdata_IDV,""))
        if(input$proteingroup=="significant in differential expression"){
            r<-listoutput[[input$selectdata_IDV]][["DEdf"]]$name}
        if(input$proteingroup=="specify below..."){
            r<-unlist(strsplit(input$proteingroupID,","))}
        r
    })
    output$IDV_table<-DT::renderDataTable({
        validate(need(rows_include(),""))
        options = list(paging=FALSE)
        react_data_collection[[input$selectdata_IDV]][rows_include(),]
    })
    output$IDV_boxplot<-renderPlot({
        validate(need(rows_include(),""))
        IDV_plot(react_data_collection[[input$selectdata_IDV]][rows_include(),])
    })
    
    output$IDV_corrplot<-renderPlot({
        validate(need(rows_include(),""))
        do.call(corrplot_customize,list(data=react_data_collection[[input$selectdata_IDV]][rows_include(),], 
                                        corr_sign=input$corr_sign, 
                                        p_threshold=input$p_threshold, 
                                        order=input$corrorder, 
                                        ncluster=input$ncluster, 
                                        colorscheme=input$colorscheme))
    })
    
    output$IDV_corrtable<-DT::renderDataTable({
        validate(need(rows_include(),""))
        options = list(paging=FALSE)
        data<-react_data_collection[[input$selectdata_IDV]][rows_include(),]
        t_data<-t(data)
        M<-cor(t_data,use="pairwise.complete.obs")
        for (i in rownames(M)){
            if (sum(!is.na(M[,i]))<2){
                M<-M[-which(rownames(M)==i),-which(rownames(M)==i)]}
        }
        M
    })
    
    output$proteinACC<-renderUI({
    validate(need((input$selectdata_IDV=="protein data")&(!is.null(data_collection()[["PTM data"]])),""))
    selectInput("proteinACC","Select the column specifying protein accession in the PTM data:",
                    choices=c(colnames(data_collection()[["PTM data"]]$other_annotation)[-1]))
    })
    output$col_position<-renderUI({
        validate(need((input$selectdata_IDV=="protein data")&(!is.null(data_collection()[["PTM data"]])),""))
        selectInput("col_position","Select the column specifying modification position in the PTM data:",
                    choices=c(colnames(data_collection()[["PTM data"]]$other_annotation)[-1]))
    })
    output$modificationType<-renderUI({
        validate(need((input$selectdata_IDV=="protein data")&(!is.null(data_collection()[["PTM data"]])),""))
        selectInput("modificationType","Select the column specifying modification type in the PTM data:",
                    choices=c(colnames(data_collection()[["PTM data"]]$other_annotation)[-1]))
    })
    
    output$IDVseperate_panel<-renderUI({
        validate(need((input$selectdata_IDV=="protein data")&(!is.null(data_collection()[["PTM data"]])),"This function is valid only when protein data selected and PTM data available."))
        available_proteins<-rows_include()
        lapply(available_proteins,function(x){
            output[[paste0("CLG_for_",x)]]<-renderUI(
                    if (length(data_collection()[["PTM data"]]$other_annotation[input$proteinACC,]==x)==0){
                        renderText("No PTM proteing groups match this protein selected")
                    }else{
                       ind<-which(data_collection()[["PTM data"]]$other_annotation[,input$proteinACC]==x)
                       renderPlot(combined_lolipop_plot(x, data_collection()[["PTM data"]],ind, input$proteinACC,input$col_position,input$modificationType)
                            )}
                )
        })
        myTabs = lapply(available_proteins, function(x){
            tabPanel(title=x, 
                     uiOutput(paste0("CLG_for_",x)))
            })
        do.call(tabsetPanel, myTabs)        
    })
    
    output$selectdata_circos<-renderUI({
        dataset<-c("termini data", "peptide data", "PTM data")[c(input$whether_termini_file,input$whether_peptide_file, input$whether_PTM_file)]
        checkboxGroupInput("selectdata_circos",
            "please select dataset for circosplot:",
            choices = dataset)
    })
    
    output$circosplot<-renderPlot({
        validate(need(input$groupproACC,""))
        dataset<-c("termini data", "peptide data", "PTM data")[input$selectdata_circos]
        proteinACC<-unlist(strsplit(input$groupproACC,","))
        names(proteinACC)<-dataset
        circosplot.fun(listoutput,data_collection(),proteinACC,dataset)
    })
    
    # report generate ----
    output$report <- downloadHandler(
        filename = "report.html",
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            #tempReport<-file.path(tempdir(), "report.rmd")
            #file.copy("report.Rmd", tempReport, overwrite = TRUE)
            # Knit the document, passing in the `params` list, and eval it in a
            src <- normalizePath("report.rmd")
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, "report.rmd", overwrite = TRUE)
            report_vals$title<-input$title
            report_vals$description<-input$description
            params0<-list(imported=report_vals)
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            out<-rmarkdown::render("report.rmd", output_file = file,params = params0)
            file.rename(out, file)
        }
    )
}

# Run the app ----
shinyApp(ui = ui, server = server)