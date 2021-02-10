library(shinydashboard)
library(shiny)
source("DataInput.R")
source("QCplots.R")
source("Preprocessing.R")
source("testingDE.R")
source("DimenReduce.R")
source("UnsupervisedLearning.R")
source("annotationtool.R")
source("EquivalenceTest.R")
source("PeptideCalculate.R")
options(shiny.maxRequestSize=20*1024^2)
# Define UI ----
ui <- dashboardPage(
    dashboardHeader(title="Proteomics Data Analysis Pipeline"),
    dashboardSidebar(sidebarMenu(
        menuItem("Input and annotation", tabName = "IA"),
        menuItem("Sample based QC", tabName = "SampleQC"),
        menuItem("Pre-processing", tabName = "PP"),
        menuItem("Statistical Inference", tabName = "StatInf"),
        menuItem("Dimensionality reduction", tabName = "DR"),
        menuItem("Clustering", tabName = "Clustering"),
        menuItem("Generate a report", tabName="GenerateReport"))),
        dashboardBody(
            tabItems(
                tabItem(tabName="IA",
                        h4("File input"),
                        fluidRow(
                        box(fileInput("raw_data", "please upload .csv/.xlsx file of raw data with quantitative intensities of proteins:",
                                   multiple = FALSE,
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv", ".xlsx")),
                         fileInput("meta_data", "please upload .csv/.xlsx file of meta data with only sample names matched with raw data, and corresponding covariates of interest:",
                                   multiple = FALSE,
                                   accept = c("text/csv",
                                              "text/comma-separated-values,text/plain",
                                              ".csv", ".xlsx"))),
                        h4("Supplemental file input"),
                        box(
                        checkboxInput("whether_termini_file", "Upload Termini file", value=FALSE),
                        conditionalPanel('input.whether_termini_file==1',
                                         fileInput("termini_data", "please upload .csv/.xlsx file of raw termini data with quantitative intensities of terminis:",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx")),
                                         h4("Annotation tools"),
                                         uiOutput("termini_idf"),
                                         uiOutput("termini_seq"),
                                         numericInput("lseqwindow", 
                                                      h5("Length of sequence window:"), 
                                                      value = 1),
                                         div(dataTableOutput("annotationpreview"), style = "font-size:60%"),
                                         downloadButton("downloadAnnotation", "Download full annotation"))),
                        box(
                            checkboxInput("whether_peptide_file","Upload Peptide file", value=FALSE),
                            conditionalPanel('input.whether_peptide_file==1',
                                             fileInput("peptide_data", "please upload .csv/.xlsx file of raw peptide data with quantitative intensities of peptides:",
                                                       multiple = FALSE,
                                                       accept = c("text/csv",
                                                                  "text/comma-separated-values,text/plain",
                                                                  ".csv", ".xlsx")),
                                             h4("Protein intensity calculator based on peptide"),
                                             selectInput("sumform", "Choose method of sum:",choices=c("sum","weighted","top3")),
                                             div(dataTableOutput("peptidepreview"), style = "font-size:60%"),
                                             checkboxInput("whether_replace_protein", "Use the dataset calculated from peptides to replace the protein data?", 
                                                           value=FALSE))
                        )
                        

                        )),
                tabItem(tabName = "SampleQC",
                        tabBox(title = "QC plots",
                               tabPanel("%CV plots for QC samples",
                                        plotOutput("graph1")),
                               tabPanel("%CV plots for QC samples stacked bar", 
                                        plotOutput("graph2")),
                               tabPanel("Distribution of identified proteins",
                                        uiOutput("group"),
                                        plotOutput("graph3")),
                               tabPanel("Overlap of identified proteins",
                                        uiOutput("group2"),
                                        plotOutput("graph4")),
                               tabPanel("Data Completeness",
                                        plotOutput("graph5")),
                               tabPanel("Correlation Plot",
                                         plotOutput("graph6"))
                )),
                tabItem(tabName = "PP",
                        h4("2. Pre-processing"),
                        # Normalization, imputation, etc.
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
                        plotOutput("viewviolin")),
                tabItem(tabName="StatInf",
                        # Equivalence test
                        h4("3.1 Inference on Equivalence"),
                        box(
                            selectInput(inputId = "dd", 
                                    label="Specify the way of testing...", 
                                    choices = c("By row, test each protein respectively",
                                                "By column, test all proteins at the same time")),
                            numericInput("lowerbound", "Lower bound:", -0.5),
                            numericInput("upperbound", "Upper bound:", 0.5),
                            uiOutput("eqFoI"),
                            uiOutput("eqFoIlevels"),
                            conditionalPanel('input.dd=="By row, test each protein respectively"',
                            plotOutput("eqtestres1")),
                            conditionalPanel('input.dd=="By column, test all proteins at the same time"',
                             plotOutput("eqtestres2"))
                        ),
                        # DE tests
                        h4("3.2 Inference on Differential Expression"),
                        box(
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
                                    selected = NULL),
                            uiOutput("FoI"),
                            h5("Volcano plot of DE testing results"),
                            plotOutput("volcanoplot"))
                        ),
                tabItem(tabName = "DR",
                        # Dimensionality Reduction
                        h4("4. Dimensionality reduction"),
                        selectInput(inputId = "DRmethod", 
                                    label="use the method of...",
                                    choices=c("PCA","t-SNE", "UMAP")),
                        uiOutput("colorfactor"),
                        conditionalPanel('input.DRmethod=="t-SNE"',
                                         numericInput("tSNEper","Specify the perplexity for t-SNE", 0)),
                        plotOutput("dimenreduction")),
                tabItem(tabName="Clustering",
                        # clustering
                        h4("5. Clustering"),
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
                                                      value = 4)),
                        plotOutput("cluster")),
                tabItem(tabName="GenerateReport",
                        downloadButton("report", "Generate report"))
            )
    )
)
# Define server logic ----
server <- function(input, output) {
    report_vals <- reactiveValues()
    termini_data<-reactive({
        req(input$termini_data)
        if(grepl(".xlsx",input$termini_data$name)){
            df <- read.xlsx(input$termini_data$datapath)}
        if (grepl(".csv",input$termini_data$name)){
            df <- read.csv(input$termini_data$datapath, header=TRUE)
        }
        df
    })
    output$termini_idf<-renderUI({
        coln<-colnames(termini_data())
        selectizeInput(inputId = "termini_idf",
                       label="Select the column of identifier for termini data",
                       choices = coln)
    })
    output$termini_seq<-renderUI({
        coln<-colnames(termini_data())
        selectizeInput(inputId = "termini_seq",
                       label="Select the column of stripped sequence for termini data",
                       choices = coln)
    })
    annotation_res<-reactive({
        annotationtool(termini_data()[,input$termini_idf], termini_data()[,input$termini_seq], input$lseqwindow)
    })
    
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
    
    peptide_data<-reactive({
        req(input$peptide_data)
        if(grepl(".xlsx",input$peptide_data$name)){
            df <- read.xlsx(input$peptide_data$datapath)}
        if (grepl(".csv",input$peptide_data$name)){
            df <- read.csv(input$peptide_data$datapath, header=TRUE)
        }
        df
    })
    
    peptide_based_fixed_data<-reactive({
        calIntPep(peptide_data(),input$sumform)
    })
    
    output$peptidepreview<-renderDataTable({
        peptide_based_fixed_data()[c(1:20),c(1:6)]
    })
    
    fixed_data<-reactive({
        req(input$raw_data)
        if(grepl(".xlsx",input$raw_data$name)){
            df <- read.xlsx(input$raw_data$datapath)}
        if (grepl(".csv",input$raw_data$name)){
            df <- read.csv(input$raw_data$datapath, header=TRUE)
        }
        if(input$whether_replace_protein==FALSE){
        inputraw(df)}else{
            peptide_based_fixed_data()
        }
    })
    
    DoE<-reactive({
        req(input$meta_data)
        if(grepl(".xlsx",input$meta_data$name)){
            df <- read.xlsx(input$meta_data$datapath)}
        if (grepl(".csv",input$meta_data$name)){
            df <- read.csv(input$meta_data$datapath, header=TRUE)
        }
        setdoe(df, fixed_data())
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
    
    output$graph1 <- renderPlot({
        report_vals$cvviolin<-cvplots(fixed_data(),DoE(),1)
        report_vals$cvviolin})
    output$graph2 <- renderPlot({
        report_vals$cvstacked<-cvplots(fixed_data(),DoE(),2)
        report_vals$cvstacked})
    output$graph3 <- renderPlot({
        report_vals$distprotein<-distIndProtein(fixed_data(),DoE(),input$group)
        report_vals$distprotein})
    output$graph4 <- renderPlot({
        report_vals$upset<-upsetplot(fixed_data(),DoE(),input$group2)
        report_vals$upset})
    output$graph5 <- renderPlot({
        report_vals$datacompleteness<-datacompleteness(fixed_data(),DoE())
        report_vals$datacompleteness})
    output$graph6 <- renderPlot({
        report_vals$corplot<-corplot(fixed_data())
        report_vals$corplot})
    
    react_res_pp<-reactive({preprocessing(fixed_data(), DoE(), input$filterlevel, input$normalization,
                                              input$whetherimpute, input$imputation)})
    react_fixed_data<-reactive({react_res_pp()[["data"]]})
    na_index<-reactive({react_res_pp()[["na.index"]]})
    
    output$viewviolin <- renderPlot({
        report_vals$violinplot<-plotviolin(react_fixed_data())
        report_vals$violinplot
    })
    
    # equivalence test FoI and FoI levels
    output$eqFoI<-renderUI({
        DoE<-DoE()
        selectInput(inputId = "eqFoI", 
                    label="test the factor of...", 
                    choices = c(colnames(DoE)[-1]),
                    selected = NULL)
    })
    output$eqFoIlevels<-renderUI({
        DoE<-DoE()
        selectInput(inputId = "eqFoIlevels", 
                    label="test the level of...", 
                    choices = unique(DoE[,input$eqFoI]),
                    selected = NULL, multiple=TRUE)
    })
    output$eqtestres1<- renderPlot({
        report_vals$eqtest1<-eqtest.row(react_fixed_data(), DoE(), input$eqFoI, 
                            input$eqFoIlevels[1], input$eqFoIlevels[2], 
                            input$lowerbound, input$upperbound)
        report_vals$eqtest1
    })
    output$eqtestres2<-renderPlot({
        report_vals$eqtest2<-eqtest.all(react_fixed_data(), DoE(), input$eqFoI, 
                            input$eqFoIlevels[1], input$eqFoIlevels[2], 
                            input$lowerbound, input$upperbound)
        report_vals$eqtest2
    })
    listoutput<-reactive({testingDE(react_fixed_data(), DoE(),input$DEtest, 
                                    input$FoI, input$whetherblock,input$blockfactor, 
                                    input$whetherweighting, input$NAweights, na_index())})

    output$volcanoplot <- renderPlot({
        report_vals$volcanoplot<-listoutput()[["graph"]]
        report_vals$volcanoplot
    })
    
    output$dimenreduction <- renderPlot({
        report_vals$dmr<-dimen.reduce(react_fixed_data(), DoE(), input$DRmethod, input$colorfactor, input$tSNEper)
        report_vals$dmr
    })
    
    output$cluster <- renderPlot({
        if (input$rows=="all"){includedrows=c(1:nrow(react_fixed_data()))}else{
            includedrows=listoutput()[["DEdf"]]$name}
        if (input$Cmethod=="Hierarchical Clustering"){
            report_vals$clustering<-fcluster(react_fixed_data(), input$Cmethod, includedrows, input$whetherlabel, 0)
        }else{
            report_vals$clustering<-fcluster(react_fixed_data(), input$Cmethod, includedrows, FALSE, input$clusternum)
        }
        report_vals$clustering
    })
    
    output$report <- downloadHandler(
        filename = "report.html",
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path("/srv/shiny-server/DA/report.Rmd")
            #tempReport<-file.path(tempdir(), "report.rmd")
            file.copy("report0.Rmd", tempReport, overwrite = TRUE)
            # Knit the document, passing in the `params` list, and eval it in a
            params0<-list(imported=report_vals)
            # child of the global environment (this isolates the code in the document
            # from the code in this app).
            rmarkdown::render(tempReport, output_file = file,
                              params = params0,
                              envir = new.env(parent = globalenv())
            )
        }
    )
}

# Run the app ----
shinyApp(ui = ui, server = server)