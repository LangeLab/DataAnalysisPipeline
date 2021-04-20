library(shinydashboard)
library(shiny)
library(shinycssloaders)
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
                        box(checkboxInput("whether_protein_file", "Upload Protein file", value=FALSE),
                        conditionalPanel('input.whether_protein_file==1',
                                         fileInput("protein_data", "please upload .csv/.xlsx file of raw data with quantitative intensities of proteins:",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx"))),
                        checkboxInput("whether_termini_file", "Upload Termini file", value=FALSE),
                        conditionalPanel('input.whether_termini_file==1',
                                         fileInput("termini_data", "please upload .csv/.xlsx file of raw termini data with quantitative intensities of terminis:",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx"))),
                        checkboxInput("whether_peptide_file","Upload Peptide file", value=FALSE),
                        conditionalPanel('input.whether_peptide_file==1',
                                         fileInput("peptide_data", "please upload .csv/.xlsx file of raw peptide data with quantitative intensities of peptides:",
                                                   multiple = FALSE,
                                                   accept = c("text/csv",
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv", ".xlsx")))),
                        
                        box(h4("Meta file input"),
                            fileInput("meta_data", "please upload .csv/.xlsx file of meta data with only sample names matched with raw data, and corresponding covariates of interest:",
                                      multiple = FALSE,
                                      accept = c("text/csv",
                                                 "text/comma-separated-values,text/plain",
                                                 ".csv", ".xlsx")),
                            checkboxInput("whether_replica", "Is there a column indicating replica in the meta data file?", value=FALSE),
                            conditionalPanel('input.whether_replica==1', 
                                             uiOutput("replicacol"),
                                             checkboxInput("whether_average_replica",
                                                           "Average the intensities across replica?",
                                                           value=FALSE)
                                             ),
                            h4("Proceed the analysis with"),
                            uiOutput("selectdata")),
                        
                        box(h4("Annotation tools"),
                            conditionalPanel('input.whether_termini_file==1',
                            uiOutput("termini_idf"),
                            uiOutput("termini_seq"),
                            numericInput("lseqwindow", 
                                         h5("Length of sequence window:"), 
                                         value = 1),
                            actionButton("annoButton", "Get annotation"),
                            div(dataTableOutput("annotationpreview")%>% withSpinner(color="#0dc5c1"), style = "font-size:60%"),
                            downloadButton("downloadAnnotation", "Download full annotation"))),
                        
                        box(h4("Protein intensity calculator based on peptide"),
                            conditionalPanel('input.whether_peptide_file==1',
                            selectInput("sumform", "Choose method of sum:",choices=c("sum","weighted","top3")),
                            actionButton("pepcalButton", "Calculate"),
                            div(dataTableOutput("peptidepreview")%>% withSpinner(color="#0dc5c1"), style = "font-size:60%"),
                            checkboxInput("whether_replace_protein", "Use the dataset calculated from peptides to replace the protein data?", 
                            value=FALSE)))
                        )),
                
                tabItem(tabName = "SampleQC",
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

                        plotOutput("viewviolin")%>% withSpinner(color="#0dc5c1"),
                        downloadButton("DownloadProcessedData", "Download")),
                tabItem(tabName="StatInf",
                        # Equivalence test
                        h4("3.1 Inference on Equivalence"),
                        box(
                            selectInput(inputId = "dd", 
                                    label="Specify the way of testing...", 
                                    choices = c("By column, test all proteins at the same time",
                                                "By row, test each protein respectively")),
                            numericInput("lowerbound", "Lower bound:", -0.5),
                            numericInput("upperbound", "Upper bound:", 0.5),
                            uiOutput("eqFoI"),
                            uiOutput("eqFoIlevels"),
                            conditionalPanel('input.dd=="By row, test each protein respectively"',
                            plotOutput("eqtestres1")%>% withSpinner(color="#0dc5c1")),
                            conditionalPanel('input.dd=="By column, test all proteins at the same time"',
                             plotOutput("eqtestres2"))
                        ),
                        # DE tests
                        h4("3.2 Inference on Differential Expression"),
                        box(
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
                                    selected = NULL),
                            uiOutput("FoI"),
                            uiOutput("FoIlevels"),
                            h5("Volcano plot of DE testing results"),
                            plotOutput("volcanoplot")%>% withSpinner(color="#0dc5c1"),
                            downloadButton("DownloadDEresults","Download significant features"))
                        ),
                tabItem(tabName = "DR",
                        # Dimensionality Reduction
                        h4("4. Dimensionality reduction"),
                        selectInput(inputId = "DRrows", 
                                    label="use the rows of...",
                                    choices=c("all","significant")),
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
    annotation_res<-eventReactive(input$annoButton,{
        annotationtool(termini_data()[,input$termini_idf], termini_data()[,input$termini_seq], input$lseqwindow)}
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
    
    peptide_data<-reactive({
        req(input$peptide_data)
        if(grepl(".xlsx",input$peptide_data$name)){
            df <- read.xlsx(input$peptide_data$datapath)}
        if (grepl(".csv",input$peptide_data$name)){
            df <- read.csv(input$peptide_data$datapath, header=TRUE, check.names=FALSE)
        }
        data.frame(df)
    })
    
    peptide_based_fixed_data<-eventReactive(
        input$pepcalButton,{
        calIntPep(peptide_data(),input$sumform)
    })
    
    output$peptidepreview<-renderDataTable({
        peptide_based_fixed_data()[c(1:20),c(1:6)]
    })
    

    output$selectdata<-renderUI({
        dataset<-c("protein data", "termini data", "peptide data")[c(input$whether_protein_file,input$whether_termini_file,input$whether_peptide_file)]
        selectInput(inputId = "selectdata", 
                    label="",
                    choices=dataset)
    })
    
    fixed_data0<-reactive({
        if(input$selectdata=="protein data"){df<-protein_data()}
        if(input$selectdata=="termini data"){df<-termini_data()}
        if(input$selectdata=="peptide data"){df<-peptide_data()}
        df<-inputraw(df)
        if(input$whether_replace_protein==TRUE){df<-peptide_based_fixed_data()}
        df
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
    
    
    DoE0<-reactive({
        req(input$meta_data)
        if(grepl(".xlsx",input$meta_data$name)){
            df <- read.xlsx(input$meta_data$datapath)}
        if (grepl(".csv",input$meta_data$name)){
            df <- read.csv(input$meta_data$datapath, header=TRUE, check.names=FALSE)
        }
        dt<-setdoe(data.frame(df), fixed_data0())
        dt
    }) 
    
    fixed_data<-reactive({
        df<-fixed_data0()
        if(input$whether_average_replica==TRUE){
            df<-averagereplica(df, replica_list())
        }
        data.frame(df)
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
    
    output$graph1 <- renderPlot({
        report_vals$cvviolin<-cvplots(fixed_data(),DoE())
        report_vals$cvviolin})
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
    
    react_res_pp<-reactive({
        dt<-preprocessing(fixed_data(), DoE(), input$filterlevel, input$normalization,
                                              input$whetherimpute, input$imputation)
        dt
        })
    
    react_fixed_data<-reactive({
        df<-react_res_pp()[["data"]]
        df
        })
    na_index<-reactive({react_res_pp()[["na.index"]]})
    
    output$DownloadProcessedData <- downloadHandler(
        filename = function() {
            paste("ProcessedData", ".csv", sep = "")
        },
        content = function(file) {
            write.csv(react_fixed_data(), file, row.names = FALSE)
        }
    )
    
    output$viewviolin <- renderPlot({
        ctitle<-paste0(c("filtered for ",input$filterlevel,"% completeness, median normalization ", input$normalization, 
                         ", imputation ", input$whetherimpute, " by ", input$imputation), collapse="")
        report_vals$violinplot<-plotviolin(react_fixed_data(),ctitle)
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
        validate(need(input$eqFoIlevels[2], "Please select two levels to conduct this test."))
        report_vals$eqtest1<-eqtest.row(react_fixed_data(), DoE(), input$eqFoI, 
                            input$eqFoIlevels[1], input$eqFoIlevels[2], 
                            input$lowerbound, input$upperbound)
        report_vals$eqtest1
    })
    output$eqtestres2<-renderPlot({
        validate(need(input$eqFoIlevels[2], "Please select two levels to conduct this test."))
        report_vals$eqtest2<-eqtest.all(react_fixed_data(), DoE(), input$eqFoI, 
                            input$eqFoIlevels[1], input$eqFoIlevels[2], 
                            input$lowerbound, input$upperbound)
        report_vals$eqtest2
    })
    
    output$FoIlevels<-renderUI({
        DoE<-DoE()
        selectInput(inputId = "FoIlevels", 
                    label="test the level of...", 
                    choices = unique(DoE[,input$FoI]),
                    selected = NULL, multiple=TRUE)
    })
    
    listoutput<-reactive({
        validate(need(input$FoIlevels[2], "Please select two levels to conduct this test."))
        testingDE(react_fixed_data(), DoE(),input$DEtest, 
                                    input$FoI, input$FoIlevels, input$whetherblock,input$blockfactor, 
                                    input$whetherweighting, input$NAweights, na_index(), input$whetherstandardDE)})
    
    output$DownloadDEresults <- downloadHandler(
        filename = function() {
            paste("DifferentialExpressionTestResults", ".csv", sep = "")
        },
        content = function(file) {
            write.csv(listoutput()[["DEdf"]], file, row.names = FALSE)
        }
    )

    output$volcanoplot <- renderPlot({
        report_vals$volcanoplot<-listoutput()[["graph"]]
        report_vals$volcanoplot
    })
    
    
    output$dimenreduction <- renderPlot({
        if (input$DRrows=="all"){includedrows=c(1:nrow(react_fixed_data()))}else{
            includedrows=listoutput()[["DEdf"]]$name}
        report_vals$dmr<-dimen.reduce(react_fixed_data(), DoE(), input$DRmethod, input$colorfactor, input$tSNEper, includedrows)
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
            #tempReport<-file.path(tempdir(), "report.rmd")
            #file.copy("report.Rmd", tempReport, overwrite = TRUE)
            # Knit the document, passing in the `params` list, and eval it in a
            src <- normalizePath("report.rmd")
            owd <- setwd(tempdir())
            on.exit(setwd(owd))
            file.copy(src, "report.rmd", overwrite = TRUE)
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