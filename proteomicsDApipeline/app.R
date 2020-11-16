library(shiny)
source("DataInput.R")
source("QCPlots.R")
source("Preprocessing.R")
source("testingDE.R")
# Define UI ----
ui <- fluidPage(
    titlePanel("Proteomics Data Analysis Pipeline"),
    sidebarLayout(
        sidebarPanel(
        h4("1. Graphs at QC stage"),
        #QC Plots
        selectInput(inputId = "options", 
        label="Quality Control graph options", 
        choices = c( "%CV plots for QC samples",
                     "%CV plots for QC samples stacked bar",
                      "Distribution of identified proteins",
                      "Overlap of identified proteins",
                      "Data Completeness"),
        selected = NULL),
        selectInput(inputId = "group",
        label="Select group column for grouped graphs",
        choices = colnames(DoE),
        selected = NULL),
        
        # Normalization, imputation, etc.
        h4("2. Processing"),
        numericInput("filterlevel", 
                     h5("Row filter level (%)"), 
                     value = 0),
        checkboxInput(inputId = "normalization", 
                      "Median Normalization", value = FALSE),
        checkboxInput(inputId = "whetherimpute", 
                      "Imputation", value = FALSE),
        selectInput(inputId = "imputation", 
                    label="Imputation Method", 
                    choices = c( "Down-shifted Normal samples",
                                 "MinProb",
                                 "knn",
                                 "min"),
                    selected = NULL),
        # tests
        h4("3. Inference on Differential Expression"),
        checkboxInput(inputId = "whetherblock", 
                      "Blocking", value = FALSE),
        selectInput(inputId = "blockfactor", 
                    label="Blocking on...(supports 2-level factors only)", 
                    choices = c(colnames(DoE)[-1]),
                    selected=NULL),
        selectInput(inputId = "DEtest", 
                    label="Statistical testing for DE", 
                    choices = c( "limma",
                                 "t-test"),
                    selected = NULL),
        selectInput(inputId = "FoI", 
                    label="test the factor of...", 
                    choices = c(colnames(DoE)[apply(DoE, 2, function(x) length(unique(x)))==2]),
                    selected = NULL),
        actionButton("update", "Update View"),
        
        # Unsupervised learning
        h4("4. Unsupervised learning"),
        checkboxInput(inputId = "whethercluster", 
                      "Clustering", value = FALSE),
        checkboxInput(inputId = "whethertsne", 
                      "T-SNE", value = FALSE),
        ),
        mainPanel(
            plotOutput("graph"),
            tableOutput("view"),
            plotOutput("volcanoplot")
            )
    )  
)

# Define server logic ----
server <- function(input, output) {
    fixed_data<-
    output$graph <- renderPlot({
        if (input$options=="%CV plots for QC samples"){
            g<-cvplots(1)}
        if (input$options=="%CV plots for QC samples stacked bar"){
            g<-cvplots(2)}
        if (input$options=="Distribution of identified proteins"){
            g<-distIndProtein(input$group)}
        if (input$options=="Overlap of identified proteins"){
            g<-upsetplot(input$group)}
        if (input$options=="Data Completeness"){
            g<-datacompleteness()}
    g
    })
    react_fixed_data<-reactive({preprocessing(input$filterlevel, input$normalization,
                                              input$whetherimpute, input$imputation)})
    output$view <- renderTable({
    sub.head<-head(react_fixed_data(), n=3)
    sub.head[,1:10]
    })
    output$volcanoplot <- renderPlot({
    g<-testingDE(react_fixed_data(), input$DEtest, input$FoI, input$whetherblock,input$blockfactor)
    g
    })
}

# Run the app ----
shinyApp(ui = ui, server = server)