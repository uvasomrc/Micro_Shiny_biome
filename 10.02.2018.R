#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(markdown)
library(dplyr)
library(tibble)
source("utils.R")

# Define UI for application that draws a histogram
ui <- navbarPage(
  "Shiny Microbiome",
  
  tabPanel("Upload Data",
           sidebarLayout(
             sidebarPanel(fileInput("countdf","Upload 'Count Data' file"),
                          fileInput("taxodf", "Upload 'Taxonomy Data' file"),
                          fileInput("sampledf", "Upload 'Sample Data' file"),
                          helpText("Select the read.table parameters below"),
                          checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
                          checkboxInput(inputId = "stringAsFactors", "stringAsFactors", TRUE),
                          radioButtons(inputId = 'sep', label = 'Separator', choices = c(Tab='\t', Comma=',',Semicolon=';', Space=''), selected = '\t'),
                          actionButton("checkbutton", label = "Check!!!")),
             mainPanel(  # Summary of file
               conditionalPanel(condition = "input.checkbutton == 0",htmlOutput("help")),
               # conditionalPanel(condition = "input.checkbutton == 1",myfunction())
               conditionalPanel(condition = "input.checkbutton == 1",htmlOutput("dimensions"))
               
             ))),
  
  tabPanel("Filtering",
           sidebarLayout(
             sidebarPanel(radioButtons("raw_or_not", "Type of data", choices = c("raw","filtered")),
                          uiOutput("ui_reads")),
             mainPanel(    
               DT::dataTableOutput("table")
             ))),
  tabPanel("visualization1"),
  tabPanel("visualization2"),
           
  navbarMenu("More",
             tabPanel("Feature1"), tabPanel("Feature2"))
  
)




server <- function(input, output) {
  
  original_data <- reactive({
    cdf <- read.table(file=input$countdf$datapath, sep=input$sep, header = input$header,
                      stringsAsFactors = input$stringAsFactors)
    tdf <- read.table(file=input$taxodf$datapath, sep=input$sep, header = input$header,
                      stringsAsFactors = input$stringAsFactors)
    sdf <- read.table(file=input$sampledf$datapath, sep=input$sep, header = input$header,
                      stringsAsFactors = input$stringAsFactors)
    
    cdf2 <- cdf[order(cdf$OTUId),]
    rownames(cdf2) <- cdf2[,1]
    cdf2[,1] <- NULL
    
    tdf2 <- tdf[order(tdf$X),]
    rownames(tdf2) <- tdf2[,1]
    tdf2[,1] <- NULL
    
    sdf2 <- sdf[order(sdf$sample),]
    rownames(sdf2) <- sdf2[,1]
    sdf2[,1] <- NULL
    
    validate(
      need(length(rownames(cdf2)) == length(rownames(tdf2))  &
             length(colnames(cdf2)) == length(rownames(sdf2)) &
             all(ifelse(rownames(cdf2) == rownames(tdf2), TRUE, FALSE)) &
             all(ifelse(colnames(cdf2) == rownames(sdf2), TRUE, FALSE)),
           "The datasets are invalid")
    )
    
    list(cdf3=cdf2, tdf3=tdf2, sdf3=sdf2)
    
  })
  
  
  
  
  
  filt1 <- reactive({
    original_data()$cdf3[colSums(original_data()$cdf3) >input$reads]
  })
  
  
  
  
  filt2 <- reactive({
    afilter <- as.data.frame(filt1() / colSums(filt1()) > input$aprop)
    
    kfil_rname <- afilter %>% 
      rownames_to_column('OTU') %>%
      filter(rowSums(afilter)/ncol(afilter) > input$sprop) %>%
      column_to_rownames('OTU') %>% 
      rownames()
    
    cdf_fil2 <- subset(filt1(), rownames(filt1()) %in% kfil_rname )
    
    list(cdf_fil2=cdf_fil2, kfil_rname=kfil_rname)
    
  })
  

  
  
  output$dimensions <- renderUI({
    # req(!is.null(original_data()$cdf3) & !is.null(original_data()$tdf3) & !is.null(original_data()$sdf3))
    if(!is.null(original_data()$cdf3) & !is.null(original_data()$tdf3) & !is.null(original_data()$sdf3)){
      str0 <- "<h2>Dimensions of Datasets:</h2>"
      str1 <- paste(" 1) The Count Data has  :::  ", ncol(original_data()$cdf3), " samples and ", nrow(original_data()$cdf3), " OTUs")
      str2 <- paste(" 2) The Taxonomy Data has  :::  ", ncol(original_data()$tdf3), " taxonomic levels")
      str3 <- paste(" 3) The Sample Data has  :::  ", ncol(original_data()$sdf3), " sample fields")
      HTML(paste(str0,str1,str2,str3, sep='<br/>'))
    } 
  })
  
  
  
  
  output$ui_reads <- renderUI({
    req(!is.null(original_data()$cdf3) & !is.null(original_data()$tdf3) & !is.null(original_data()$sdf3))
    if(input$raw_or_not == "filtered") {
      list(
        hr(),
        helpText("Filter by Minimum Reads"),
        numericInput("reads", "Minimum Reads:", 0, min = 0, max = max(colSums(original_data()$cdf3))),
        hr(),
        helpText("Filter by Taxanomy Prevalence: K over A"),
        helpText("Proportion of Taxamonmy Abundance(0 to 1) : A"),
        numericInput("aprop", "Minimum Taxanomy Abundance:", 0, min = 0, max = 1),
        helpText("Proportion of Samples (0 to 1) : K "),
        numericInput("sprop", "Minimum Proportion of Samples:", 0, min = 0, max = 1))
      
      # actionButton("readsfilter", "Filter by Reads"))
      
    }
    
  })
  
  
  
  
  output$table <- DT::renderDataTable({ 
    if(input$raw_or_not == "raw") {
      original_data()$cdf3 
      
    } else if (input$raw_or_not == "filtered") {
      req(!is.null(input$reads) & !is.null(input$aprop) & !is.null(input$sprop))
      # filt1() & 
      filt2()$cdf_fil2
    }
    
  })
  
  
  
  output$help <- renderText({
    "<h2>To Check Validity of the Datasets:</h2>
    <ol>
    <li>Check the relevant parameter and separator setting</li>
    <li>Browse and upload the datasets</li>
    <li>Press Check!!!</li>
    </ol>"
  })
  
}




# Run the application 
shinyApp(ui = ui, server = server)

