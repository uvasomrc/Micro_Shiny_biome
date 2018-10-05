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
library(stringi)
library(reshape2)
library(ggplot2)
library(scales)
library(plotly)

# source("utils.R")


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
                          uiOutput("ui_filter")),
             mainPanel(    
               DT::dataTableOutput("table")
             ))),
  tabPanel("Stacked Bar",
           sidebarLayout(
             sidebarPanel(radioButtons("raw_or_not2", "Type of data", choices = c("raw","filtered")),
                          uiOutput("ui_bar")),
             mainPanel(    
               plotlyOutput("stackedbar")
             ))),
  tabPanel("taxo_visual_data",DT::dataTableOutput("table_bar")),
  tabPanel("tdf table",DT::dataTableOutput("table_tdf")),
  tabPanel("sdf table",DT::dataTableOutput("table_sdf")),
  
  navbarMenu("More",
             tabPanel("Feature1"), tabPanel("Feature2"))
  
)




server <- function(input, output) {
  
  ## Import and store uploaded original datasets
  original_data <- reactive({
    req(!is.null(input$countdf$datapath) & !is.null(input$taxodf$datapath) & !is.null(input$sampledf$datapath))
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
    
    list(cdf_orig=cdf2, tdf_orig=tdf2, sdf_orig=sdf2)
  })
  
  
  ## Filtering process on the original datasets
  filtered_data <- reactive({
    
    reads_filter <- original_data()$cdf_orig[colSums(original_data()$cdf_orig) >input$reads]
    a_filter <- as.data.frame(reads_filter / colSums(reads_filter) > input$aprop)
    
    kfil_rname <- a_filter %>% 
      rownames_to_column('OTU') %>%
      filter(rowSums(a_filter)/ncol(a_filter) > input$sprop) %>%
      column_to_rownames('OTU') %>% 
      rownames()
    
    cdf_fil2 <- subset(reads_filter, rownames(reads_filter) %in% kfil_rname )
    tdf_fil <- subset(original_data()$tdf_orig, rownames(original_data()$tdf_orig) %in% kfil_rname )
    sdf_fil <- subset(original_data()$sdf_orig, rownames(original_data()$sdf_orig) %in% colnames(a_filter))
    
    
    list(cdf_filt=cdf_fil2, tdf_filt=tdf_fil, sdf_filt=sdf_fil) 
  })
  
  
  ## Instruction on how to upload datasets
  output$help <- renderText({
    "<h3>To Check Validity of the Datasets:</h3>
    <ol>
    <li>Check the relevant parameter and separator setting</li>
    <li>Browse and upload the datasets</li>
    <li>Press Check!!!</li>
    </ol>"
  })
  
  
  ## Showing dimensions of datasets
  output$dimensions <- renderUI({   
    str0 <- "<h2>Dimensions of Datasets:</h2>"
    str1 <- paste(" 1) The Count Data has  :::  ", ncol(original_data()$cdf_orig), " samples and ", nrow(original_data()$cdf_orig), " OTUs")
    str2 <- paste(" 2) The Taxonomy Data has  :::  ", ncol(original_data()$tdf_orig), " taxonomic levels")
    str3 <- paste(" 3) The Sample Data has  :::  ", ncol(original_data()$sdf_orig), " sample fields")
    HTML(paste(str0,str1,str2,str3, sep='<br/>'))
    
  })
  
  
  ## Filter tab UI
  output$ui_filter <- renderUI({
    if(input$raw_or_not == "filtered") {
      list(
        hr(),
        helpText("Filter by Minimum Reads"),
        numericInput("reads", "Minimum Reads:", 0, min = 0, max = max(colSums(original_data()$cdf_orig))),
        hr(),
        helpText("Filter by Taxanomy Prevalence: K over A"),
        helpText("Proportion of Taxamonmy Abundance(0 to 1) : A"),
        numericInput("aprop", "Minimum Taxanomy Abundance:", 0, min = 0, max = 1),
        helpText("Proportion of Samples (0 to 1) : K "),
        numericInput("sprop", "Minimum Proportion of Samples:", 0, min = 0, max = 1))
    }
    
  })
  
  ## Show original and filtered OTU table
  output$table <- DT::renderDataTable({ 
    if(input$raw_or_not == "raw") {
      original_data()$cdf_orig 
      
    } else if (input$raw_or_not == "filtered") {
      req(!is.null(input$reads) & !is.null(input$aprop) & !is.null(input$sprop))
      filtered_data()$cdf_filt
    }
    
  })
  
  
  ## show original and filtered taxo table
  output$table_tdf <- DT::renderDataTable({ 
    if(input$raw_or_not == "raw") {
      original_data()$tdf_orig 
      
    } else if (input$raw_or_not == "filtered") {
      req(!is.null(input$reads) & !is.null(input$aprop) & !is.null(input$sprop))
      filtered_data()$tdf_filt
    }
  })
  
  
  ## show original and filtered sample table
  output$table_sdf <- DT::renderDataTable({ 
    
    if(input$raw_or_not == "raw") {
      original_data()$sdf_orig
      
    } else if (input$raw_or_not == "filtered") {
      req(!is.null(input$reads) & !is.null(input$aprop) & !is.null(input$sprop))
      filtered_data()$sdf_filt
    }
  })
  
  
  ## stacked bar tab UI
  output$ui_bar <- renderUI({
    #req(!is.null(original_data()$cdf3) & !is.null(original_data()$tdf3) & !is.null(original_data()$sdf3))
    #if(input$raw_or_not2 == "filtered") {
    list(
      hr(),
      selectInput("taxo_level", label = h5("Select Taxonomic level of interest"),
                  choices = colnames(original_data()$tdf_orig)),
      hr())
  })
  
  
  ## Data wrangling for stacked bar
  visual_data <- reactive({
    
    if(input$raw_or_not2 == "raw") {
      
      taxo_raw_vis <- as.data.frame(original_data()$tdf_orig[input$taxo_level])
      colnames(taxo_raw_vis) <- c('taxo_class')
      taxo_raw_vis <- as.data.frame(lapply(taxo_raw_vis , function(x) {gsub(".*unclassified.*", "Unclassified", x)}))  # all factor
      cdf_raw_vis <- as.data.frame(original_data()$cdf_orig)
      cdf_raw_vis <- as.data.frame(lapply(cdf_raw_vis, as.numeric))
      df_raw <- cbind(cdf_raw_vis, taxo_raw_vis)
      
      df_raw2 <- df_raw %>%
        group_by(taxo_class) %>%
        summarise_all(funs(sum)) %>%
        as.data.frame()
      
      rownames(df_raw2) <- df_raw2[,1]
      df_raw2[,1] <- NULL
      df_raw2
      
    } else if(input$raw_or_not2 == "filtered") {
      req(!is.null(input$reads) & !is.null(input$aprop) & !is.null(input$sprop))
      
      taxo_filt_vis <- as.data.frame(filtered_data()$tdf_filt[input$taxo_level])
      colnames(taxo_filt_vis) <- c('taxo_class')
      taxo_filt_vis <- as.data.frame(lapply(taxo_filt_vis , function(x) {gsub(".*unclassified.*", "Unclassified", x)}))  # all factor
      cdf_filt_vis <- as.data.frame(filtered_data()$cdf_filt)
      cdf_filt_vis <- as.data.frame(lapply(cdf_filt_vis, as.numeric))
      df_filt <- cbind(cdf_filt_vis, taxo_filt_vis)
      
      df_filt2 <- df_filt %>%
        group_by(taxo_class) %>%
        summarise_all(funs(sum)) %>%
        as.data.frame()
      
      rownames(df_filt2) <- df_filt2[,1]
      df_filt2[,1] <- NULL
      df_filt2
    }
  })
  
  
  ## Show stacked bar table just for reference
  output$table_bar <- DT::renderDataTable({
    visual_data()
  })
  
  
  
  ## Render stacked bar graph
  output$stackedbar <- renderPlotly({
    
    visual_data_melted <- melt(cbind(visual_data(), ind=rownames(visual_data())), id.vars = c('ind'))
    
    print(
      ggplotly(
        ggplot(visual_data_melted,aes(x = variable, y = value, fill = ind)) + 
          geom_bar(position = "fill",stat = "identity") +
          scale_y_continuous(labels = percent_format())
      )
    )
    
  })
  
  }




# Run the application 
shinyApp(ui = ui, server = server)

