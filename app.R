#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)



# read.table("xxx.csv", sep = "\t", dec = ".", header = TRUE,
#            encoding = "UTF-8", stringsAsFactors = FALSE, quote = "")

# Define UI for application that draws a histogram
ui <- fluidPage(
  # App title
  titlePanel("Shiny Microbiome"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("countdf","Upload 'Count Data' file"), # fileinput() function is used to get the file upload contorl option
      fileInput("taxodf", "Upload 'Taxonomy Table' file"),
      fileInput("sampledf", "Upload 'Sample Data' file"),
      helpText("Select the read.table parameters below"),
      checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
      checkboxInput(inputId = "stringAsFactors", "stringAsFactors", TRUE),
      radioButtons(inputId = 'sep', label = 'Separator', choices = c(Tab='\t', Comma=',',Semicolon=';', Space=''), selected = '\t'),
      actionButton("checkbutton", label = "Check!!!")
      ),
    mainPanel(
      conditionalPanel(condition = "input.checkbutton == 0",
                       htmlOutput("help")),
      htmlOutput("checktxt")
    )
      
    ))


  

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  validation <- eventReactive(input$checkbutton, {
    if(nrow(input$countdf)=nrow(input$taxodf) & ncol(input$countdf)=nrow(input$sampledf)) {
      "The datasets are valid"} 
    else {
      "The datasets are invalid"}
  })
  
  output$help <- renderText({
    
    "<h3>To Check Validity of the Datasets:</h3>
    <ol>
    <li>Check the relevant parameter and separator setting</li>
    <li>Browse and upload the datasets</li>
    <li>Press Check!!!</li>
    </ol>"
    
  })
  
  output$checktxt <- renderText({
    validation()
  })
  
  # output$checkbutton <- renderPrint({   #renderText
  #   if(nrow(input$countdf)==nrow(input$taxodf) & ncol(input$countdf)==nrow(input$sampledf)) {
  #     "The datasets are valid"
  #   } else {
  #     "The datasets are invalid"}
  #   })
}

   
# Run the application 
shinyApp(ui = ui, server = server)

