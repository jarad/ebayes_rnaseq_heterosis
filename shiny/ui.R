library(shiny)

d = read.csv("gene_names.csv", stringsAsFactors=FALSE)


shinyUI(fluidPage(

  titlePanel("Heterosis simulations"),
  
  sidebarLayout(
  
    sidebarPanel(
      selectizeInput('geneid', 'Gene ID', 1:40000)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data", dataTableOutput('data') )
      )
    )
  )
))
