library(shiny)

d = read.csv("gene_names.csv")


shinyUI(fluidPage(

  titlePanel("Heterosis simulations"),
  
  sidebarLayout(
  
    sidebarPanel(
      selectInput('geneid', 'Gene ID', d$GeneID)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data", dataTableOutput(data) )
      )
    )
  )
))
