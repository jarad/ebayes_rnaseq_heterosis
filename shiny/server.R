library(shiny)
library(plyr)
library(ggplot2)

# d = read.table("rna_seq.txt", header=TRUE)

shinyServer(function(input,output) {
  
#   d_sub = reactive({
#     d[d$GeneID == input$geneid,]
#   })
#   
  output$data = renderDataTable({ data.frame(a=1) })
  
})
