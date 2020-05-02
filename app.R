library(shiny)
library(readxl)
library(magrittr)
library(tidyverse)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("airway")
# library(airway)
# BiocManager::install("DESeq2")
# library(DESeq2)

# data(airway)
# airway$dex %<>% relevel('untrt')
# 
# dds <- DESeqDataSet(airway, design = ~ cell + dex)
# dds <- DESeq(dds, betaPrior = FALSE)
# res1 <- results(dds, contrast = c('dex','trt','untrt'))
# res1 <- as.data.frame(lfcShrink(dds, contrast = c('dex','trt','untrt'), res = res1, type = 'normal'))
# res1 <- res1[, c(2,6)]
# res2 <- results(dds, contrast = c('cell', 'N061011', 'N61311'))
# res2 <- as.data.frame(lfcShrink(dds, contrast = c('cell', 'N061011', 'N61311'), res = res2, type = 'normal'))
# res2 <- res2[, c(2,6)]
# exp.data.frame <- cbind(res1, res2)
# rownames(exp.data.frame) <- NULL
# exp.data.frame$GeneName <- as.character(paste0("GeneName", 1:64102))
# exp.data.frame$Alt.GeneName <- as.character(paste0("Alt.GeneName", 1:64102))
# exp.data.frame <- exp.data.frame[,c(5,6,1:4)]
# colnames(exp.data.frame) <- c("GeneName", "Alt.GeneName", "Control_log2FC", "Control_pval.adj",
#                               "Treated_log2FC", "Treated_pval.adj")

exp.data.frame <- readRDS("exp.data.frame.RDS")
data <- exp.data.frame

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("MyDEAr - My Differential Expression Analysis in R"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Selector for choosing FC value from the uploaded dataset ----
      selectInput(inputId = "GeneID",
                  label = "Choose the column with gene id:",
                  choices = NULL),
      
      # Input: Selector for choosing FC value from the uploaded dataset ----
      selectInput(inputId = "FoldChange",
                  label = "Choose the column with fold change values:",
                  choices = NULL),
      
      # Input: Selector for choosing  p-values from the uploaded dataset ----
      selectInput(inputId = "p.value",
                  label = "Choose the column with p-value values:",
                  choices = NULL),
      
      # Horizontal line ----
      tags$hr(),
      
      numericInput(inputId = "UP",
                   label = "Fold change upregulation:",
                   value = 0.6,
                   step = 0.1),
      
      numericInput(inputId = "DOWN",
                   label = "Fold change downregulation:",
                   value = -0.6,
                   step = 0.1),
      
      numericInput(inputId = "pval",
                   label = "p-value:",
                   value = 0.05),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select number of rows to display ----
      radioButtons(inputId = "disp", 
                   label = "Arrange data based on:",
                   choices = c(`p-value` = "pvalue",
                               `Descending fold change` = "desc.fc",
                               `Ascending fold change` = "fc"),
                   selected = "pvalue"),
      
      # Horizontal line ----
      tags$hr()
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Content", tableOutput("contents")),
                  tabPanel("Volcano Plot", plotOutput("volcano_plot")))
      
    )
  )
)


server <- function(input, output, session) {
  
  info <- eventReactive(list(input$contents), {
    
    df <- na.omit(data)
    
    ref <- c()
    
    for (i in names(data)){
      name <- colnames(data[i])
      if (class(data[,i]) == "character") {
        ref <- c(ref, i) } 
      else {
        next()
      }}
    
    vars1 <- ref
    vars2 <- grep("FC", colnames(df), value = TRUE)
    vars3 <- grep("val", colnames(df), value = TRUE)
    
    # Update select input immediately after clicking on the action button.
    updateSelectInput(session, "GeneID", "Choose the column with gene id:", choices = vars1)
    updateSelectInput(session, "FoldChange", "Choose the column with fold change values:", choices = vars2)
    updateSelectInput(session, "p.value", "Choose the column with p-value values:", choices = vars3)
    
    df
    
  })
  
  output$contents <- renderTable({
    
    df <- info()
    
    df <- select(df, input$GeneID, input$FoldChange, input$p.value)
    
    # check if there is need to convert fold change data
    
    df$Threshold <- ifelse(test = df[ ,input$FoldChange] >= input$UP & 
                             df[ ,input$p.value] < input$pval, 
                           yes = "Upregulated", 
                           ifelse(df[ ,input$FoldChange] <= input$DOWN & 
                                    df[ ,input$p.value] < input$pval,
                                  yes = "Downregulated", no = "Not significant"))
    
    # try isolate number of observations to show
    
    if (input$disp == "pvalue") {
      return(df %>% 
               dplyr::arrange(!!rlang::sym(input$p.value)))
    } else if (input$disp == "desc.fc") {
      return(df %>% 
               dplyr::arrange(dplyr::desc(!!rlang::sym(input$FoldChange))))
    } else {
      return(df %>% 
               dplyr::arrange(!!rlang::sym(input$FoldChange)))
    }
    
  }, digits = 7)
  
  output$volcano_plot <- renderPlot({
    
    df <- info()
    
    df <- select(df, input$GeneID, input$FoldChange, input$p.value)
    
    df$Threshold <- ifelse(test = df[ ,input$FoldChange] >= input$UP & 
                             df[ ,input$p.value] < input$pval, 
                           yes = "Upregulated", 
                           ifelse(df[ ,input$FoldChange] <= input$DOWN & 
                                    df[ ,input$p.value] < input$pval,
                                  yes = "Downregulated", no = "Not significant"))
    
    vp <- ggplot(data = df,
                 mapping = aes_string(x = df[,input$FoldChange],
                                      y = -log10(df[,input$p.value]),
                                      color = "Threshold")) +
      geom_point(alpha = 0.3, size = 1) +
      scale_color_manual(values = c("dodgerblue", "gold", "deeppink2")) +
      ggtitle(paste0("Volcano plot showing \n ", input$FoldChange, " and ", input$p.value)) +
      labs(color = "Expression pattern") +
      xlab("log2FC") +
      ylab("-log10(p-value)") +
      theme_bw() +
      theme(axis.text = element_text(size = 14, face = "bold", color = "black"),
            axis.title = element_text(size = 14, face = "bold", color = "black"), 
            plot.title = element_text(size = 21, face = "bold", color = "black", hjust = 0.5),
            legend.title = element_text(size = 14, face = "bold", colour = "black"),
            legend.text = element_text(size = 14, face = "bold", colour = "black"),
            legend.position = "bottom") +
      guides(color = guide_legend(override.aes = list(size = 2)))
    
    print(vp)
    
  })
  
}


# Create Shiny app ----
shinyApp(ui, server)