


library(shiny)
library(DT)
library(tidyverse)
library(yardstick)
library(patchwork)
library(shinycssloaders)
library(glue)
# library(rstanarm)
library(rtemis)
library(valr)
library(parsnip)

# ------------------------------------------------------------------------------
# GENERAL FUNCTIONS
# ------------------------------------------------------------------------------

theme_roc <- function(major.breaks = c(0, .1, .25, .5, .75, .9, 1), 
                      minor.breaks = c(seq(0, .1, by = .01), seq(.9, 1, by = .01)), 
                      guide = TRUE, xlab = "False positive fraction", 
                      ylab = "True positive fraction", theme = theme_bw){
  
  
  res <- list(scale_x_continuous(xlab, breaks = major.breaks, minor_breaks = minor.breaks),
              scale_y_continuous(ylab, breaks = major.breaks, minor_breaks = minor.breaks), 
              theme())
  
  if(guide){
    
    pcol <- theme()$panel.grid.major$colour
    if(is.null(pcol)) pcol <- "white"
    res <- append(res, geom_abline(slope = 1, intercept = 0, color = pcol))
    
  }
  
  res
  
}


theme_pr <- function(major.breaks = c(0, .1, .25, .5, .75, .9, 1), 
                     minor.breaks = c(seq(0, .1, by = .01), seq(.9, 1, by = .01)), 
                     guide = TRUE, xlab = "Recall", 
                     ylab = "Precision", theme = theme_bw){
  
  
  res <- list(scale_x_continuous(xlab, breaks = major.breaks, minor_breaks = minor.breaks),
              scale_y_continuous(ylab, breaks = major.breaks, minor_breaks = minor.breaks), 
              theme())
  
  if(guide){
    
    pcol <- theme()$panel.grid.major$colour
    if(is.null(pcol)) pcol <- "white"
    res <- append(res, geom_abline(slope = 1, intercept = 0, color = pcol))
    
  }
  
  res
  
}

# ------------------------------------------------------------------------------
# PLOT TEST
# ------------------------------------------------------------------------------


plot_test <- function(x, roc = TRUE) {
  
  
  if (isTRUE(roc)) {
    
    p1 <- x %>% 
      mutate(model = as.factor(model)) %>%
      group_by(model) %>%
      roc_curve(clinical, .pred_pathogenic) %>%
      ggplot(aes(1-specificity, sensitivity)) +
      geom_path(aes(group = model, color = model),  show.legend = TRUE,  size = 2) +
      theme_roc() +
      theme(plot.title = element_text(size=20, face="bold"),
            axis.title.x = element_text(size=17, face="bold"),
            axis.title.y = element_text(size=17, face="bold"),
            axis.text.x = element_text(size=14, face="bold"),
            axis.text.y = element_text(size=14, face="bold"),
            legend.text = element_text(size=12)) +
      labs(title = 'ROC curve')
    
    p1
    
  } else {
    
    p2 <- x %>%
      mutate(model = as.factor(model)) %>%
      group_by(model) %>%
      pr_curve(clinical, .pred_pathogenic) %>%
      ggplot(aes(recall, precision)) +
      geom_path(aes(group = model, color = model),  show.legend = TRUE,  size = 2) +
      theme_pr() +
      theme(plot.title = element_text(size=20, face="bold"),
            axis.title.x = element_text(size=17, face="bold"),
            axis.title.y = element_text(size=17, face="bold"),
            axis.text.x = element_text(size=14, face="bold"),
            axis.text.y = element_text(size=14, face="bold"),
            legend.text = element_text(size=12)) +
      labs(title = 'Precision-Recall curve')
    
    
    p2
  }
}

# setwd('/data-cbl/frequena_data/cnvscore/bancco_app')

# setwd('bancco_app')

# ------------------------------------------------------------------------------
# EXAMPLE FILES
# ------------------------------------------------------------------------------

# write_tsv(output_clinvar_deletion %>% select(chrom, start, end, type_variant, clinical), 'test_clinvar_del.tsv',col_names = FALSE)
# write_tsv(output_clinvar_duplication %>% select(chrom, start, end, type_variant, clinical), 'test_clinvar_dup.tsv',col_names = FALSE)


# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------
# save(clinvar_ens_del3, file = 'bancco_app/models_ens/clinvar_ens_del3.RData')
# save(decipher_ens_del3, file = 'bancco_app/models_ens/decipher_ens_del3.RData')

source('load_data.R')
# 
# # Ensemble models
load('models_ens/logistic_del_length.RData')
load('models_ens/logistic_del_length.RData')

load('models_ens/clinvar_ens_del.RData')
load('models_ens/clinvar_ens_dup.RData')

load('models_ens/decipher_ens_del.RData')
load('models_ens/decipher_ens_dup.RData')

load('models_ens/logistic_del_length.RData')
load('models_ens/logistic_dup_length.RData')

load('models_ens/logistic_del_n_genes.RData')
load('models_ens/logistic_dup_n_genes.RData')

load('models_ens/logistic_del_omim.RData')
load('models_ens/logistic_dup_omim.RData')

load('models_ens/clinvar_ens_del3.RData')
load('models_ens/decipher_ens_del3.RData')


# ------------------------------------------------------------------------------
# APPLICATION
# ------------------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("CNVscore - Bancco dataset"),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(width = 2,
                 fileInput("file1", "Upload .tsv file here",
                           multiple = FALSE,
                           accept = c(".txt",
                                      ".tsv")),
                 downloadLink('download_file_1', label = "Download file example [#1]"),
                 downloadLink('download_file_2', label = "Download results file")
                 
                 
                 
                 # Button
                 # downloadButton("downloadData", "Download the Predictions")
    ),
    # Show the table with the predictions
    mainPanel(
      # textOutput('check_path')
      plotOutput("roc_curve_plot") %>% withSpinner(5)
      # DT::dataTableOutput("cnvs_annot") %>% withSpinner(5)
    )
  )
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  read_user_input <- reactive({
    
    req(input$file1)
    
    user_input <- read_tsv(input$file1$datapath, 
                           col_names = c('chrom', 'start', 'end', 'variant_class', 'clinical'),
                           col_types = list(chrom = col_character()))

    
    user_input <- user_input %>% 
      # mutate(clinical = as.factor(clinical)) %>%
      mutate(id_tmp = row_number()) %>%
      mutate(source = 'Bancco') %>%
      mutate(length_cnv = end - start + 1)
    
    
    test141 <<- user_input
    
  })
  
  reactive_annotated <- reactive({
    
    test123 <<- check_cnv_v2(read_user_input())
    
    tmp_df <- check_cnv_v2(read_user_input())
    
    list('df_del' = tmp_df %>% filter(type_variant == 'deletion'),
         'df_dup' = tmp_df %>% filter(type_variant == 'duplication'))
    
  })
  
  predicted_user_input <- reactive({
    
    test160 <<- reactive_annotated()
    
    tmp_del <- reactive_annotated()[[1]]
    tmp_dup <- reactive_annotated()[[2]]

    bind_rows(
    predict_ensemble(clinvar_ens_del, tmp_del) %>% mutate(model = 'ClinVar - DEL')  %>% rename(.pred_pathogenic = mean_score),
    predict_ensemble(clinvar_ens_dup, tmp_dup) %>% mutate(model = 'ClinVar - DUP')  %>% rename(.pred_pathogenic = mean_score),
    
    predict_ensemble(clinvar_ens_del3, tmp_del) %>% mutate(model = 'ClinVar - DEL - 1:3')  %>% rename(.pred_pathogenic = mean_score),
    predict_ensemble(decipher_ens_del3, tmp_del) %>% mutate(model = 'DECIPHER - DEL - 1:3')  %>% rename(.pred_pathogenic = mean_score),
    
    predict_ensemble(decipher_ens_del, tmp_del) %>% mutate(model = 'DECIPHER - DEL')  %>% rename(.pred_pathogenic = mean_score),
    predict_ensemble(decipher_ens_dup, tmp_dup) %>% mutate(model = 'DECIPHER - DUP')  %>% rename(.pred_pathogenic = mean_score),
    
    predict_chrom_aware(logistic_clinvar_del_length, tmp_del) %>% mutate(model = 'Naive - length - DEL'),
    predict_chrom_aware(logistic_clinvar_dup_length, tmp_dup) %>% mutate(model = 'Naive - length - DUP'),
    
    predict_chrom_aware(logistic_clinvar_del_n_genes, tmp_del) %>% mutate(model = 'Naive - n_genes - DEL'),
    predict_chrom_aware(logistic_clinvar_dup_n_genes, tmp_dup) %>% mutate(model = 'Naive - n_genes - DUP'),
    
    predict_chrom_aware(logistic_clinvar_del_omim, tmp_del) %>% mutate(model = 'Naive - omim - DEL'),
    predict_chrom_aware(logistic_clinvar_dup_omim, tmp_dup) %>% mutate(model = 'Naive - omim - DUP')
    )
    
  })
  
  
  output$check_path <- renderText({
    
    paste(getwd())
    
  })
  
  
  output$roc_curve_plot <- renderPlot({
    
    test141 <<- predicted_user_input()
    
  
    p1_roc <- plot_test(predicted_user_input())
    p1_pr <- plot_test(predicted_user_input(), roc = FALSE)

    
    p1_roc + p1_pr

  })
  
  output$cnvs_annot <- renderDataTable({
    
    reactive_annotated()[[1]]
    
  })
  
  
  down_file1 <- reactive({
    
    bind_rows(
      read_tsv('test_clinvar_del.tsv', col_names = FALSE),
      read_tsv('test_clinvar_dup.tsv', col_names = FALSE)
    )
    
  })
  
  
  output$download_file_1 <- downloadHandler(
    filename = 'file_example.tsv',
    content = function(file) {
      
      write_tsv(down_file1(), file, col_names = FALSE)
    }
  )
  

  
  output$download_file_2 <- downloadHandler(
    filename = 'file_results.tsv',
    content = function(file) {
      
      write_tsv(predicted_user_input(), file, col_names = TRUE)
    }
  )
  
  
  
}
# Run the application 
shinyApp(ui = ui, server = server)

