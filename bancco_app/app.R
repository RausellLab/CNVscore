library(shiny)
library(DT)
library(tidyverse)
library(yardstick)
library(patchwork)
library(shinycssloaders)
library(glue)
library(rstanarm)
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


plot_test <- function(..., roc = TRUE) {
  
  
  if (isTRUE(roc)) {
    
    p1 <- bind_rows(...) %>% 
      mutate(tag = as.factor(tag)) %>%
      ggplot(aes(1-specificity, sensitivity)) +
      geom_path(aes(group = tag, color = tag),  show.legend = TRUE,  size = 2) +
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
    
    p2 <- bind_rows(...) %>%
      mutate(tag = as.factor(tag)) %>%
      ggplot(aes(recall, precision)) +
      geom_path(aes(group = tag, color = tag),  show.legend = TRUE,  size = 2) +
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




# ------------------------------------------------------------------------------
# WRITE CLINVAR TEST DATA
# ------------------------------------------------------------------------------
# 
# write_tsv(output_clinvar_deletion %>% select(chrom, start, end, type_variant, clinical), 'test_clinvar_del.tsv',col_names = FALSE)
# write_tsv(output_clinvar_duplication %>% select(chrom, start, end, type_variant, clinical), 'test_clinvar_dup.tsv',col_names = FALSE)

# ------------------------------------------------------------------------------
# READ MODELS
# ------------------------------------------------------------------------------

# setwd('/data-cbl/frequena_data/cnvscore/bancco_app')

# setwd('bancco_app')


# print(paste0(getwd()))

bancco_bayesian_clinvar_del_nohuman <- list()
bancco_bayesian_clinvar_del_human  <- list()
#
bancco_bayesian_clinvar_dup_nohuman  <- list()
bancco_bayesian_clinvar_dup_human  <- list()
#
bancco_logistic_clinvar_del_length  <- list()
bancco_logistic_clinvar_del_n_genes  <- list()
bancco_logistic_clinvar_del_omim  <- list()
#
bancco_logistic_clinvar_dup_length  <- list()
bancco_logistic_clinvar_dup_n_genes  <- list()
bancco_logistic_clinvar_dup_omim  <- list()
#
bancco_bayesian_decipher_del_nohuman <- list()
bancco_bayesian_decipher_dup_nohuman <- list()

for (i in 1:23) {
  print(i)
  bancco_bayesian_clinvar_del_nohuman[[i]] <- readRDS(glue('models/bayesian_clinvar_del_nohuman/bayesian_clinvar_del_nohuman_{i}.RData'))
  bancco_bayesian_clinvar_del_human[[i]] <- readRDS(glue('models/bayesian_clinvar_del_human/bayesian_clinvar_del_human_{i}.RData'))
  bancco_bayesian_clinvar_dup_nohuman[[i]] <- readRDS(glue('models/bayesian_clinvar_dup_nohuman/bayesian_clinvar_dup_nohuman_{i}.RData'))
  bancco_bayesian_clinvar_dup_human[[i]] <- readRDS(glue('models/bayesian_clinvar_dup_human/bayesian_clinvar_dup_human_{i}.RData'))
  bancco_logistic_clinvar_del_length[[i]] <- readRDS(glue('models/logistic_clinvar_del_length/logistic_clinvar_del_length_{i}.RData'))
  bancco_logistic_clinvar_del_n_genes[[i]] <- readRDS(glue('models/logistic_clinvar_del_n_genes/logistic_clinvar_del_n_genes_{i}.RData'))
  bancco_logistic_clinvar_del_omim[[i]] <- readRDS(glue('models/logistic_clinvar_del_omim/logistic_clinvar_del_omim_{i}.RData'))
  bancco_logistic_clinvar_dup_length[[i]] <- readRDS(glue('models/logistic_clinvar_dup_length/logistic_clinvar_dup_length_{i}.RData'))
  bancco_logistic_clinvar_dup_n_genes[[i]] <- readRDS(glue('models/logistic_clinvar_dup_n_genes/logistic_clinvar_dup_n_genes_{i}.RData'))
  bancco_logistic_clinvar_dup_omim[[i]] <- readRDS(glue('models/logistic_clinvar_dup_omim/logistic_clinvar_dup_omim_{i}.RData'))

  bancco_bayesian_decipher_del_nohuman[[i]] <- readRDS(glue('models/bayesian_decipher_del_nohuman/bayesian_decipher_del_nohuman_{i}.RData'))
  bancco_bayesian_decipher_dup_nohuman[[i]] <- readRDS(glue('models/bayesian_decipher_dup_nohuman/bayesian_decipher_dup_nohuman_{i}.RData'))
}

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

source('load_data.R') 

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
        
        test140 <<- user_input
        
        
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
        
        # if (nrow(reactive_annotated()[[1]]) > 0 & nrow(reactive_annotated()[[2]]) > 0) {
          
          pre_bayesian_clinvar_del_nohuman <- predict_chrom_aware_rtemis(bancco_bayesian_clinvar_del_nohuman, reactive_annotated()[[1]] , 'deletion', 'unbiased')
          pre_bayesian_clinvar_del_human <-  predict_chrom_aware_rtemis(bancco_bayesian_clinvar_del_human, reactive_annotated()[[1]] , 'deletion', 'knowledge-based')
          pre_bayesian_clinvar_dup_nohuman <-  predict_chrom_aware_rtemis(bancco_bayesian_clinvar_dup_nohuman, reactive_annotated()[[2]] , 'duplication', 'unbiased')
          pre_bayesian_clinvar_dup_human <-  predict_chrom_aware_rtemis(bancco_bayesian_clinvar_dup_human, reactive_annotated()[[2]] , 'duplication', 'knowledge-based')

          
          
          pre_bayesian_decipher_del_nohuman <- predict_chrom_aware_rtemis(bancco_bayesian_decipher_del_nohuman, reactive_annotated()[[1]] , 'deletion - decipher', 'unbiased')
          pre_bayesian_decipher_dup_nohuman <- predict_chrom_aware_rtemis(bancco_bayesian_decipher_dup_nohuman, reactive_annotated()[[2]] , 'duplication - decipher', 'unbiased')
          
          
          pre_logistic_clinvar_del_length <- predict_chrom_aware(bancco_logistic_clinvar_del_length, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic - clinical ~ length')
          pre_logistic_clinvar_del_n_genes <- predict_chrom_aware(bancco_logistic_clinvar_del_n_genes, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic -  clinical ~ nº genes')
          pre_logistic_clinvar_del_omim <- predict_chrom_aware(bancco_logistic_clinvar_del_omim, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic - clinical ~ OMIM genes')

          pre_logistic_clinvar_dup_length <- predict_chrom_aware(bancco_logistic_clinvar_dup_length, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ length')
          pre_logistic_clinvar_dup_n_genes <- predict_chrom_aware(bancco_logistic_clinvar_dup_n_genes, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ nº genes')
          pre_logistic_clinvar_dup_omim <- predict_chrom_aware(bancco_logistic_clinvar_dup_omim, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ OMIM genes')

          pre_bayesian_clinvar_del_nohuman[[3]] <- pre_bayesian_clinvar_del_nohuman[[3]] %>% select(-c(map_estimate, mad, chrom))
          pre_bayesian_clinvar_del_human[[3]] <- pre_bayesian_clinvar_del_human[[3]] %>% select(-c(map_estimate, mad, chrom))
          pre_bayesian_clinvar_dup_nohuman[[3]] <- pre_bayesian_clinvar_dup_nohuman[[3]] %>% select(-c(map_estimate, mad, chrom))
          pre_bayesian_clinvar_dup_human[[3]] <- pre_bayesian_clinvar_dup_human[[3]] %>% select(-c(map_estimate, mad, chrom))
          
          pre_bayesian_decipher_del_nohuman[[3]] <- pre_bayesian_decipher_del_nohuman[[3]] %>% select(-c(map_estimate, mad, chrom))
          pre_bayesian_decipher_dup_nohuman[[3]] <- pre_bayesian_decipher_dup_nohuman[[3]] %>% select(-c(map_estimate, mad, chrom))
          
          
          pre_logistic_clinvar_del_length[[3]] <- pre_logistic_clinvar_del_length[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
          pre_logistic_clinvar_del_n_genes[[3]] <- pre_logistic_clinvar_del_n_genes[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
          pre_logistic_clinvar_del_omim[[3]] <- pre_logistic_clinvar_del_omim[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
            
          pre_logistic_clinvar_dup_length[[3]] <- pre_logistic_clinvar_dup_length[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
          pre_logistic_clinvar_dup_n_genes[[3]] <- pre_logistic_clinvar_dup_n_genes[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
          pre_logistic_clinvar_dup_omim[[3]] <- pre_logistic_clinvar_dup_omim[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
          
          
          
          
          list(
          pre_bayesian_clinvar_del_nohuman,
          pre_bayesian_clinvar_del_human,
          pre_bayesian_clinvar_dup_nohuman,
          pre_bayesian_clinvar_dup_human,
          pre_logistic_clinvar_del_length ,
          pre_logistic_clinvar_del_n_genes,
          pre_logistic_clinvar_del_omim,
          pre_logistic_clinvar_dup_length,
          pre_logistic_clinvar_dup_n_genes,
          pre_logistic_clinvar_dup_omim,
          pre_bayesian_decipher_del_nohuman,
          pre_bayesian_decipher_dup_nohuman
          )
          

        

    })
    
    
    output$check_path <- renderText({
      
      paste(getwd())
      
    })
    

    output$roc_curve_plot <- renderPlot({
      
      test141 <<- predicted_user_input()
      



      p1_roc <- plot_test(
        predicted_user_input()[[1]][[1]],
        predicted_user_input()[[2]][[1]],
        predicted_user_input()[[3]][[1]],
        predicted_user_input()[[4]][[1]],
        predicted_user_input()[[5]][[1]],
        predicted_user_input()[[6]][[1]],
        predicted_user_input()[[7]][[1]],
        predicted_user_input()[[8]][[1]],
        predicted_user_input()[[9]][[1]],
        predicted_user_input()[[10]][[1]],
        predicted_user_input()[[11]][[1]],
        predicted_user_input()[[12]][[1]]
        
        )

      p1_pr <- plot_test(
        predicted_user_input()[[1]][[2]],
        predicted_user_input()[[2]][[2]],
        predicted_user_input()[[3]][[2]],
        predicted_user_input()[[4]][[2]],
        predicted_user_input()[[5]][[2]],
        predicted_user_input()[[6]][[2]],
        predicted_user_input()[[7]][[2]],
        predicted_user_input()[[8]][[2]],
        predicted_user_input()[[9]][[2]],
        predicted_user_input()[[10]][[2]],
        predicted_user_input()[[11]][[2]],
        predicted_user_input()[[12]][[2]],
        
        roc = FALSE)

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
        
    down_file2 <- reactive({
      
      test199 <<- predicted_user_input()
      
      # if (nrow(reactive_annotated()[[1]]) > 0 & nrow(reactive_annotated()[[2]]) > 0) {
        
        
        bind_rows(predicted_user_input()[[1]][[3]],
                  predicted_user_input()[[2]][[3]],
                  predicted_user_input()[[3]][[3]],
                  predicted_user_input()[[4]][[3]],
                  predicted_user_input()[[5]][[3]],
                  predicted_user_input()[[6]][[3]],
                  predicted_user_input()[[7]][[3]],
                  predicted_user_input()[[8]][[3]],
                  predicted_user_input()[[9]][[3]],
                  predicted_user_input()[[10]][[3]],
                  predicted_user_input()[[11]][[3]],
                  predicted_user_input()[[12]][[3]]
        )

      
    })
    
    output$download_file_2 <- downloadHandler(
      filename = 'file_results.tsv',
      content = function(file) {
        
        write_tsv(down_file2(), file, col_names = TRUE)
      }
    )

    
    
}
# Run the application 
shinyApp(ui = ui, server = server)

