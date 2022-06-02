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
          pre_bayesian_decipher_dup_nohuman <- predict_chrom_aware_rtemis(bancco_bayesian_decipher_dup_nohuman, reactive_annotated()[[1]] , 'duplication - decipher', 'unbiased')
          
          
          pre_logistic_clinvar_del_length <- predict_chrom_aware(bancco_logistic_clinvar_del_length, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic - clinical ~ length')
          pre_logistic_clinvar_del_n_genes <- predict_chrom_aware(bancco_logistic_clinvar_del_n_genes, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic -  clinical ~ nº genes')
          pre_logistic_clinvar_del_omim <- predict_chrom_aware(bancco_logistic_clinvar_del_omim, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic - clinical ~ OMIM genes')

          pre_logistic_clinvar_dup_length <- predict_chrom_aware(bancco_logistic_clinvar_dup_length, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ length')
          pre_logistic_clinvar_dup_n_genes <- predict_chrom_aware(bancco_logistic_clinvar_dup_n_genes, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ nº genes')
          pre_logistic_clinvar_dup_omim <- predict_chrom_aware(bancco_logistic_clinvar_dup_omim, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ OMIM genes')

          pre_bayesian_clinvar_del_nohuman[[3]] <- pre_bayesian_clinvar_del_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
          pre_bayesian_clinvar_del_human[[3]] <- pre_bayesian_clinvar_del_human[[3]] %>% select(-c(id, map_estimate, mad, chrom))
          pre_bayesian_clinvar_dup_nohuman[[3]] <- pre_bayesian_clinvar_dup_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
          pre_bayesian_clinvar_dup_human[[3]] <- pre_bayesian_clinvar_dup_human[[3]] %>% select(-c(id, map_estimate, mad, chrom))
          
          pre_bayesian_decipher_del_nohuman <- pre_bayesian_decipher_del_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
          pre_bayesian_decipher_dup_nohuman <- pre_bayesian_decipher_dup_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
          
          
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
          
          
    #     } else if (reactive_annotated()[[1]] > 0) {
    #       
    #       pre_bayesian_clinvar_del_nohuman <- predict_chrom_aware_rtemis(bancco_bayesian_clinvar_del_nohuman, reactive_annotated()[[1]] , 'deletion', 'unbiased')
    #       pre_bayesian_clinvar_del_human <-  predict_chrom_aware_rtemis(bancco_bayesian_clinvar_del_human, reactive_annotated()[[1]] , 'deletion', 'knowledge-based')
    #       
    #       pre_logistic_clinvar_del_length <- predict_chrom_aware(bancco_logistic_clinvar_del_length, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic - clinical ~ length')
    #       pre_logistic_clinvar_del_n_genes <- predict_chrom_aware(bancco_logistic_clinvar_del_n_genes, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic - clinical ~ nº genes')
    #       pre_logistic_clinvar_del_omim <- predict_chrom_aware(bancco_logistic_clinvar_del_omim, reactive_annotated()[[1]], is_bancco = TRUE, tag_bancco = 'Deletion - logistic - clinical ~ OMIM genes')
    #       
    #       
    #       pre_bayesian_clinvar_del_nohuman[[3]] <- pre_bayesian_clinvar_del_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #       pre_bayesian_clinvar_del_human[[3]] <- pre_bayesian_clinvar_del_human[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #       pre_bayesian_clinvar_dup_nohuman[[3]] <- pre_bayesian_clinvar_dup_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #       pre_bayesian_clinvar_dup_human[[3]] <- pre_bayesian_clinvar_dup_human[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #       
    #       pre_logistic_clinvar_del_length[[3]] <- pre_logistic_clinvar_del_length[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
    #       pre_logistic_clinvar_del_n_genes[[3]] <- pre_logistic_clinvar_del_n_genes[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
    #       pre_logistic_clinvar_del_omim[[3]] <- pre_logistic_clinvar_del_omim[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
    #       
    #       
    #       list(
    #         pre_bayesian_clinvar_del_nohuman,
    #         pre_bayesian_clinvar_del_human,
    #         pre_logistic_clinvar_del_length ,
    #         pre_logistic_clinvar_del_n_genes,
    #         pre_logistic_clinvar_del_omim
    #       )
    #   
    #     
    #       
    # } else if (reactive_annotated()[[2]] > 0) {
    #   
    #   pre_bayesian_clinvar_dup_nohuman <-  predict_chrom_aware_rtemis(bancco_bayesian_clinvar_dup_nohuman, reactive_annotated()[[2]] , 'duplication', 'unbiased')
    #   pre_bayesian_clinvar_dup_human <-  predict_chrom_aware_rtemis(bancco_bayesian_clinvar_dup_human, reactive_annotated()[[2]] , 'duplication', 'knowledge-based')
    #   
    #   pre_logistic_clinvar_dup_length <- predict_chrom_aware(bancco_logistic_clinvar_dup_length, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ length')
    #   pre_logistic_clinvar_dup_n_genes <- predict_chrom_aware(bancco_logistic_clinvar_dup_n_genes, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ nº genes')
    #   pre_logistic_clinvar_dup_omim <- predict_chrom_aware(bancco_logistic_clinvar_dup_omim, reactive_annotated()[[2]], is_bancco = TRUE, tag_bancco = 'Duplication - logistic - clinical ~ OMIM genes')
    # 
    #   
    #   pre_bayesian_clinvar_del_nohuman[[3]] <- pre_bayesian_clinvar_del_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #   pre_bayesian_clinvar_del_human[[3]] <- pre_bayesian_clinvar_del_human[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #   pre_bayesian_clinvar_dup_nohuman[[3]] <- pre_bayesian_clinvar_dup_nohuman[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #   pre_bayesian_clinvar_dup_human[[3]] <- pre_bayesian_clinvar_dup_human[[3]] %>% select(-c(id, map_estimate, mad, chrom))
    #   
    #   pre_logistic_clinvar_dup_length[[3]] <- pre_logistic_clinvar_dup_length[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
    #   pre_logistic_clinvar_dup_n_genes[[3]] <- pre_logistic_clinvar_dup_n_genes[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
    #   pre_logistic_clinvar_dup_omim[[3]] <- pre_logistic_clinvar_dup_omim[[3]] %>% select(.pred_pathogenic, clinical, tag) %>% mutate(sd = NA)
    #   
    #   
    #   list(
    #     pre_bayesian_clinvar_dup_nohuman,
    #     pre_bayesian_clinvar_dup_human,
    #     pre_logistic_clinvar_dup_length,
    #     pre_logistic_clinvar_dup_n_genes,
    #     pre_logistic_clinvar_dup_omim
    #   )
    # 
    # }
        

    })
    
    
    output$check_path <- renderText({
      
      paste(getwd())
      
    })
    

    output$roc_curve_plot <- renderPlot({
      
      test141 <<- predicted_user_input()
      

      # if (nrow(reactive_annotated()[[1]]) > 0 & nrow(reactive_annotated()[[2]]) > 0) {


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

      # } else if (reactive_annotated()[[1]] > 0) {
      # 
      #   p1_roc <- plot_test(
      #     predicted_user_input()[[1]][[1]],
      #     predicted_user_input()[[2]][[1]],
      #     predicted_user_input()[[3]][[1]],
      #     predicted_user_input()[[4]][[1]],
      #     predicted_user_input()[[5]][[1]]
      #     )
      # 
      #   p1_pr <- plot_test(
      #     predicted_user_input()[[1]][[2]],
      #     predicted_user_input()[[2]][[2]],
      #     predicted_user_input()[[3]][[2]],
      #     predicted_user_input()[[4]][[2]],
      #     predicted_user_input()[[5]][[2]],
      #     roc = FALSE)
      # 
      #   p1_roc + p1_pr
      # 
      # 
      # } else if (reactive_annotated()[[2]] > 0) {
      # 
      #   p1_roc <- plot_test(
      #     predicted_user_input()[[1]][[1]],
      #     predicted_user_input()[[2]][[1]],
      #     predicted_user_input()[[3]][[1]],
      #     predicted_user_input()[[4]][[1]],
      #     predicted_user_input()[[5]][[1]]
      #   )
      # 
      #   p1_pr <- plot_test(
      #     predicted_user_input()[[1]][[2]],
      #     predicted_user_input()[[2]][[2]],
      #     predicted_user_input()[[3]][[2]],
      #     predicted_user_input()[[4]][[2]],
      #     predicted_user_input()[[5]][[2]],
      #     roc = FALSE)
      # 
      # 
      #     p1_roc + p1_pr
      # }
      
  
     
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
        
        
        bind_rows(p1_roc %>% rename(p1 = specificity, p2 = sensitivity) %>% mutate(tag3 = 'roc'), 
                  p1_pr %>% rename(p1 = recall, p2 = precision) %>% mutate(tag3 = 'pr'))
        
        # p1_roc <- bind_rows(
        #   predicted_user_input()[[1]][[1]], 
        #   predicted_user_input()[[2]][[1]],
        #   predicted_user_input()[[3]][[1]], 
        #   predicted_user_input()[[4]][[1]],
        #   predicted_user_input()[[5]][[1]], 
        #   predicted_user_input()[[6]][[1]],
        #   predicted_user_input()[[7]][[1]], 
        #   predicted_user_input()[[8]][[1]]
        #   
        #   )
        # 
        # p1_pr <- bind_rows(
        #   predicted_user_input()[[1]][[2]], 
        #   predicted_user_input()[[2]][[2]],
        #   predicted_user_input()[[3]][[2]], 
        #   predicted_user_input()[[4]][[2]],
        #   predicted_user_input()[[5]][[2]], 
        #   predicted_user_input()[[6]][[2]],
        #   predicted_user_input()[[7]][[2]], 
        #   predicted_user_input()[[8]][[2]])
        # 
        # bind_rows(p1_roc %>% rename(p1 = specificity, p2 = sensitivity) %>% mutate(tag3 = 'roc'), 
        #           p1_pr %>% rename(p1 = recall, p2 = precision) %>% mutate(tag3 = 'pr'))
        # 
      # } else if (reactive_annotated()[[1]] > 0) {
      #   
      #   
      #   bind_rows(predicted_user_input()[[1]][[3]],
      #             predicted_user_input()[[2]][[3]],
      #             predicted_user_input()[[3]][[3]],
      #             predicted_user_input()[[4]][[3]]
      #   )
        
        # p1_roc <- bind_rows(
        #   predicted_user_input()[[1]][[1]],
        #   predicted_user_input()[[2]][[1]],
        #   predicted_user_input()[[3]][[1]],
        #   predicted_user_input()[[4]][[1]]
        # )
        # 
        # p1_pr <- bind_rows(
        #   predicted_user_input()[[1]][[2]],
        #   predicted_user_input()[[2]][[2]],
        #   predicted_user_input()[[3]][[2]],
        #   predicted_user_input()[[4]][[2]])
        # 
        # bind_rows(p1_roc %>% rename(p1 = specificity, p2 = sensitivity) %>% mutate(tag3 = 'roc'), 
        #           p1_pr %>% rename(p1 = recall, p2 = precision) %>% mutate(tag3 = 'pr'))
        # 
        
      # } else if (reactive_annotated()[[2]] > 0) {
      #   
      #   
      #   bind_rows(predicted_user_input()[[1]][[3]],
      #             predicted_user_input()[[2]][[3]],
      #             predicted_user_input()[[3]][[3]],
      #             predicted_user_input()[[4]][[3]]
      #     )
      #   
        # p1_roc <- bind_rows(
        #   predicted_user_input()[[1]][[1]],
        #   predicted_user_input()[[2]][[1]],
        #   predicted_user_input()[[3]][[1]],
        #   predicted_user_input()[[4]][[1]]
        # )
        # 
        # p1_pr <- bind_rows(
        #   predicted_user_input()[[1]][[2]],
        #   predicted_user_input()[[2]][[2]],
        #   predicted_user_input()[[3]][[2]],
        #   predicted_user_input()[[4]][[2]])
        
 
      # }
      
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


# 
# bannco_predict_chrom_aware <- function(list_models, df_to_predict, only_table = FALSE, is_bancco = TRUE, tag_bancco) {
# 
# 
#   # list_models <- bancco_rf_del_nohuman
#   # df_to_predict <- test913141
#   # is_bancco <- TRUE
#   # tag_bancco <- 'jeje'
# 
#   
#   if (isTRUE(is_bancco)) {
#     
#     tag_used <- tag_bancco
#     
#   }
#   
#   vector_chrom <- unlist(map(list_models, function(x) x[['chrom_target']]))
#   
#   for (i in 1:length(vector_chrom)) {
#     
#     tmp_chrom <- vector_chrom[i]
#     
#     test_split <- df_to_predict %>% filter(chrom == tmp_chrom)
#     tmp_model <- list_models[[i]][["model_trained"]]
#     
#     if (nrow(test_split) == 0) {
#       
#       tmp_predicted <- tibble()
#       
#     } else {
#       
#       tmp_predicted <- stats::predict(tmp_model, test_split, type = 'prob')
#       
#       tmp_predicted <- tmp_predicted %>%
#         bind_cols(test_split %>% select(clinical, id)) 
# 
#     }  
#     
#     if (i == 1) {
#       
#       prob_predicted <- tmp_predicted
#       
#     } else {
#       
#       prob_predicted <- prob_predicted %>% bind_rows(tmp_predicted)
#       
#     }
#   }
#   
#   prob_predicted <- prob_predicted %>% mutate(tag = tag_used)
#   
#   if (nrow(prob_predicted %>% filter(is.na(.pred_pathogenic))) > 0) {
#     
#     print(glue('There is {nrow(prob_predicted %>% filter(is.na(.pred_pathogenic)))} NA values'))
#     
#   }
#   
#   if (isTRUE(only_table)) return(prob_predicted)
#   
#   auc_result <- prob_predicted %>%
#     roc_auc(clinical, .pred_pathogenic)  %>%
#     pull(.estimate) %>%
#     round(3)
#   
#   pr_auc_result <- prob_predicted %>%
#     pr_auc(clinical, .pred_pathogenic)  %>%
#     pull(.estimate) %>%
#     round(3)
#   
#   roc_result <- prob_predicted %>%
#     roc_curve(clinical, .pred_pathogenic) %>%
#     mutate(tag = paste(tag_used, '-', auc_result))
#   
#   pr_result <- prob_predicted %>%
#     pr_curve(clinical, .pred_pathogenic) %>%
#     mutate(tag = paste(tag_used, '-', pr_auc_result))
#   
#   return(list(roc_result, pr_result,  prob_predicted))
#   
# }
# ------------------------------------------------------------------------------
# FUNCTION - PREDICT ACROSS CHROMOSOMES
# ------------------------------------------------------------------------------
# 
# predict_chrom_aware_rtemis <- function(list_models, test_split,
#                                        tag_variant = 'deletion',
#                                        tag_features, only_table = FALSE) {
#   # list_models <- rtemis_lasso
#   # test_split <- output_df_deletion_test
#   
#   # list_models <- bancco_bayesian_del_nohuman
#   # test_split <- output_gw_df
#   # tag_variant <- 'deletion'
#   # tag_features <- 'human_control'
#   
#   tag_used <- paste(tag_variant, 'bayesian', tag_features, sep = ' - ')
#   
#   vector_chrom <- unlist(map(list_models, function(x) x[['chrom_target']]))
#   
#   iter_vector <- 1:length(vector_chrom)
#   
#   tmp_predicted <- iter_vector %>%
#     map_dfr(~ predict_rtemis(list_models, .x, vector_chrom, test_split))
#   
#   if (isTRUE(only_table)) return(tmp_predicted)
#   
#   
#   ss_auc <- tmp_predicted %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
#   pr_auc <- tmp_predicted %>% pr_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
#   
#   
#   tmp_predicted <- tmp_predicted %>% mutate(tag = paste(tag_used, '-', ss_auc))
#   
#   tmp_roc_curve <- tmp_predicted %>% roc_curve(clinical, .pred_pathogenic)  %>% mutate(tag = paste(tag_used, '-', ss_auc))
#   tmp_pr_curve <- tmp_predicted %>% pr_curve(clinical, .pred_pathogenic)  %>% mutate(tag = paste(tag_used, '-', pr_auc))
#   
#   result <- list('tmp_roc_curve' = tmp_roc_curve, 'tmp_pr_curve' = tmp_pr_curve, 'tmp_predicted' = tmp_predicted)
#   
#   return(result)
#   
# }


# ------------------------------------------------------------------------------
# FUNCTION - PREDICT BAY_RTEMIS CHROMOSOMES
# ------------------------------------------------------------------------------
# 
# predict_rtemis <- function(x, chrom_tmp, vector_chrom, df_predict) {
#   
#   # x <- list_models
#   # chrom_tmp <- 15
#   # vector_chrom <- vector_chrom
#   # df_predict <- test_split
#   
#   tmp_rulefit <- x[[chrom_tmp]][["set_rules"]]
#   tmp_model <- x[[chrom_tmp]][["model_trained"]]
#   
#   df_predict <- df_predict %>% filter(chrom == vector_chrom[chrom_tmp])
#   # df_predict <- df_predict %>% filter(chrom == 22)
#   
#   
#   if (nrow(df_predict) == 0) {
#     
#     return(tibble())
#     
#   }
#   
#   testing_set_annotated <- generate_rtemis(df_predict, tmp_rulefit)
#   
#   tmp_posterior <- posterior_linpred(tmp_model, newdata = testing_set_annotated, transform = TRUE) %>%
#     as_tibble()
#   
#   prob_tmp <- tmp_posterior %>% 
#     map_dbl(~ median(.x)) %>%
#     enframe(name = NULL) %>%
#     rename(.pred_pathogenic = value) %>%
#     mutate(.pred_pathogenic = 1 - .pred_pathogenic)
#   
#   
#   sd_tmp <- tmp_posterior %>% map_dbl(~ sd(.x)) %>%
#     enframe(name = NULL) %>%
#     rename(sd = value)
#   
#   tmp_predicted <- prob_tmp %>% 
#     bind_cols(sd_tmp, df_predict %>% select(clinical, id)) %>%
#     mutate(chrom = vector_chrom[chrom_tmp])
#   
#   # print(glue('Chrom {vector_chrom[chrom_tmp]} DONE!'))
#   
#   return(tmp_predicted)
# }


# ------------------------------------------------------------------------------
# FUNCTION - GENERATE PARALLEL MATCHING RULES
# ------------------------------------------------------------------------------

# generate_rtemis <- function(newdata, rules_set) {
#   
#   # newdata <- df_predict
#   # rules_set <- tmp_rulefit
#   
#   matrix_zeros <- matrix(0, nrow(newdata), nrow(rules_set))
#   
#   newdata <- newdata %>% mutate(id = row_number())
#   
#   for (i in seq(nrow(rules_set))) {
#     
#     match <- newdata %>% filter_(rules_set$rule[i]) %>% pull(id)
#     matrix_zeros[match, i] <- 1
#   }
#   df_zeros <- matrix_zeros %>% as_tibble() %>% bind_cols(newdata %>% select(clinical))
#   return(df_zeros)
# }
# 
# # ------------------------------------------------------------------------------
# # CNVs ANNOTATION
# # ------------------------------------------------------------------------------
# 
# 
# check_cnv_v2 <- function(input_df, mode_reg = FALSE, factor_clinical = TRUE) {
#   
#   # input_df <- clinvar_match_deletion %>% filter(id_tmp == 561957)
#   
#   get_systems <- function(x) {
#     
#     hpo_from_gene <- hpo_genes %>% filter(gene %in% x) %>% pull(hp)
#     
#     result_n_systems <- unlist(map(hpo_from_gene, function(x) get_ancestors(hpo_dbs, x))) %>% 
#       enframe() %>%
#       filter(value %in% anato_df$name) %>%
#       count(value) %>% 
#       nrow()
#     
#     return(result_n_systems)
#     
#   }
#   
#   
#   # MODE_REG = ON
#   
#   if (isTRUE(mode_reg)) {
#     
#     mode_reg_enhancers <- bed_intersect(df_enhancers %>%
#                                           filter(phast100 >= 0.20 |
#                                                    phast46pla >= 0.20 |
#                                                    phast46pri >= 0.20), input_df) %>%
#       # filter(.overlap > threshold_30_tmp) %>%
#       rename(id_tmp = id_tmp.y, gene = gene.x) %>%
#       select(id_tmp, gene) %>% 
#       distinct()
#     
#     mode_reg_on <- mode_reg_enhancers
#   } else {
#     
#     mode_reg_on <- tibble()
#   }
#   
#   
#   
#   number_of_genes <- input_df %>%
#     bed_intersect(hgcn_genes %>% select(chrom, start, end, gene)) %>%
#     rename(id_tmp = id_tmp.x, gene = gene.y) %>%
#     select(id_tmp, gene) %>%
#     count(id_tmp) %>%
#     rename(n_genes = n)
#   
#   # Blacklist regions
#   result_df_blacklist <- bed_coverage(input_df, blacklist_encode) %>% rename(blacklist = .frac) %>% select(id_tmp, blacklist)
#   # Recombination rate
#   result_recomb <- bed_closest(input_df, recomb, suffix = c('', '.y')) %>% group_by(id_tmp) %>% 
#     filter(cm_mb.y == max(cm_mb.y)) %>% slice(1) %>% ungroup() %>% rename(recombination = cm_mb.y) %>% select(id_tmp, recombination)
#   # Centromeric distance
#   dist_cent <- bed_closest(input_df, region_gaps %>% filter(type == 'centromere')) %>% 
#     rename(id_tmp = id_tmp.x, dist_cent = .dist) %>% mutate(dist_cent = abs(dist_cent)/1e6) %>% select(id_tmp, dist_cent)
#   # Telomeric distance
#   dist_tel <- bed_closest(input_df, region_gaps %>% filter(type == 'telomere')) %>% rename(id_tmp = id_tmp.x, dist_tel = .dist) %>% 
#     group_by(id_tmp) %>% mutate(dist_tel = abs(min(dist_tel))/1e6) %>% ungroup() %>% select(id_tmp, dist_tel) %>% distinct()
#   
#   # SV hotspots
#   n_hotspot <- bed_intersect(input_df, hotspot, suffix = c('', '.y')) %>% select(id_tmp) %>% distinct() %>% mutate(hotspot = 1)
#   # HARs
#   n_hars <- bed_intersect(input_df, hars, suffix = c('', '.y')) %>% select(id_tmp) %>% distinct() %>% mutate(hars = 1)
#   # LADs
#   n_lads <- bed_intersect(input_df, lads, suffix = c('', '.y')) %>% select(id_tmp) %>% distinct() %>% mutate(lads = 1)
#   # Gene density
#   gene_density <- bed_intersect(gene_density_tbl, input_df) %>% group_by(id_tmp.y) %>%
#     filter(gene_density.x == max(gene_density.x)) %>% slice(1) %>% ungroup() %>% 
#     rename(id_tmp = id_tmp.y, gene_density = gene_density.x) %>% select(id_tmp, gene_density)
#   # PubMed - Deletions and duplications
#   pubmed_total <- bed_intersect(input_df, pubmed_df, suffix = c('', '.y')) %>% group_by(id_tmp) %>% 
#     filter(hits_del.y == max(hits_del.y)) %>% slice(1) %>% ungroup() %>% 
#     rename(hits_del = hits_del.y, hits_dup = hits_dup.y) %>% select(id_tmp, hits_del, hits_dup)
#   # Ensembl - CTCF
#   result_ctcf <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'CTCF_binding_site')) %>% select(id_tmp, .frac) %>%
#     rename(ctcf = .frac)
#   # Ensembl - Enhancers
#   result_enhancer <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'enhancer')) %>% select(id_tmp, .frac) %>%
#     rename(enhancer = .frac)
#   # Ensembl - Open chromatin
#   result_open <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'open_chromatin_region')) %>% select(id_tmp, .frac) %>%
#     rename(open = .frac)
#   # Ensembl - Promoter
#   result_promoter <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'promoter')) %>% select(id_tmp, .frac) %>%
#     rename(promoter = .frac)
#   # Ensembl - Promoter flank
#   result_promoterflank <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'promoter_flanking_region')) %>% 
#     select(id_tmp, .frac) %>% rename(promoterflank = .frac)
#   
#   # Ensembl - TFBS
#   result_tfbs <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'TF_binding_site')) %>% select(id_tmp, .frac) %>%
#     rename(tfbs = .frac)
#   # UCNE regions
#   result_ucne <- bed_coverage(input_df, ucne) %>% select(id_tmp, .frac) %>% rename(ucne = .frac)
#   result_ucne <- result_ucne %>% mutate(ucne = if_else(ucne == 0, 0, 1))
#   # Clinvar regions
#   result_df_clinvar <- input_df %>% 
#     bed_intersect(clinvar_variants %>% 
#                     filter(variant_class %in% c('indel', 'single nucleotide variant') & clinical == "pathogenic"),
#                   suffix = c('', 'y')) %>%
#     select(id_tmp) %>%
#     distinct() %>%
#     mutate(clinvar = 1)
#   
#   # GWAS variants
#   result_df_gwas <- input_df %>% bed_intersect(gwas_variants %>% filter(INTERGENIC == "No"), suffix = c('', 'y')) %>%
#     select(id_tmp) %>%
#     distinct() %>%
#     mutate(gwas = 1)
#   
#   # TADs
#   
#   result_tads <- input_df %>%
#     bed_intersect(tad %>%
#                     pivot_longer(-c(id, chrom), names_to = 'coord', values_to = 'start') %>%
#                     mutate(end = start)) %>%
#     count(id_tmp.x, id.y) %>%
#     filter(n == 1) %>%
#     select(-id.y) %>%
#     distinct() %>%
#     rename(id_tmp = id_tmp.x, tads = n)
#   
#   # CADD maximum
#   
#   result_cadd <- input_df %>%
#     bed_intersect(cadd_max %>% mutate(chrom = as.character(chrom))) %>%
#     select(id_tmp.x, max_cadd.y) %>%
#     rename(id_tmp = id_tmp.x, max_cadd = max_cadd.y) %>%
#     group_by(id_tmp) %>%
#     summarise(max_cadd = max(max_cadd))
#   
#   # GERP maximum
#   
#   
#   
#   result_gerp <- input_df %>%
#     bed_intersect(gerp_max) %>%
#     select(id_tmp.x, max_gerp.y) %>%
#     rename(id_tmp = id_tmp.x, max_gerp = max_gerp.y) %>%
#     group_by(id_tmp) %>%
#     summarise(max_gerp = max(max_gerp))
#   
#   # Remot-GW
#   
#   result_obs_exp <- input_df %>%
#     bed_intersect(remot_cnvscore) %>%
#     select(id_tmp.x, obs_exp.y) %>%
#     rename(id_tmp = id_tmp.x, obs_exp = obs_exp.y) %>%
#     group_by(id_tmp) %>%
#     summarise(max_obs_exp = max(obs_exp))
#   
#   region_level <- input_df %>%
#     left_join(result_df_blacklist, by = 'id_tmp') %>%
#     left_join(result_recomb, by = 'id_tmp') %>%
#     left_join(dist_cent, by = 'id_tmp') %>%
#     left_join(dist_tel, by = 'id_tmp') %>%
#     left_join(n_hotspot, by = 'id_tmp') %>%
#     left_join(n_hars, by = 'id_tmp') %>%
#     left_join(n_lads, by = 'id_tmp') %>%
#     left_join(gene_density, by = 'id_tmp') %>%
#     left_join(pubmed_total, by = 'id_tmp') %>%
#     left_join(result_ctcf, by = 'id_tmp') %>%
#     left_join(result_enhancer, by = 'id_tmp') %>%
#     left_join(result_open, by = 'id_tmp') %>%
#     left_join(result_promoter, by = 'id_tmp') %>%
#     left_join(result_promoterflank, by = 'id_tmp') %>%
#     left_join(result_tfbs, by = 'id_tmp') %>%
#     left_join(result_ucne, by = 'id_tmp') %>%
#     left_join(result_df_clinvar, by = 'id_tmp') %>%
#     left_join(result_df_gwas, by = 'id_tmp') %>%
#     left_join(result_tads, by = 'id_tmp') %>%
#     left_join(result_cadd, by = 'id_tmp') %>%
#     left_join(result_gerp, by = 'id_tmp') %>%
#     left_join(result_obs_exp, by = 'id_tmp') %>%
#     left_join(number_of_genes, by = 'id_tmp') %>%
#     # mutate(gene_density = ifelse(is.na(gene_density), 2, gene_density)) %>%
#     # mutate(dist_tel = ifelse(is.na(dist_tel), 27, dist_tel)) %>%
#     # mutate(tads = ifelse(is.na(tads), 0, tads)) %>%
#     # mutate(max_cadd = ifelse(is.na(max_cadd), 0, max_cadd)) %>%
#     # mutate(max_gerp = ifelse(is.na(max_gerp), 0, max_gerp)) %>%
#     # mutate(max_obs_exp = ifelse(is.na(max_obs_exp), 0, max_obs_exp)) %>%
#     # mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes)) %>%
#     # mutate(lads = ifelse(is.na(lads), 0, lads)) %>%
#     # mutate(hars = ifelse(is.na(hars), 0, hars)) %>%
#     # mutate(hotspot = ifelse(is.na(hotspot), 0, hotspot)) %>%
#     # mutate(clinvar = if_else(is.na(clinvar), 0, clinvar)) %>%
#   # mutate(gwas = if_else(is.na(gwas), 0, gwas)) %>%
#   mutate(across(everything(), ~replace_na(.x, 0))) %>%
#     select(-c('chrom', 'start', 'end', 'variant_class', 'source', 'clinical', 'length_cnv'))
#   
#   # GENE-LEVEL RESULTS
#   
#   gene_level <- input_df %>%
#     bed_intersect(hgcn_genes %>% select(chrom, start, end, gene)) %>%
#     rename(id_tmp = id_tmp.x, gene = gene.y) %>%
#     select(id_tmp, gene) %>%
#     bind_rows(mode_reg_on) %>%
#     distinct() %>%
#     mutate(disease = if_else(gene %in% (hgcn_genes %>% filter(disease == 'Yes') %>% pull(gene)), 1, 0)) %>% 
#     mutate(omim = if_else(gene %in% (hgcn_genes %>% filter(omim == 'Yes') %>% pull(gene)), 1, 0)) %>% 
#     mutate(haplo = if_else(gene %in% (haplo_triplo_genes %>% filter(haplo == 'yes') %>% pull(gene)), 1, 0)) %>% 
#     mutate(triplo = if_else(gene %in% (haplo_triplo_genes %>% filter(triplo == 'yes') %>% pull(gene)), 1, 0)) %>%
#     mutate(mouse_embryo = if_else(gene %in% (mgi %>% filter(str_detect(pheno, 'MP:0010768|MP:0005380')) %>% pull(gene)), 1, 0)) %>%
#     left_join(hgcn_genes %>% select(gene, pLI) %>% rename(pli = pLI), by = 'gene') %>%
#     left_join(loeuf_score, by = 'gene') %>%
#     left_join(eds, by = 'gene') %>%
#     left_join(pnull, by = 'gene') %>%
#     left_join(hgcn_genes %>% select(gene, ccr), by = 'gene') %>%
#     left_join(hgcn_genes %>% select(gene, hi), by = 'gene') %>%
#     left_join(expression_features, by = 'gene') %>%
#     left_join(genes_promoter, by = 'gene') %>%
#     left_join(crispr_score_df, by = 'gene') %>%
#     left_join(para_genes %>% rename(paralogous_genes = n), by = 'gene') %>%
#     left_join(string_db, by = 'gene') %>%
#     mutate(gene_hpo = if_else(gene %in% hpo_genes$gene, 1, 0)) %>% 
#     mutate(prot_complex = if_else(gene %in% prot_complex$gene, 1, 0)) %>% 
#     mutate(prot_complex_nohuman = if_else(gene %in% hu_map$gene, 1, 0)) %>% 
#     mutate(ohnolog = if_else(gene %in% hgcn_genes[hgcn_genes$ohnolog == 'Yes',]$gene, 1, 0)) %>%
#     mutate(essent_cl = if_else(gene %in% (hgcn_genes %>% filter(str_detect(fusil, 'CL')) %>% pull(gene)), 1, 0)) %>%
#     mutate(essent_dl = if_else(gene %in% (hgcn_genes %>% filter(str_detect(fusil, 'DL')) %>% pull(gene)), 1, 0)) %>%
#     mutate(tf = if_else(gene %in% tf_genes, 1, 0)) %>%
#     select(-gene) %>%
#     # mutate(disease = if_else(is.na(disease), 0, disease)) %>%
#     # mutate(haplo = if_else(is.na(haplo), 0, haplo)) %>%
#     # mutate(triplo = if_else(is.na(triplo), 0, triplo)) %>%
#     # mutate(mouse_embryo = if_else(is.na(mouse_embryo), 0, mouse_embryo)) %>%
#     # mutate(tf = ifelse(is.na(tf), 0, tf)) %>%
#     # mutate(gene_hpo = ifelse(is.na(gene_hpo), 0, gene_hpo)) %>%
#     # mutate(prot_complex = ifelse(is.na(prot_complex), 0, prot_complex)) %>%
#     # mutate(ohnolog = ifelse(is.na(ohnolog), 0, ohnolog)) %>%
#     # mutate(paralogous_genes = ifelse(is.na(paralogous_genes), 0, paralogous_genes)) %>%
#     # mutate(essent_cl = ifelse(is.na(essent_cl), 0, essent_cl)) %>%
#     # mutate(essent_dl = ifelse(is.na(essent_dl), 0, essent_dl)) %>%
#   # mutate(pli = ifelse(is.na(pli), 0, pli)) %>%
#   # mutate(loeuf = ifelse(is.na(loeuf), 0, loeuf)) %>%
#   # mutate(eds = ifelse(is.na(eds), 0, eds)) %>%
#   # mutate(pnull = ifelse(is.na(pnull), 0, pnull)) %>%
#   # mutate(ccr = ifelse(is.na(ccr), 0, ccr)) %>%
#   # mutate(hi = ifelse(is.na(hi), 0, hi)) %>%
#   # mutate(mean_phast = ifelse(is.na(mean_phast), 0, mean_phast)) %>%
#   # mutate(cpg_density = ifelse(is.na(cpg_density), 0, cpg_density)) %>%
#   # mutate(crispr_score = ifelse(is.na(crispr_score), 0, crispr_score)) %>%
#   # mutate(mean_expression = ifelse(is.na(mean_expression), 0, mean_expression)) %>%
#   # mutate(min_expression = ifelse(is.na(min_expression), 0, min_expression)) %>%
#   # mutate(degree_to_triplo = ifelse(is.na(degree_to_triplo), 0, degree_to_triplo)) %>%
#   # mutate(degree_to_haplo = ifelse(is.na(degree_to_haplo), 0, degree_to_haplo)) %>%
#   # mutate(degree = ifelse(is.na(degree), 0, degree)) %>%
#   # mutate(page_rank = ifelse(is.na(page_rank), 0, page_rank)) %>%
#   mutate(haplo_short = ifelse(is.na(haplo_short), 4, haplo_short)) %>%
#     mutate(triplo_short = ifelse(is.na(triplo_short), 4, triplo_short)) %>%
#     mutate(across(everything(), ~replace_na(.x, 0))) %>%
#     # Aggregation step
#     group_by(id_tmp) %>%
#     mutate(across(where(is.numeric) & !c(haplo_short, triplo_short), max)) %>%
#     mutate(across(c(haplo_short, triplo_short), min)) %>%
#     ungroup() %>%
#     distinct()
#   
#   
#   # 2 imputated features
#   # haplo_short -> 4 (median)
#   # triplo_short -> 4 (median)
#   # REMOVED dist_tel -> 27 (median)
#   # REMOVED gene_density -> 2 (median)
#   
#   cnvs_annotated <- input_df %>%
#     left_join(gene_level, by = 'id_tmp') %>%
#     left_join(region_level, by = 'id_tmp') %>%
#     mutate(haplo_short = ifelse(is.na(haplo_short), 4, haplo_short)) %>%
#     mutate(triplo_short = ifelse(is.na(triplo_short), 4, triplo_short)) %>%
#     mutate(across(everything(), ~replace_na(.x, 0))) %>%
#     mutate(clinical = as.factor(clinical)) %>%
#     rename(type_variant = variant_class, id = id_tmp)
#   
#   # mutate(clinical = as.factor(clinical)) %>%
#   # mutate(clinical = fct_relevel(clinical, 'pathogenic', 'benign'))
#   # map(~ sum(is.na(.x)))
#   
#   
#   if (isTRUE(factor_clinical)) {
#     
#     cnvs_annotated <- cnvs_annotated %>%
#       mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign')))
#     
#   }
#   
#   return(cnvs_annotated)
#   
# }