# plumber.R

#* @apiTitle CNVscore API
#* @apiDescription CNVscore API


library(plumber)
library(tidyverse)
library(glue)
library(rtemis)
library(rstanarm)
library(valr)

setwd('/cnvscore')


# api_bayesian_clinvar_del_nohuman <- bayesian_clinvar_del_nohuman

api_bayesian_clinvar_del_nohuman <- list()
for (i in 1:23) {
  api_bayesian_clinvar_del_nohuman[[i]] <- readRDS(glue('bayesian_clinvar_del_nohuman_{i}.RData'))
}

source('load_data.R')

#* @param input_chrom Chromosome - genomic interval
#* @param input_start Start - genomic interval
#* @param input_end End - genomic interval 
#* @param input_type Deletion or Duplication
#* @post /classifier
function(input_chrom, input_start, input_end, input_type){
  

  user_input <- tibble('chrom' = as.character(input_chrom),
                       'start' = as.double(input_start),
                       'end' = as.double(input_end),
                       'variant_class' = tolower(input_type),
                       'clinical' = 'benign',
                       'length_cnv' = 100,
                       'source' = 'provided_by_user') %>%
    mutate(id_tmp = row_number())
  
  # user_input <- tibble('chrom' = input_chrom, 'start' = input_start, 'end' = input_end, 'variant_class' = input_type) %>%
  #   mutate(chrom = as.character(chrom), start = as.integer(start), 
  #          end = as.integer(end), variant_class = tolower(variant_class)) %>%
  #   mutate(id_tmp = row_number(), source = 'provided_by_user', clinical = 'benign', length_cnv = 100)

  test1 <<- user_input
  
  tmp_df <- check_cnv_v2(test1) 
  
  test2 <<- tmp_df
  
  tmp_predicted <- predict_chrom_aware_rtemis(api_bayesian_clinvar_del_nohuman,
                                              tmp_df, 'deletion', 'unbiased approach',
                             only_table = TRUE)
  
  test654 <<- tmp_predicted
  

  tmp_predicted <- tmp_predicted %>% mutate(chrom = input_chrom,
                                            start = input_start,
                                            end = input_end)
  
  test16 <<- tmp_predicted
  
  
  tmp_predicted

}


# model_chosen <- tibble(chrom = api_bayesian_clinvar_del_nohuman %>% map_chr(~ .x$chrom_target)) %>% 
#   mutate(id = row_number()) %>% filter(chrom == input_chrom) %>% pull(id)
# 
# 
# which_model_chosen <- api_bayesian_clinvar_del_nohuman[[model_chosen]][['set_rules']]
# 
# test4 <<- tmp_df
# test5 <<- which_model_chosen
# 
# testing_set_annotated <- generate_rtemis(tmp_df, which_model_chosen)
# 
# vector_chrom <- unlist(map(api_bayesian_clinvar_del_nohuman, function(x) x[['chrom_target']]))
# 
# 
# iter_vector <- 1:length(vector_chrom)
# 
# 
# tmp_predicted <- iter_vector %>%
#   map_dfr(~ predict_rtemis(api_bayesian_clinvar_del_nohuman, .x, vector_chrom, tmp_df))
# 
# rules_with_bay_coeff <- api_bayesian_clinvar_del_nohuman[[model_chosen]]$model_trained$coefficients %>% 
#   as_tibble(rownames = 'rule_id')
# 
# if (tmp_predicted$.pred_pathogenic > 0.5) {
#   
#   # IMPORTANT -> value < 0 
#   rules_with_bay_coeff <- rules_with_bay_coeff %>% filter(value < 0)
#   
# } else {
#   
#   rules_with_bay_coeff <- rules_with_bay_coeff %>% filter(value > 0)
#   
# }
# 
# test6 <<- rules_with_bay_coeff
# 
# rules_with_bay_coeff <- rules_with_bay_coeff %>% 
#   mutate(value = abs(value)) %>% 
#   filter(rule_id != '(Intercept)') %>% 
#   arrange(desc(value)) %>%
#   slice_head(n = 25)
# 
# # # IMPORTANT -> yes_or_no == 1
# met_rules <- testing_set_annotated %>% select(-clinical) %>%
#   pivot_longer(everything(), names_to = 'rule_id', values_to = 'yes_or_no') %>%
#   left_join(which_model_chosen %>% select(-coefficient), by = 'rule_id') %>%
#   inner_join(rules_with_bay_coeff, by = 'rule_id') %>%
#   filter(yes_or_no == 1)
# 
# tmp_predicted <- tmp_predicted %>% mutate(chrom = input_chrom, 
#                                           start = input_start, 
#                                           end = input_end) %>%
#   select(-clinical) %>%
#   mutate(rules = paste(met_rules$rule, collapse = ', '))
