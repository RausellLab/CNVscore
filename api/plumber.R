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
# api_bayesian_clinvar_dup_nohuman <- bayesian_clinvar_dup_nohuman

api_bayesian_clinvar_del_nohuman <- list()
api_bayesian_clinvar_dup_nohuman <- list()

for (i in 1:23) {
  api_bayesian_clinvar_del_nohuman[[i]] <- readRDS(glue('bayesian_clinvar_del_nohuman_{i}.RData'))
  api_bayesian_clinvar_dup_nohuman[[i]] <- readRDS(glue('bayesian_clinvar_dup_nohuman_{i}.RData'))

  }

source('load_data.R')


#* @filter checkAuth
function(req, res){
  if (is.null(req$input_chrom)){
    
    res$status <- 401 # Unauthorized
    return(list(error="Missing chromosome"))
    
  } else if (is.null(req$input_start)) {
    
    res$status <- 401 # Unauthorized
    return(list(error="Missing genomic interval - start"))
    
  } else if (is.null(req$input_end)) {
    
    res$status <- 401 # Unauthorized
    return(list(error="Missing genomic interval - end"))
    
  } else if (is.null(req$input_type)) {

    res$status <- 401 # Unauthorized
    return(list(error="Missing variant type"))
  
  } else {
    plumber::forward()
  }
}

#* @param input_type deletion or duplication
#* @param input_end End - genomic interval 
#* @param input_start Start - genomic interval
#* @param input_chrom Chromosome - genomic interval
#* @post /classifier
function(input_chrom, input_start, input_end, input_type){
  
    input_mod_chrom <- as.character(input_chrom)
    input_mod_start <- as.double(str_remove_all(input_start, ','))
    input_mod_end <- as.double(str_remove_all(input_end, ','))
    input_mod_type <- tolower(input_type)
    

  user_input <- tibble('chrom' = input_mod_chrom,
                       'start' = input_mod_start,
                       'end' = input_mod_end,
                       'variant_class' = input_mod_type,
                       'clinical' = 'benign',
                       'length_cnv' = 100,
                       'source' = 'provided_by_user') %>%
    mutate(id_tmp = row_number())
  
  # user_input <- tibble('chrom' = input_chrom, 'start' = input_start, 'end' = input_end, 'variant_class' = input_type) %>%
  #   mutate(chrom = as.character(chrom), start = as.integer(start), 
  #          end = as.integer(end), variant_class = tolower(variant_class)) %>%
  #   mutate(id_tmp = row_number(), source = 'provided_by_user', clinical = 'benign', length_cnv = 100)

  test1 <<- user_input
  
  tmp_df <- check_cnv_v2(user_input) 
  
  test2 <<- tmp_df
  
  if (input_mod_type == 'deletion') {
    
    tmp_predicted <- predict_chrom_aware_rtemis(api_bayesian_clinvar_del_nohuman,
                                                tmp_df, 'deletion', 'unbiased approach',
                                                only_table = TRUE)
    
    tmp_predicted <- get_reliability_score_mid(ref_quantiles, tmp_predicted) %>%
      rename(uncertainty_level = reliability_score)
    
  } else {
    
    tmp_predicted <- predict_chrom_aware_rtemis(api_bayesian_clinvar_dup_nohuman,
                                                tmp_df, 'duplication', 'unbiased approach',
                                                only_table = TRUE)
    
    tmp_predicted <- get_reliability_score_mid(ref_quantiles_dup, tmp_predicted) %>%
      rename(uncertainty_level = reliability_score)
    
  }
  

  
  
  test654 <<- tmp_predicted
  
  tmp_rules <- tmp_predicted %>%
    select(chrom, id, rules) %>%
    separate_rows(rules, sep = ';') %>%
    rename(rule_id = rules) %>%
    left_join(risk_support_del_clinvar_unbiased, by = c('rule_id', 'chrom')) %>%
    mutate(risk = round(risk, 2)) %>%
    filter(!str_detect(rule, 'max_obs_exp')) %>%
    filter(coefficient != 0) %>%
    mutate(sign_coefficient = ifelse(coefficient > 0, 'pos', 'neg')) %>%
    filter(sign_coefficient == ifelse(test654$.pred_pathogenic > 0.5, 'pos', 'neg')) %>%
    filter(risk > 0.8) %>%
    arrange(desc(coefficient)) %>%
    slice_head(n = 5) %>%
    mutate(rule_def = glue('{rule} (risk:{risk} - support:{support})')) %>%
    group_by(id) %>%
    summarise(rules = str_c(rule_def, collapse =";"))
  
  
  tmp_predicted <- tmp_predicted %>% select(-rules) %>% left_join(tmp_rules, by = 'id')
  

  tmp_predicted <- tmp_predicted %>% mutate(chrom = input_chrom,
                                            start = input_start,
                                            end = input_end,
                                            variant_class = input_mod_type)
  
  test16 <<- tmp_predicted
  
  
  tmp_predicted <- tmp_predicted %>% 
    rename(cnvscore = .pred_pathogenic) %>% 
    select(chrom, start, end, variant_class, cnvscore, uncertainty_level, rules)
  
  tmp_predicted

}



  
