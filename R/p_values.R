
a <- bind_rows(input_check_cnv_del_benign, 
                               input_check_cnv_del_pathogenic_clinvar,
                               input_check_cnv_dup_pathogenic_clinvar,
                               input_check_cnv_del_pathogenic_decipher,
                               input_check_cnv_dup_pathogenic_decipher,
                               input_check_cnv_dup_benign) %>%
  select(source, length_cnv, variant_class)

a %>% 
  filter(source %in% c('audano_et_al', 'beyter_et_al', 'chaisson_et_al', 'gnomad_v2.1')) %>%
 pull(length_cnv) %>% median()

a %>% 
  filter(source == 'decipher') %>% pull(length_cnv) %>% median()

a %>% 
  filter(source == 'clinvar') %>% pull(length_cnv) %>% median()

a %>% 
  filter(source %in% c('dgv', 'decipher_control')) %>% pull(length_cnv) %>% median()


a %>% 
  filter(source %in% 
           c('audano_et_al', 'beyter_et_al', 
             'chaisson_et_al', 'gnomad_v2.1', 'decipher')) %>%
  mutate(source = if_else(source == 'decipher', 'pathogenic_decipher', 'benign_wgs')) %>%
  select(-variant_class) %>%
  wilcox_test(length_cnv ~ source)



just_test <- a %>% 
  filter(source %in% 
           c('audano_et_al', 'beyter_et_al', 
             'chaisson_et_al', 'gnomad_v2.1', 'decipher')) %>%
  mutate(source = if_else(source == 'decipher', 'pathogenic_decipher', 'benign_wgs')) %>%
  select(-variant_class)


just_test1 <- just_test[just_test$source == 'pathogenic_decipher',]$length_cnv
just_test2 <- just_test[just_test$source == 'benign_wgs',]$length_cnv

wilcox.test(just_test1, just_test2, alternative = 'greater')

a %>% 
  filter(source %in% 
           c('decipher_control', 'clinvar')) %>%
  mutate(source = if_else(source == 'clinvar', 'pathogenic_decipher', 'decipher_control')) %>%
  select(-variant_class) %>%
  wilcox_test(length_cnv ~ source, alternative = 'greater')

a %>% 
  filter(source %in% 
           c('dgv', 'clinvar')) %>%
  mutate(source = if_else(source == 'clinvar', 'pathogenic_decipher', 'dgv')) %>%
  select(-variant_class) %>%
  wilcox_test(length_cnv ~ source, alternative = 'greater')


#---------------------

cor_test1 <- result_clinvar_del %>%
  filter(tag %in% c('bayesian_unbiased', 'strvctvre', 'tada', 'xcnv')) %>%
  pivot_wider(id_cols = id, names_from = tag, values_from = .pred_pathogenic) %>%
  select(-id)
  

cor_test1 %>% cor_test(bayesian_unbiased, strvctvre, method = 'spearman')
cor_test1 %>% cor_test(bayesian_unbiased, tada, method = 'spearman')
cor_test1 %>% cor_test(bayesian_unbiased, xcnv, method = 'spearman')

cor.test(cor_test1$bayesian_unbiased, cor_test1$strvctvre, method = 'spearman')
cor.test(cor_test1$bayesian_unbiased, cor_test1$tada, method = 'spearman')
cor.test(cor_test1$bayesian_unbiased, cor_test1$xcnv, method = 'spearman')

#---------------------


result_clinvar_del %>%
  left_join(res_df %>% select(id, reliability_score, clinical), by = 'id') %>%
  filter(tag == 'bayesian_unbiased') %>%
  cor_test(.pred_pathogenic, reliability_score, method = 'spearman')


res_df_dup %>%
  cor_test(.pred_pathogenic, reliability_score, method = 'spearman')


#---------------------

result_clinvar_del %>%
  left_join(res_df %>% select(id, reliability_score, clinical), by = 'id') %>%
  group_by(tag) %>%
  cor_test(.pred_pathogenic, reliability_score, method = 'spearman') %>%
  arrange(desc(p))
  

  

