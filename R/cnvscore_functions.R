# ------------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------------


report_split_cnvs <- function(x) {


from_hg19_to_hg38 = import.chain('/data-cbl/liftover/hg19ToHg38.over.chain')


granges_total_df <- x %>% select(id_tmp, chrom, start, end) %>% as.data.frame() %>% GRanges()

seqlevelsStyle(granges_total_df) = "UCSC"  # necessary

total_liftover = liftOver(granges_total_df, from_hg19_to_hg38)

good_after_liftover <- total_liftover %>% 
  as_tibble() %>%
  count(id_tmp) %>%
  filter(n == 1) %>%
  pull(id_tmp)

return(good_after_liftover)

}


# ------------------------------------------------------------------------------
# DIFFERENCE BETWEEN CODING MAPPING AND + NON-CODING TARGET REGIONS NO MAPPING
# ------------------------------------------------------------------------------

diff_coding_noncoding <- function(input_id, input_chrom, input_start, input_end, 
                                  input_rule_model = rulefit_model_deletion,
                                  input_model = bayesian_model_deletion) {
  
  
  # test_deletion_diff_coding_noncoding %>% filter(ref == 'more') %>% arrange(desc(dist)
  
  # tmp1 <- tibble(input_id = test_tbl_deletion %>% filter(id == 241024) %>% pull(id),
  #                input_clinical = 'pensionista',
  #                input_variant = 'deletion',
  #                input_chrom = test_tbl_deletion %>% filter(id == 241024) %>% pull(chrom),
  #                input_start = test_tbl_deletion %>% filter(id == 241024) %>% pull(start),
  #                input_end = test_tbl_deletion %>% filter(id == 241024) %>% pull(end))
  
  # input_model <- bayesian_model_deletion
  # input_rule_model <- rulefit_model_deletion
  
  tmp1 <- tibble(input_id = input_id,
                 input_clinical = 'pensionista',
                 input_variant = 'deletion',
                 input_chrom = input_chrom,
                 input_start = input_start,
                 input_end = input_end)
  
  
  
  prev_tmp1 <- check_cnv(tmp1$input_id, tmp1$input_clinical, tmp1$input_variant, tmp1$input_chrom,
                         tmp1$input_start, tmp1$input_end)
  
  prev_tmp2 <- check_cnv(tmp1$input_id, tmp1$input_clinical, tmp1$input_variant, tmp1$input_chrom,
                         tmp1$input_start, tmp1$input_end, mode_reg = TRUE)
  
  
  prev_tmp1 <- generate_rules(input_rule_model, prev_tmp1)
  prev_tmp1 <- prev_tmp1 %>% mutate_at(vars(!contains(c('max', 'dist', 'pubmed', 'n_system', 
                                                        'min_', 'density'))), as.factor)
  
  prev_tmp2 <- generate_rules(input_rule_model, prev_tmp2)
  prev_tmp2 <- prev_tmp2 %>% mutate_at(vars(!contains(c('max', 'dist', 'pubmed', 'n_system', 
                                                        'min_', 'density'))), as.factor)
  
  
  # a <- prev_tmp1 %>% select(where(is.integer) | where(is.double)) %>%
  #   pivot_longer(everything(), names_to = 'name', values_to = 'value1')
  #
  # b <- prev_tmp2 %>% select(where(is.integer) | where(is.double)) %>%
  #   pivot_longer(everything(), names_to = 'name', values_to = 'value2')
  #
  # a %>% left_join(b) %>% mutate(diff = value1 - value2) %>% filter(diff != 0) %>%
  #   left_join(coef(input_rule_model), by = c('name' = 'term')) %>% View()
  
  
  result_tmp1 <- posterior_epred(input_model, newdata = prev_tmp1) %>%
    as_tibble() %>% rename(no_reg = `1`)
  
  result_tmp2 <- posterior_epred(input_model, newdata = prev_tmp2) %>%
    as_tibble() %>% rename(reg = `1`)
  
  if (median(1 - result_tmp1$no_reg) == median(1 - result_tmp2$reg)) {
    
    result_ref <- 'equal'
  } else {
    
    result_ref <- if_else(median(1 - result_tmp1$no_reg) < median(1 - result_tmp2$reg), 'more', 'less')
  }
  
  # result_tmp1 %>%
  #   bind_cols(result_tmp2) %>%
  #   pivot_longer(everything(), names_to = 'type', values_to = 'value') %>%
  #   mutate(value = 1 - value) %>%
  #   ggplot(aes(value)) +
  #   geom_density(aes(fill = type), alpha = 0.3) +
  #   # scale_x_log10() +
  #   theme_minimal()
  
  x <- overlap(list('pre_reg' = result_tmp1$no_reg, 'post_reg' = result_tmp2$reg))
  
  return(tibble('id' = input_id,'dist' = 1 - x$OV, 'ref' = result_ref))
  
}



# ------------------------------------------------------------------------------
# PLOT DIFFERENCES CODING -> NON-CODING
# ------------------------------------------------------------------------------

plot_diff <- function(input_tbl) {
  
  # input_tbl <- result_deletion_diff_coding_noncoding
  
  tbl_tmp <-  input_tbl %>%
    mutate(more_0 = if_else(dist > 0, 'yes', 'no')) %>%
    group_by(more_0) %>%
    count(more_0, ref) %>%
    mutate(ref = factor(ref, levels = c('more', 'less'))) %>%
    mutate(perc = round(n / sum(n)*100, 2))
  
  total_n <- tbl_tmp$n %>% sum()
  
  p1_results <- input_tbl %>%
    ggplot(aes(dist)) +
    geom_histogram(aes(fill = ref), color = 'black') +
    theme_minimal()
  
  p2_results <- tbl_tmp %>%
    summarise(total = sum(n)) %>%
    mutate(perc = total / sum(total)) %>%
    mutate(more_0 = factor(more_0, levels = c('yes', 'no'))) %>%
    ggplot(aes('Total', perc)) +
    geom_col(aes(fill = more_0), color = 'black') +
    scale_y_continuous(label = percent) +
    geom_text(aes(label = paste0(100*round(perc, 2), '% ', '(', total,'/',total_n, ')' )),
              size = 5, position = position_stack(vjust = 0.5)) +
    theme_minimal() +
    labs(fill = 'Difference?')
  
  p3_results <- tbl_tmp %>%
    filter(more_0 == 'yes') %>%
    mutate(ref = factor(ref, levels = c('more', 'less'))) %>%
    ggplot(aes('Yes', perc)) +
    geom_col(aes(fill = ref), color = 'black') +
    geom_text(aes(label = paste0(round(perc, 2), '% ', '(', n,'/',total_n, ')' )),
              size = 5, position = position_stack(vjust = 0.5)) +
    theme_minimal()
  
  p1_results | (p2_results / p3_results)
}


# ------------------------------------------------------------------------------
# MATCHING CNVs BY LENGTH
# ------------------------------------------------------------------------------

matching_length <- function(bin_length = 100, tbl_input) {
  
  
  ranking_tech <- tibble(ranking = 1:6, source = c('beyter_et_al', 'audano_et_al', 'chaisson_et_al',
                                   'gnomad_v2.1', 'dgv', 'decipher_control'))
  
  # bin_length <- 100
  # tbl_input <- vus_decipher_before_match
  
  
  tbl_bins <- tibble('start' = seq(1, 1e7, bin_length), 'end' = seq(bin_length, 1e7, bin_length)) %>%
    mutate(chrom = 'chr1') %>%
    mutate(id = paste(start, end, sep = '-'))
  
  tmp_input_check_cnv <- tbl_input %>% mutate(chrom = 'chr1',
                                              start = length_cnv,
                                              end = length_cnv) %>%
    select(chrom, start, end, clinical, source, id_tmp)
  
  tmp_input_check_cnv <- tmp_input_check_cnv %>%
    bed_intersect(tbl_bins)
  
  # Select genomic intervals with control + decipher CNVs
  bins_with_decipher_and_control <- tmp_input_check_cnv %>%
    # count(id.y, source.x) %>%
    count(id.y, clinical.x) %>%
    count(id.y) %>%
    filter(n > 1) %>%
    pull(id.y)
  
  # Select ClinVar CNVs mapping the genomic intervals
  # previously selected
  selected_decipher_cnvs <- tmp_input_check_cnv %>%
    filter(clinical.x == 'pathogenic') %>%
    filter(id.y %in% bins_with_decipher_and_control) %>%
    pull(id_tmp.x)
  
  # only 35 intervals with > 1 decipher cnv
  number_decipher_cnvs <- tmp_input_check_cnv %>%
    filter(id.y %in% bins_with_decipher_and_control) %>%
    filter(clinical.x == 'pathogenic') %>%
    count(id.y) %>%
    rename(n_decipher = n)
  
  final_keep_tmp_ids <- c()
  
  for (i in 1:nrow(number_decipher_cnvs)) {
    
    print(glue('{i} / {nrow(number_decipher_cnvs)} genomic intervals'))
    
    tmp_control_id_tmp <- tmp_input_check_cnv %>%
      filter(id.y %in% bins_with_decipher_and_control) %>%
      filter(clinical.x != 'pathogenic') %>%
      filter(id.y == number_decipher_cnvs$id.y[i]) %>%
      left_join(ranking_tech, by = c('source.x' = 'source')) %>%
      # if we do not sample randomly, most benign CNVs are mapping chr1
      sample_n(size = nrow(.)) %>%
      arrange(ranking) %>%
      slice(1:number_decipher_cnvs$n_decipher[i]) %>%
      pull(id_tmp.x)
    
    # for (i in 1:)
    
    tmp_decipher_id_tmp <-  tmp_input_check_cnv %>%
      filter(clinical.x == 'pathogenic') %>%
      filter(id.y %in% number_decipher_cnvs$id.y[i]) %>%
      pull(id_tmp.x)
    
    # filter out extra DECIPHER CNVs
    if (number_decipher_cnvs$n_decipher[i] > length(tmp_control_id_tmp)) {
      
      n_extra <- length(tmp_control_id_tmp)
      
      tmp_decipher_id_tmp <- tmp_input_check_cnv %>%
        filter(clinical.x == 'pathogenic') %>%
        filter(id.y == number_decipher_cnvs$id.y[i]) %>%
        slice_sample(n = n_extra) %>%
        pull(id_tmp.x)
    }
    
    tmp_id_keeped <- c(tmp_control_id_tmp, tmp_decipher_id_tmp)
    
    final_keep_tmp_ids <- c(final_keep_tmp_ids, tmp_id_keeped)
    
  }
  
  tbl_input <- tbl_input %>%
    filter(id_tmp %in% final_keep_tmp_ids)
  
  return(tbl_input)
  
  
}

# ------------------------------------------------------------------------------
# COMPARISON Nº CNVS - CYTOBANDS LENGTH
# ------------------------------------------------------------------------------


cytolength_count <- function(cnv_tbl, tag = 'Deletion') {
  
  # cnv_tbl <- input_check_cnv_deletion
  
  
  cnv_tbl1 <- cnv_tbl %>% filter(clinical == 'benign')
  cnv_tbl2 <- cnv_tbl %>% filter(clinical != 'benign')
  
  
  tmp_distribution1 <- coord_cytobands %>%
    rowwise() %>%
    mutate(n_cnvs = valr::bed_intersect(tibble('chrom' = chrom, 'start' = start, 'end' = end), cnv_tbl1) %>% nrow()) %>%
    ungroup()
  
  
  tmp_distribution2 <- coord_cytobands %>%
    rowwise() %>%
    mutate(n_cnvs = valr::bed_intersect(tibble('chrom' = chrom, 'start' = start, 'end' = end), cnv_tbl2) %>% nrow()) %>%
    ungroup()
  
  p1 <- tmp_distribution1 %>%
    mutate(length_cytoband = end - start + 1) %>%
    filter(n_cnvs != 0) %>%
    select(Name, length_cytoband, n_cnvs) %>%
    mutate(ratio = n_cnvs / length_cytoband) %>%
    ggplot(aes(length_cytoband, n_cnvs)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    scale_x_log10() +
    stat_cor(label.y = 20) +
    theme_minimal() +
    labs(x = 'Log10(cytobands length)', y = 'CNVs count',
         title = glue('Nº CNVs per cytobands length - Benign - {tag}'))
  
  p2 <- tmp_distribution2 %>%
    mutate(length_cytoband = end - start + 1) %>%
    filter(n_cnvs != 0) %>%
    select(Name, length_cytoband, n_cnvs) %>%
    mutate(ratio = n_cnvs / length_cytoband) %>%
    ggplot(aes(length_cytoband, n_cnvs)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    scale_x_log10() +
    stat_cor(label.y = 20) +
    theme_minimal() +
    labs(x = 'Log10(cytobands length)', y = 'CNVs count',
         title = glue('Nº CNVs per cytobands length - Pathogenic - {tag}'))
  
  
  # cytoband with CNVs (yes/no)
  
  
  
  tmp_tmp_1 <- tmp_distribution1 %>%
    count(n_cnvs) %>%
    mutate(is_cnvs = if_else(n_cnvs == 0, 'no', 'yes')) %>%
    group_by(is_cnvs) %>%
    mutate(total = sum(n)) %>%
    select(is_cnvs, total) %>%
    ungroup() %>%
    slice(1:2) %>%
    mutate(perc = total / sum(total)) %>%
    mutate(type = 'Benign')
  
  
  
  tmp_tmp_2 <- tmp_distribution2 %>%
    count(n_cnvs) %>%
    mutate(is_cnvs = if_else(n_cnvs == 0, 'no', 'yes')) %>%
    group_by(is_cnvs) %>%
    mutate(total = sum(n)) %>%
    select(is_cnvs, total) %>%
    ungroup() %>%
    slice(1:2) %>%
    mutate(perc = total / sum(total)) %>%
    mutate(type = 'Pathogenic')
  
  p3 <- tmp_tmp_1 %>%
    bind_rows(tmp_tmp_2) %>%
    ggplot(aes(type, perc)) +
    geom_col(aes(fill = is_cnvs), color = 'black') +
    scale_y_continuous(label = percent) +
    theme_minimal() +
    labs(x = glue('{tag}'), y = 'Percentage cytobands', fill = 'Cytoband with CNV(s)')
  
  
  
  p5 <- tmp_distribution1 %>% 
    left_join(tmp_distribution2 %>% select(chrom, Name, n_cnvs), by = c('chrom', 'Name')) %>%
    mutate(type = case_when(
      n_cnvs.x > 0 & n_cnvs.y > 0 ~ 'both',
      n_cnvs.x > 0 & n_cnvs.y == 0 ~ 'benign',
      n_cnvs.x == 0 & n_cnvs.y > 0 ~ 'pathogenic',
      TRUE ~ 'none'
    )) %>%
    count(type) %>%
    mutate(perc = n / sum(n)) %>%
    mutate(type = factor(type, levels = c('pathogenic', 'none', 'both', 'benign'))) %>%
    ggplot(aes('Total', perc))+
    geom_col(aes(fill = type), color = 'black') +
    theme_minimal() +
    scale_y_continuous(label = percent) +
    geom_text(aes(label = paste0(100*round(perc, 2), '% ')),
              size = 5, position = position_stack(vjust = 0.5)) +
    labs(x = glue('{tag}'), y = 'Percentage cytobands', fill = 'Category')
  
  
  
  # Cumulative percentage
  
  cum_tmp1 <- tmp_distribution1 %>%
    select(Name, n_cnvs) %>%
    filter(n_cnvs != 0) %>%
    arrange(desc(n_cnvs)) %>%
    mutate(perc = n_cnvs / sum(n_cnvs)) %>%
    mutate(cum_perc = cumsum(perc)) %>%
    mutate(position = row_number()) %>%
    select(-Name) %>%
    mutate(type = 'Benign')
  
  cum_tmp2 <- tmp_distribution2 %>%
    select(Name, n_cnvs) %>%
    filter(n_cnvs != 0) %>%
    arrange(desc(n_cnvs)) %>%
    mutate(perc = n_cnvs / sum(n_cnvs)) %>%
    mutate(cum_perc = cumsum(perc)) %>%
    mutate(position = row_number()) %>%
    select(-Name) %>%
    mutate(type = 'Pathogenic')
  
  p4 <- cum_tmp1 %>%
    bind_rows(cum_tmp2) %>%
    mutate(type = factor(type, levels = c('Pathogenic', 'Benign'))) %>%
    ggplot(aes(position, cum_perc)) +
    geom_point(aes(fill = type), shape = 21) +
    scale_y_continuous(label = percent) +
    geom_path(aes(group = type, color = type), show.legend = FALSE) +
    theme_minimal() +
    labs(x = 'Rank cytobands', y = 'Cumulative percentage',
         fill = 'Clinical', title = glue('{tag}'))
  
  return(list(p1 + p2, p3, p4, p5))
  
}

# ------------------------------------------------------------------------------
# PLOT DISTRIBUTION ACROSS THE GENOME
# ------------------------------------------------------------------------------

plot_distribution_genome <- function(coord_cytobands, cnv_tbl, tag = 'Deletions') {
  
  # tag <- 'Deletions'
  # cnv_tbl <- input_check_cnv_deletion
  
  tmp_distribution <- coord_cytobands %>%
    rowwise() %>%
    mutate(n_cnvs = valr::bed_intersect(tibble('chrom' = chrom, 'start' = start, 'end' = end), cnv_tbl) %>% nrow()) %>%
    ungroup()
  
  
  # test_ratio <-  cnv_tbl %>%
  #    bed_intersect(coord_cytobands %>% mutate(Name = paste0(chrom, Name))) %>%
  #    group_by(id_tmp.x) %>%
  #    # there is not 2 cytobands 100% overlap
  #    filter(.overlap == max(.overlap)) %>%
  #    ungroup() %>%
  #    count(Name.y) %>%
  #    rename(n_cnvs = n, name_two = Name.y) %>%
  #    left_join(coord_cytobands %>% mutate(name_two = paste0(chrom, Name)),
  #              by = 'name_two') %>%
  #    mutate(length_cyto = (end - start + 1) / 1e6) %>%
  #    mutate(ratio = n_cnvs / length_cyto) %>%
  #    filter(n_cnvs > 0) %>%
  #    mutate(Name = paste0(chrom, Name)) 
  #  
  #  global_ratio <- (sum(test_ratio$n_cnvs))/(sum(test_ratio$length_cyto))
  # 
  # 
  #  test_ratio %>%
  #    mutate(chrom = factor(chrom, levels = as.character(seq(1,22)))) %>%
  #    ggplot(aes(reorder(Name, -ratio), ratio)) +
  #      geom_col() +
  #    geom_hline(yintercept = global_ratio ) +
  #    facet_wrap(vars(chrom), scales = 'free') +
  #    labs(y = 'Nº CNVs / Mb') +
  #    theme_minimal()
  
  
  
  
  
  p_distribution1 <- tmp_distribution %>%
    group_by(chrom) %>%
    summarise(total_cnvs = sum(n_cnvs), .groups = 'keep') %>%
    ggplot(aes(reorder(chrom, -total_cnvs), total_cnvs)) +
    geom_col() +
    theme_minimal() +
    xlab('Chromosome') +
    ylab('Count CNVs') +
    ggtitle(paste('Number of CNVs per chromosome -', tag))
  
  total_cnvs <- tmp_distribution %>% pull(n_cnvs) %>% sum()
  
  p_distribution2 <- tmp_distribution %>%
    mutate(perc = n_cnvs / total_cnvs) %>%
    group_by(chrom) %>%
    arrange(desc(perc)) %>%
    ungroup() %>%
    ggplot(aes(reorder(Name, -perc), perc)) +
    geom_col(aes(fill = chrom), color = 'black', show.legend = FALSE) +
    scale_y_continuous(label = percent) +
    facet_wrap(~ chrom, scales = 'free') +
    theme_minimal() +
    ggtitle(paste('Distribution of CNVs across genome -', tag)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  
  
  p_distribution3 <- tmp_distribution %>%
    mutate(perc = n_cnvs / total_cnvs) %>%
    group_by(chrom) %>%
    mutate(perc = n_cnvs / sum(n_cnvs)) %>%
    arrange(desc(perc)) %>%
    ungroup() %>%
    ggplot(aes(reorder(Name, -perc), perc)) +
    geom_col(aes(fill = chrom), color = 'black', show.legend = FALSE) +
    scale_y_continuous(label = percent) +
    facet_wrap(~ chrom, scales = 'free') +
    theme_minimal() +
    ggtitle(paste('Distribution of CNVs across genome (percentages by chromosome) -', tag)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  
  p_distribution4 <- tmp_distribution %>%
    ggplot(aes(reorder(Name, -n_cnvs), n_cnvs)) +
    geom_col(aes(fill = chrom), color = 'black', show.legend = FALSE) +
    facet_wrap(~ chrom, scales = 'free') +
    theme_minimal() +
    ggtitle(paste('Distribution of CNVs across genome (counts by chromosome) -', tag)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
  
  
  p_distribution5 <- cnv_tbl %>%
    mutate(clinical = if_else(clinical == 'benign', 'benign', 'pathogenic')) %>%
    count(clinical, chrom) %>%
    left_join(coord_chrom_hg19 %>% rename(length_chrom = length), by = 'chrom') %>%
    ggplot(aes(length_chrom, n)) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x) +
    geom_label(aes(label = chrom)) +
    ggtitle(paste('Nº CNVs and chromosome length -', tag)) +
    labs(x = 'Chromosome length', y = 'Number of CNVs - {tag}') +
    theme_minimal() +
    stat_cor(label.y = 40) +
    facet_wrap(~ clinical)
  
  
  tmp_plot <- cnv_tbl %>%
    mutate(clinical = if_else(clinical == 'benign', 'benign', 'pathogenic')) %>%
    count(clinical, chrom) %>%
    group_by(clinical) %>%
    mutate(perc = n / sum(n)) %>%
    ungroup() %>%
    mutate(chrom = fct_relevel(chrom, c(as.character(1:22), 'X')))
  
  mean_general <- tmp_plot  %>% pull(perc) %>% mean()
  
  p_distribution6 <- tmp_plot %>%
    ggplot(aes(chrom, perc)) +
    geom_col(aes(fill = clinical), color = 'black', position = 'dodge') +
    geom_hline(aes(yintercept = mean_general), linetype="dashed") +
    scale_y_continuous(label = percent) +
    labs(x = 'Chromosome', y = 'Proportion of CNVs',
         title = paste('Proportion of CNVs across chromosomes -', tag)) +
    theme_minimal()
  
  
  p_distribution7 <- cnv_tbl %>%
    filter(chrom != '17') %>%
    
    mutate(dist_cent = bed_closest(tibble('chrom' = chrom, 'start' = start, 'end' = end),
                                   region_gaps %>%
                                     filter(type == 'centromere')) %>% pull(.dist) %>% abs()) %>%
    mutate(clinical = if_else(clinical == 'benign', 'benign', 'pathogenic')) %>%
    mutate(dist_cent_cont = dist_cent / 1e6) %>%
    select(id_tmp, dist_cent_cont, clinical) %>%
    group_by(clinical) %>%
    mutate(dist_cent_disc = cut(dist_cent_cont, seq(1, 150, 1))) %>%
    count(dist_cent_disc) %>%
    mutate(n_tile = row_number()) %>%
    ggplot(aes(n_tile, n)) +
    geom_point(aes(fill = clinical), color = 'black', shape = 21) +
    facet_wrap(~ clinical) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    labs(y = 'Cnvs per 1Mbp', x = 'centromeric distance (Mbp)',
         title = paste('Nº CNVs and centromeric distance - ', tag))
  
  p_distribution8 <- cnv_tbl %>%
    filter(chrom != '17') %>%
    rowwise() %>%
    mutate(clinical = if_else(clinical == 'benign', 'benign', 'pathogenic')) %>%
    mutate(dist_cent = bed_closest(tibble('chrom' = chrom, 'start' = start, 'end' = end),
                                   region_gaps %>%
                                     filter(type == 'telomere')) %>% pull(.dist) %>% abs() %>% min()) %>%
    mutate(dist_cent_cont = dist_cent / 1e6) %>%
    select(id_tmp, dist_cent_cont, clinical) %>%
    group_by(clinical) %>%
    mutate(dist_cent_disc = cut(dist_cent_cont, seq(1, 150, 1))) %>%
    count(dist_cent_disc) %>%
    mutate(n_tile = row_number()) %>%
    ggplot(aes(n_tile, n)) +
    geom_point(aes(fill = clinical), color = 'black', shape = 21) +
    facet_wrap(~ clinical) +
    geom_smooth(method='lm', formula= y~x) +
    theme_minimal() +
    labs(y = 'Cnvs per 1Mbp', x = 'Telomere distance (Mbp)',
         title = paste('Nº CNVs and telomere distance - ', tag))
  
  
  
  p9_tmp <-   cnv_tbl %>%
    mutate(clinical = if_else(clinical == 'benign', 'benign', 'pathogenic')) %>%
    count(clinical, chrom) %>%
    left_join(coord_chrom_hg19 %>% rename(length_chrom = length), by = 'chrom') %>%
    mutate(length_chrom = length_chrom / 1e6) %>%
    mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign')))
  
  
  average_ratio_pathogenic <- p9_tmp %>% 
    filter(clinical == 'pathogenic') %>% summarise(sum_n = sum(n), sum_length = sum(length_chrom)) %>%
    mutate(ratio = sum_n / sum_length) %>% pull(ratio)
  
  average_ratio_benign <- p9_tmp %>% 
    filter(clinical != 'pathogenic') %>% summarise(sum_n = sum(n), sum_length = sum(length_chrom)) %>%
    mutate(ratio = sum_n / sum_length) %>% pull(ratio)
  
  p_distribution9 <- p9_tmp %>%
    mutate(ratio = (n / length_chrom)) %>% 
    ggplot(aes(reorder(chrom, -ratio), ratio)) + 
    geom_col(color = 'black', aes(fill = clinical), position = 'dodge') +
    geom_hline(yintercept = average_ratio_pathogenic, linetype = 'dashed', color = 'red') +
    geom_hline(yintercept = average_ratio_benign, linetype = 'dashed', color = 'steelblue') +
    theme_minimal() +
    labs(y = "Nº CNVs / chromosome length (Mbs)", x = 'Chromosome', 
         title = paste0('Nº CNVs / chromosome length (Mbs) - ', tag))
  
  
  return(list(p_distribution1,
              p_distribution2,
              p_distribution3,
              p_distribution4,
              p_distribution5,
              p_distribution6,
              p_distribution7,
              p_distribution8,
              p_distribution9))
}


# ------------------------------------------------------------------------------
# RULEFIT - EXPLORATION
# ------------------------------------------------------------------------------

# SUMMARY CV
# print(rulefit_model_deletion$glm$model)
# # Df -> nº of nonzero coefficients
# print(rulefit_model_deletion$glm$model$glmnet.fit)
#
# # CROSS-VALIDATION GRID LAMBDA PLOT
# plot(rulefit_model_deletion$glm$model)


rules_rtemis <- function(rule_model, bayesian_model, training_tbl) {
  
  
  # rule_model <- which_model_chosen$set_rules
  # bayesian_model <- which_model_chosen$model_trained
  # training_tbl <- which_database_chosen
  
  binary_rules <- generate_rtemis(training_tbl, rule_model)
  
  count_pathogenic <- binary_rules %>% 
    filter(clinical == 'pathogenic') %>% 
    select(-clinical) %>% 
    # removed a *2
    map_dfc(~ sum(.x)) %>%
    pivot_longer(everything(), names_to = 'rule', values_to = 'count_patho')
  
  count_total <- binary_rules %>% 
    select(-clinical) %>% 
    # removed a *2
    map_dfc(~ sum(.x)) %>%
    pivot_longer(everything(), names_to = 'rule', values_to = 'count_total')
  
  count_result <- count_pathogenic %>% 
    left_join(count_total, by = 'rule') %>%
    mutate(risk = count_patho / count_total) %>%
    select(rule, risk)

  
  support_rules <- binary_rules %>% 
    select(-clinical) %>% 
    map_dfc(~ sum(.x)) %>%
    pivot_longer(everything(), names_to = 'rule', values_to = 'support')
  
  
  output_df <- bayesian_model$coefficients %>% enframe(name = 'rule') %>% 
    filter(rule != '(Intercept)') %>%
    left_join(count_result, by = 'rule') %>%
    left_join(support_rules, by = 'rule') %>%
    left_join(rule_model %>% select(-coefficient) %>% rename(description = rule), by = c('rule' = 'rule_id')) %>%
    rename(estimate = value)
  
  return(output_df)
}


plot_rules <- function(rule_model, bayesian_model, training_tbl) {
  
  # rule_model <- rulefit_model_deletion plot_rules(rulefit_model_deletion, bayesian_model_deletion)
  
  
  
  rule_model <- remove_numbers(rule_model)
  
  
  rules_tbl <- rule_model %>%
    select(-coefficient_lambda.min)
  
  bay_tbl <- bayesian_model$coefficients %>% enframe() %>%
    mutate(name = str_sub(name, end=-2))
  
  bay_tbl$name[1] <- '(Intercept)'
  
  bay_tbl <- bay_tbl %>%
    rename(term = name) %>%
    rename(coefficient = value) %>%
    left_join(rules_tbl, by = 'term') %>%
    filter(str_detect(term, '^[r]'))
  
  bay_tbl <- bay_tbl %>%
    rowwise() %>%
    mutate(tmp_column = get_rules_information(training_tbl, rule)) %>%
    separate(col = tmp_column, into = c('support', 'risk'), sep = ' ') %>%
    mutate(support = as.integer(support), risk = as.numeric(risk)) %>%
    relocate(term, rule, coefficient, support, risk) %>%
    arrange(desc(risk))
  
  p1_plot <- bay_tbl %>%
    mutate(color_support = if_else(support <= 5, 'yes', 'no')) %>%
    ggplot(aes(coefficient, risk)) +
    geom_point(shape = 21, color = 'black', aes(fill = support), size = 4) +
    scale_fill_viridis_c() +
    geom_hline(yintercept = 0.5, linetype="dashed") +
    geom_vline(xintercept = 0, linetype="dashed") +
    theme_minimal()
  
  p2_plot <- bay_tbl %>%
    ggplot(aes(support, risk)) +
    geom_point() +
    theme_minimal()
  
  p3_plot <- bay_tbl %>%
    select(rule) %>%
    separate_rows(rule, sep = ' & ') %>%
    mutate(rule = str_remove(rule, '\\<.*')) %>%
    mutate(rule = str_remove(rule, '\\>.*')) %>%
    count(rule) %>%
    na.omit() %>%
    ggplot(aes(reorder(rule, n), n)) +
    geom_col(aes(fill = n), color = 'black') +
    coord_flip() +
    scale_fill_viridis_c() +
    theme_minimal()
  
  return(list(p1_plot + p2_plot + p3_plot, rules_tbl))
  
  
}

# ------------------------------------------------------------------------------
# EXPLORE RANGE THRESHOLDS STANDARD DEVIATION (SD)
# ------------------------------------------------------------------------------

range_thresholds <- function(input_tbl, tag = NULL) {
  
  # input_tbl <- enriched_test_deletion
  
  max_sd <- input_tbl %>% arrange(desc(sd)) %>% slice(1) %>% pull(sd)
  range_sd <- seq(0, max_sd - 0.01, 0.01)
  
  result_dbl_sd <- c()
  
  for (i in 1:length(range_sd)) {
    # print(range_sd[i])
    
    tmp_result <- input_tbl %>%
      filter(sd <= range_sd[i]) %>%
      roc_auc(clinical, pred_target) %>%
      pull(.estimate)
    
    result_dbl_sd <- c(tmp_result, result_dbl_sd)
    
  }
  
  
  tibble('auc' = result_dbl_sd, 'sd_threshold' = range_sd) %>%
    ggplot(aes(range_sd, auc)) +
    geom_line(size = 2) +
    theme_minimal() +
    labs(xlab = 'Uncertainty (sd) threshold', title = glue('{tag}'))
  
}

# ------------------------------------------------------------------------------
# RETAINED DATA
# ------------------------------------------------------------------------------

retained_data <- function(input_tbl, tag = NULL) {
  
  # input_tbl <- enriched_test_deletion
  range_retained <- seq(0, 1, 0.01)
  
  result_dbl_sd <- c()
  result_dbl_mad <- c()
  result_dbl_random <- c()
  
  
  for (i in 1:length(range_retained)) {
    
    tmp_result <- input_tbl %>%
      arrange(sd) %>%
      slice_head(prop = range_retained[i]) %>%
      roc_auc(clinical, pred_target)  %>%
      pull(.estimate)
    
    tmp_result2 <- input_tbl %>%
      arrange(mad) %>%
      slice_head(prop = range_retained[i]) %>%
      roc_auc(clinical, pred_target)  %>%
      pull(.estimate)
    
    tmp_result3 <- input_tbl %>%
      slice_sample(prop = range_retained[i]) %>%
      roc_auc(clinical, pred_target)  %>%
      pull(.estimate)
    
    
    result_dbl_sd <- c(result_dbl_sd, tmp_result)
    result_dbl_mad <- c(result_dbl_mad, tmp_result2)
    result_dbl_random <- c(result_dbl_random, tmp_result3)
  }
  
  
  tibble('auc_sd' = result_dbl_sd,
         'auc_mad' = result_dbl_mad,
         'auc_random' = result_dbl_random,
         'data_retained' = range_retained) %>%
    pivot_longer(-data_retained, names_to = 'metric', values_to = 'value') %>%
    separate(metric, c('remove', 'metric'), sep = '_') %>%
    select(-remove) %>%
    ggplot(aes(data_retained, value)) +
    geom_line(size = 2, aes(color = metric)) +
    theme_minimal() +
    scale_x_continuous(label = percent) +
    labs(xlab = 'Uncertainty (sd) threshold', title = glue('{tag}'))
  
}


# ------------------------------------------------------------------------------
# RULEFIT - UNCERTAINTY
# laplacian problem -> in order to have a lot of zeros in the signal,
# you are also forcing the non-zero elements to be very small
# BAYESIAN APPROACH (laplacian -> horsehoe -> finnish horseshoe)
# ------------------------------------------------------------------------------
# plot(bayesian_model, "areas", prob = 0.95, prob_outer = 1)

library(rstan)
library(rstanarm)

# nº iterations per chain (iter)
# nº chains (chains)
# algorithm (algorithm = "sampling", "optimizing", "meanfield", "fullrank")
# QR
# posterior distribution -> point estimate (mean, median)
#
#
#
#
#
# bayesian_model <- rstanarm::stan_glm(formula = formule_models,
#                                        family = 'binomial',
#                                        data = rulefit_model_deletion$full_data,
#                                        cores = 4,
#                                        iter = 2000,
#                                        chains = 4,
#                                        algorithm = 'sampling', # variational inference algorithms
#                                        QR = TRUE,
#                                        prior = hs(),
#                                        prior_intercept = hs())
# 
# library(bayesplot)
# print(bayesian_model)
# summary(bayesian_model)
# 
# as_tibble(bayesian_model_duplication)[2:4] %>%
#   
#   mcmc_areas(prob = 0.8)



# ------------------------------------------------------------------------------
# GET APPLICABILITY RESULTS
# https://www.marlycormar.com/presentations/R-Pharma-2019/presentation.html#1
# ------------------------------------------------------------------------------

get_distances <- function(training_tbl, formule_input, input_tbl, threshold = 0.25, only_origin = FALSE) {

# 
  # training_tbl <- output_clinvar_deletion
  # input_tbl <- output_decipher_deletion
  # formule_input <- human_no_control
  # threshold <- 0.95
  
  
  tmp_result <- input_tbl %>% select(id)
  
  input_tbl <-  input_tbl %>%
    select(-c(type_variant, source, clinical, id)) %>%
    mutate_if(is.factor, as.character) %>%
    mutate_if(is.character, as.double)
  
  # Pathogenic CNVs
  
  recipe_applicable <-
    recipe( ~ ., data = training_tbl %>% filter(clinical == 'pathogenic') %>%
              select(any_of(all.vars(formule_input)[-1]))) %>%
    step_dummy(all_nominal()) %>%
    # Remove variables that have the same value for every data point.
    step_zv(all_predictors()) %>%
    # Transform variables to be distributed as Gaussian-like as possible.
    step_YeoJohnson(all_numeric()) %>%
    # Normalize numeric data to have a mean of zero and
    # standard deviation of one.
    step_normalize(all_numeric())
  
  applicability_pca <- apd_pca(recipe_applicable, data = training_tbl %>% 
                                 filter(clinical == 'pathogenic'), threshold = threshold)
  
  origin_pathogenic <- applicable::score(applicability_pca, training_tbl %>% 
                                           filter(clinical == 'pathogenic'))
  
  origin_pathogenic <- origin_pathogenic %>% 
    rename(pc1_origin = PC01, pc2_origin = PC02) %>%
    select(pc1_origin, pc2_origin)
  
  
  tmp_distance_pathogenic <- applicable::score(applicability_pca, input_tbl) %>%
    rename(dist_pathogenic = distance_pctl, pc1_pathogenic = PC01, pc2_pathogenic = PC02) %>%
    select(dist_pathogenic, pc1_pathogenic, pc2_pathogenic)
    
  
  ## Benign CNVs
  
  recipe_applicable <-
    recipe( ~ ., data = training_tbl %>% filter(clinical == 'benign') %>%
              select(any_of(all.vars(formule_input)[-1]))) %>%
    step_dummy(all_nominal()) %>%
    # Remove variables that have the same value for every data point.
    step_zv(all_predictors()) %>%
    # Transform variables to be distributed as Gaussian-like as possible.
    step_YeoJohnson(all_numeric()) %>%
    # Normalize numeric data to have a mean of zero and
    # standard deviation of one.
    step_normalize(all_numeric())
  
  applicability_pca <- apd_pca(recipe_applicable, data = training_tbl %>% 
                                 filter(clinical == 'benign'), threshold = threshold)
  
  origin_benign <- applicable::score(applicability_pca, training_tbl %>% filter(clinical == 'benign'))
  
  origin_benign <- origin_benign %>% 
    rename(pc1_origin = PC01, pc2_origin = PC02) %>%
    select(pc1_origin, pc2_origin)
  
  
  
  if (isTRUE(only_origin)) {
    
    tmp_origin_pos <- training_tbl %>% filter(clinical == 'pathogenic') %>% select(id, clinical)
    origin_pathogenic <- tmp_origin_pos %>% bind_cols(origin_pathogenic)
    
    tmp_origin_neg <- training_tbl %>% filter(clinical == 'benign') %>% select(id, clinical)
    origin_benign <- tmp_origin_neg %>% bind_cols(origin_benign)

    tmp_result <- origin_pathogenic %>% bind_rows(origin_benign)
    
    return(tmp_result)
    
    }
  
  tmp_distance_benign <- applicable::score(applicability_pca, input_tbl) %>%
    rename(dist_benign = distance_pctl, pc1_benign = PC01, pc2_benign = PC02) %>%
    select(dist_benign, pc1_benign, pc2_benign)
  
  
  # Distance total 
  
  
  recipe_applicable <-
    recipe( ~ ., data = training_tbl %>%
              select(any_of(all.vars(formule_input)[-1]))) %>%
    step_dummy(all_nominal()) %>%
    # Remove variables that have the same value for every data point.
    step_zv(all_predictors()) %>%
    # Transform variables to be distributed as Gaussian-like as possible.
    step_YeoJohnson(all_numeric()) %>%
    # Normalize numeric data to have a mean of zero and
    # standard deviation of one.
    step_normalize(all_numeric())
  
  applicability_pca <- apd_pca(recipe_applicable, data = training_tbl, threshold = threshold)
  
  tmp_distance_total <- applicable::score(applicability_pca, input_tbl) %>%
    rename(dist_total = distance_pctl, pc1_total = PC01, pc2_total = PC02) %>%
    select(dist_total, pc1_total, pc2_total)
  
  
  
  # Results
  
  tmp_result <- tmp_result %>% 
    bind_cols(tmp_distance_pathogenic, tmp_distance_benign, tmp_distance_total)
  
  
  return(tmp_result)
}

# ------------------------------------------------------------------------------
# PLOT APPLICABILITY RESULTS
# https://www.marlycormar.com/presentations/R-Pharma-2019/presentation.html#1
# ------------------------------------------------------------------------------

plot_applicability <- function(pca_input, dataframe_input, formule_input) {
  
  
  # pca_score <- score(pca_input, dataframe_input)
  
  dataframe_input <- total_df_deletion
  
  pca_output <- dataframe_input %>%
    select(any_of(all.vars(formule_input)[-1])) %>%
    prcomp(scale = TRUE, center = TRUE)
  
  
  pca_output$x %>%
    as_tibble() %>%
    select(PC1, PC2) %>%
    bind_cols(dataframe_input %>% select(clinical)) %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(fill = clinical), shape = 21, color = 'black') +
    theme_minimal()
  
}

# ------------------------------------------------------------------------------
# ENRICH RESULTS
# ------------------------------------------------------------------------------


enrich_prediction <- function(input_model, input_test, applicability_input) {
  
  # input_test <- df_manolo
  # applicability_input <- applicability_pca_deletion
  
  # tmp_result <- posterior_epred(input_model, newdata = input_test) %>%
  #   as_tibble()
  # 
  # # point estimate
  # tmp_point_estimate <- tmp_result %>%
  #   map_dbl(~ median(.x)) %>%
  #   enframe(name = NULL) %>%
  #   rename(pred_target = value) %>%
  #   mutate(pred_target = 1 - pred_target)
  
  # uncertainty
  # tmp_uncertainty <- tmp_result %>% map_dbl(~ sd(.x)) %>% enframe(name = NULL) %>%
  #   rename(sd = value)
  # 
  # tmp_mad <- tmp_result %>% map_dbl(~ mad(.)) %>% enframe(name = NULL) %>%
  #   rename(mad = value)
  # 
  
  # applicability (larger distance -> lower rank)
  tmp_applicability <- applicable::score(applicability_input, input_test %>%
                                           mutate_if(is.factor, as.character) %>%
                                           mutate_if(is.character, as.double) %>%
                                           select(-c(type_variant, source, clinical, id))) %>%
    select(distance_pctl)
  # mutate(distance_pctl = 100 - distance_pctl)
  
  # final_result <- input_test %>%
  #   select(clinical) %>%
  #   bind_cols(tmp_point_estimate,
  #             tmp_uncertainty,
  #             tmp_mad,
  #             tmp_applicability)
  
  
  return(final_result)
}


# ------------------------------------------------------------------------------
# PCA - ENRICHED RESULTS
# ------------------------------------------------------------------------------


pca_enrich <- function(enrich_tbl, training_tbl, query_tbl, tag = 'Deletions' ) {
  
  # enrich_tbl <- enriched_test_deletion
  # training_tbl <- training_tbl_deletion
  # query_tbl <- test_tbl_deletion
  
  
  pca_plot1 <- enrich_tbl %>%
    mutate(cat = case_when(
      distance_pctl >= 90 & sd >= 0.2 ~ 'outliers',
      distance_pctl <= 10 & sd >= 0.2 ~ 'between',
      sd < 0.2 ~ 'ok',
      TRUE ~ 'no_clear'
    )) %>%
    ggplot(aes(sd, distance_pctl)) +
    geom_point(aes(fill = cat), shape = 21, size = 2) +
    geom_segment(x = 0.2, xend = 0.4, y = 10, yend = 10, linetype = 'dashed') +
    geom_segment(x = 0.2, xend = 0.4, y = 90, yend = 90, linetype = 'dashed') +
    coord_cartesian(ylim = c(0,100)) +
    geom_segment(x = 0.2, xend = 0.2, y = 0, yend = 100, linetype = 'dashed') +
    theme_minimal() +
    labs(title = glue('{tag}'))
  
  
  pre_tbl <- query_tbl %>%
    bind_cols(enrich_tbl %>%     mutate(cat = case_when(
      distance_pctl >= 90 & sd >= 0.2 ~ 'outliers',
      distance_pctl <= 10 & sd >= 0.2 ~ 'between',
      sd < 0.2 ~ 'ok',
      TRUE ~ 'no_clear'
    )) %>% select(-clinical)) %>%
    mutate(type = 'test') %>%
    bind_rows(training_tbl %>% mutate(type = 'training'))
  
  
  pca_plot2 <- pre_tbl %>%
    select(any_of(all.vars(formule_models)[-1])) %>%
    prcomp(scale = TRUE, center = TRUE) %>%
    pluck('x') %>%
    as_tibble() %>%
    select(PC1, PC2) %>%
    bind_cols(pre_tbl %>% select(type, cat, clinical)) %>%
    filter(type == 'training') %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(fill = factor(clinical)), shape = 21, color = 'black', size = 2) +
    stat_ellipse(aes(color = clinical)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_minimal()
  
  
  pca_plot3 <- pre_tbl %>%
    select(any_of(all.vars(formule_models)[-1])) %>%
    prcomp(scale = TRUE, center = TRUE) %>%
    pluck('x') %>%
    as_tibble() %>%
    select(PC1, PC2) %>%
    bind_cols(pre_tbl %>% select(type, type, cat, clinical)) %>%
    filter(type == 'test') %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(fill = factor(clinical)), shape = 21, color = 'black', size = 2) +
    stat_ellipse(aes(color = clinical)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_minimal()
  
  pca_plot4 <- pre_tbl %>%
    select(any_of(all.vars(formule_models)[-1])) %>%
    prcomp(scale = TRUE, center = TRUE) %>%
    pluck('x') %>%
    as_tibble() %>%
    select(PC1, PC2) %>%
    bind_cols(pre_tbl %>% select(type, type, cat, clinical)) %>%
    filter(type == 'test') %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(fill = factor(cat)), shape = 21, color = 'black', size = 2) +
    stat_ellipse(aes(color = cat)) +
    geom_vline(xintercept = 0, linetype="dashed") +
    geom_hline(yintercept = 0, linetype="dashed") +
    theme_minimal()
  
  return(pca_plot1 + pca_plot2 + pca_plot3 + pca_plot4)
}

# ------------------------------------------------------------------------------
# PLOT CORRELATIONS ENRICHED RESULTS
# ------------------------------------------------------------------------------


plot_enriched <- function(enrich_tbl) {
  # enrich_tbl <- enriched_test_deletion
  
  
  p_corr1 <- enrich_tbl %>%
    select(-clinical) %>%
    correlate(method = 'spearman') %>%
    stretch(remove.dups = TRUE) %>%
    na.omit() %>%
    mutate(id = paste(x, y, sep =' - ')) %>%
    mutate(positive = if_else(r >= 0, 'yes', 'no')) %>%
    filter(x != 'mad') %>%
    filter(y != 'mad') %>%
    filter(id != 'pred_target - sd') %>%
    ggplot(aes(id, r)) +
    geom_col(aes(fill = positive), color = 'black') +
    coord_flip() +
    labs(x = 'Association', y = 'Spearman correlation', title = 'Correlation') +
    theme_minimal()
  
  p_corr2 <- enrich_tbl %>%
    ggplot(aes(sd, distance_pctl)) +
    geom_point() +
    theme_minimal()
  
  p_corr3 <- enrich_tbl %>%
    ggplot(aes(pred_target, sd)) +
    geom_point() +
    theme_minimal()
  
  p_corr4 <- enrich_tbl %>%
    ggplot(aes(pred_target, distance_pctl)) +
    geom_point() +
    theme_minimal()
  
  p_corr5 <- enrich_tbl %>%
    mutate(cat = case_when(
      pred_target <= 0.05  ~ 'less than 0.05',
      pred_target >= 0.95 ~ 'more than 0.95',
      TRUE ~ 'between 0.05 and 0.95'
    )) %>%
    mutate(cat = fct_relevel(cat, 'less than 0.05', 'between 0.05 and 0.95', 'more than 0.95')) %>%
    ggplot(aes(cat, distance_pctl)) +
    geom_boxplot(aes(fill = cat)) +
    theme_minimal() +
    labs(y = 'Distance percentiles (based on PCA)', x = 'Category', fill = 'Legend')
  
  p_corr1 + p_corr2 + p_corr3 + p_corr4 + p_corr5 + plot_layout(ncol = 2)
  
}


# ------------------------------------------------------------------------------
# GENERATE RULES - ITERATE OVER THE ROWS
# ------------------------------------------------------------------------------


per_row <- function(input_id, input_obs, input_model) {
  
  # input_id <- input_obs$id[1]
  # input_model <- model
  
  # tmp_df <- input_model %>%
  #   separate_rows(rule, sep = '&') %>%
  #   mutate(rule = str_remove_all(rule, ' ')) %>%
  #   mutate(rule = case_when(
  #     str_detect(rule, '>=') ~ str_replace(rule, '>=', ' >= '),
  #     str_detect(rule, '<=') ~ str_replace(rule, '<=', ' <= '),
  #     str_detect(rule, '>') ~ str_replace(rule, '>', ' > '),
  #     str_detect(rule, '<') ~ str_replace(rule, '<', ' < ')
  #   )) %>%
  #   group_by(term) %>%
  #   summarise(rule = paste(rule, collapse = ' & '), .groups = 'drop')
  # 
  # input_model <- input_model %>% 
  #   select(coefficient_lambda.min, term) %>%
  #   left_join(tmp_df, by = 'term')
  
  # pivot_wider(id_cols = -rule, names_from = rule, values_from = rule, values_fn = list(rule = str_c)) %>%
  # unite(-c(coefficient_lambda.min,term), col = 'rule', sep = ' & ', na.rm = TRUE)
  
  
  obs <- input_obs %>% filter(id %in% input_id)
  
  tmp_obs <- input_obs %>%
    filter(id %in% input_id) %>%
    bind_cols(
      input_model %>%
        rowwise() %>%
        mutate(yes_no = if_else(nrow(obs %>% filter_(rule)) == 1, 1, 0)) %>%
        # mutate(yes_no = if_else(nrow(obs %>% filter(eval_tidy(parse_expr(rule)))) == 1, 1, 0)) %>%
        # mutate(yes_no = if_else(nrow(obs %>% filter(rlang::eval_tidy(rlang::parse_expr(rule)))) == 1, 1, 0)) %>%
        ungroup() %>%
        select(term, yes_no) %>%
        pivot_wider(values_from = yes_no, names_from = term))
  
  return(tmp_obs)
}

# ------------------------------------------------------------------------------
# CONVERT VARIABLES DBL TO FACTORS
# ------------------------------------------------------------------------------

# to_factor <- function(input_tbl) {
#
#   test_tbl_deletion_annotated %>% mutate_at(across())
#
#
#   input_tbl <- test_tbl_deletion_annotated
#
#   input_tbl %>%
#     mutate()
#
#
# }

# ------------------------------------------------------------------------------
# ADD 1 - RULEFIT MODEL
# ------------------------------------------------------------------------------


# add_number_one <- function(testing_tbl) {
#   
#   # testing_tbl <- test_tbl_deletion
#   
#   result <- testing_tbl %>%
#     mutate(cnv_syndromes1 = cnv_syndromes) %>%
#     mutate(embryo_mouse1 = embryo_mouse) %>%
#     mutate(n_target_drugs1 = n_target_drugs) %>%
#     mutate(patho_cnv1 = patho_cnv) %>%
#     mutate(essent_cl1 = essent_cl) %>%
#     mutate(n_prot_complex1 = n_prot_complex) %>%
#     mutate(n_genes_hpo1 = n_genes_hpo)
#   
#   return(result)
#   
# }


remove_numbers <- function(rulefit_model) {
  
  result <- coef(rulefit_model) %>%
    as_tibble() %>%
    mutate(rule = str_replace_all(rule, 'cnv_syndromes1', 'cnv_syndromes')) %>%
    mutate(rule = str_replace_all(rule, 'embryo_mouse1', 'embryo_mouse')) %>%
    mutate(rule = str_replace_all(rule, 'n_target_drugs1', 'n_target_drugs')) %>%
    mutate(rule = str_replace_all(rule, 'patho_cnv1', 'patho_cnv')) %>%
    mutate(rule = str_replace_all(rule, 'essent_cl1', 'essent_cl')) %>%
    mutate(rule = str_replace_all(rule, 'n_prot_complex1', 'n_prot_complex')) %>%
    mutate(rule = str_replace_all(rule, 'n_genes_hpo1', 'n_genes_hpo'))
  
  
  return(result)
  
}

# ------------------------------------------------------------------------------
# GENERATE RULES FOR THE DATASET
# ------------------------------------------------------------------------------


generate_rules <- function(model, input_obs) {
  
  # model <- rulefit_model_deletion
  # input_obs <- test_tbl_deletion
  
  
  rulefit_model <- remove_numbers(model)
  
  
  model <- rulefit_model %>% filter(str_detect(term, '^[r]'))
  
  tmp_result <- tibble()
  
  for (i in 1:nrow(input_obs)) {
    
    print(glue('{i}/{nrow(input_obs)}'))
    
    obs <- input_obs[i,]
    
    tmp_obs <- input_obs[i,] %>%
      bind_cols(
        model %>%
          rowwise() %>%
          mutate(yes_no = if_else(nrow(obs %>% filter_(rule)) == 1, 1, 0)) %>%
          ungroup() %>%
          select(term, yes_no) %>%
          pivot_wider(values_from = yes_no, names_from = term))
    
    tmp_result <- tmp_result %>% bind_rows(tmp_obs)
  }
  
  return(tmp_result)
  
}
# ------------------------------------------------------------------------------
# API CNVSCORE
# https://medium.com/@JB_Pleynet/how-to-do-an-efficient-r-api-81e168562731
# ------------------------------------------------------------------------------

# library(plumber)


# ------------------------------------------------------------------------------
# FUNCTION - LENGTH DISTRIBUTION - PLOT
# ------------------------------------------------------------------------------

plot_length_distribution <- function(data) {
  
  # data <- input_check_cnv_deletion
  
  p_distribution0 <- data %>%
    ggplot(aes(length_cnv)) +
    geom_histogram(aes(fill = class_plot),color = 'black', alpha = 0.4, bins = 30) +
    scale_x_log10() +
    # scale_y_log10() +
    facet_wrap(~ class_plot) +
    theme_minimal() +
    theme(axis.title.y = element_text(size=15, face="bold"),
          axis.title.x = element_text(size=15, face="bold"))
  
  
  
  p_distribution0_2 <- data %>%
    group_by(class_plot) %>%
    mutate(p_length_cnv = ntile(length_cnv, 100)) %>%
    group_by(class_plot, p_length_cnv) %>%
    mutate(median_length_cnv = median(length_cnv)) %>%
    ungroup() %>%
    ggplot(aes(p_length_cnv, median_length_cnv)) +
    geom_point(aes(fill = class_plot), shape = 21) +
    scale_y_log10() +
    theme_minimal() +
    theme(axis.title.y = element_text(size=15, face="bold"),
          axis.title.x = element_text(size=15, face="bold"))
  
  
  p_distribution0_3 <- data %>%
    ggplot(aes(length_cnv)) +
    geom_density(aes(fill = class_plot), alpha = 0.4) +
    scale_x_log10() +
    facet_wrap(~ class_plot) +
    theme_minimal() +
    theme(axis.title.y = element_text(size=15, face="bold"),
          axis.title.x = element_text(size=15, face="bold"))
  
  p_distribution0_4 <- data %>%
    ggplot(aes(length_cnv)) +
    geom_density(aes(fill = clinical), alpha = 0.4) +
    scale_x_log10() +
    theme_minimal(base_size = 20) +
    theme(axis.title.y = element_text(size=15, face="bold"),
          axis.title.x = element_text(size=15, face="bold"))
  
  # p_distribution0 + p_distribution0_2 + p_distribution0_3 + p_distribution0_4 + plot_layout(nrow = 2)
  (p_distribution0  + p_distribution0_3) / p_distribution0_4 + plot_layout(nrow = 2)
  
}

# ------------------------------------------------------------------------------
# FUNCTION - CUMULATIVE PLOT
# ------------------------------------------------------------------------------

cum_plot <- function(by_percentage = TRUE, label = NULL, target_n, total_n, benign = FALSE,  select_filter, ...) {
  
  
  # by_percentage <- TRUE
  # target_n <- 146
  # total_n <- clinvar_ind_deletion %>% nrow()
  # select_filter <- 'benign'
  # input_total <- logistic_clinvar_cum_dist_deletion %>% bind_rows(xgboost_clinvar_cum_dist_deletion)
  select_filter <- if_else(benign, 'benign', 'pathogenic')
  input_total <- bind_rows(...)
  
  random_perc <- (target_n /  total_n) / target_n
  random_tbl <- tibble('model' = 'random', cum_line = rep(random_perc, total_n))
  random_tbl <- random_tbl %>%
    mutate(perc_line = cumsum(cum_line)) %>%
    mutate(number_cnv = row_number())
  
  if (by_percentage & !benign) {
    
    input_total %>%
      group_by(model) %>%
      arrange(desc(pred_target)) %>% # desc() is the only change
      mutate(cum_line1 = if_else(clinical == select_filter, 1, 0)) %>%
      mutate(cum_line = cumsum(cum_line1)) %>%
      mutate(perc_line = cum_line / target_n) %>%
      mutate(number_cnv = row_number()) %>%
      ggplot(aes(number_cnv, perc_line)) +
      geom_path(aes(group = model, color = model), size = 2) +
      geom_path(data = random_tbl, aes(number_cnv, perc_line, group = model, color = model), size = 2, linetype = 'dashed') +
      scale_y_continuous(label = scales::percent) +
      xlab('Total number of CNVs') +
      ylab(glue('Cumulative frequency')) +
      ggtitle(glue('Cumulative frequency {select_filter} CNVs - {label}')) +
      theme_minimal()
    
  } else {
    
    
    input_total %>%
      group_by(model) %>%
      arrange(pred_target) %>%
      mutate(cum_line1 = if_else(clinical == select_filter, 1, 0)) %>%
      mutate(cum_line = cumsum(cum_line1)) %>%
      mutate(perc_line = cum_line / target_n) %>%
      mutate(number_cnv = row_number()) %>%
      ggplot(aes(number_cnv, perc_line)) +
      geom_path(aes(group = model, color = model), size = 2) +
      geom_path(data = random_tbl, aes(number_cnv, perc_line, group = model, color = model), size = 2, linetype = 'dashed') +
      scale_y_continuous(label = scales::percent) +
      xlab('Total number of CNVs') +
      ylab(glue('Cumulative frequency')) +
      ggtitle(glue('Cumulative frequency {select_filter} CNVs - {label}')) +
      theme_minimal()
  }
  # } else {
  
  # input_total %>%
  #   group_by(model) %>%
  #   arrange(desc(pred_target)) %>%
  #   mutate(cum_line1 = if_else(clinical == select_filter, 1, 0)) %>%
  #   mutate(cum_line = cumsum(cum_line1)) %>%
  #   mutate(number_cnv = row_number()) %>%
  #   ggplot(aes(number_cnv, cum_line)) +
  #   geom_path(aes(group = model, color = model), size = 2) +
  #   xlab('Total number of CNVs') +
  #   ylab(glue('Cumulative frequency')) +
  #   ggtitle(glue('Cumulative frequency {select_filter} CNVs - {label}')) +
  #   geom_hline(yintercept = target_n, linetype = 'dashed') +
  #   theme_minimal()
  
  
  # }
  
}

# ------------------------------------------------------------------------------
# FUNCTION - RECIPROCAL OVERLAP
# ------------------------------------------------------------------------------

reciprocal_overlap <- function(main_tbl) {
  
  # main_tbl <- input_check_cnv_del_benign
  
  ##
  
  main_tbl <- main_tbl %>% mutate(length_cnv = end - start + 1)
  
  main_intersect <- main_tbl %>% 
    bed_intersect(main_tbl) %>%
    mutate(length_cnv = end.x - start.x + 1) %>%
    mutate(overlap = (.overlap + 1)/ length_cnv) %>%
    filter(overlap >= 0.9) %>%
    rename(ref = id_tmp.x, sub = id_tmp.y) %>%
    select(ref, sub)
  
  main_intersect <- main_intersect %>% 
    left_join(main_intersect, by = c('sub' = 'ref')) %>%
    filter(ref == sub.y) %>%
    select(ref, sub) %>% 
    filter(ref != sub)
  
  if (nrow(main_intersect) == 0) stop('No reciprocal overlap found.')
  
  
  highest_hits <- main_intersect %>% count(ref) %>% arrange(desc(n)) %>%
    pull(ref)
  
  discard_cnvs <- c()
  
  for (i in 1:length(highest_hits)) {
    
    print(glue('{i} / {length(highest_hits)}'))
    
    if (highest_hits[i] %in% discard_cnvs) next
    
    overlap_vector <- main_intersect %>% filter(ref == highest_hits[i]) %>%
      pull(sub)
    
    discard_cnvs_tmp <- main_tbl %>%
      filter(id_tmp %in% c(highest_hits[i],overlap_vector )) %>%
      arrange(length_cnv) %>%
      select(id_tmp, length_cnv) %>%
      slice(-1) %>%
      pull(id_tmp)
    
    discard_cnvs <- c(discard_cnvs, discard_cnvs_tmp)
  }
  return(discard_cnvs)
}

# ------------------------------------------------------------------------------
#  FUNCTION - RECIPROCAL OVERLAP - INDEPENDENT
# ------------------------------------------------------------------------------

reciprocal_overlap_ind <- function(main_tbl, secondary_tbl) {
  
  # main_tbl <- output_clinvar_deletion
  # secondary_tbl <- output_df_deletion
  # #
  
  main_tbl <- main_tbl %>% mutate(length_cnv = end - start + 1)
  secondary_tbl <- secondary_tbl %>% mutate(length_cnv = end - start + 1)
  
  main_intersect <- main_tbl %>% 
    bed_intersect(secondary_tbl) %>%
    mutate(overlap = (.overlap + 1)/ length_cnv.x) %>%
    filter(overlap >= 0.9) %>%
    rename(ref = id_tmp.x, sub = id.y) %>%
    select(ref, sub)
  
  second_intersect <-  secondary_tbl %>% 
    bed_intersect(main_tbl) %>%
    mutate(overlap = (.overlap + 1)/ length_cnv.x) %>%
    filter(overlap >= 0.9) %>%
    rename(ref = id.x, sub = id_tmp.y) %>%
    select(ref, sub)
  
  main_intersect <- main_intersect %>%
    left_join(second_intersect, by = c('sub' = 'ref')) %>%
    filter(ref == sub.y) %>%
    select(ref, sub)
  
  
  if (nrow(main_intersect) == 0) stop('No reciprocal overlap found.')
  
  
  discard_cnvs <- main_intersect %>% count(ref) %>% arrange(desc(n)) %>%
    pull(ref)
  
  
  return(discard_cnvs)
}

# ------------------------------------------------------------------------------
#  FUNCTION - GENERATE RANDOM CNVs FOR EACH DECIPHER CNV
# ------------------------------------------------------------------------------


gen_random <- function(input_chrom, input_start, input_end, input_tmp_id, n_rep = 9, input_n_try) {
  
  result_df <- tibble()
  
  tmp_chrom <- input_chrom
  tmp_start <- input_start
  tmp_end <- input_end
  tmp_id <- input_tmp_id
  
  tmp_length <- tmp_end - tmp_start + 1
  
  target_cnv <- tibble('chrom' = tmp_chrom,'start' = tmp_start, 'end' = tmp_end)
  
  
  tmp_cyto <- bed_intersect(coord_cytobands, target_cnv) %>%
    pull(Name.x)
  
  
  tmp_n_genes <- hgcn_genes %>% bed_intersect(target_cnv) %>% nrow()
  
  result_tmp_cnv <- tibble()
  
  n_try <- 0
  
  while (nrow(result_tmp_cnv) < n_rep) {
    
    n_try <- n_try + 1
    
    if (n_try < input_n_try) {
      
      # print(n_try)
      
      
      selected_cytoband <- coord_cytobands %>% filter(Name %in% tmp_cyto & chrom == tmp_chrom)
      random_from <- selected_cytoband %>% slice_head() %>% pull(start)
      random_to <-  selected_cytoband %>% slice_tail() %>% pull(end) # - tmp_length
      random_start <- sample(random_from:random_to, 1)
      random_end <- random_start + tmp_length
      
      random_tbl <- tibble('chrom' = tmp_chrom,
                           'start' = random_start,
                           'end' = random_end)
      
      # Filter 1: Remove random CNVs mapping the original CNV
      overlap_original_cnv <- random_tbl %>% bed_intersect(target_cnv)
      if (nrow(overlap_original_cnv) > 0) next
      
      # Filter 2: eliminate cnvs with no genes if the original cnv has genes
      random_genes <- hgcn_genes %>% bed_intersect(random_tbl) %>% nrow()
      if (tmp_n_genes > 0 & random_genes == 0) next
      
      result_tmp_cnv <- result_tmp_cnv %>% bind_rows(random_tbl)
      
    } else {
      
      result_tmp_cnv <- tibble('chrom' = 'error', 'start' = 0, 'end' = 0)
      break
    }
    
  }
  result_df <- result_df %>% bind_rows(result_tmp_cnv %>% mutate('id_tmp' = tmp_id))
  # }
  
  return(result_df)
  
}

# ------------------------------------------------------------------------------
# FUNCTION - EVALUATE BIAS ASSOCIATED WITH CHROMOSOMES
# ------------------------------------------------------------------------------

plot_bias <- function(main_tbl, tag = 'Deletion') {
  
  # main_tbl <- output_df_deletion
  
  first_split <- initial_split(main_tbl, strata = clinical)
  
  training_tbl <- first_split %>% training()
  test_tbl <- first_split %>% testing()
  
  
  logistic_bias_chrom <- logistic_reg() %>%
    set_mode('classification') %>%
    set_engine('glm') %>%
    fit(clinical ~ chrom , data = training_tbl)
  
  logistic_bias_length <- logistic_reg() %>%
    set_mode('classification') %>%
    set_engine('glm') %>%
    fit(clinical ~ length_cnv , data = training_tbl)
  
  logistic_bias_n_genes <- logistic_reg() %>%
    set_mode('classification') %>%
    set_engine('glm') %>%
    fit(clinical ~ n_genes , data = training_tbl)
  
  logistic_bias_disease <- logistic_reg() %>%
    set_mode('classification') %>%
    set_engine('glm') %>%
    fit(clinical ~ disease , data = training_tbl)
  
  logistic_bias_result_chrom <- evaluate_model(logistic_bias_chrom, test_object = test_tbl,
                                               model_name = 'Clinical ~ chrom', cv = FALSE,
                                               cum_dist = FALSE)
  
  logistic_bias_result_length <- evaluate_model(logistic_bias_length, test_object = test_tbl,
                                                model_name = 'Clinical ~ length', cv = FALSE,
                                                cum_dist = FALSE)
  
  logistic_bias_result_n_genes <- evaluate_model(logistic_bias_n_genes, test_object = test_tbl,
                                                 model_name = 'Clinical ~ nº genes', cv = FALSE,
                                                 cum_dist = FALSE)
  
  logistic_bias_disease <- evaluate_model(logistic_bias_disease, test_object = test_tbl,
                                          model_name = 'Clinical ~ Disease genes (Yes/No)', cv = FALSE,
                                          cum_dist = FALSE)
  
  
  
  logistic_bias_result_chrom[[1]] %>%
    bind_rows(logistic_bias_result_length[[1]],
              logistic_bias_result_n_genes[[1]],
              logistic_bias_disease[[1]]) %>%
    ggplot(aes(1-specificity, sensitivity)) +
    geom_path(aes(group = model, color = model), size = 2.5, show.legend = TRUE) +
    theme_minimal() +
    theme(plot.title = element_text(size=10, face="bold"),
          axis.title.x = element_text(size=15, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    ggtitle(glue('{tag} - logistic model (Clinical ~ chrom) -  = {logistic_bias_result_chrom[[3]]},
                 logistic model (Clinical ~ length) -  = , {logistic_bias_result_length[[3]]},
                 logistic model (Clinical ~ n_genes) = {logistic_bias_result_n_genes[[3]]},
                 logistic model (Clinical ~ disease genes (yes/no)) = {logistic_bias_disease[[3]]}'))
  
  
}


# ------------------------------------------------------------------------------
# FUNCTION - ANNOTATE IN PARALLEL
# ------------------------------------------------------------------------------

annot_parallel <- function(input_tbl, n_cores = 20, reg = FALSE) {
  
  # input_tbl <- tmp_merge_clinvar
  
  test_ok <- check_cnv(input_tbl$id_tmp[1],
                       input_tbl$clinical[1],
                       input_tbl$variant_class[1],
                       input_tbl$chrom[1],
                       input_tbl$start[1],
                       input_tbl$end[1])
  
  if (nrow(test_ok) == 1) {
    
    print('OK - REG FALSE')
  } else {
    stop('The test did not work!')
  }
  
  test_ok <- check_cnv(input_tbl$id_tmp[1],
                       input_tbl$clinical[1],
                       input_tbl$variant_class[1],
                       input_tbl$chrom[1],
                       input_tbl$start[1],
                       input_tbl$end[1], mode_reg = TRUE)
  
  if (nrow(test_ok) == 1) {
    
    print('OK - REG TRUE')
  } else {
    stop('The test did not work!')
  }
  
  
  plan("multiprocess", workers = n_cores)
  
  
  tic()
  
  if (isFALSE(reg)) {
    
    output_df <- future_pmap_dfr(list(input_tbl$id_tmp, 
                                      input_tbl$clinical, 
                                      input_tbl$variant_class,
                                      input_tbl$chrom, 
                                      input_tbl$start, 
                                      input_tbl$end), 
                                 check_cnv,
                                 .progress = TRUE)
  } else {
    
    output_df <- future_pmap_dfr(list(input_tbl$id_tmp, 
                                      input_tbl$clinical, 
                                      input_tbl$variant_class,
                                      input_tbl$chrom, 
                                      input_tbl$start, 
                                      input_tbl$end), 
                                 check_cnv(mode_reg = TRUE),
                                 .progress = TRUE)
  }
  
  toc()
  
  return(output_df)
}

# ------------------------------------------------------------------------------
# FUNCTION - RUN XGBOOST IN PARALLEL
# ------------------------------------------------------------------------------


parallel_xgboost <- function(training_set, chrom_to_run, formule_used, list_hyper) {
  
  # training_set <- output_df_deletion
  # split <- FALSE
  # chrom_to_run <- '1'
  # formule_used <- formule_total
  
  # if (isFALSE(split)) {
  
  train_split <- training_set %>% filter(chrom != chrom_to_run)
  test_split <- training_set %>% filter(chrom == chrom_to_run)
  #   
  # } else {
  # 
  # train_split <- training_set %>% filter(chrom != chrom_to_run)
  # test_split <- test_set %>% filter(chrom == chrom_to_run)
  # }
  
  # tmp_model <- boost_tree(trees = list_hyper$trees, 
  #                         tree_depth = list_hyper$tree_depth, 
  #                         learn_rate = list_hyper$learn_rate) %>%
  #   set_mode("classification") %>%
  #   set_engine("xgboost") %>%
  #   fit(formule_used, data = train_split)
  
  tmp_model <- boost_tree() %>%
    set_mode("classification") %>%
    set_engine("xgboost") %>%
    fit(formule_used, data = train_split)
  
  
  prob_predicted <- predict(tmp_model, test_split, type = 'prob') %>%
    bind_cols(test_split %>% select(clinical, id))
  
  
  # return(list('prob_df' = prob_predicted, 'model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  return(list('model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
}

# ------------------------------------------------------------------------------
# FUNCTION - RUN RF IN PARALLEL
# ------------------------------------------------------------------------------


parallel_rf <- function(training_set, chrom_to_run, formule_used, list_hyper) {
  
  # training_set <- total_df_deletion
  # chrom_to_run <- '1'
  # formule_used <- formule_length
  
  
  
  train_split <- training_set %>% filter(chrom != chrom_to_run)
  test_split <- training_set %>% filter(chrom == chrom_to_run)
  
  
  if (is.null(list_hyper)) {
    
    tmp_model <-
      rand_forest() %>%
      set_engine('randomForest') %>%
      set_mode('classification') %>%
      fit(formule_used, data = train_split)
    
  } else {
    
    
    tmp_model <-
      rand_forest(trees = list_hyper$trees, 
                  min_n = list_hyper$min_n) %>%
      set_engine('randomForest') %>%
      set_mode('classification') %>%
      fit(formule_used, data = train_split)
    
  }
  
  prob_predicted <- predict(tmp_model, test_split, type = 'prob') %>%
    bind_cols(test_split %>% select(clinical, id))
  
  
  # return(list('prob_df' = prob_predicted, 'model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  return(list('model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
}


# ------------------------------------------------------------------------------
# FUNCTION - RUN GBM IN PARALLEL
# ------------------------------------------------------------------------------


parallel_gbm <- function(training_set, chrom_to_run, formule_used, list_hyper) {
  
  # training_set <- total_df_deletion
  # chrom_to_run <- '1'
  # formule_used <- formule_length
  
  
  
  train_split <- training_set %>% 
    filter(chrom != chrom_to_run) %>%
    mutate(clinical = if_else(
      clinical == 'benign', 0, 1))
  
  
  test_split <- training_set %>% filter(chrom == chrom_to_run)
  
  
  if (is.null(list_hyper)) {
    
    tmp_model <- gbm.fit <- gbm(
      formula = formule_used,
      data = train_split,
      n.trees = 500,
      distribution = 'bernoulli',
      bag.fraction = 0.5,
      shrinkage = 0.001,
      n.minobsinnode = 1,
      interaction.depth = 2
    )  
    
  } else {
    
    
    # tmp_model <-
    #   rand_forest(trees = list_hyper$trees, 
    #               min_n = list_hyper$min_n) %>%
    #   set_engine('randomForest') %>%
    #   set_mode('classification') %>%
    #   fit(formule_used, data = train_split)
    
  }
  
  prob_predicted <- predict(tmp_model, newdata = test_split, type = 'response') %>%
    as_tibble() %>% rename(.pred_pathogenic = value) %>%
    bind_cols(test_split %>% select(clinical, id))
  
  
  
  # return(list('prob_df' = prob_predicted, 'model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  return(list('model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  
}

# ------------------------------------------------------------------------------
# FUNCTION - RUN LOGISTIC FUNCTION IN PARALLEL
# ------------------------------------------------------------------------------


parallel_logistic <- function(training_set, chrom_to_run, formule_used, list_hyper) {
  
  # training_set <- total_df_deletion
  # chrom_to_run <- '1'
  # formule_used <- formule_length
  
  
  
  train_split <- training_set %>% filter(chrom != chrom_to_run)
  test_split <- training_set %>% filter(chrom == chrom_to_run)
  
  
  if (is.null(list_hyper)) {
    
    tmp_model <-
      logistic_reg() %>%
      set_engine('glm') %>%
      set_mode('classification') %>%
      fit(formule_used, data = train_split)
    
  } else {
    
    
    # tmp_model <-
    #   rand_forest(trees = list_hyper$trees, 
    #               min_n = list_hyper$min_n) %>%
    #   set_engine('randomForest') %>%
    #   set_mode('classification') %>%
    #   fit(formule_used, data = train_split)
    
  }
  
  prob_predicted <- predict(tmp_model, test_split, type = 'prob') %>%
    bind_cols(test_split %>% select(clinical, id))
  
  
  # return(list('prob_df' = prob_predicted, 'model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  return(list('model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  
}

# ------------------------------------------------------------------------------
# FUNCTION - RUN RULEFIT MODELS IN PARALLEL
# ------------------------------------------------------------------------------

parallel_rulefit <- function(training_set, chrom_to_run, formule_used, 
                             list_hyper = NULL) {
  
  
  # training_set <- output_df_deletion_train
  # chrom_to_run <- '1'
  # formule_used <- human_no_control
  # list_hyper
  
  print(chrom_to_run)
  train_split <- training_set %>% filter(chrom != chrom_to_run)
  
  
  if (!is.null(list_hyper)) {
    
    tmp_model <-
      rules::rule_fit(learn_rate = list_hyper$learning_rate, 
                      trees = list_hyper$nrounds, tree_depth = list_hyper$max_dept) %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("xrf") %>%
      fit(formule_used, data = train_split )
  } else {
    
    tmp_model <-
      rules::rule_fit() %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("xrf") %>%
      fit(formule_used, data = train_split )
    
    
  }
  
  
  
  # tmp_model <- xrftest::xrf(formule_used,
  #                           data = train_split,
  #                           xgb_control = list_hyper,
  #                           family = "binomial")
  
  
  
  return(list('model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  
  # return(list('model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
  
  # return(list('prob_df' = prob_predicted, 'model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
}

# ------------------------------------------------------------------------------
# FUNCTION - GENERATE PARALLEL MATCHING RULES
# ------------------------------------------------------------------------------

generate_rtemis <- function(newdata, rules_set) {
  
  # newdata <- df_predict
  # rules_set <- tmp_rulefit
  
  matrix_zeros <- matrix(0, nrow(newdata), nrow(rules_set))
  
  newdata <- newdata %>% mutate(id = row_number())
  
  for (i in seq(nrow(rules_set))) {
    
    match <- newdata %>% filter_(rules_set$rule[i]) %>% pull(id)
    matrix_zeros[match, i] <- 1
  }
  df_zeros <- matrix_zeros %>% as_tibble() %>% bind_cols(newdata %>% select(clinical))
  return(df_zeros)
}

# ------------------------------------------------------------------------------
# FUNCTION - RUN BAYESIAN MODELS IN PARALLEL
# ------------------------------------------------------------------------------
# 
# parallel_bayesian <- function(training_set, chrom_to_run, 
#                               list_hyper = NULL,
#                               formule_used, split = FALSE, prior = 'normal') {
#   
#   print(chrom_to_run)
#   # training_set <- output_df_deletion
#   # tag_variant <- 'deletion'
#   # nc <- 10
#   # list_hyper <- tibble(nrounds = 10, learning_rate = 0.2, max_dept = 2, min_child_weight = 1)
#   # prior <- 'normal'
#   # tag_formule <- 'human_no_control'
#   # model_name <- 'bayesian'
#   # formule_used <- human_no_control
#   # split <- FALSE
#   # chrom_to_run <- '1'
#   
#   
#   if (is.null(list_hyper)) {
#     list_hyper <- list(nrounds = 20, scale_pos_weight = 9)
#   } else {
#     list_hyper <- list(nrounds = list_hyper$nrounds,
#                        learning_rate = list_hyper$learning_rate, 
#                        max_dept = list_hyper$max_dept, 
#                        min_child_weight = list_hyper$min_child_weight)
#   }
#   
#   if (isFALSE(split)) {
#     
#     train_split <- training_set %>% filter(chrom != chrom_to_run)
#     test_split <- training_set %>% filter(chrom == chrom_to_run)
#     
#   } else {
#     
#     train_split <- training_set %>% filter(chrom != chrom_to_run)
#     test_split <- test_set %>% filter(chrom == chrom_to_run)
#   }
#   
#   if (prior == 'normal') {
#     
#     prior_selected <- normal()
#     
#   } else if (prior == 'hs') {
#     
#     prior_selected <- hs()
#     
#   } else if (prior == 'laplace') {
#     
#     prior_selected <- laplace()
#     
#   } else if (prior == 'hs_plus') {
#     
#     prior_selected <- hs_plus()
#     
#   } else if (prior == 'lasso') {
#     
#     prior_selected <- lasso()
#   }
#   
#   
#   tmp_rulefit <- xrftest::xrf(formule_used,
#                               data = train_split,
#                               xgb_control = list_hyper,
#                               family = "binomial")
#   
#   set_rules <- coef(tmp_rulefit)$term
#   set_rules <- set_rules[str_detect(set_rules, '^[r]')]
#   set_rules <- set_rules[!str_detect(set_rules, 'recomb_rate')]
#   
#   
#   
#   formule_used <- reformulate(termlabels = c(all.vars(formule_used[-2]
#   ), set_rules), response = 'clinical')
#   # tic()
#   tmp_model <- rstanarm::stan_glm(formula = formule_used,
#                                   family = 'binomial',
#                                   data = tmp_rulefit$full_data,
#                                   # select_if(~  sum(.x == 1) > 2),
#                                   cores = 2,
#                                   iter = 1500,
#                                   chains = 2,
#                                   algorithm = 'sampling', # variational inference algorithms
#                                   QR = FALSE,
#                                   prior = prior_selected,
#                                   prior_intercept = student_t())
#   
#   
#   return(list('prob_df' = prob_predicted, 'model_trained_rulefit' = tmp_rulefit,
#               'model_trained' = tmp_model, 'chrom_target' = chrom_to_run))
#   
# }

# ------------------------------------------------------------------------------
# FUNCTION - RTEMIS 1
# ------------------------------------------------------------------------------



rtemis_step1 <- function(input_tbl, tag_variant, vector_features, tag_features, input_prior, nc, input_hyper, input_iter = 3000) {
  
  
  # input_tbl = output_clinvar_duplication
  # tag_variant = 'duplication'
  # vector_features = c('length_cnv')
  # tag_features = 'test1231'
  # input_prior = 'hs'
  # nc = 23
  # input_hyper = tibble(trees = 500, depth = 2, min_n = 1)
  # input_iter <- 500
  
  if (input_prior == 'normal') {
    
    prior_selected <- normal()
    
  } else if (input_prior == 'hs') {
    
    prior_selected <- hs()
    
  } else if (input_prior == 'laplace') {
    
    prior_selected <- laplace()
    
  } else if (input_prior == 'hs_plus') {
    
    prior_selected <- hs_plus()
    
  } else if (input_prior == 'lasso') {
    
    prior_selected <- lasso()
  }
  
  plan("multiprocess", workers = nc)
  
  tag_used <- paste(tag_variant, 'bayesian', input_prior, tag_features, sep = ' - ')
  
  vector_chrom_input_tbl <- input_tbl %>% count(chrom) %>% pull(chrom)
  
  output_list <- vector_chrom_input_tbl %>%
    future_map(~ rtemis_step2(input_tbl, .,
                              vector_features = vector_features,
                              prior = input_prior,
                              input_hyper = input_hyper,
                              input_iter = input_iter))
  return(output_list)
}



# ------------------------------------------------------------------------------
# FUNCTION - RUN RTEMIS + BAYESIAN MODELS IN PARALLEL
# ------------------------------------------------------------------------------

rtemis_step2 <- function(training_set, chrom_to_run, 
                         vector_features,
                         formule_used, prior = NULL, input_hyper, input_iter = 3000) {
  
  # chrom_to_run = 8
  # training_set <- input_tbl
  # vector_features = vector_features
  # prior = input_prior
  # input_hyper = input_hyper
  
  
  train_split <- training_set %>% filter(chrom != chrom_to_run) %>% select(-chrom)
  
  park.rf <- s.RULEFEAT(x = train_split %>% select(all_of(vector_features)),
                        y = train_split %>% select(clinical) %>% 
                          mutate(clinical = as.factor(if_else(clinical == 'pathogenic', 1, 0))) %>%
                          mutate(clinical = factor(clinical, levels = c('1', '0'))) %>%
                          pull(),
                        n.trees = input_hyper$trees,
                        gbm.params = list(interaction.depth = input_hyper$depth, n.minobsinnode = input_hyper$min_n))
  
  tmp_rulefit <- park.rf$mod$rule.coefs %>% as_tibble(rownames = 'rule_id') %>%
    rename(rule = Rule, coefficient = Coefficient) %>%
    mutate(rule = as.character(rule))
  
  train_split_annotated <- generate_rtemis(train_split, tmp_rulefit)
  
  if (prior == 'normal') {
    
    prior_selected <- normal()
    
  } else if (prior == 'hs') {
    
    prior_selected <- hs()
    
  } else if (prior == 'laplace') {
    
    prior_selected <- laplace()
    
  } else if (prior == 'hs_plus') {
    
    prior_selected <- hs_plus()
    
  } else if (prior == 'lasso') {
    
    prior_selected <- lasso()
  }
  
  tmp_model <- rstanarm::stan_glm(formula = clinical ~ .,
                                  family = 'binomial',
                                  data = train_split_annotated,
                                  cores = 2,
                                  iter = input_iter,
                                  chains = 2,
                                  algorithm = 'sampling', 
                                  QR = FALSE,
                                  prior = prior_selected)
  
  print(glue('Chromosome {chrom_to_run} DONE!'))
  
  return(list('model_trained_rulefit' = park.rf,
              'model_trained' = tmp_model, 'chrom_target' = chrom_to_run, 'set_rules' = tmp_rulefit))
}


# ------------------------------------------------------------------------------
# FUNCTION - CHROMOSOME-AWARE DISTRIBUTION
# ------------------------------------------------------------------------------
# 
chrom_aware <- function(input_tbl, tag_variant = 'deletion', tag_formule = 'human_control',
                        model_name = NULL,
                        formule_model = NULL,
                        input_prior = 'normal',
                        list_hyper = NULL,
                        nc = 23) {

  # input_tbl <- output_df_deletion_train
  # tag_variant <- 'deletion'
  # model_name <- 'random_forest'
  # tag_formule <- 'formule_total'
  # formule_model <- formule_total
  # input_prior <- 'normal'
  # nc <- 12
  # list_hyper <- tibble(nrounds = 10, learning_rate = 0.2, max_dept = 2, min_child_weight = 1)
  # list_hyper <- NULL



  # plan("multiprocess", workers = nc)


  vector_chrom_input_tbl <- input_tbl %>% count(chrom) %>% pull(chrom)

  p3 <- input_tbl %>%
    count(chrom, clinical) %>%
    mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign'))) %>%
    mutate(chrom = factor(chrom, levels = c(seq(1,22), 'X'))) %>%
    ggplot(aes(chrom, n)) +
    geom_col(aes(fill = clinical), position = 'dodge', color = 'black') +
    theme_minimal()

  # Specify model

  if (model_name == 'random_forest') {

    output_list <- vector_chrom_input_tbl %>%
      future_map(~ parallel_rf(input_tbl, ., formule_model, list_hyper))
    
  } else if (model_name == 'logistic') {
    
    output_list <- vector_chrom_input_tbl %>%
      future_map(~ parallel_logistic(input_tbl, ., formule_model, list_hyper))
    
  } else if (model_name == 'gbm') {
    
    output_list <- vector_chrom_input_tbl %>%
      future_map(~ parallel_gbm(input_tbl, ., formule_model, list_hyper))
    


  } else if (model_name == 'xgboost') {

    output_list <- vector_chrom_input_tbl %>%
      future_map(~ parallel_xgboost(input_tbl, ., formule_model, list_hyper))

  } else if (model_name == 'bayesian') {

    output_list <- vector_chrom_input_tbl %>%
      future_map(~ parallel_bayesian(input_tbl, .,
                                     formule_used = formule_model,
                                     prior = input_prior,
                                     list_hyper = list_hyper))

  } else if (model_name == 'rulefit') {

    output_list <- vector_chrom_input_tbl %>%
      future_map(~ parallel_rulefit(input_tbl, ., formule_model, list_hyper))

  } else if (model_name == 'bay_rtemis') {

    output_list <- vector_chrom_input_tbl %>%
      future_map(~ parallel_bay_rtemis(input_tbl, .,
                                       formule_used = formule_model,
                                       prior = input_prior,
                                       list_hyper = list_hyper))
  }

  if (model_name == 'bayesian' | model_name == 'bay_rtemis') {

    tag_used <- paste(tag_variant, model_name, input_prior, tag_formule, sep = ' - ')

  } else {
    tag_used <- paste(tag_variant, model_name, tag_formule, sep = ' - ')
  }

  result_list <- list(output_list, 'tag' = tag_used)

  return(result_list)
}

# ------------------------------------------------------------------------------
# FUNCTION - PREDICT BAY_RTEMIS CHROMOSOMES
# ------------------------------------------------------------------------------

predict_rtemis <- function(x, chrom_tmp, vector_chrom, df_predict, 
                           only_pos_dist = FALSE) {
  
  # x <- bayesian_clinvar_del_nohuman
  # chrom_tmp <- 16
  # vector_chrom <- 3
  # df_predict <- output_decipher_deletion
  
  tmp_rulefit <- x[[chrom_tmp]][["set_rules"]]
  tmp_model <- x[[chrom_tmp]][["model_trained"]]
  
  df_predict <- df_predict %>% filter(chrom == vector_chrom[chrom_tmp])
  # df_predict <- df_predict %>% filter(chrom == 3)
  
  
  if (nrow(df_predict) == 0) {
    
    return(tibble())
    
  }
  
  testing_set_annotated <- generate_rtemis(df_predict, tmp_rulefit)
  
  rules_tmp <- testing_set_annotated %>%
    select(-clinical) %>%
    bind_cols(df_predict %>% select(id)) %>%
    pivot_longer(-id, names_to = 'rule_id', values_to = 'value') %>%
    filter(value == 1) %>%
    left_join(tmp_rulefit, by = 'rule_id') %>%
    left_join(tmp_model$coefficients %>% as_tibble(rownames = 'rule_id') %>%
                rename(coeff_bayesian = value), by = 'rule_id') %>%
    # filter(coefficient != 0) %>%
    # mutate(sign_coefficient = ifelse(coefficient > 0, 'pos', 'neg')) %>%
    # # filter(sign_coefficient == 'pos') %>%
    # filter(risk > 0.8) %>%
    # group_by(sign_coefficient) %>%
    # arrange(desc(coefficient)) %>%
    # slice_head(n = 5)
    # mutate(coeff_bayesian = abs(coeff_bayesian)) %>%
    # filter(coeff_bayesian >= 0.015 | coeff_bayesian <= -0.015) %>%
    # group_by(id) %>%
    # arrange(desc(coeff_bayesian)) %>%
    # slice_head(n = 5) %>%
# 
#     mutate(rule = paste0(rule, ' - risk:', risk, ' - support:', support,
#                          ' - coeff:', coeff_bayesian)) %>%
    select(id, rule_id) %>%
    group_by(id) %>%
    summarise(rules = str_c(rule_id, collapse =";"))
    
  tmp_posterior <- posterior_epred(tmp_model, newdata = testing_set_annotated) %>%
    as_tibble()
  
  tmp2_posterior <- tmp_posterior
  
  colnames(tmp2_posterior) <- df_predict$id
  
  if (isTRUE(only_pos_dist)) return(tmp2_posterior)
  
  
  
  avoid_too_sparse_error <- function(x){
    tryCatch(
      # This is what I want to do...
      {
        return(map_estimate(x))
      },
      # ... but if an error occurs, tell me what happened: 
      error=function(error_message) {
        return(1)
      }
    )
  }

  
  prob_tmp <- tmp_posterior %>% 
    map_dbl(~ median(.x)) %>%
    enframe(name = NULL) %>%
    rename(.pred_pathogenic = value) %>%
    mutate(.pred_pathogenic = 1 - .pred_pathogenic)
  
  prob2_tmp <- tmp_posterior %>% 
    map_dbl(~ avoid_too_sparse_error(.x)) %>%
    enframe(name = NULL) %>%
    rename(map_estimate = value) %>%
    mutate(map_estimate = 1 - map_estimate)
  
  
  sd_tmp <- tmp_posterior %>% map_dbl(~ sd(.x)) %>%
    enframe(name = NULL) %>%
    rename(sd = value)
  
  mad_tmp <- tmp_posterior %>% map_dbl(~ mad(.x)) %>%
    enframe(name = NULL) %>%
    rename(mad = value)
  
  tmp_predicted <- prob_tmp %>% 
    bind_cols(sd_tmp, df_predict %>% select(clinical, id)) %>%
    bind_cols(prob2_tmp) %>%
    bind_cols(mad_tmp) %>%
    left_join(rules_tmp, by = 'id') %>%
    mutate(chrom = vector_chrom[chrom_tmp])
    # mutate(vector_posterior = paste(tmp_posterior, collapse = ', '))
  
  return(tmp_predicted)
}


# ------------------------------------------------------------------------------
# FUNCTION - PREDICT ACROSS CHROMOSOMES
# ------------------------------------------------------------------------------

predict_chrom_aware_rtemis <- function(list_models, test_split,
                                       tag_variant = 'deletion',
                                       tag_features, only_table = FALSE) {
  
  # list_models <- bayesian_clinvar_del_nohuman
  # test_split <- just_test
  # 
  # tag_variant <- 'deletion'
  # tag_features <- 'unbiased_approach'
  
  tag_used <- paste(tag_variant, 'bayesian', tag_features, sep = ' - ')
  
  vector_chrom <- unlist(map(list_models, function(x) x[['chrom_target']]))
  
  iter_vector <- 1:length(vector_chrom)
  
  tmp_predicted <- iter_vector %>%
    map_dfr(~ predict_rtemis(list_models, .x, vector_chrom, test_split))
  
  # tmp_post <- iter_vector %>%
  #   map(~ predict_rtemis(list_models, .x, vector_chrom, test_split, only_pos_dist = TRUE))
  
  if (isTRUE(only_table)) return(tmp_predicted)
  
  
  ss_auc <- tmp_predicted %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
  pr_auc <- tmp_predicted %>% pr_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
  
  
  tmp_predicted <- tmp_predicted %>% mutate(tag = paste(tag_used, '-', ss_auc))
  
  tmp_roc_curve <- tmp_predicted %>% roc_curve(clinical, .pred_pathogenic)  %>% mutate(tag = paste(tag_used, '-', ss_auc))
  tmp_pr_curve <- tmp_predicted %>% pr_curve(clinical, .pred_pathogenic)  %>% mutate(tag = paste(tag_used, '-', pr_auc))
  
  # result <- list('tmp_roc_curve' = tmp_roc_curve, 'tmp_pr_curve' = tmp_pr_curve, 'tmp_predicted' = tmp_predicted, 'post_dist' = tmp_post)
  
  result <- list('tmp_roc_curve' = tmp_roc_curve, 'tmp_pr_curve' = tmp_pr_curve, 'tmp_predicted' = tmp_predicted)
  
  return(result)
  
}

# ------------------------------------------------------------------------------
# FUNCTION - PREDICT ACROSS CHROMOSOMES
# ------------------------------------------------------------------------------

predict_chrom_aware <- function(list_models, df_to_predict, only_table = FALSE,
                                is_bancco = FALSE, tag_bancco = NULL) {
  
  # list_models <- gbm_clinvar_del_nohuman
  # df_to_predict <- decipher_setting_3
  
  if (is_bancco) {
    
    tag_used <- tag_bancco
    vector_chrom <- unlist(map(list_models, function(x) x[['chrom_target']]))

  } else {
  
  tag_used <- list_models[["tag"]]
  vector_chrom <- unlist(map(list_models[[1]], function(x) x[['chrom_target']]))
  
  }
  
  for (i in 1:length(vector_chrom)) {
    
    tmp_chrom <- vector_chrom[i]
    
    test_split <- df_to_predict %>% filter(chrom == tmp_chrom)
    
    if (is_bancco) {
      
      tmp_model <- list_models[[i]][["model_trained"]]
      
    } else {
      tmp_model <- list_models[[1]][[i]][["model_trained"]]
      
    }
    
    
    if (nrow(test_split) == 0) {
      
      
      tmp_predicted <- tibble()
      
    } else if (str_detect(tag_used, 'random_forest')) {
      
      tmp_predicted <- predict(tmp_model, test_split, type = 'prob')
      
      tmp_predicted <- tmp_predicted %>%
        bind_cols(test_split %>% select(clinical, id)) 
      
    } else if (str_detect(tag_used, 'xgboost')) {
      
      tmp_predicted <- predict(tmp_model, test_split, type = 'prob')
      
      tmp_predicted <- tmp_predicted %>%
        bind_cols(test_split %>% select(clinical, id)) 
      
    } else if (str_detect(tag_used, 'gbm')) {
      
      tmp_predicted <- predict(tmp_model, test_split, type = 'response')
      
      tmp_predicted <- tmp_predicted %>%
        as_tibble() %>% rename(.pred_pathogenic = value) %>%
        bind_cols(test_split %>% select(clinical, id)) 
      
    } else if (str_detect(tag_used, 'logistic')) {
      
      tmp_predicted <- predict(tmp_model, test_split, type = 'prob')
      
      tmp_predicted <- tmp_predicted %>%
        bind_cols(test_split %>% select(clinical, id)) 
      
      
    } else if (str_detect(tag_used, 'rulefit')) {
      
      if (nrow(test_split) == 1) {
        
        tmp_predicted <- tibble('.pred_pathogenic' = NA)
      } else {
        
        tmp_predicted <- predict(tmp_model, test_split, type = 'prob')
        
        tmp_predicted <- tmp_predicted %>%
          bind_cols(test_split %>% select(clinical, id)) 
        
      }
      
    }  
    
    if (i == 1) {
      
      prob_predicted <- tmp_predicted
      
    } else {
      
      prob_predicted <- prob_predicted %>% bind_rows(tmp_predicted)
      
    }
  }
  
  prob_predicted <- prob_predicted %>% mutate(tag = tag_used)
  
  if (nrow(prob_predicted %>% filter(is.na(.pred_pathogenic))) > 0) {
    
    print(glue('There is {nrow(prob_predicted %>% filter(is.na(.pred_pathogenic)))} NA values'))
    
  }
  
  if (isTRUE(only_table)) return(prob_predicted)
  
  auc_result <- prob_predicted %>%
    roc_auc(clinical, .pred_pathogenic)  %>%
    pull(.estimate) %>%
    round(3)
  
  pr_auc_result <- prob_predicted %>%
    pr_auc(clinical, .pred_pathogenic)  %>%
    pull(.estimate) %>%
    round(3)
  
  roc_result <- prob_predicted %>%
    roc_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = paste(tag_used, '-', auc_result))
  
  pr_result <- prob_predicted %>%
    pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = paste(tag_used, '-', pr_auc_result))
  
  prob_predicted <- prob_predicted %>% mutate(tag = paste(tag, '-',  auc_result ))
  
  return(list('tmp_roc_curve' = roc_result, 'tmp_pr_curve' = pr_result, 'tmp_predicted' = prob_predicted))
  
}


# ------------------------------------------------------------------------------
# FUNCTION - EVALUATE BAYESIAN MODEL
# ------------------------------------------------------------------------------

evaluate_baymodel <- function(model_list) {
  
  # https://jrnold.github.io/bayesian_notes/mcmc-diagnostics.html
  # ESS (effective sample size) - amount of independent information
  # RHAT (Potential Scale Reduction)
  # MCSE (Monte Carlo Standard Error (MCSE)) - how big the estimation noise is sd / effect size
  
  # model_list <- model_vs9[[1]]
  # model_list <- model_list[[1]]
  # model_list <- bayesian_clinvar_del_nohuman
  
  
  vector_rows <- c()
  
  for ( i in 1:length(model_list)) {
    
    tmp_model <- model_list[[i]]$model_trained
    tmp_table <- describe_posterior(tmp_model, test = NULL, centrality = NULL, diagnostic = 'all')
    
    tmp_n_row <-   tmp_table %>% 
      as_tibble() %>%
      filter(Rhat >= 1.1 | ESS < 1000) %>% nrow()
    
    # vector_rows <- c(vector_rows, tmp_n_row)
    
    print(glue('Model {i} - {tmp_n_row} features'))
    
  }
  
  
}


# ------------------------------------------------------------------------------
# FUNCTION - CALCULATE ODDS RATIO BETWEEN TWO CATEGORIES OF GENES (BACKGROUND = HGNC_GENES)
# ------------------------------------------------------------------------------

run_odds <- function(x, y, tag, tag2) {
  
  # x <- genes_pyu
  # y <- hgcn_genes %>% filter(fusil == 'CL (Cellular lethal)') %>% pull(gene)
  
  a <- sum(x %in% y)
  b <- length(y) - a
  c <- length(x) - a
  d <- hgcn_genes %>% filter(! gene %in% c(x,y)) %>% pull(gene) %>% length()
  
  result <- fisher.test( matrix(c(a,b,c,d), nrow = 2, byrow = TRUE)) %>% 
    glance() %>%
    select(-method, -alternative) %>%
    # pivot_longer(everything(), names_to = 'param', values_to = 'value') %>%
    mutate(tag = tag) %>%
    mutate(tag2 = tag2)
  
  return(result)
  
}


# ------------------------------------------------------------------------------
# FUNCTION - PLOT  RESULTS (ROC-AUC + PR-AUC + CUMULATIVE FREQUENCY)
# ------------------------------------------------------------------------------



evaluate_model <- function(model, model_name,
                           test_object = NULL,
                           rulefit = FALSE,
                           bay = FALSE,
                           cv = FALSE,
                           cum_dist = FALSE) {
  
  
  
  # model <- bayesian_model_duplication
  # test_object <- clinvar_tbl_duplication_annotated
  # bay = TRUE
  # model_name = 'Bayesian'
  # cv = FALSE
  # cum_dist = FALSE
  
  if (isTRUE(cv) & isFALSE(cum_dist)) {
    
    # ROC AUC
    tmp_mean <- model %>%
      unnest(.predictions) %>%
      group_by(id) %>%
      roc_auc(clinical, target_class) %>%
      summarise(mean(.estimate)) %>%
      round(2)
    
    tmp_err <- model %>%
      unnest(.predictions) %>%
      group_by(id) %>%
      roc_auc(clinical, target_class) %>%
      summarise(sd(.estimate) / sqrt(10)) %>%
      round(3)
    
    p1 <- model %>%
      unnest(.predictions) %>%
      # filter(id == 'Fold01') %>%
      group_by(id) %>%
      # ggplot(aes(.pred_Pathogenic)) + geom_density(aes(fill = clinical), alpha = 0.3) + facet_wrap(~ id, scales = 'free')
      roc_curve(clinical, target_class) %>%
      ggplot(aes(1-specificity, sensitivity)) +
      geom_path(aes(group = id, color = id),  show.legend = T) +
      theme_bw() +
      theme(plot.title = element_text(size=19, face="bold"),
            axis.title.x = element_text(size=19, face="bold"),
            axis.title.y = element_text(size=19, face="bold"),
            axis.text.x = element_text(size=14, face="bold"),
            axis.text.y = element_text(size=14, face="bold")) +
      
      labs(title = glue('{model_name}  (Mean AUC: {tmp_mean} ± {tmp_err})'))
    
    # PR AUC
    tmp_mean <- model %>%
      unnest(.predictions) %>%
      group_by(id) %>%
      pr_auc(clinical, target_class) %>%
      summarise(mean(.estimate)) %>%
      round(2)
    
    tmp_err <-model %>%
      unnest(.predictions) %>%
      group_by(id) %>%
      pr_auc(clinical, target_class) %>%
      summarise(sd(.estimate) / sqrt(10)) %>%
      round(3)
    
    p2 <- model %>%
      unnest(.predictions) %>%
      group_by(id) %>%
      pr_curve(clinical, target_class) %>%
      ggplot(aes(recall, precision)) +
      geom_path(aes(group = id, color = id), show.legend = FALSE) +
      theme_bw() +
      theme(plot.title = element_text(size=19, face="bold"),
            axis.title.x = element_text(size=19, face="bold"),
            axis.title.y = element_text(size=19, face="bold"),
            axis.text.x = element_text(size=14, face="bold"),
            axis.text.y = element_text(size=14, face="bold")) +
      labs(title = glue('{model_name} (Mean AUCpr: {tmp_mean} ± {tmp_err})'))
    
    return(list(p1, p2))
    
    
  } else if (isTRUE(cv) & isTRUE(cum_dist)) {
    
    
    if (isFALSE(rulefit)) {
      perc_10_cv_tmp <- model %>%
        unnest(.predictions) %>%
        group_by(id) %>%
        rename(p_patho = target_class) %>%
        mutate(p_patho = ntile(-p_patho, 10)) %>%
        select(id, p_patho, clinical) %>%
        group_by(id, p_patho) %>%
        count(clinical) %>%
        filter(clinical == select_filter) %>%
        group_by(id) %>%
        mutate(percentage = n / sum(n)) %>%
        group_by(id) %>%
        mutate(cum_perc = cumsum(percentage)) %>%
        mutate(model = model_name)
      
    } else {
      
      perc_10_cv_tmp <- model %>%
        select(id, prob_predicted, clinical) %>%
        group_by(id) %>%
        mutate(prob_predicted = ntile(-prob_predicted, 10)) %>%
        select(id, prob_predicted, clinical) %>%
        group_by(id, prob_predicted) %>%
        count(clinical) %>%
        filter(clinical == select_filter) %>%
        group_by(id) %>%
        mutate(percentage = n / sum(n)) %>%
        mutate(cum_perc = cumsum(percentage)) %>%
        mutate(model = model_name) %>%
        rename(p_patho = prob_predicted)
    }
    
    return(perc_10_cv_tmp)
    
  } else if (isFALSE(cv) & isTRUE(cum_dist)) {
    
    if (isFALSE(rulefit) & isFALSE(bay)) {
      
      tmp_object <- predict(model, test_object, type = 'prob') %>%
        rename(pred_target = target_class) %>%
        select(pred_target) %>%
        bind_cols(test_object %>% select(clinical)) %>% mutate(model = model_name)
      
    } else if (isTRUE(bay)) {
      
      
      tmp_object <-  posterior_epred(model, newdata = test_object) %>%
        as_tibble() %>%
        map_dbl(~ median(.x)) %>% # mean, median, map_estimate
        enframe(name = NULL) %>%
        bind_cols(test_object %>% select(clinical)) %>%
        rename(pred_target = value) %>%
        mutate(pred_target = 1 - pred_target) %>%
        mutate(model = model_name)
      
      
    } else {
      
      tmp_object <-  tibble(target_class = 1 - as.vector(predict(model, test_object, type = 'response'))) %>%
        rename(pred_target = target_class) %>%
        select(pred_target) %>%
        bind_cols(test_object %>% select(clinical))  %>%
        mutate(model = model_name)
      
      
    }
    
    return(tmp_object)
    
  } else if (isFALSE(cv) & isFALSE(cum_dist)) {
    
    if (isFALSE(rulefit) & isFALSE(bay)) {
      
      tmp_object <- predict(model, test_object, type = 'prob') %>%
        rename(pred_target = target_class) %>%
        select(pred_target) %>%
        bind_cols(test_object %>% select(clinical)) %>% mutate(model = model_name)
      
    } else if (isTRUE(bay)) {
      
      tmp_object <- posterior_epred(model, newdata = test_object) %>%
        as_tibble() %>%
        map_dbl(~ median(.x)) %>% # mean, median, map_estimate
        enframe(name = NULL) %>%
        bind_cols(test_object %>% select(clinical)) %>%
        rename(pred_target = value) %>%
        mutate(pred_target = 1 - pred_target) %>%
        mutate(model = model_name)
      
      
    } else {
      
      tmp_object <-  tibble(target_class = 1 - as.vector(predict(model, test_object, type = 'response'))) %>%
        rename(pred_target = target_class) %>%
        select(pred_target) %>%
        bind_cols(test_object %>% select(clinical))  %>%
        mutate(model = model_name)
    }
    
    tmp_roc <- tmp_object %>%
      roc_curve(clinical, pred_target) %>%
      mutate(model = model_name)
    
    tmp_pr <- tmp_object %>%
      pr_curve(clinical, pred_target) %>%
      mutate(model = model_name)
    
    tmp_auc <- tmp_object %>%
      roc_auc(clinical, pred_target) %>%
      mutate(model = model_name) %>%
      pull(.estimate) %>%
      round(2)
    
    tmp_pr_auc <- tmp_object %>%
      pr_auc(clinical, pred_target) %>%
      mutate(model = model_name) %>%
      pull(.estimate) %>%
      round(2)
    
    
    
    return (list(tmp_roc, tmp_pr, tmp_auc, tmp_pr_auc))
    
  }
  
}


# ------------------------------------------------------------------------------
# FUNCTION (ONLY FOR RULEFIT) - PLOT CROSS-VALIDATION RESULTS (ROC-AUC + PR-AUC)
# ------------------------------------------------------------------------------


from_cv_rulefit <- function(model, model_name = 'RuleFit') {
  
  # model <- rulefit_rs
  
  model <- model %>% mutate(prob_predicted = 1 - prob_predicted)
  
  # ROC AUC
  tmp_mean <- model %>%
    group_by(id) %>%
    roc_auc(clinical, prob_predicted) %>%
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <- model %>%
    group_by(id) %>%
    roc_auc(clinical, prob_predicted) %>%
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p1 <- model %>%
    group_by(id) %>%
    roc_curve(clinical, prob_predicted) %>%
    ggplot(aes(1-specificity, sensitivity)) +
    geom_path(aes(group = id, color = id),  show.legend = T) +
    theme_bw() +
    theme(plot.title = element_text(size=19, face="bold"),
          axis.title.x = element_text(size=19, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    
    labs(title = glue('{model_name}  (Mean AUC: {tmp_mean} ± {tmp_err})'))
  
  # PR AUC
  tmp_mean <- model %>%
    group_by(id) %>%
    pr_auc(clinical, prob_predicted) %>%
    summarise(mean(.estimate)) %>%
    round(2)
  
  tmp_err <- model %>%
    group_by(id) %>%
    pr_auc(clinical, prob_predicted) %>%
    summarise(sd(.estimate) / sqrt(10)) %>%
    round(3)
  
  p2 <- model %>%
    group_by(id) %>%
    pr_curve(clinical, prob_predicted) %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = id, color = id), show.legend = FALSE) +
    theme_bw() +
    theme(plot.title = element_text(size=19, face="bold"),
          axis.title.x = element_text(size=19, face="bold"),
          axis.title.y = element_text(size=19, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold")) +
    labs(title = glue('{model_name} (Mean AUCpr: {tmp_mean} ± {tmp_err})'))
  
  return(list(p1, p2))
  
}

# ------------------------------------------------------------------------------
# FUNCTION - BAYESIAN RULEFIT - CROSS-VALIDATION
# ------------------------------------------------------------------------------


bay_rulefit_cv <- function(cv_object) {
  
  bayesian_result <- tibble()
  
  for (i in 1:10) {
    
    print(glue('Rulefit - split nº {i}/10 '))
    
    
    tmp_training <- cv_object$splits[[i]] %>% training()
    tmp_testing <- cv_object$splits[[i]] %>% testing()
    tmp_name_fold <- cv_object %>% slice(i) %>% pull(id)
    
    
    data_pre_lasso <- xrftest::xrf(formule_models,
                                   data = tmp_training,
                                   # xgb_control = list(nrounds = 100, max_depth = 5, min_child_weight = 3),
                                   family = "binomial")
    
    print(glue('Bayesian - split nº {i}/10 '))
    
    
    # stan(file = 'stan_model.stan' , iter= 2000, data = data_pre_lasso$full_data,
    #      chains= 4, seed=194838)
    
    # special_columns <-    data_pre_lasso %>%
    #   coef(s = "lambda.min") %>%
    #   as_tibble() %>% filter(coefficient_lambda.min >= 1 | coefficient_lambda.min <= -1) %>%
    #   na.omit() %>%
    #   pull(term)
    
    bayesian_model <- rstanarm::stan_glm(formule_models,
                                         family = 'binomial',
                                         data = data_pre_lasso$full_data,
                                         # data = data_pre_lasso$full_data[,colnames(data_pre_lasso$full_data) %in% special_columns],
                                         # cores = 4,
                                         iter = 2000,
                                         # chains = 4,
                                         algorithm = 'meanfield', # variational inference algorithms
                                         QR = TRUE,
                                         prior = normal())
    # QR = TRUE,
    # prior = laplace())
    
    pred_bayesian  <- posterior_epred(bayesian_model, newdata = tmp_testing) %>%
      as_tibble() %>%
      map_dbl(~ mean(.x))
    
    tmp_probs <- pred_bayesian %>% as_tibble() %>% rename(prob_predicted = value)
    
    tmp_result <- tmp_probs %>% bind_cols(tmp_testing %>% select(clinical)) %>% mutate(id = tmp_name_fold)
    
    if (i == 1) {
      
      bayesian_result <- tmp_result
      
    } else {
      
      bayesian_result <- bayesian_result %>% bind_rows(tmp_result)
      
    }
    
  }
  return(bayesian_result)
}

# ------------------------------------------------------------------------------
# FUNCTION - RULEFIT MODEL - CROSS-VALIDATION
# ------------------------------------------------------------------------------


rulefit_cv <- function(cv_object) {
  
  tbl_result <- tibble()
  
  for (i in 1:10) {
    print(glue('Split nº {i}/10'))
    
    tmp_training <- cv_object$splits[[i]] %>% training()
    tmp_testing <- cv_object$splits[[i]] %>% testing()
    tmp_name_fold <- cv_object %>% slice(i) %>% pull(id)
    
    
    tmp_xrf <- xrf(formule_models,
                   data = tmp_training,
                   xgb_control = list(nrounds = 20, scale_pos_weight = 9),
                   family = "binomial")
    
    tmp_probs <- predict(tmp_xrf, tmp_testing) %>% as_tibble() %>% rename(prob_predicted = `1`)
    
    tmp_result <- tmp_probs %>% bind_cols(tmp_testing %>% select(clinical)) %>% mutate(id = tmp_name_fold)
    
    if (i == 1) {
      
      tbl_result <- tmp_result
    } else {
      tbl_result <- tbl_result %>% bind_rows(tmp_result)
      
      
    }
  }
  return(tbl_result)
}


# ------------------------------------------------------------------------------
# FUNCTION - EXTRACT SUPPORT AND RISK
# ------------------------------------------------------------------------------

get_rules_information <- function(data, rule) {
  
  
  tmp_filter <-   data %>%
    filter_(rule)
  
  result_support <- tmp_filter %>% nrow()
  result_risk <- tmp_filter %>%
    filter(clinical == 'pathogenic') %>%
    nrow() / result_support
  
  return(paste(result_support, result_risk))
  
}


# ------------------------------------------------------------------------------
# GET PREDICTIONS FROM STRUCTURE
# ------------------------------------------------------------------------------

get_structure <- function(input_df, input_type = c('DEL', 'DUP')) {
  
  # Setting liftover
  library(rtracklayer)
  library(liftOver)
  
  from_hg19_to_hg38 = import.chain('/data-cbl/liftover/hg19ToHg38.over.chain')
  
  input_df_liftover <- clinvar_ind_deletion %>%
    mutate(tmp_id = row_number()) %>%
    select(chrom, start, end, tmp_id) %>%
    GRanges()
  
  seqlevelsStyle(input_df_liftover) = "UCSC"
  
  # WHY ALMOST X2
  input_df_liftover <- input_df_liftover %>%
    liftOver(from_hg19_to_hg38) %>%
    as_tibble()
  
  main_path <- '/data-cbl/frequena_data/rival_cnvscore/structure/StrVCTVRE-v.1.6/'
  
  use_python('/home/frequena/.conda/envs/py3.6/bin/python', required = TRUE)
  use_condaenv(condaenv = "py3.6", conda = "/usr/local/miniconda3/condabin/conda", required = TRUE)
  
  
  input_df <- input_df_liftover %>%
    mutate(start = start - 1) %>%
    mutate(type = input_type)
  # 1-based -> 0-based
  
  write_tsv(input_df, paste0(main_path, 'test_frequena/test.bed'),
            col_names = FALSE)
  
  setwd(main_path)
  
  system(glue('python {main_path}StrVCTVRE.py \\
           -i {main_path}test_frequena/test.bed \\
           -o {main_path}test_frequena/test_annotated.bed \\
           -p {main_path}data/hg38.phyloP100way.bw \\
           -f bed'))
  
  setwd('/data-cbl/frequena_data/cnvxplorer')
  
  output_df <- read_tsv('/data-cbl/frequena_data/rival_cnvscore/structure/StrVCTVRE-v.1.6/test_frequena/test_annotated.bed',
                        col_names = c('chrom', 'start', 'end', 'variant_type', 'prob')) %>%
    mutate(start = start + 1) # 0-based -> 1-based
  
  return(output_df)
  
}


# ------------------------------------------------------------------------------
# ROC THEME
# Source: https://github.com/sachsmc/plotROC/blob/master/R/style_roc.R
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
# JUST GET AUC (BOTH)
# ------------------------------------------------------------------------------

just_auc <- function(x) {
  
  # x <- result_annotsv_clinvar_deletion
  
  tmp_roc_auc <- x %>% 
    roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
  
  tmp_pr_auc <- x %>% 
    pr_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
  
  
  c(tmp_roc_auc, tmp_pr_auc)
  
}


# ------------------------------------------------------------------------------
# CALCULATE RISK AND SUPPORT RULES FROM TRAIN ANNOTATION DATASET
# ------------------------------------------------------------------------------


find_risk <- function(train_annotated) {
  
  # train_annotated <- train_part_annotated
  vector_risk <- c()
  vector_support <- c()
  
  for (i in 1:(ncol(train_annotated)-1)) {
    
    print(i)
    
    tmp_filter <- paste(colnames(train_annotated)[i], '== 1')
    
    tmp_support <- train_annotated %>%
      count(across(colnames(train_annotated)[i])) %>%
      filter_(tmp_filter) %>%
      pull(n) 
    
    tmp_risk <- train_annotated %>% 
      count(across(colnames(train_annotated)[i]), clinical) %>% 
      filter_(tmp_filter) %>% 
      mutate(perc = n / sum(n)) %>%
      filter(clinical == 'pathogenic') %>%
      pull(perc) %>%
      round(3)
    
    if (identical(tmp_risk, numeric(0))) tmp_risk <- 0
    if (identical(tmp_support, numeric(0))) tmp_support <- 0
    
    
    vector_support <- c(vector_support, tmp_support)
    
    vector_risk <- c(vector_risk, tmp_risk)
    
  }
  
  result_df <- tibble(term = colnames(train_annotated[-ncol(train_annotated)]),
                      risk = vector_risk, support = vector_support
  )
  
  return(result_df)
  
}


# ------------------------------------------------------------------------------
# MATCH PATHOGENIC AND BENIGN CNVs BY LENGTH
# ------------------------------------------------------------------------------


match_patho_benign <- function(df, max_length) {
  
  
  # df <- clinvar_ind_deletion
  # max_length <- max(clinvar_ind_deletion$length_cnv)
  
  tbl_bins <- tibble('start' = seq(1, 1e7, 100), 'end' = seq(100, 1e7, 100)) %>%
    mutate(chrom = 'chr1') %>%
    mutate(id = paste(start, end, sep = '-'))
  
  tmp_benign <- df %>% filter(clinical == 'benign')
  tmp_pathogenic <- df %>% filter(clinical == 'pathogenic')
  
  tmp_tbl_bins <- tbl_bins %>% filter(end <= max_length + 100)
  
  df_match <- tibble()
  rubbish_vector <- c()
  for (i in 1:nrow(tmp_pathogenic)) {
    
    print(glue('{i}/{nrow(tmp_pathogenic)}'))
    
    tmp_lenght_cnv <- tmp_pathogenic %>% slice(i) %>% pull(length_cnv)
    tmp_tmp_id <- tmp_pathogenic %>% slice(i) %>% pull(id)
    
    tmp_bin <- tmp_tbl_bins %>%
      rowwise() %>%
      filter(between(tmp_lenght_cnv, start, end))
    
    tmp_benign_selected <- tmp_benign %>%
      rowwise() %>%
      filter(between(length_cnv, tmp_bin$start, tmp_bin$end)) %>%
      ungroup()
    
    if (nrow(tmp_benign_selected) > 1) {
      
      tmp_id_benign_selected <- tmp_benign_selected %>% 
        filter(!id %in% rubbish_vector) %>%
        slice_sample(n = 1) %>% 
        pull(id)
    } else {
      tmp_id_benign_selected <- tmp_benign_selected %>% pull(id)
    }
    
    rubbish_vector <- c(rubbish_vector, tmp_id_benign_selected)
    
    tmp_df_match <- tibble('id' = tmp_tmp_id, 'id_match' = tmp_id_benign_selected)
    df_match <- df_match %>% bind_rows(tmp_df_match)
  }
  return(df_match)
}


# ------------------------------------------------------------------------------
# RUN CROSS-VALIDATION
# ------------------------------------------------------------------------------

cross_validation <- function(df_1, df_2, tag_variant = 'deletion', 
                             tag_formule = 'formule_total',
                             model_name = 'random_forest', 
                             input_hyper = NULL,
                             formule_model = formule_total,
                             vector_features = NULL, nc = 23) {
  
  # df_1 = output_df_deletion
  # df_2 = output_df_deletion_train
  # tag_variant = 'deletion'
  # tag_formule = 'human_control'
  # model_name = 'random_forest'
  # formule_model = human_control
  # vector_features = vector_yes_human_control
  
  
  # plan("multiprocess", workers = nc)
  
  tmp_df <- df_2 %>%
    mutate(cross_id = sample(rep(1:10, length.out = nrow(.))))
  
  cross_roc_auc <- c()
  cross_raw <- tibble()
  cross_ss <- tibble()
  cross_pr_auc <- c()
  cross_pr <- c()
  
  qc1 <- tibble()
  qc2 <- tibble()
  
  for (i in 1:10) {
    
    print(glue('{i} / 10 - cross-validation'))
    # Training set
    tmp_train <- tmp_df %>% filter(cross_id != i)
    
    
    tmp_train <- tmp_train %>% 
      # filter(id %in% keep_patho_standard) %>%
      bind_rows(df_1 %>% filter(id %in% tmp_train$id_match))
    
    
    keep_patho_standard <- NoiseFiltersR::IPF(reformulate(vector_features, response = 'clinical'), 
                                              tmp_train, k = input_hyper$k)$cleanData %>%
      filter(clinical == 'pathogenic') %>% 
      pull(id)
    
    tmp_train <- tmp_df %>% filter(id %in% keep_patho_standard) %>%
      bind_rows(df_1 %>% filter(id %in% tmp_train$id_match))
    
    # Testing set
    tmp_test <- tmp_df %>% filter(cross_id == i)
    tmp_test <- tmp_test %>% bind_rows(df_1 %>% filter(id %in% tmp_test$id_match))
    
    qc1 <- qc1 %>% bind_rows(tmp_train %>% count(cross_id, clinical) %>% mutate(tag = 'train', tag2 = i),
                             tmp_test %>% count(cross_id, clinical) %>% mutate(tag = 'test', tag2 = i))
    
    qc2 <- qc2 %>% bind_rows(tmp_train %>% select(length_cnv, clinical) %>% mutate(tag = 'train', tag2 = i),
                             tmp_test %>% select(length_cnv, clinical) %>% mutate(tag = 'test', tag2 = i))
    
    if (!str_detect(model_name, 'bayesian')) {
      
      
      rf_model <- chrom_aware(tmp_train, 
                              tag_variant = tag_variant, 
                              tag_formule = tag_formule,
                              model_name = model_name, 
                              formule_model = formule_model,
                              list_hyper = input_hyper
      )
      
      model_result <- predict_chrom_aware(rf_model, tmp_test)
      
    } else {
      
      bayesian_model <- rtemis_step1(input_tbl = tmp_train,
                                     tag_variant = tag_variant,
                                     vector_features = vector_features,
                                     tag_features = tag_formule,
                                     input_prior = 'hs',
                                     nc = 23,
                                     input_hyper = input_hyper)
      
      model_result <- predict_chrom_aware_rtemis(bayesian_model, tmp_test, tag_variant, tag_formule)
      
    }
    
    cross_raw <- cross_raw %>% bind_rows(model_result[[3]] %>% 
                                           mutate(tag = paste(tag, '-', 'cross-val.', i)))
    
    # ROC
    
    cross_tmp_roc_auc <- round(model_result[[3]] %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate), 3)
    
    cross_roc_auc <- c(cross_roc_auc, cross_tmp_roc_auc)
    
    cross_ss <- cross_ss %>% bind_rows(model_result[[3]] %>% 
                                         roc_curve(clinical, .pred_pathogenic) %>% 
                                         mutate(tag = paste('cross-val', i, '-', cross_tmp_roc_auc)))
    
    # PR
    cross_tmp_pr_auc <- round(model_result[[3]] %>% pr_auc(clinical, .pred_pathogenic) %>% pull(.estimate), 3)
    
    cross_pr_auc <- c(cross_pr_auc, cross_tmp_pr_auc)
    
    cross_pr <- cross_pr %>% bind_rows(model_result[[3]] %>% 
                                         pr_curve(clinical, .pred_pathogenic) %>% 
                                         mutate(tag = paste('cross-val', i, '-', cross_tmp_pr_auc)))
  }
  
  cross_mean_roc_auc <- round(mean(cross_roc_auc),3)
  cross_se_roc_auc <-  round(sd(cross_roc_auc)/sqrt(length(cross_roc_auc)),3)
  
  cross_mean_pr_auc <- round(mean(cross_pr_auc),3)
  cross_se_pr_auc <-  round(sd(cross_pr_auc)/sqrt(length(cross_pr_auc)),3)
  
  merged_roc_auc <- cross_raw %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
  merged_pr_auc <- cross_raw %>% pr_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
  
  merged_roc_curve <- cross_raw %>% roc_curve(clinical, .pred_pathogenic)
  merged_pr_curve <- cross_raw %>% pr_curve(clinical, .pred_pathogenic)
  
  
  tag_title_ss <- paste(paste(tag_variant,tag_formule,model_name, sep = ' - '), cross_mean_roc_auc, '±', cross_se_roc_auc)
  tag_title_pr <- paste(paste(tag_variant,tag_formule,model_name, sep = ' - '), cross_mean_pr_auc, '±', cross_se_pr_auc)
  
  result <- list('cross_mean_roc_auc' = cross_mean_roc_auc, 'cross_se_roc_auc' = cross_se_roc_auc,
                 'cross_mean_pr_auc' = cross_mean_pr_auc, 'cross_se_pr_auc' = cross_se_pr_auc,
                 'cross_raw' = cross_raw, 'cross_ss' = cross_ss, 'cross_pr' = cross_pr,
                 'tag_title_ss' = tag_title_ss, 'tag_title_pr' = tag_title_pr, 'qc1' = qc1, 'qc2' = qc2,
                 'merged_roc_auc' = merged_roc_auc, 'merged_pr_auc' = merged_pr_auc, 
                 'merged_roc_curve' = merged_roc_curve, 'merged_pr_curve' = merged_pr_curve)
  
  return(result)
}



# ------------------------------------------------------------------------------
# PLOT CROSS-VALIDATION
# ------------------------------------------------------------------------------


plot_cross <- function(x, qc = FALSE, single_one = FALSE) {
  
  # x <- rf_del_human
  
  
  p1 <- x[[6]] %>%
    ggplot(aes(1 - specificity, sensitivity)) +
    geom_path(aes(group = tag, color = tag),  show.legend = TRUE) +
    ggtitle(paste0(str_remove(x[[8]], pattern = "(\\d{1}).*"), ' (merged:', x[[12]], ')', '')) +
    theme_roc() +
    theme(legend.title = element_blank())
  
  p2 <- x[[7]] %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = tag, color = tag),  show.legend = TRUE) +
    ggtitle(paste0(str_remove(x[[9]], pattern = "(\\d{1}).*"), ' (merged:', x[[13]], ')')) +
    theme_pr() +
    theme(legend.title = element_blank())
  
  
  
  p3 <-  x[[10]] %>%
    ggplot(aes(factor(tag), n)) +
    geom_col(aes(fill = clinical)) + 
    facet_wrap(vars(tag2)) +
    theme_minimal() 
  
  
  p4 <- x[[11]] %>%
    ggplot(aes(tag,length_cnv)) +
    geom_boxplot(aes(color = clinical)) +
    facet_wrap(vars(tag2)) +
    theme_minimal()
  
  p5 <- x[[14]] %>%
    mutate(tag = 'merged') %>%
    ggplot(aes(1 - specificity, sensitivity)) +
    geom_path(aes(group = tag, color = tag),  show.legend = FALSE) +
    ggtitle(paste0(str_remove(x[[8]], pattern = "(\\d{1}).*"), ' (merged:', x[[12]], '±',x[[2]], ')')) +
    theme_roc() +
    theme(legend.title = element_blank())
  
  p6 <- x[[15]] %>%
    mutate(tag = 'merged') %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = tag, color = tag),  show.legend = FALSE) +
    ggtitle(paste0(str_remove(x[[9]], pattern = "(\\d{1}).*"), ' (merged:', x[[13]], '±',x[[4]], ')')) +
    theme_pr() +
    theme(legend.title = element_blank())
  
  
  
  if (isTRUE(qc)) {
    p1 + p2 + p3 + p4
  } else if (isFALSE(single_one)) {
    p1 + p2
  } else {
    
    p5 + p6
    
  }
  
}

# ------------------------------------------------------------------------------
# PLOT TEST
# ------------------------------------------------------------------------------


plot_test <- function(..., roc = TRUE) {
  
  
  if (isTRUE(roc)) {
    
    p0 <<- bind_rows(...)
    
    p1 <- bind_rows(...) %>% 
      mutate(tag = as.factor(tag))
      ggplot(aes(1-specificity, sensitivity)) +
      geom_path(aes(group = tag, color = tag2),  show.legend = TRUE,  size = 1) +
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
      geom_path(aes(group = tag, color = tag),  show.legend = TRUE,  size = 1) +
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



plot_test_raw <- function(x, tag) {
  library(ggpubr)
  
  
  # x <- merged_df_deletion_test
  
  p3 <-  x %>%
    ggplot(aes(length_cnv)) +
    geom_density(aes(fill = clinical), alpha = 0.4, show.legend = FALSE) +
    theme_minimal() +
    ggtitle(glue('{nrow(x)} {tag}'))
  
  p4 <-   x %>%
    count(clinical) %>%
    ggplot(aes(clinical, n)) +
    geom_col(aes(fill = clinical), color = 'black', show.legend = FALSE) +
    labs(y = 'Frequency') +
    theme_minimal()
  
  p5 <-  x %>%
    ggplot(aes(length_cnv)) +
    geom_histogram(aes(fill = clinical), color = 'black', bins = 30, show.legend = FALSE) +
    facet_wrap(vars(clinical)) +
    theme_minimal()
  
  
  p6 <- x %>%
    summarise(length_cnv = quantile(length_cnv)) %>%
    mutate(name = c('Min', 'Q1', 'Median (Q2)', 'Q3', 'Max')) %>%
    select(name, length_cnv) %>%
    ggtexttable(rows = NULL)
  
  
  p3 + p4 + p5
}


# ------------------------------------------------------------------------------
# CNVs ANNOTATION
# ------------------------------------------------------------------------------

check_cnv_v2 <- function(input_df, mode_reg = FALSE, factor_clinical = TRUE) {
  
  # input_df <- just_test
  
  # get_systems <- function(x) {
  #   
  #   hpo_from_gene <- hpo_genes %>% filter(gene %in% x) %>% pull(hp)
  #   
  #   result_n_systems <- unlist(map(hpo_from_gene, function(x) get_ancestors(hpo_dbs, x))) %>% 
  #     enframe() %>%
  #     filter(value %in% anato_df$name) %>%
  #     count(value) %>% 
  #     nrow()
  #   
  #   return(result_n_systems)
  #   
  # }
  
  
  # MODE_REG = ON
  
  if (isTRUE(mode_reg)) {
    
    mode_reg_enhancers <- bed_intersect(df_enhancers %>%
                                          filter(phast100 >= 0.20 |
                                                   phast46pla >= 0.20 |
                                                   phast46pri >= 0.20), input_df) %>%
      # filter(.overlap > threshold_30_tmp) %>%
      rename(id_tmp = id_tmp.y, gene = gene.x) %>%
      select(id_tmp, gene) %>% 
      distinct()
    
    mode_reg_on <- mode_reg_enhancers
  } else {
    
    mode_reg_on <- tibble()
  }
  
  
  
  number_of_genes <- input_df %>%
    bed_intersect(hgcn_genes %>% select(chrom, start, end, gene)) %>%
    rename(id_tmp = id_tmp.x, gene = gene.y) %>%
    select(id_tmp, gene) %>%
    dplyr::count(id_tmp) %>%
    rename(n_genes = n)
  
  # Blacklist regions
  result_df_blacklist <- bed_coverage(input_df, blacklist_encode) %>% rename(blacklist = .frac) %>% select(id_tmp, blacklist)
  # Recombination rate
  result_recomb <- bed_closest(input_df, recomb, suffix = c('', '.y')) %>% group_by(id_tmp) %>% 
    filter(cm_mb.y == max(cm_mb.y)) %>% slice(1) %>% ungroup() %>% rename(recombination = cm_mb.y) %>% select(id_tmp, recombination)
  # Centromeric distance
  dist_cent <- bed_closest(input_df, region_gaps %>% filter(type == 'centromere')) %>% 
    rename(id_tmp = id_tmp.x, dist_cent = .dist) %>% mutate(dist_cent = abs(dist_cent)/1e6) %>% select(id_tmp, dist_cent)
  # Telomeric distance
  dist_tel <- bed_closest(input_df, region_gaps %>% filter(type == 'telomere')) %>% rename(id_tmp = id_tmp.x, dist_tel = .dist) %>% 
    group_by(id_tmp) %>% mutate(dist_tel = abs(min(dist_tel))/1e6) %>% ungroup() %>% select(id_tmp, dist_tel) %>% distinct()
  
  # SV hotspots
  n_hotspot <- bed_intersect(input_df, hotspot, suffix = c('', '.y')) %>% select(id_tmp) %>% distinct() %>% mutate(hotspot = 1)
  # HARs
  n_hars <- bed_intersect(input_df, hars, suffix = c('', '.y')) %>% select(id_tmp) %>% distinct() %>% mutate(hars = 1)
  # LADs
  n_lads <- bed_intersect(input_df, lads, suffix = c('', '.y')) %>% select(id_tmp) %>% distinct() %>% mutate(lads = 1)
  # Gene density
  gene_density <- bed_intersect(gene_density_tbl, input_df) %>% group_by(id_tmp.y) %>%
    filter(gene_density.x == max(gene_density.x)) %>% slice(1) %>% ungroup() %>% 
    rename(id_tmp = id_tmp.y, gene_density = gene_density.x) %>% select(id_tmp, gene_density)
  # PubMed - Deletions and duplications
  pubmed_total <- bed_intersect(input_df, pubmed_df, suffix = c('', '.y')) %>% group_by(id_tmp) %>% 
    filter(hits_del.y == max(hits_del.y)) %>% slice(1) %>% ungroup() %>% 
    rename(hits_del = hits_del.y, hits_dup = hits_dup.y) %>% select(id_tmp, hits_del, hits_dup)
  # Ensembl - CTCF
  result_ctcf <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'CTCF_binding_site')) %>% select(id_tmp, .frac) %>%
    rename(ctcf = .frac)
  # Ensembl - Enhancers
  result_enhancer <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'enhancer')) %>% select(id_tmp, .frac) %>%
    rename(enhancer = .frac)
  # Ensembl - Open chromatin
  result_open <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'open_chromatin_region')) %>% select(id_tmp, .frac) %>%
    rename(open = .frac)
  # Ensembl - Promoter
  result_promoter <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'promoter')) %>% select(id_tmp, .frac) %>%
    rename(promoter = .frac)
  # Ensembl - Promoter flank
  result_promoterflank <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'promoter_flanking_region')) %>% 
    select(id_tmp, .frac) %>% rename(promoterflank = .frac)
  
  # Ensembl - TFBS
  result_tfbs <- bed_coverage(input_df, ensembl_reg %>% filter(type == 'TF_binding_site')) %>% select(id_tmp, .frac) %>%
    rename(tfbs = .frac)
  # UCNE regions
  result_ucne <- bed_coverage(input_df, ucne) %>% select(id_tmp, .frac) %>% rename(ucne = .frac)
  result_ucne <- result_ucne %>% mutate(ucne = if_else(ucne == 0, 0, 1))
  # Clinvar regions
  result_df_clinvar <- input_df %>% 
    bed_intersect(clinvar_variants %>% 
                    filter(variant_class %in% c('indel', 'single nucleotide variant') & clinical == "pathogenic"),
                  suffix = c('', 'y')) %>%
    select(id_tmp) %>%
    distinct() %>%
    mutate(clinvar = 1)
  
  # GWAS variants
  result_df_gwas <- input_df %>% bed_intersect(gwas_variants %>% filter(INTERGENIC == "No"), suffix = c('', 'y')) %>%
    select(id_tmp) %>%
    distinct() %>%
    mutate(gwas = 1)
  
  # TADs
  
  result_tads <- input_df %>%
    bed_intersect(tad %>%
                    pivot_longer(-c(id, chrom), names_to = 'coord', values_to = 'start') %>%
                    mutate(end = start)) %>%
    count(id_tmp.x, id.y) %>%
    filter(n == 1) %>%
    select(-id.y) %>%
    distinct() %>%
    rename(id_tmp = id_tmp.x, tads = n)
  
  # CADD maximum
  
  result_cadd <- input_df %>%
    bed_intersect(cadd_max %>% mutate(chrom = as.character(chrom))) %>%
    select(id_tmp.x, max_cadd.y) %>%
    rename(id_tmp = id_tmp.x, max_cadd = max_cadd.y) %>%
    group_by(id_tmp) %>%
    summarise(max_cadd = max(max_cadd))
  
  # GERP maximum
  
  result_gerp <- input_df %>%
    bed_intersect(gerp_max) %>%
    select(id_tmp.x, max_gerp.y) %>%
    rename(id_tmp = id_tmp.x, max_gerp = max_gerp.y) %>%
    group_by(id_tmp) %>%
    summarise(max_gerp = max(max_gerp))
  
  # Remot-GW
  
  result_obs_exp <- input_df %>%
    bed_intersect(remot_cnvscore) %>%
    select(id_tmp.x, obs_exp.y) %>%
    rename(id_tmp = id_tmp.x, obs_exp = obs_exp.y) %>%
    group_by(id_tmp) %>%
    summarise(max_obs_exp = max(obs_exp))

  region_level <- input_df %>%
    select(id_tmp) %>%
    left_join(result_df_blacklist, by = 'id_tmp') %>%
    left_join(result_recomb, by = 'id_tmp') %>%
    left_join(dist_cent, by = 'id_tmp') %>%
    left_join(dist_tel, by = 'id_tmp') %>%
    left_join(n_hotspot, by = 'id_tmp') %>%
    left_join(n_hars, by = 'id_tmp') %>%
    left_join(n_lads, by = 'id_tmp') %>%
    left_join(gene_density, by = 'id_tmp') %>%
    left_join(pubmed_total, by = 'id_tmp') %>%
    left_join(result_ctcf, by = 'id_tmp') %>%
    left_join(result_enhancer, by = 'id_tmp') %>%
    left_join(result_open, by = 'id_tmp') %>%
    left_join(result_promoter, by = 'id_tmp') %>%
    left_join(result_promoterflank, by = 'id_tmp') %>%
    left_join(result_tfbs, by = 'id_tmp') %>%
    left_join(result_ucne, by = 'id_tmp') %>%
    left_join(result_df_clinvar, by = 'id_tmp') %>%
    left_join(result_df_gwas, by = 'id_tmp') %>%
    left_join(result_tads, by = 'id_tmp') %>%
    left_join(result_cadd, by = 'id_tmp') %>%
    left_join(result_gerp, by = 'id_tmp') %>%
    left_join(result_obs_exp, by = 'id_tmp') %>%
    left_join(number_of_genes, by = 'id_tmp') %>%
    # mutate(gene_density = ifelse(is.na(gene_density), 2, gene_density)) %>%
    # mutate(dist_tel = ifelse(is.na(dist_tel), 27, dist_tel)) %>%
    # mutate(tads = ifelse(is.na(tads), 0, tads)) %>%
    # mutate(max_cadd = ifelse(is.na(max_cadd), 0, max_cadd)) %>%
    # mutate(max_gerp = ifelse(is.na(max_gerp), 0, max_gerp)) %>%
    # mutate(max_obs_exp = ifelse(is.na(max_obs_exp), 0, max_obs_exp)) %>%
    # mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes)) %>%
    # mutate(lads = ifelse(is.na(lads), 0, lads)) %>%
    # mutate(hars = ifelse(is.na(hars), 0, hars)) %>%
    # mutate(hotspot = ifelse(is.na(hotspot), 0, hotspot)) %>%
    # mutate(clinvar = if_else(is.na(clinvar), 0, clinvar)) %>%
    # mutate(gwas = if_else(is.na(gwas), 0, gwas)) %>%
    # mutate(clinical = as.character(clinical)) %>%
    mutate(across(everything(), ~replace_na(.x, 0)))
    # select(-c('chrom', 'start', 'end', 'variant_class', 'source', 'clinical', 'length_cnv'))

  # GENE-LEVEL RESULTS
  
  gene_level <- input_df %>%
    bed_intersect(hgcn_genes %>% select(chrom, start, end, gene)) %>%
    rename(id_tmp = id_tmp.x, gene = gene.y) %>%
    select(id_tmp, gene) %>%
    bind_rows(mode_reg_on) %>%
    distinct() %>%
    mutate(disease = if_else(gene %in% (hgcn_genes %>% filter(disease == 'Yes') %>% pull(gene)), 1, 0)) %>% 
    mutate(omim = if_else(gene %in% (hgcn_genes %>% filter(omim == 'Yes') %>% pull(gene)), 1, 0)) %>% 
    mutate(haplo = if_else(gene %in% (haplo_triplo_genes %>% filter(haplo == 'yes') %>% pull(gene)), 1, 0)) %>% 
    mutate(triplo = if_else(gene %in% (haplo_triplo_genes %>% filter(triplo == 'yes') %>% pull(gene)), 1, 0)) %>%
    mutate(mouse_embryo = if_else(gene %in% (mgi %>% filter(str_detect(pheno, 'MP:0010768|MP:0005380')) %>% pull(gene)), 1, 0)) %>%
    left_join(hgcn_genes %>% select(gene, pLI) %>% rename(pli = pLI), by = 'gene') %>%
    left_join(loeuf_score, by = 'gene') %>%
    left_join(eds, by = 'gene') %>%
    left_join(pnull, by = 'gene') %>%
    left_join(hgcn_genes %>% select(gene, ccr), by = 'gene') %>%
    left_join(hgcn_genes %>% select(gene, hi), by = 'gene') %>%
    left_join(expression_features, by = 'gene') %>%
    left_join(genes_promoter, by = 'gene') %>%
    left_join(crispr_score_df, by = 'gene') %>%
    left_join(para_genes %>% rename(paralogous_genes = n), by = 'gene') %>%
    left_join(string_db, by = 'gene') %>%
    mutate(gene_hpo = if_else(gene %in% hpo_genes$gene, 1, 0)) %>% 
    mutate(prot_complex = if_else(gene %in% prot_complex$gene, 1, 0)) %>% 
    mutate(prot_complex_nohuman = if_else(gene %in% hu_map$gene, 1, 0)) %>% 
    mutate(ohnolog = if_else(gene %in% hgcn_genes[hgcn_genes$ohnolog == 'Yes',]$gene, 1, 0)) %>%
    mutate(essent_cl = if_else(gene %in% (hgcn_genes %>% filter(str_detect(fusil, 'CL')) %>% pull(gene)), 1, 0)) %>%
    mutate(essent_dl = if_else(gene %in% (hgcn_genes %>% filter(str_detect(fusil, 'DL')) %>% pull(gene)), 1, 0)) %>%
    mutate(tf = if_else(gene %in% tf_genes, 1, 0)) %>%
    select(-gene) %>%
    # mutate(disease = if_else(is.na(disease), 0, disease)) %>%
    # mutate(haplo = if_else(is.na(haplo), 0, haplo)) %>%
    # mutate(triplo = if_else(is.na(triplo), 0, triplo)) %>%
    # mutate(mouse_embryo = if_else(is.na(mouse_embryo), 0, mouse_embryo)) %>%
    # mutate(tf = ifelse(is.na(tf), 0, tf)) %>%
    # mutate(gene_hpo = ifelse(is.na(gene_hpo), 0, gene_hpo)) %>%
    # mutate(prot_complex = ifelse(is.na(prot_complex), 0, prot_complex)) %>%
    # mutate(ohnolog = ifelse(is.na(ohnolog), 0, ohnolog)) %>%
    # mutate(paralogous_genes = ifelse(is.na(paralogous_genes), 0, paralogous_genes)) %>%
    # mutate(essent_cl = ifelse(is.na(essent_cl), 0, essent_cl)) %>%
    # mutate(essent_dl = ifelse(is.na(essent_dl), 0, essent_dl)) %>%
    # mutate(pli = ifelse(is.na(pli), 0, pli)) %>%
    # mutate(loeuf = ifelse(is.na(loeuf), 0, loeuf)) %>%
    # mutate(eds = ifelse(is.na(eds), 0, eds)) %>%
    # mutate(pnull = ifelse(is.na(pnull), 0, pnull)) %>%
    # mutate(ccr = ifelse(is.na(ccr), 0, ccr)) %>%
    # mutate(hi = ifelse(is.na(hi), 0, hi)) %>%
    # mutate(mean_phast = ifelse(is.na(mean_phast), 0, mean_phast)) %>%
    # mutate(cpg_density = ifelse(is.na(cpg_density), 0, cpg_density)) %>%
    # mutate(crispr_score = ifelse(is.na(crispr_score), 0, crispr_score)) %>%
    # mutate(mean_expression = ifelse(is.na(mean_expression), 0, mean_expression)) %>%
    # mutate(min_expression = ifelse(is.na(min_expression), 0, min_expression)) %>%
    # mutate(degree_to_triplo = ifelse(is.na(degree_to_triplo), 0, degree_to_triplo)) %>%
    # mutate(degree_to_haplo = ifelse(is.na(degree_to_haplo), 0, degree_to_haplo)) %>%
    # mutate(degree = ifelse(is.na(degree), 0, degree)) %>%
    # mutate(page_rank = ifelse(is.na(page_rank), 0, page_rank)) %>%
    mutate(haplo_short = ifelse(is.na(haplo_short), 4, haplo_short)) %>%
    mutate(triplo_short = ifelse(is.na(triplo_short), 4, triplo_short)) %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    # Aggregation step
    group_by(id_tmp) %>%
    mutate(across(where(is.numeric) & !c(haplo_short, triplo_short), max)) %>%
    mutate(across(c(haplo_short, triplo_short), min)) %>%
    ungroup() %>%
    distinct()
  
  
  # 2 imputated features
  # haplo_short -> 4 (median)
  # triplo_short -> 4 (median)
  # REMOVED dist_tel -> 27 (median)
  # REMOVED gene_density -> 2 (median)
  
  
  cnvs_annotated <- input_df %>%
    mutate(clinical = as.character(clinical)) %>%
    left_join(gene_level, by = 'id_tmp') %>%
    left_join(region_level, by = 'id_tmp') %>%
    mutate(haplo_short = ifelse(is.na(haplo_short), 4, haplo_short)) %>%
    mutate(triplo_short = ifelse(is.na(triplo_short), 4, triplo_short))
  
  # skip_chrom <- colnames(cnvs_annotated)[-1]
  
  cnvs_annotated <- cnvs_annotated %>%
    mutate(across(!where(is.character), ~replace_na(.x, 0))) %>%
    mutate(clinical = as.factor(clinical)) %>%
    rename(type_variant = variant_class, id = id_tmp)
  
  # cnvs_annotated <- cnvs_annotated %>%
  #   mutate(across(all_of(skip_chrom), ~replace_na(.x, 0))) %>%
  #   mutate(clinical = as.factor(clinical)) %>%
  #   rename(type_variant = variant_class, id = id_tmp)

  # mutate(clinical = as.factor(clinical)) %>%
  # mutate(clinical = fct_relevel(clinical, 'pathogenic', 'benign'))
  # map(~ sum(is.na(.x)))
  
  
  if (isTRUE(factor_clinical)) {
    
    cnvs_annotated <- cnvs_annotated %>%
      mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign')))

  }

  return(cnvs_annotated)
  
}


# # ------------------------------------------------------------------------------
# # REMOVE NOISY OBSERVATIONS
# # ------------------------------------------------------------------------------
# 
# 
# test_removed_ids_def <- tibble()
# for (i in 1:20) {
#   print(i)
#   test_kept_ids <- NoiseFiltersR::IPF(human_no_control, output_df_deletion, k = 3)$cleanData %>%
#     filter(clinical == 'pathogenic') %>% pull(id)
#   
#   test_removed_ids <- output_df_deletion %>% filter(clinical == 'pathogenic' & (!id %in% test_kept_ids)) %>% 
#     mutate(iter_n = i) %>%
#     select(iter_n, id)
#   
#   test_removed_ids_def <- test_removed_ids_def %>% bind_rows(test_removed_ids)
# }


# ------------------------------------------------------------------------------
# CALCULATE PATHOGENICITY DIFF. SCORES
# ------------------------------------------------------------------------------

get_mode_noncoding <- function(x, model_used, input_tag) {
  
  # x <- output_clinvar_deletion
  
  input_mode_reg_off <- check_cnv_v2(x)
  
  input_mode_reg_on <- check_cnv_v2(x, mode_reg = TRUE)
  
  result_mode_reg_off <- predict_chrom_aware_rtemis(model_used, input_mode_reg_off,
                                                    'deletion', 'human_control')
  result_mode_reg_on <-predict_chrom_aware_rtemis(model_used, input_mode_reg_on,
                                                  'deletion', 'human_control')
  
  tmp_p1 <- result_mode_reg_off[[3]] %>%
    rename(patho_off = .pred_pathogenic) %>% 
    select(patho_off, id) %>%
    left_join(result_mode_reg_on[[3]] %>% rename(patho_on = .pred_pathogenic) %>% select(patho_on, id), by = 'id') %>%
    mutate(diff = patho_on - patho_off) %>%
    mutate(tag = case_when(
      diff == 0 ~ 'no difference',
      diff > 0 ~ 'higher pathogenicity score',
      diff < 0 ~ 'lower pathogenicity score'
    )) %>%
    count(tag) %>%
    mutate(perc = n / sum(n)) %>%
    ggplot(aes(reorder(tag, -perc), perc)) +
    geom_col(aes(fill = tag), color = 'black') +
    scale_y_continuous(label = percent) +
    geom_text(aes(label = paste0(100*round(perc, 4), '%')), size = 5,vjust = 0) +
    theme_minimal(base_size = 15) +
    labs(x = 'Category', y = 'Percentage (%)') +
    labs(title = input_tag)
  
  
  tmp_p2 <- result_mode_reg_off[[3]] %>%
    rename(patho_off = .pred_pathogenic) %>% 
    select(patho_off, id) %>%
    left_join(result_mode_reg_on[[3]] %>% rename(patho_on = .pred_pathogenic) %>% select(patho_on, id), by = 'id') %>%
    mutate(diff = patho_on - patho_off) %>%
    mutate(tag = case_when(
      patho_off < 0.5 & patho_on > 0.5  ~ 'change to pathogenic',
      patho_off < 0.5 & patho_on < 0.5 ~ 'no change',
      patho_off > 0.5 & patho_on > 0.5 ~ 'no change',
      patho_off > 0.5 & patho_on < 0.5 ~ 'change to non-pathogenic'
    )) %>%
    count(tag) %>%
    mutate(perc = n / sum(n)) %>%
    ggplot(aes(reorder(tag, -perc), perc)) +
    geom_col(aes(fill = tag), color = 'black') +
    scale_y_continuous(label = percent) +
    geom_text(aes(label = paste0(100*round(perc, 4), '%')), size = 5,vjust = 0) +
    theme_minimal(base_size = 15) +
    labs(x = 'Category', y = 'Percentage (%)') +
    labs(title = input_tag)
  
  tmp_p1 + tmp_p2
  
}


# ------------------------------------------------------------------------------
# Simple function for the calculation of the AUC
# ------------------------------------------------------------------------------

just_auc <- function(x) {
  
  tmp_auc <- x %>% 
    roc_auc(clinical, .pred_pathogenic) %>%
    pull(.estimate) %>%
    round(2)
  
  tmp_pr <- x %>% 
    pr_auc(clinical, .pred_pathogenic) %>%
    pull(.estimate) %>%
    round(2)
  
  return(c(tmp_auc, tmp_pr))
}


# ------------------------------------------------------------------------------
# WRITE COMPETITOR'S RESULTS
# ------------------------------------------------------------------------------


write_competitors <- function(to_write, which_class = 'DEL') {
  
  from_hg19_to_hg38 = import.chain('/data-cbl/liftover/hg19ToHg38.over.chain')
  
  to_write_tmp <- to_write %>% 
    select(chrom, start, end, id) %>% 
    GRanges()
  
  seqlevelsStyle(to_write_tmp) = "UCSC"  # necessary
  to_write_hg38 = liftOver(to_write_tmp, from_hg19_to_hg38) %>% 
    as_tibble() %>%
    rename(chrom = seqnames) %>%
    mutate(type = which_class)

# ------------------------------------------------------------------------------
# CADD-SV
# ------------------------------------------------------------------------------

  
write_tsv(to_write_hg38 %>%
            mutate(start = start - 1) %>%
            mutate(type = which_class) %>%
            select(chrom, start, end, type), 'rival_cnvscore/CADD-SV/input/id_input.bed', col_names = FALSE)
  



# ------------------------------------------------------------------------------
# AnnotSV
# ------------------------------------------------------------------------------


write_tsv(to_write %>%
            mutate(start = start - 1) %>%
            mutate(type = which_class) %>%
            select(chrom, start, end, type), 'rival_cnvscore/annot_sv/input.bed', col_names = FALSE)

# ------------------------------------------------------------------------------
# TADA
# conda activate tada
# cd /data-cbl/frequena_data/cnvscore/rival_cnvscore/TADA/
# predict_variants -c config_del_comparison.yml -o comp_folder/del/ && predict_variants -c config_dup_comparison.yml -o comp_folder/dup/

# ------------------------------------------------------------------------------


write_tsv(to_write %>%
            mutate(start = start - 1) %>%
            select(chrom, start, end), 'rival_cnvscore/TADA/comp_folder/both/input.bed', col_names = FALSE)




# ------------------------------------------------------------------------------
# STRVCTRUE
# conda create --prefix=/data-cbl/frequena_data/conda_frequena/structure python=3.7.2
# conda install scikit-learn=0.21.3
# conda activate structure
# cd /data-cbl/frequena_data/cnvscore/rival_cnvscore/structure/StrVCTVRE-v.1.6/
# python StrVCTVRE.py -i ../input_del.bed -o ../output_del.bed -f bed && python StrVCTVRE.py -i ../input_dup.bed -o ../output_dup.bed -f bed
# 
# ------------------------------------------------------------------------------

# 0-BASED - 1-BASED




  to_write_hg38 %>%
  mutate(start = start - 1) %>%
  select(chrom, start, end, type) %>%
  write_tsv('rival_cnvscore/structure/input.bed', col_names = FALSE)




# ------------------------------------------------------------------------------
# TAD FUSION
# conda activate tad_fusion
# cd /data-cbl/frequena_data/cnvscore/rival_cnvscore/tad_fusion/
# src/cal_tad_fusion_score -md Model/GM_Rao_5kb -f input.tsv -mnl 10000 -mxl 5000000 -w 100 -d 0.06 -o output_del.tsv
# ------------------------------------------------------------------------------

write_tsv(to_write %>%
            mutate(chrom = paste0('chr', chrom)) %>%
            select(chrom, start, end), 'rival_cnvscore/tad_fusion/input.tsv', col_names = FALSE)

# ------------------------------------------------------------------------------
# X-CNV
# git clone https://github.com/kbvstmd/XCNV.git
# cd XCNV
# sh Install.sh
# changed bin/XCNV to #!/data-cbl/frequena_data/conda_frequena/xcnv/bin/Rscript
# conda activate xcnv
# cd /data-cbl/frequena_data/cnvscore/rival_cnvscore/x_cnv
# ./bin/XCNV input_del.bed && ./bin/XCNV input_dup.bed
# ------------------------------------------------------------------------------

write_tsv(to_write %>%
            mutate(start = start - 1) %>%
            mutate(type = if_else(which_class == 'DEL','loss', 'gain')) %>%
            select(chrom, start, end, type), 'rival_cnvscore/x_cnv/input.bed', col_names = FALSE)


# ------------------------------------------------------------------------------
# ClassifyCNV
# conda activate classifycnv
# cd /data-cbl/frequena_data/cnvscore/rival_cnvscore/classifycnv
# rm -rf results_del/Intermediate_files
# python3 ClassifyCNV.py --infile input_del.bed  --GenomeBuild hg19 --precise --outdir result_del && python3 ClassifyCNV.py --infile input_dup.bed  --GenomeBuild hg19 --precise --outdir result_dup
# 
# ------------------------------------------------------------------------------

write_tsv(to_write %>%
            mutate(start = start - 1) %>%
            mutate(type = which_class) %>%
            mutate(chrom = paste0('chr', chrom)) %>%
            select(chrom, start, end, type), 'rival_cnvscore/classifycnv/input.bed', col_names = FALSE)

# ------------------------------------------------------------------------------
# JARVIS
# find . -name 'jarvis*' | parallel -j 10 tabix -p bed {}
# conda activate jarvis
# rm output_*.bed
# parallel --jobs 22 'bedtools intersect -sorted -wb -a input_{}_del.bed -b /data-cbl/jarvis/jarvis.chr{}.both-features.sorted.bed.bgz > output_{}_del.bed ' ::: {1..22}
# parallel --jobs 22 'bedtools intersect -sorted -wb -a input_{}_dup.bed -b /data-cbl/jarvis/jarvis.chr{}.both-features.sorted.bed.bgz > output_{}_dup.bed ' ::: {1..22}
# THERE ARE MISSING CNVs - it seems they do not appear on the scores - ask antonio
# Some Jarvis values are Inf - I had to remove them
# ------------------------------------------------------------------------------



1:22 %>% map(function(x) {
  
  to_write %>% 
    filter(chrom == x) %>%
    arrange(start) %>%
    mutate(chrom = paste0('chr', x)) %>%
    select(chrom, start, end, id) %>%
    write_tsv(glue('rival_cnvscore/jarvis/input_{x}.bed'), col_names = FALSE)
  
})

# 1:22 %>% map(function(x) {
#   
#   output_clinvar_duplication %>% 
#     filter(chrom == x) %>%
#     arrange(start) %>%
#     mutate(chrom = paste0('chr', x)) %>%
#     select(chrom, start, end, id) %>%
#     write_tsv(glue('rival_cnvscore/jarvis/input_{x}_dup.bed'), col_names = FALSE)
#   
# })


# ------------------------------------------------------------------------------
# gwRVIS
# conda activate jarvis
# parallel --jobs 22 'bedtools intersect -sorted -wb -a input_{}_del.bed -b /data-cbl/gwrvis/gwrvis_single_nt.chr{}.bed.gz > output_{}_del.bed ' ::: {1..22}
# parallel --jobs 22 'bedtools intersect -sorted -wb -a input_{}_dup.bed -b /data-cbl/gwrvis/gwrvis_single_nt.chr{}.bed.gz > output_{}_dup.bed ' ::: {1..22}
# ------------------------------------------------------------------------------

1:22 %>% map(function(x) {
  
  to_write %>% 
    filter(chrom == x) %>%
    arrange(start) %>%
    mutate(chrom = paste0('chr', x)) %>%
    select(chrom, start, end, id) %>%
    write_tsv(glue('rival_cnvscore/gwrvis/input_{x}.bed'), col_names = FALSE)
  
})

}


# ------------------------------------------------------------------------------
# READ COMPETITOR'S RESULTS
# ------------------------------------------------------------------------------

run_competitors <- function(x) {
  
  file.remove('/data-cbl/frequena_data/cnvscore/run_rivals_check.txt')
  
  system(glue('/data-cbl/frequena_data/cnvscore/run_rivals.sh {x}'))
  
  read_tsv('/data-cbl/frequena_data/cnvscore/run_rivals_check.txt')
  
  
}

# ------------------------------------------------------------------------------
# READ COMPETITOR'S RESULTS
# ------------------------------------------------------------------------------

read_competitors2 <- function(label_title,x, list_models) {
  
  
  # x <- input_df
  # list_models <- list_results
  # label_title <- 'prueba - borrar'
  

  
  # Convert from hg38 to hg19
  
  from_hg19_to_hg38 = import.chain('/data-cbl/liftover/hg19ToHg38.over.chain')
  
  to_write_tmp <- x %>% 
    select(chrom, start, end, id) %>% 
    GRanges()

  seqlevelsStyle(to_write_tmp) = "UCSC"  # necessary
  to_write_hg38 = liftOver(to_write_tmp, from_hg19_to_hg38) %>% 
    as_tibble() %>%
    rename(chrom = seqnames) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    select(chrom, start, end, id) %>%
    left_join(x %>% select(id, clinical), by = 'id')
  
  ## STRUCTRUE ------------------------
  
  
  result_structure_manolo_deletion <-  read_tsv('rival_cnvscore/structure/output.bed', col_names = c('chrom', 'start', 'end'))
  
  result_structure_manolo_deletion <- result_structure_manolo_deletion %>% 
    mutate(start = start + 1) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    rename(.pred_pathogenic = X5) %>%
    mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
    mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
    left_join(to_write_hg38, by = c('chrom', 'start', 'end')) %>%
    select(.pred_pathogenic, id)
  
  
  ## CADD-SV ------------------------
  
  
  result_cadd_manolo_deletion <- read_tsv('rival_cnvscore/CADD-SV/output/input2_score.bed', skip = 1, col_names = c('chrom', 'start', 'end', 'class','name', '.pred_pathogenic', letters))
  
  result_cadd_manolo_deletion <- result_cadd_manolo_deletion %>% 
    mutate(chrom = as.character(chrom)) %>%
    mutate(start = start + 1) %>%
    select(chrom, start, end, .pred_pathogenic) %>%
    left_join(to_write_hg38 , by = c('chrom', 'start', 'end')) %>%
    na.omit() %>%
    select(id, .pred_pathogenic)
  
  ## TADA------------------------
  
  result_tada_manolo_deletion <-
    read_tsv('rival_cnvscore/TADA/comp_folder/both/Annotated_Predicted_PATHOGENIC.csv') %>%
    rename(chrom = CHR, start = START, end = END, .pred_pathogenic = `Pathogenicity Score`) %>%
    mutate(chrom = str_remove(chrom, 'chr'), start = start + 1) %>%
    # select(chrom, start, end, .pred_pathogenic) %>%
    left_join(x %>% select(chrom, start, end, clinical, id), by = c('chrom', 'start', 'end')) %>%
    na.omit() %>%
    select(id, .pred_pathogenic)
  
  ## TAD FUSION ------------------------
  
  result_tadfusion_manolo_deletion <-  read_tsv('rival_cnvscore/tad_fusion/output.tsv', col_names = c('chrom', 'start', 'end','score')) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    mutate(score = (score - min(score)) / (max(score) - min(score))) %>%
    right_join(x %>% select(chrom, start, end, clinical, id), by = c('chrom', 'start', 'end')) %>%
    mutate(score = if_else(is.na(score), 0, score)) %>%
    rename(.pred_pathogenic = score) %>%
    select(id, .pred_pathogenic)
  
  ## Classify CNV -----------------
  
  
  result_classifycnv_manolo_deletion <- read_tsv('rival_cnvscore/classifycnv/result/Scoresheet.txt') %>%
    select(Chromosome, Start, End, `Total score`) %>%
    rename(score = `Total score`, chrom = Chromosome, start = Start, end = End) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    mutate(score = (score - min(score)) / (max(score) - min(score))) %>%
    rename(.pred_pathogenic = score) %>%
    mutate(start = start + 1) %>%
    left_join(x %>% select(chrom, start, end, clinical, id), by = c('chrom', 'start', 'end')) %>%
    select(id, .pred_pathogenic)

  
  ## X-CNV------------------------
  
  result_xcnv_manolo_deletion <- read_csv('rival_cnvscore/x_cnv/input.output.csv') %>%
    rename(chrom = Chr, start = Start, end = End, .pred_pathogenic = MVP_score) %>%
    mutate(start = start + 1) %>%
    select(chrom, start, end, .pred_pathogenic) %>%
    mutate(chrom = as.character(chrom)) %>%
    left_join(x %>% select(chrom, start, end, clinical, id), by = c('chrom', 'start', 'end')) %>%
    select(id, .pred_pathogenic)
  
  ## JARVIS------------------------
  
  
  result_jarvis_manolo_deletion <- 1:22 %>% map_dfr(function(x) {
    
    tmp_df <-  read_tsv(glue('rival_cnvscore/jarvis/output_{x}.bed'), col_names = FALSE)
    
    if (nrow(tmp_df) == 0) {
      
      return(tibble())
      
    } else {
      
      tmp_df <-  tmp_df %>%
        rename(id = X4, score = X8) %>%
        select(id, score) %>%
        mutate(id = as.character(id))
      
    }
    
    return(tmp_df)
    
  })
  
  
  result_jarvis_manolo_deletion<- result_jarvis_manolo_deletion %>%
    filter(!is.infinite(score)) %>%
    na.omit() %>%
    group_by(id) %>%
    summarise(median_score = median(score),
              mean_score = mean(score),
              q1_score = quantile(score, 0.25),
              q3_score = quantile(score, 0.75)
    ) %>%
    pivot_longer(cols = -id) %>%
    group_by(id) %>%
    summarise(.pred_pathogenic = mean(value)) %>%
    na.omit() 
  
  ## GWRVIS------------------------
  
  
  result_gwrvis_manolo_deletion <- 1:22 %>% map_dfr(function(x) {
    
    print(x)
    
    tmp_df <-  read_tsv(glue('rival_cnvscore/gwrvis/output_{x}.bed'), col_names = FALSE)
    
    if (nrow(tmp_df) != 0) {
      
      tmp_df <-  tmp_df %>%
        rename(id = X4, score = X8) %>%
        select(id, score) %>%
        mutate(id = as.character(id))
      
    } else {
      
      return(tibble())
      
    }
    
    return(tmp_df)
    
    
  })
  
  
  result_gwrvis_manolo_deletion<- result_gwrvis_manolo_deletion %>%
    filter(!is.infinite(score)) %>%
    na.omit() %>%
    group_by(id) %>%
    summarise(median_score = median(score),
              mean_score = mean(score),
              q1_score = quantile(score, 0.25),
              q3_score = quantile(score, 0.75)
    ) %>%
    pivot_longer(cols = -id) %>%
    group_by(id) %>%
    summarise(.pred_pathogenic = 1 - mean(value)) %>%
    na.omit()
  
  predictions_tables <- bind_rows(
    result_cadd_manolo_deletion %>% mutate(tag = 'cadd'),
    result_tada_manolo_deletion %>% mutate(tag = 'tada'),
    result_structure_manolo_deletion %>% mutate(tag = 'strvctvre'),
    result_jarvis_manolo_deletion %>% mutate(tag = 'jarvis') %>% mutate(id = as.numeric(id)),
    result_gwrvis_manolo_deletion %>% mutate(tag = 'gwrvis') %>% mutate(id = as.numeric(id)),
    result_classifycnv_manolo_deletion %>% mutate(tag = 'classifycnv') %>% mutate(id = as.numeric(id)),
    result_xcnv_manolo_deletion %>% mutate(tag = 'xcnv') %>% mutate(id = as.numeric(id))
  )
  
  report_own_and_competitors <- bind_rows(list_models, 
                                          predictions_tables)
  
  return(report_own_and_competitors)

}
  


# ------------------------------------------------------------------------------
# READ COMPETITOR'S RESULTS
# ------------------------------------------------------------------------------

read_competitors <- function(label_title,x, list_models) {
  
  # x <- output_decipher_duplication
  # 
  # list_models <- list(result_bay_x_bias_dup_length, result_bay_x_bias_dup_n_genes,
  #                                                                 result_bay_x_bias_dup_disease_omim, result_bay_x_dup_nohuman, 
  #                                                                 result_bay_x_dup_human,
  #                                                                 result_rf_x_dup_nohuman,
  #                                                                 result_rf_x_dup_human)
  
  roc_curves_table <- tibble()
  pr_curves_table <- tibble()
  result_main_table <- tibble()

  if (length(list_models) > 0) {
    
    for(i in 1:length(list_models)) {
      
      roc_curves_table <- roc_curves_table %>% bind_rows(list_models[[i]]$tmp_roc_curve %>% mutate(tool = paste0('CNVscore', i)))
      pr_curves_table <- pr_curves_table %>% bind_rows(list_models[[i]]$tmp_pr_curve  %>% mutate(tool = paste0('CNVscore', i)))
      
      
      result_main_table <- result_main_table %>% bind_rows(tibble(tool = paste0('CNVscore', i), 
                                                                  read_rows = nrow(list_models[[i]]$tmp_predicted),
                                                                  roc_auc = list_models[[i]]$tmp_predicted  %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3),
                                                                  pr_auc = list_models[[i]]$tmp_predicted  %>% pr_auc(clinical, .pred_pathogenic) %>% pull(.estimate) %>% round(3)
      ))
    }
  }
  
  
  # For CADD-SV and STRVCTURE (hg38)

  
  from_hg19_to_hg38 = import.chain('/data-cbl/liftover/hg19ToHg38.over.chain')
  
  to_write_tmp <- x %>% 
    select(chrom, start, end, id) %>% 
    GRanges()
  
  seqlevelsStyle(to_write_tmp) = "UCSC"  # necessary
  to_write_hg38 = liftOver(to_write_tmp, from_hg19_to_hg38) %>% 
    as_tibble() %>%
    rename(chrom = seqnames) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    select(chrom, start, end, id) %>%
    left_join(x %>% select(id, clinical), by = 'id') %>%
    select(-id)
    
  
  result_cadd_manolo_deletion <- read_tsv('rival_cnvscore/CADD-SV/output/input2_score.bed', skip = 1, col_names = c('chrom', 'start', 'end', 'class','name', '.pred_pathogenic', letters))
  
  result_cadd_manolo_deletion <- result_cadd_manolo_deletion %>% 
    mutate(chrom = as.character(chrom)) %>%
    mutate(start = start + 1) %>%
    select(chrom, start, end, .pred_pathogenic) %>%
    left_join(to_write_hg38 , by = c('chrom', 'start', 'end')) %>%
    na.omit()
  
  cadd_del_manolo_auc <- just_auc(result_cadd_manolo_deletion)
  
  
  
  ## TADA------------------------
  
  result_tada_manolo_deletion <-
    read_tsv('rival_cnvscore/TADA/comp_folder/both/Annotated_Predicted_PATHOGENIC.csv') %>%
    rename(chrom = CHR, start = START, end = END, .pred_pathogenic = `Pathogenicity Score`) %>%
    mutate(chrom = str_remove(chrom, 'chr'), start = start + 1) %>%
    select(chrom, start, end, .pred_pathogenic) %>%
    left_join(x %>% select(chrom, start, end, clinical), by = c('chrom', 'start', 'end')) %>%
    na.omit()
  
  tada_del_manolo_auc <- just_auc(result_tada_manolo_deletion)
  
  
  ## TAD FUSION ------------------------
  
  result_tadfusion_manolo_deletion <-  read_tsv('rival_cnvscore/tad_fusion/output.tsv', col_names = c('chrom', 'start', 'end','score')) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    mutate(score = (score - min(score)) / (max(score) - min(score))) %>%
    right_join(x %>% select(chrom, start, end, clinical), by = c('chrom', 'start', 'end')) %>%
    mutate(score = if_else(is.na(score), 0, score)) %>%
    rename(.pred_pathogenic = score)
  
  tadfusion_del_manolo_auc <- just_auc(result_tadfusion_manolo_deletion)
  
  
  ## STRUCTRUE ------------------------
  
  result_structure_manolo_deletion <-  read_tsv('rival_cnvscore/structure/output.bed', col_names = c('chrom', 'start', 'end'))
  
  
  result_structure_manolo_deletion <- result_structure_manolo_deletion %>% 
    mutate(start = start + 1) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%

    rename(.pred_pathogenic = X5) %>%
    mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
    mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
    left_join(to_write_hg38, by = c('chrom', 'start', 'end')) %>%
    select(.pred_pathogenic, clinical)
  
  structure_del_manolo_auc <- just_auc(result_structure_manolo_deletion)
  
  
  ## JARVIS  ------------------------
  
  
  result_jarvis_manolo_deletion <- 1:22 %>% map_dfr(function(x) {
    
    tmp_df <-  read_tsv(glue('rival_cnvscore/jarvis/output_{x}.bed'), col_names = FALSE)
    
    if (nrow(tmp_df) == 0) {
      
      return(tibble())
      
    } else {
      
      tmp_df <-  tmp_df %>%
        rename(id = X4, score = X8) %>%
        select(id, score) %>%
        mutate(id = as.character(id))
      
    }
    
    return(tmp_df)
    
  })
  
  
  result_jarvis_manolo_deletion<- result_jarvis_manolo_deletion %>%
    filter(!is.infinite(score)) %>%
    na.omit() %>%
    group_by(id) %>%
    summarise(median_score = median(score),
              mean_score = mean(score),
              q1_score = quantile(score, 0.25),
              q3_score = quantile(score, 0.75)
    ) %>%
    pivot_longer(cols = -id) %>%
    group_by(id) %>%
    summarise(.pred_pathogenic = mean(value)) %>%
    left_join(x %>% select(id, clinical) %>% mutate(id = as.character(id))) %>%
    na.omit() 
  
  jarvis_del_manolo_auc <- just_auc(result_jarvis_manolo_deletion)
  
  
  ## gwRVIS  ------------------------
  
  result_gwrvis_manolo_deletion <- 1:22 %>% map_dfr(function(x) {
    
    print(x)
    
    tmp_df <-  read_tsv(glue('rival_cnvscore/gwrvis/output_{x}.bed'), col_names = FALSE)
    
    if (nrow(tmp_df) != 0) {
      
      tmp_df <-  tmp_df %>%
        rename(id = X4, score = X8) %>%
        select(id, score) %>%
        mutate(id = as.character(id))

    } else {
      
      return(tibble())
      
    }
    
    return(tmp_df)
    
    
  })
  
  
  result_gwrvis_manolo_deletion<- result_gwrvis_manolo_deletion %>%
    filter(!is.infinite(score)) %>%
    na.omit() %>%
    group_by(id) %>%
    summarise(median_score = median(score),
              mean_score = mean(score),
              q1_score = quantile(score, 0.25),
              q3_score = quantile(score, 0.75)
    ) %>%
    pivot_longer(cols = -id) %>%
    group_by(id) %>%
    summarise(.pred_pathogenic = 1 - mean(value)) %>%
    left_join(x %>% select(id, clinical) %>% mutate(id = as.character(id))) %>%
    na.omit() # ASK ANTONIO
  
  gwrvis_del_manolo_auc <- just_auc(result_gwrvis_manolo_deletion)
  
  ## Classify CNV -----------------
  
  
  result_classifycnv_manolo_deletion <- read_tsv('rival_cnvscore/classifycnv/result/Scoresheet.txt') %>%
    select(Chromosome, Start, End, `Total score`) %>%
    rename(score = `Total score`, chrom = Chromosome, start = Start, end = End) %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    mutate(score = (score - min(score)) / (max(score) - min(score))) %>%
    rename(.pred_pathogenic = score) %>%
    mutate(start = start + 1) %>%
    left_join(x %>% select(chrom, start, end, clinical), by = c('chrom', 'start', 'end'))
  
  classifycnv_del_manolo_auc <- just_auc(result_classifycnv_manolo_deletion)
  
  
  ## X-CNV------------------------
  
  result_xcnv_manolo_deletion <- read_csv('rival_cnvscore/x_cnv/input.output.csv') %>%
    rename(chrom = Chr, start = Start, end = End, .pred_pathogenic = MVP_score) %>%
    mutate(start = start + 1) %>%
    select(chrom, start, end, .pred_pathogenic) %>%
    mutate(chrom = as.character(chrom)) %>%
    left_join(x %>% select(chrom, start, end, clinical), by = c('chrom', 'start', 'end'))
  
  xcnv_del_manolo_auc <- just_auc(result_xcnv_manolo_deletion)
  
  
  
  ## AnnotSV------------------------
  
  
  # result_annotsv_manolo_deletion <- read_tsv('rival_cnvscore/annot_sv/output.bed') %>%
  #   filter(Annotation_mode == 'full') %>%
  #   rename(chrom = SV_chrom, start = SV_start, end = SV_end, .pred_pathogenic = AnnotSV_ranking_score) %>%
  #   mutate(start = start + 1) %>%
  #   select(chrom, start, end, .pred_pathogenic) %>%
  #   mutate(chrom = as.character(chrom)) %>%
  #   left_join(x %>% select(chrom, start, end, clinical), by = c('chrom', 'start', 'end')) %>%
  #   select(.pred_pathogenic, clinical)
  # 
  # 
  # annotsv_del_manolo_auc <- just_auc(result_annotsv_manolo_deletion)
  
  annotsv_del_manolo_auc <- c(0, 0)
  result_annotsv_manolo_deletion <- tibble()
  
  
  result_main_table <- result_main_table %>% bind_rows(tibble('tool' = c('CADD-SV', 'TADA', 'TADFusion', 'STRVCTURE', 'JARVIS', 'GWRVIS', 'ClassifyCNV', 'X-CNV', 'AnnotSV'),
                              'read_rows' = c(nrow(result_cadd_manolo_deletion), nrow(result_tada_manolo_deletion), nrow(result_tadfusion_manolo_deletion), nrow(result_structure_manolo_deletion), nrow(result_jarvis_manolo_deletion), nrow(result_gwrvis_manolo_deletion), nrow(result_classifycnv_manolo_deletion), nrow(result_xcnv_manolo_deletion), nrow(result_annotsv_manolo_deletion)),
                              'roc_auc' = c(cadd_del_manolo_auc[1], tada_del_manolo_auc[1], tadfusion_del_manolo_auc[1], structure_del_manolo_auc[1], jarvis_del_manolo_auc[1], gwrvis_del_manolo_auc[1], classifycnv_del_manolo_auc[1], xcnv_del_manolo_auc[1], annotsv_del_manolo_auc[1]),
                              'pr_auc' = c(cadd_del_manolo_auc[2], tada_del_manolo_auc[2], tadfusion_del_manolo_auc[2], structure_del_manolo_auc[2], jarvis_del_manolo_auc[2], gwrvis_del_manolo_auc[2], classifycnv_del_manolo_auc[2], xcnv_del_manolo_auc[2], annotsv_del_manolo_auc[2])
  ))
  
  
  result_main_table <- result_main_table %>% filter(!tool %in% c('TADFusion', 'AnnotSV'))
    
  
  result_main_table <-  result_main_table  %>%
    mutate(tag2 = case_when(
      str_detect(tool, 'CNVscore') ~ 'CNVscore',
      tool %in% c('ClassifyCNV', 'TADA') ~ 'Knowledge-based',
      tool %in% c('TADFusion','CADD-SV', 'STRVCTURE', 'JARVIS', 'GWRVIS', 'X-CNV') ~ 'Unbiased approach',
    ))
  
  
  predictions_tables <- bind_rows(
    result_cadd_manolo_deletion %>% select(.pred_pathogenic, id) %>% rename(tool = 'cadd'),
    result_tada_manolo_deletion %>% select(.pred_pathogenic, id) %>% rename(tool = 'tada'),
    result_structure_manolo_deletion  %>% select(.pred_pathogenic, id) %>% rename(tool = 'strvctvre'),
    result_jarvis_manolo_deletion  %>% select(.pred_pathogenic, id) %>% rename(tool = 'jarvis'),
    result_gwrvis_manolo_deletion  %>% select(.pred_pathogenic, id) %>% rename(tool = 'gwrvis'),
    result_classifycnv_manolo_deletion  %>% select(.pred_pathogenic, id) %>% rename(tool = 'classifycnv'),
    result_xcnv_manolo_deletion  %>% select(.pred_pathogenic, id) %>% rename(tool = 'xcnv')
    # result_annotsv_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    #   mutate(tag = glue('AnnotSV - {annotsv_del_manolo_auc[1]}'), tool = 'AnnotSV')
  )
  
  
  
  roc_curves_table <- roc_curves_table %>% bind_rows(
    result_cadd_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>%
      mutate(tag = glue('CADD-SV - {cadd_del_manolo_auc[1]}'), tool = 'CADD-SV'),
    result_tada_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('TADA (*) - {tada_del_manolo_auc[1]}'), tool = 'TADA'),
    # result_tadfusion_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    #   mutate(tag = glue('TADFusion - {tadfusion_del_manolo_auc[1]}'), tool = 'TADFusion'),
    result_structure_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('STRVCTURE  - {structure_del_manolo_auc[1]}'), tool = 'STRVCTURE'),
    result_jarvis_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('JARVIS - {jarvis_del_manolo_auc[1]}'), tool = 'JARVIS'),
    result_gwrvis_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('GWRVIS - {gwrvis_del_manolo_auc[1]}'), tool = 'GWRVIS'),
    result_classifycnv_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('ClassifyCNV - {classifycnv_del_manolo_auc[1]}'), tool = 'ClassifyCNV'),
    result_xcnv_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('X-CNV - {xcnv_del_manolo_auc[1]}'), tool = 'X-CNV')
    # result_annotsv_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    #   mutate(tag = glue('AnnotSV - {annotsv_del_manolo_auc[1]}'), tool = 'AnnotSV')
  )
  
  
  
  pr_curves_table <- pr_curves_table %>% bind_rows(
    result_cadd_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
      mutate(tag = glue('CADD-SV - {cadd_del_manolo_auc[2]}'), tool = 'CADD-SV'),
    result_tada_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('TADA  (*) - {tada_del_manolo_auc[2]}'), tool = 'TADA'),
    # result_tadfusion_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
    #   mutate(tag = glue('TADFusion - {tadfusion_del_manolo_auc[2]}'), tool = 'TADFusion'),
    result_structure_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('STRVCTURE - {structure_del_manolo_auc[2]}'), tool = 'STRVCTURE'),
    result_jarvis_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('JARVIS - {jarvis_del_manolo_auc[2]}'), tool = 'JARVIS'),
    result_gwrvis_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('GWRVIS - {gwrvis_del_manolo_auc[2]}'), tool = 'GWRVIS'),
    result_classifycnv_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('ClassifyCNV - {classifycnv_del_manolo_auc[2]}'), tool = 'ClassifyCNV'),
    result_xcnv_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
      mutate(tag = glue('X-CNV - {xcnv_del_manolo_auc[2]}'), tool = 'X-CNV')
    # result_annotsv_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
    #   mutate(tag = glue('AnnotSV - {annotsv_del_manolo_auc[2]}'), tool = 'AnnotSV')
  )

  tmp_result_main_table <- result_main_table  %>%
    left_join(roc_curves_table) %>%
    mutate(tag = fct_reorder(tag, -roc_auc))
  
  tmp2_result_main_table <- tibble(tag = levels(tmp_result_main_table$tag)) %>%
    mutate(tag2 = case_when(
      str_detect(tag, 'dummy') ~ 'red',
      str_detect(tag, 'bayesian|random_forest') ~ 'blue',
      str_detect(tag, 'ClassifyCNV |TADA') ~ 'yellow',
      str_detect(tag, 'TADFusion|CADD-SV|STRVCTURE|JARVIS|GWRVIS|X-CNV') ~ 'green'
      
    ))
  
  tmp3_result_main_table <- result_main_table  %>%
    left_join(pr_curves_table) %>%
    mutate(tag = fct_reorder(tag, -pr_auc))
  
  tmp4_result_main_table <- tibble(tag = levels(tmp3_result_main_table$tag)) %>%
    mutate(tag2 = case_when(
      str_detect(tag, 'dummy') ~ 'red',
      str_detect(tag, 'bayesian|random_forest') ~ 'blue',
      str_detect(tag, 'ClassifyCNV |TADA') ~ 'yellow',
      str_detect(tag, 'TADFusion|CADD-SV|STRVCTURE|JARVIS|GWRVIS|X-CNV') ~ 'green',
      
    ))
  
  
  
  result_main_plot1 <- tmp_result_main_table  %>%
    ggplot(aes(1-specificity, sensitivity)) +
    geom_path(aes(group = tag, color = reorder(tag, -roc_auc)),  show.legend = TRUE,  size = 1) +
    scale_color_manual(values = tmp2_result_main_table$tag2) +
    theme_roc() +
    theme(plot.title = element_text(size=18, face="bold"),
          axis.title.x = element_text(size=17, face="bold"),
          axis.title.y = element_text(size=17, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold"),
          legend.text = element_text(size=12)) +
    labs(title = paste(label_title), color = 'Tools')
  
  result_main_plot2 <- tmp3_result_main_table %>%
    ggplot(aes(recall, precision)) +
    geom_path(aes(group = tag, color = reorder(tag, -pr_auc)),  show.legend = TRUE,  size = 1) +
    scale_color_manual(values = tmp4_result_main_table$tag2) +
    theme_roc() +
    theme(plot.title = element_text(size=18, face="bold"),
          axis.title.x = element_text(size=17, face="bold"),
          axis.title.y = element_text(size=17, face="bold"),
          axis.text.x = element_text(size=14, face="bold"),
          axis.text.y = element_text(size=14, face="bold"),
          legend.text = element_text(size=12)) +
  labs(title = paste(label_title), color = 'Tools')
    
  
  return(list(result_main_table, roc_curves_table, pr_curves_table, predictions_tables))

}

# ------------------------------------------------------------------------------
# REMOVE NOISY OBSERVATIONS
# ------------------------------------------------------------------------------

remove_noisy_observations <- function(x) {
  
  # x <- clinvar_ind_deletion
  
  df_match_nomatter <- match_patho_benign(x, max(x$length_cnv))
  
  test_removed_ids_del <- tibble()
  for (i in 1:50) {
    print(i)
    test_kept_ids <- NoiseFiltersR::HARF(human_no_control, x, nfolds = 5, ntrees = 100)$cleanData %>%
      filter(clinical == 'pathogenic') %>% pull(id)
    
    test_removed_ids <- x %>% filter(clinical == 'pathogenic' & (!id %in% test_kept_ids)) %>% 
      mutate(iter_n = i) %>%
      select(iter_n, id)
    
    test_removed_ids_del <- test_removed_ids_del %>% bind_rows(test_removed_ids)
  }
  
  
  test_removed_ids_del <- test_removed_ids_del %>% count(id) %>% filter(n >= 35) %>% pull(id)
  
  x <- x %>%
    filter(clinical == 'pathogenic') %>%
    filter(!id %in% test_removed_ids_del) %>%
    bind_rows(x %>%
                filter(id %in% (df_match_nomatter %>% 
                                  filter(!id %in% test_removed_ids_del) %>% .$id_match)))
  
  return(x)
  
}


# ------------------------------------------------------------------------------
# REMOVE CONFLICTIVE OBSERVATIONS
# ------------------------------------------------------------------------------


overlap_benign_pathogenic <- function(df_patho, df_benign) {
  
  tmp_ref_sub1 <- df_patho %>%
    bed_intersect(df_benign) %>%
    mutate(overlap = (.overlap + 1)/ length_cnv.x) %>%
    filter(overlap >= 0.9) %>%
    rename(ref = id_tmp.x, sub = id_tmp.y) %>%
    select(ref, sub)
  
  tmp_ref_sub2 <- df_benign %>%
    bed_intersect(df_patho) %>%
    mutate(overlap = (.overlap + 1)/ length_cnv.x) %>%
    filter(overlap >= 0.9) %>%
    rename(ref = id_tmp.x, sub = id_tmp.y) %>%
    select(ref, sub)
  
  tmp_ref_sub3 <- tmp_ref_sub1 %>%
    left_join(tmp_ref_sub2, by = c('sub' = 'ref')) %>%
    filter(ref == sub.y) %>%
    select(ref, sub)
  
  return(list(id_tmp_patho = tmp_ref_sub3$ref %>% unique(), id_tmp_benign = tmp_ref_sub3$sub %>% unique()))
}

# ------------------------------------------------------------------------------
# DRAW FLOW
# ------------------------------------------------------------------------------

fill_flow <- function(x, tag) {
  
  n_del_patho_clinvar <- x %>% filter(source == 'clinvar' & variant_class == 'deletion' & clinical == 'pathogenic') %>% nrow()
  n_del_patho_decipher <- x %>% filter(source == 'decipher' & variant_class == 'deletion' & clinical == 'pathogenic') %>% nrow()
  n_dup_patho_clinvar <- x %>% filter(source == 'clinvar' & variant_class == 'duplication' & clinical == 'pathogenic') %>% nrow()
  n_dup_patho_decipher <- x %>% filter(source == 'decipher' & variant_class == 'duplication' & clinical == 'pathogenic') %>% nrow()
  
  n_del_benign <- x %>% filter(variant_class == 'deletion' & clinical == 'benign') %>% nrow()
  n_dup_benign <- x %>% filter(variant_class == 'duplication' & clinical == 'benign') %>% nrow()
  
  tibble('tag' = tag, 
         'n_del_patho_clinvar' = n_del_patho_clinvar,
         'n_del_patho_decipher' = n_del_patho_decipher,
         'n_dup_patho_clinvar' = n_dup_patho_clinvar,
         'n_dup_patho_decipher' = n_dup_patho_decipher,
         'n_del_benign' = n_del_benign,
         'n_dup_benign' = n_dup_benign)
  
}

# ------------------------------------------------------------------------------
# GET RESULTS
# ------------------------------------------------------------------------------


get_results <- function(input_df, type_variant = 'DEL', tag_title, do_all = TRUE) {
  
  # input_df <- output_decipher_duplication
  
  n_pos <- input_df %>% filter(clinical == 'pathogenic') %>% nrow()
  n_neg <- input_df %>% filter(clinical == 'benign') %>% nrow()
  
  if (type_variant == 'DEL') {
    
    result_naive_length <- predict_chrom_aware(logistic_clinvar_del_length, input_df)
    result_naive_n_genes <- predict_chrom_aware(logistic_clinvar_del_n_genes, input_df)
    result_naive_omim <- predict_chrom_aware(logistic_clinvar_del_omim, input_df)
    
    result_gbm_nohuman <- predict_chrom_aware(gbm_clinvar_del_nohuman, input_df)
    result_gbm_human <- predict_chrom_aware(gbm_clinvar_del_human, input_df)
    
    result_bay_nohuman  <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, input_df, 'deletion', 'unbiased approach')
    result_bay_human <- predict_chrom_aware_rtemis(bayesian_clinvar_del_human, input_df, 'deletion', 'knowledge-based')
    
    # result_bay_omim <- predict_chrom_aware_rtemis(bayesian_clinvar_del_omim, input_df, 'deletion', 'OMIM gene(s) both patho and non-patho')
    
    result_both <- predict_chrom_aware_rtemis(bayesian_clinvar_del_both, input_df, 'deletion', 'both')
    
    if (isTRUE(do_all)) {
      

    write_competitors(input_df)
    run_competitors('DEL')
    
    }
    
    list_results <- bind_rows(result_naive_length[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'naive_model_length'),
                         result_naive_n_genes[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'naive_model_n_genes'),
                         result_naive_omim[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'naive_model_omim'), 
                         result_gbm_nohuman[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'gbm_unbiased'),
                         result_gbm_human[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'gbm_knowledge_based'),
                         result_bay_nohuman[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_unbiased'),
                         result_bay_human[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_knowledge_based'),
                         # result_bay_omim[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_special_omim'),
                         result_both[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_both'))
    
    final_results_x_del <- read_competitors2(label_title = 
                                              glue('DEL - {tag_title} dataset ({n_pos} (+) - {n_neg} (-))'), 
                                            x = input_df,
                                            list_results
    )
    
    return(final_results_x_del)
    
  } else {
    
    result_naive_length <- predict_chrom_aware(logistic_clinvar_dup_length, input_df)
    result_naive_n_genes <- predict_chrom_aware(logistic_clinvar_dup_n_genes, input_df)
    result_naive_omim <- predict_chrom_aware(logistic_clinvar_dup_omim, input_df)
    
    result_gbm_nohuman <- predict_chrom_aware(gbm_clinvar_dup_nohuman, input_df)
    result_gbm_human <- predict_chrom_aware(gbm_clinvar_dup_human, input_df)
    
    result_bay_nohuman  <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_nohuman, input_df, 'duplication', 'unbiased approach')
    result_bay_human <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_human, input_df, 'duplication', 'knowledge-based')
    
    result_both <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_both, input_df, 'duplication', 'both')
    
    
     if (isTRUE(do_all)) {
      
      
      write_competitors(input_df)
      run_competitors('DUP')
      
    }
    
    # list_results <- list(result_naive_length,
    #                      result_naive_n_genes,
    #                      result_naive_omim,
    #                      result_gbm_nohuman,
    #                      result_gbm_human,
    #                      result_bay_nohuman,
    #                      result_bay_human,
    #                      result_both)
    
    
    
    list_results <- bind_rows(result_naive_length[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'naive_model_length'),
                              result_naive_n_genes[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'naive_model_n_genes'),
                              result_naive_omim[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'naive_model_omim'), 
                              result_gbm_nohuman[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'gbm_unbiased'),
                              result_gbm_human[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'gbm_knowledge_based'),
                              result_bay_nohuman[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_unbiased'),
                              result_bay_human[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_knowledge_based'),
                              # result_bay_omim[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_special_omim'),
                              result_both[[3]] %>% select(id, .pred_pathogenic) %>% mutate(tag = 'bayesian_both'))
    
    
    final_results_x_dup <- read_competitors2(
      label_title = glue('DUP - {tag_title} dataset ({n_pos} (+) - {n_neg} (-))'), 
      x = input_df,
      list_results
    )
    
    return(final_results_x_dup)
    
  }
  
}


# ------------------------------------------------------------------------------
# PALETTE COLORS
# ------------------------------------------------------------------------------



branded_colors <- list(
  "color1"   = "#f94144",
  "color2"    = "#f3722c",
  "color3" = "#f8961e",
  "color4"  = "#f9c74f",
  "color5"   = "#90be6d", 
  "color6"   = "#43aa8b",
  "color7" = '#577590'
)


branded_pal <- function(
  primary = "color1", 
  other = "grey", 
  direction = 1
) {
  stopifnot(primary %in% names(branded_colors))
  
  function(n) {
    if (n > 6) warning("Branded Color Palette only has 6 colors.")
    
    if (n == 2) {
      other <- if (!other %in% names(branded_colors)) {
        other
      } else {
        branded_colors[other]
      }
      color_list <- c(other, branded_colors[primary])
    } else {
      color_list <- branded_colors[1:n]
    }
    
    color_list <- unname(unlist(color_list))
    if (direction >= 0) color_list else rev(color_list)
  }
}

scale_color_branded <- function(
  primary = "color1", 
  other = "grey", 
  direction = 1, 
  ...
) {
  ggplot2::discrete_scale(
    "colour", "branded", 
    branded_pal(primary, other, direction), 
    ...
  )
}

scale_fill_branded <- function(
  primary = "color1", 
  other = "grey", 
  direction = 1, 
  ...
) {
  ggplot2::discrete_scale(
    "fill", "branded", 
    branded_pal(primary, other, direction), 
    ...
  )
}


# ------------------------------------------------------------------------------
# TO GOOGLESHEET
# ------------------------------------------------------------------------------

to_googlesheet <- function(x) {

  result_df <- x[[2]] %>%
    select(tag) %>%
    distinct() %>%
    mutate(roc_auc = str_extract(tag, '[0-9]+(\\.[0-9]{1,3})?')) %>%
    mutate(tag = str_remove(tag, '[0-9]+(\\.[0-9]{1,3})?')) %>%
    mutate(tag = str_remove(tag, 'deletion - ')) %>%
    mutate(tag = if_else(tag == 'TADA (*) - ', 'TADA  (*) - ', tag)) %>%
    mutate(tag = if_else(tag == 'STRVCTURE  - ', 'STRVCTURE - ', tag)) %>%
    left_join(
      x[[3]] %>%
        select(tag) %>%
        distinct() %>%
        mutate(pr_auc = str_extract(tag, '[0-9]+(\\.[0-9]{1,3})?')) %>%
        mutate(tag = str_remove(tag, '[0-9]+(\\.[0-9]{1,3})?')) %>%
        mutate(tag = str_remove(tag, 'deletion - '))

    , by = 'tag') %>%
    mutate(roc_auc = as.numeric(roc_auc)) %>%
    mutate(pr_auc = as.numeric(pr_auc)) %>%
    arrange(desc(roc_auc))
  return(result_df)
}



# ------------------------------------------------------------------------------
# PLOT PRECISION AND NPV IMPROVED IV
# ------------------------------------------------------------------------------


plot_rates4 <- function(x, 
                        y,
                        # z,
                        tag, tag2, 
                        input_patho = 0.5, 
                        input_benign = 0.5, 
                        fixed_interval = 'fixed-interval',
                        limit_y = 0.75,
                        n_split_score = NULL,
                        n_split = 20, 
                        tmp_sub) {
  
  
  # x <- almost_clinvar
  # y <- res_df_both
  # z <- res_df_combined
  # tag <- 'ClinVar'
  # tag2 <- 'Model (unbiased features)'
  # fixed_interval <- 'equally-sized'
  # n_split <- 3
  # n_split_score <- 5
  # what_is_pathogenic <- 0.5
  # what_is_benign <- 0.5
  # limit_y <- 0.6
  # tmp_sub <- 'test'
  
  what_is_pathogenic <- input_patho
  what_is_benign <- input_benign
  
  conditional_interval <- function(fixed_interval, a) {
    
    if (fixed_interval == 'fixed-interval') {
      
      # cut_interval(a, n = 10)
      cut_interval(a, n = n_split)
      
      
    } else {
      
      ntile(a, n = n_split)
      
    }
    
  }
  
  
  # tmp_precision_combined <- z %>%
  #   mutate(result_cnvscore_combined = case_when(
  #     cnvscore_score_combined >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
  #     cnvscore_score_combined >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
  #   select(reliability_score, cnvscore_score_combined, contains('result')) %>%
  #   pivot_longer(-c(reliability_score, cnvscore_score_combined),
  #                names_to = 'tool', values_to = 'result') %>%
  #   filter(result %in% c('TP', 'FP')) %>%
  #   group_by(reliability_score) %>%
  #   mutate(score_interval =
  #            factor(ntile(cnvscore_score_combined, n = n_split_score))) %>%
  #   ungroup() %>%
  #   count(reliability_score, score_interval, tool, result) %>%
  #   pivot_wider(id_cols = c(reliability_score, score_interval, tool),
  #               names_from = 'result', values_from = 'n') %>%
  #   group_by(reliability_score, score_interval, tool) %>%
  #   mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
  #   mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
  #   mutate(precision = TP / (TP + FP)) %>%
  #   ungroup() %>%
  #   mutate(reliability_score = ntile(reliability_score, 3)) %>%
  #   mutate(reliability_score = factor(reliability_score))
  
  
  
  tmp_precision_both <- y %>%
    mutate(result_cnvscore_both = case_when(
      cnvscore_score_both >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
      cnvscore_score_both >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
    select(reliability_score, cnvscore_score_both, contains('result')) %>%
    pivot_longer(-c(reliability_score, cnvscore_score_both),
                 names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TP', 'FP')) %>%
    group_by(reliability_score) %>%
    mutate(score_interval =
             factor(ntile(cnvscore_score_both, n = n_split_score))) %>%
    ungroup() %>%
    count(reliability_score, score_interval, tool, result) %>%
    pivot_wider(id_cols = c(reliability_score, score_interval, tool),
                names_from = 'result', values_from = 'n') %>%
    group_by(reliability_score, score_interval, tool) %>%
    mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
    mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
    mutate(precision = TP / (TP + FP)) %>%
    ungroup() %>%
    mutate(reliability_score = factor(reliability_score))

  
  tmp_precision <- x %>%
      mutate(result_cnvscore = case_when(
        cnvscore_score >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
        cnvscore_score >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
    mutate(result_structure = case_when(
      structure_score >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
      structure_score >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
    select(reliability_score, cnvscore_score, structure_score, contains('result')) %>%
    pivot_longer(-c(reliability_score, cnvscore_score, structure_score), 
                 names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TP', 'FP')) %>%
    group_by(reliability_score) %>%
    mutate(score_interval = factor(ntile(cnvscore_score, n = n_split_score))) %>%
    ungroup() %>%
    # filter(tool == 'result_cnvscore') %>%
    # ggplot(aes(score_interval, cnvscore_score)) +
    # geom_boxplot(aes(color = factor(reliability_score))) + theme_minimal()
    count(reliability_score, score_interval, tool, result) %>%
    pivot_wider(id_cols = c(reliability_score, score_interval, tool),
                names_from = 'result', values_from = 'n') %>%
    # mutate(score_interval = factor(ntile(cnvscore_score, n = n_split_score))) %>%
    group_by(reliability_score, score_interval, tool) %>%
    mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
    mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
    mutate(precision = TP / (TP + FP)) %>%
    ungroup() %>%
    mutate(reliability_score = factor(reliability_score)) %>%
    bind_rows(tmp_precision_both)
  
  
  
  
  # tmp_precision_both %>%
  #   filter(tool == 'result_cnvscore_both') %>%
  #   mutate(total_total = FP + TP)
  # 
  # tmp_precision_both %>%
  #   filter(tool == 'result_cnvscore_both') %>%
  #   group_by(reliability_score) %>%
  #   # group_by(reliability_score, score_interval) %>%
  #   summarise(total_tp = sum(TP),
  #             total_fp = sum(FP)) %>%
  #   mutate(precision = total_tp / (total_tp + total_fp))
    
  
  
  pre_reli_sd1 <- tmp_precision %>%
    ggplot(aes(score_interval, precision)) +
    geom_line(aes(group = reliability_score, 
                    color = reliability_score)) +
    geom_point(aes(fill = reliability_score), shape = 21, color = 'black') +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(limit_y,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (Precision) - ', tag2),
         subtitle = tmp_sub,
         x = 'Equally-sized bins (CNVscore score))', y = 'Precision (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    facet_wrap(vars(tool))
  
  
  
  pre_reli_sd2 <- tmp_precision %>%
    mutate(total_obs = TP + FP) %>%
    group_by(score_interval) %>%
    mutate(perc = total_obs / sum(total_obs)) %>%
    # arrange(score_interval)
    ggplot(aes(as.factor(score_interval), perc)) +
    geom_col(aes(group = reliability_score, fill = reliability_score, 
                 color = reliability_score),
             color = 'black') +
    scale_y_continuous(label = percent) +
    geom_text(aes(label = paste0(100*round(perc, 4), '%', paste0('(', total_obs,')')),
                  group = reliability_score), 
              position = position_stack(vjust = 0.3),
              size = 3,vjust = 0) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (Precision) - ', tag2),
         subtitle = tmp_sub,
         x = 'Equally-sized bins (CNVscore score))', y = 'Precision (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  
  tmp_npv_both <- y %>%
    mutate(result_cnvscore_both = case_when(
      cnvscore_score_both < what_is_pathogenic & clinical == 'benign' ~ 'TN',
      cnvscore_score_both < what_is_pathogenic & clinical == 'pathogenic' ~ 'FN')) %>%
    select(reliability_score, cnvscore_score_both, contains('result')) %>%
    pivot_longer(-c(reliability_score, cnvscore_score_both),
                 names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TN', 'FN')) %>%
    group_by(reliability_score) %>%
    mutate(score_interval =
             factor(ntile(cnvscore_score_both, n = n_split_score))) %>%
    ungroup() %>%
    count(reliability_score, score_interval, tool, result) %>%
    pivot_wider(id_cols = c(reliability_score, score_interval, tool),
                names_from = 'result', values_from = 'n') %>%
    group_by(reliability_score, score_interval, tool) %>%
    mutate(FN = ifelse(is.na(FN), 0, FN)) %>%
    mutate(TN = ifelse(is.na(TN), 0, TN)) %>%
    mutate(npv = TN / (TN + FN)) %>%
    ungroup() %>%
    mutate(reliability_score = factor(reliability_score))

  tmp_npv <- x %>%
    mutate(result_cnvscore = case_when(
      cnvscore_score < what_is_pathogenic & clinical == 'benign' ~ 'TN',
      cnvscore_score < what_is_pathogenic & clinical == 'pathogenic' ~ 'FN')) %>%
    mutate(result_structure = case_when(
      structure_score < what_is_pathogenic & clinical == 'benign' ~ 'TN',
      structure_score < what_is_pathogenic & clinical == 'pathogenic' ~ 'FN')) %>%
    select(reliability_score, cnvscore_score, structure_score, contains('result')) %>%
    pivot_longer(-c(reliability_score, cnvscore_score, structure_score), 
                 names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TN', 'FN')) %>%
    group_by(reliability_score) %>%
    mutate(score_interval = factor(ntile(cnvscore_score, n = n_split_score))) %>%
    ungroup() %>%
    count(reliability_score, score_interval, tool, result) %>%
    pivot_wider(id_cols = c(reliability_score, score_interval, tool),
                names_from = 'result', values_from = 'n') %>%
    group_by(reliability_score, score_interval, tool) %>%
    mutate(FN = ifelse(is.na(FN), 0, FN)) %>%
    mutate(TN = ifelse(is.na(TN), 0, TN)) %>%
    mutate(npv = TN / (TN + FN)) %>%
    ungroup() %>%
    mutate(reliability_score = factor(reliability_score)) %>%
    bind_rows(tmp_npv_both)

  
  pre_reli_sd3 <- tmp_npv %>%
    ggplot(aes(score_interval, npv)) +
    geom_line(aes(group = reliability_score, 
                    color = reliability_score)) +
    geom_point(aes(fill = reliability_score), shape = 21, color = 'black') +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(limit_y,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (NPV) - ', tag2),
         subtitle = tmp_sub,
         x = 'Equally-sized bins (CNVscore score))', y = 'NPV (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    facet_wrap(vars(tool))
  
  
  pre_reli_sd4 <- tmp_npv %>%
    mutate(total_obs = FN + TN) %>%
    group_by(score_interval) %>%
    mutate(perc = total_obs / sum(total_obs)) %>%
    # arrange(score_interval)
    ggplot(aes(as.factor(score_interval), perc)) +
    geom_col(aes(group = reliability_score, fill = reliability_score, color = reliability_score),
             color = 'black') +
    scale_y_continuous(label = percent) +
    geom_text(aes(label = paste0(100*round(perc, 4), '%', paste0('(', total_obs,')')),
                  group = reliability_score), 
              position = position_stack(vjust = 0.3),
              size = 3,vjust = 0) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (NPV) - ', tag2),
         subtitle = tmp_sub,
         x = 'Equally-sized bins (CNVscore score))', y = 'NPV (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  
  pre_reli_sd5 <- x %>%
    select(contains('score'), clinical) %>%
    pivot_longer(cols = -c(reliability_score, clinical), names_to = 'tool',
                 values_to = 'score') %>%
    bind_rows(y %>% rename(score = cnvscore_score_both) %>%
                mutate(tool = 'cnvscore_score_both') %>%
                select(-id)) %>%
    group_by(reliability_score, tool) %>%
    roc_auc(clinical, score) %>%
    ggplot(aes(factor(reliability_score), .estimate)) +
    geom_smooth(aes(group = tool, fill = tool, 
                    color = tool)) +
      geom_point() +
    labs(x = 'Reliability_score', y = 'AUROC') +
    theme_minimal()

  pre_reli_sd1 + 
    pre_reli_sd2 +
    pre_reli_sd3 + 
    pre_reli_sd4 +
    pre_reli_sd5 + plot_layout(nrow = 2)
  
  
}

# ------------------------------------------------------------------------------
# PLOT PRECISION AND NPV IMPROVED III
# ------------------------------------------------------------------------------

plot_rates3 <- function(x, y, z, tag, tag2, 
                        input_patho = 0.5, 
                        input_benign = 0.5, 
                        fixed_interval = TRUE,
                        limit_y = 0.75,
                        n_split = 20, 
                        tmp_sub) {
  

  # x <- almost_decipher
  # fixed_interval <- 'fixed-interval'
  # n_split <- 3
  # y <- res_df_both
  # z <- res_df_knowledge_based
  # tag <- 'DECIPHER'
  # tag2 <- 'Model (unbiased features)'
  # input_patho <- 0.5
  # input_benign <- 0.5
  # limit_y <- 0.6
  # tmp_sub <- 'test'
  
  what_is_pathogenic <- input_patho
  what_is_benign <- input_benign
  
  total_precision_sd <- tibble()
  total_npv_sd <- tibble()
  total_auc <- tibble()
  
  for (i in 1:100) {
    
    # tmp_precision_knowledge_based <- z %>%
    #   group_by(clinical) %>%
    #   sample_frac(0.8) %>%
    #   ungroup() %>%
    #   mutate(result_cnvscore_knowledge_based = case_when(
    #     cnvscore_score_knowledge_based >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
    #     cnvscore_score_knowledge_based >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
    #   select(reliability_score, contains('result')) %>%
    #   pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
    #   filter(result %in% c('TP', 'FP')) %>%
    #   mutate(tile_general = reliability_score) %>%
    #   count(tile_general, tool, result) %>%
    #   pivot_wider(id_cols = c(tile_general, tool),
    #               names_from = 'result', values_from = 'n')
    
    
    # if (! 'TP' %in% colnames(tmp_precision_knowledge_based)) {
    #   
    #   tmp_precision_knowledge_based <- tmp_precision_knowledge_based %>% mutate(FN = 0)
    # } 
    # 
    # if (! 'FP' %in% colnames(tmp_precision_knowledge_based)) {
    #   
    #   tmp_precision_knowledge_based <- tmp_precision_knowledge_based %>% mutate(FP = 0)
    # }
    
    # tmp_precision_knowledge_based <- tmp_precision_knowledge_based %>%
    #   group_by(tile_general, tool) %>%
    #   mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
    #   mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
    #   mutate(precision = TP / (TP + FP)) %>%
    #   ungroup() %>%
    #   select(tile_general, tool, precision) %>%
    #   mutate(n_iter = i)
    
    
      # tmp_precision_both <- y %>%
      # group_by(clinical) %>%
      # sample_frac(0.8) %>%
      # ungroup() %>%
      # mutate(result_cnvscore_both = case_when(
      #   cnvscore_score_both >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
      #   cnvscore_score_both >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
      # select(reliability_score, contains('result')) %>%
      # pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
      # filter(result %in% c('TP', 'FP')) %>%
      # mutate(tile_general = reliability_score) %>%
      # count(tile_general, tool, result) %>%
      # pivot_wider(id_cols = c(tile_general, tool),
      #             names_from = 'result', values_from = 'n')
      # 
    
    # if (! 'TP' %in% colnames(tmp_precision_both)) {
    #   
    #   tmp_precision_both <- tmp_precision_both %>% mutate(FN = 0)
    # } 
    # 
    # if (! 'FP' %in% colnames(tmp_precision_both)) {
    #   
    #   tmp_precision_both <- tmp_precision_both %>% mutate(FP = 0)
    # } 
    
    # tmp_precision_both <- tmp_precision_both %>%
    #   group_by(tile_general, tool) %>%
    #   mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
    #   mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
    #   mutate(precision = TP / (TP + FP)) %>%
    #   ungroup() %>%
    #   select(tile_general, tool, precision) %>%
    #   mutate(n_iter = i)

  tmp_precision_sd <- x %>%
    group_by(clinical) %>%
    sample_frac(0.8) %>%
    ungroup() %>%
    mutate(result_cnvscore = case_when(
      cnvscore_score >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
      cnvscore_score >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
    mutate(result_structure = case_when(
      structure_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
      structure_score >= 0.5 & clinical == 'benign' ~ 'FP')) %>%
    select(reliability_score, contains('result')) %>%
    pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TP', 'FP')) %>%
    mutate(tile_general = reliability_score) %>%
    count(tile_general, tool, result) %>%
    pivot_wider(id_cols = c(tile_general, tool),
                names_from = 'result', values_from = 'n') %>%
    group_by(tile_general, tool) %>%
    mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
    mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
    mutate(precision = TP / (TP + FP)) %>%
    ungroup() %>%
    select(tile_general, tool, precision) %>%
    mutate(n_iter = i)
  
  
  total_precision_sd <- total_precision_sd %>%
    bind_rows(tmp_precision_sd)
    # bind_rows(tmp_precision_both) %>%
    # bind_rows(tmp_precision_knowledge_based)
  
  
  # tmp_npv_knowledge_based <- z %>%
  #   group_by(clinical) %>%
  #   sample_frac(0.8) %>%
  #   ungroup() %>%
  #   mutate(result_cnvscore_knowledge_based = case_when(
  #     cnvscore_score_knowledge_based < what_is_benign & clinical == 'benign' ~ 'TN',
  #     cnvscore_score_knowledge_based < what_is_benign & clinical == 'pathogenic' ~ 'FN')) %>%
  #   select(reliability_score, contains('result')) %>%
  #   pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
  #   filter(result %in% c('TN', 'FN')) %>%
  #   mutate(tile_general = reliability_score) %>%
  #   count(tile_general, tool, result) %>%
  #   pivot_wider(id_cols = c(tile_general, tool), 
  #               names_from = 'result', values_from = 'n') %>%
  #   group_by(tile_general, tool)
  # 
  # if (! 'FN' %in% colnames(tmp_npv_knowledge_based)) {
  #   
  #   tmp_npv_knowledge_based <- tmp_npv_knowledge_based %>% mutate(FN = 0)
  # } 
  # 
  # if (! 'TN' %in% colnames(tmp_npv_knowledge_based)) {
  #   
  #   tmp_npv_knowledge_based <- tmp_npv_knowledge_based %>% mutate(TN = 0)
  # } 
  # tmp_npv_knowledge_based <- tmp_npv_knowledge_based %>%
  #   mutate(FN = ifelse(is.na(FN), 0, FN)) %>%
  #   mutate(TN = ifelse(is.na(TN), 0, TN)) %>%
  #   mutate(npv = TN / (TN + FN)) %>%
  #   ungroup() %>%
  #   select(tile_general, tool, npv) %>%
  #   mutate(n_iter = i)
  
  
  # tmp_npv_both <- y %>%
  #   group_by(clinical) %>%
  #   sample_frac(0.8) %>%
  #   ungroup() %>%
  #   mutate(result_cnvscore_both = case_when(
  #     cnvscore_score_both < what_is_benign & clinical == 'benign' ~ 'TN',
  #     cnvscore_score_both < what_is_benign & clinical == 'pathogenic' ~ 'FN')) %>%
  #   select(reliability_score, contains('result')) %>%
  #   pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
  #   filter(result %in% c('TN', 'FN')) %>%
  #   mutate(tile_general = reliability_score) %>%
  #   count(tile_general, tool, result) %>%
  #   pivot_wider(id_cols = c(tile_general, tool), 
  #               names_from = 'result', values_from = 'n') %>%
  #   group_by(tile_general, tool)
  #   
  #   if (! 'FN' %in% colnames(tmp_npv_both)) {
  #     
  #     tmp_npv_both <- tmp_npv_both %>% mutate(FN = 0)
  #   } 
  # 
  # if (! 'TN' %in% colnames(tmp_npv_both)) {
  #   
  #   tmp_npv_both <- tmp_npv_both %>% mutate(TN = 0)
  # } 
  #   tmp_npv_both <- tmp_npv_both %>%
  #   mutate(FN = ifelse(is.na(FN), 0, FN)) %>%
  #   mutate(TN = ifelse(is.na(TN), 0, TN)) %>%
  #   mutate(npv = TN / (TN + FN)) %>%
  #   ungroup() %>%
  #   select(tile_general, tool, npv) %>%
  #   mutate(n_iter = i)
  
 
  tmp_npv_sd <- x %>%
    group_by(clinical) %>%
    sample_frac(0.8) %>%
    ungroup() %>%
    mutate(result_cnvscore = case_when(
      cnvscore_score < what_is_benign & clinical == 'benign' ~ 'TN',
      cnvscore_score < what_is_benign & clinical == 'pathogenic' ~ 'FN')) %>%
    mutate(result_structure = case_when(
      structure_score < 0.5 & clinical == 'benign' ~ 'TN',
      structure_score < 0.5 & clinical == 'pathogenic' ~ 'FN')) %>%
    select(reliability_score, contains('result')) %>%
    pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TN', 'FN')) %>%
    mutate(tile_general = reliability_score) %>%
    count(tile_general, tool, result) %>%
    pivot_wider(id_cols = c(tile_general, tool), 
                names_from = 'result', values_from = 'n') %>%
    group_by(tile_general, tool) %>%
    mutate(FN = ifelse(is.na(FN), 0, FN)) %>%
    mutate(TN = ifelse(is.na(TN), 0, TN)) %>%
    mutate(npv = TN / (TN + FN)) %>%
    ungroup() %>%
    select(tile_general, tool, npv) %>%
    mutate(n_iter = i)
  
  total_npv_sd <- total_npv_sd %>%
    bind_rows(tmp_npv_sd)
    # bind_rows(tmp_npv_both) %>%
    # bind_rows(tmp_npv_knowledge_based)
  
  
  
  # total_auc_both <- y %>%
  #       group_by(clinical) %>%
  #       sample_frac(0.8) %>%
  #       ungroup() %>%
  #       group_by(reliability_score) %>%
  #       roc_auc(clinical, cnvscore_score_both) %>%
  #       mutate(tool = 'cnvscore_score_both') %>%
  #       mutate(n_iter = i)
  
  # total_auc_knowledge_based <- z %>%
  #   group_by(clinical) %>%
  #   sample_frac(0.8) %>%
  #   ungroup() %>%
  #   group_by(reliability_score) %>%
  #   roc_auc(clinical, cnvscore_score_knowledge_based) %>%
  #   mutate(tool = 'cnvscore_score_knowledge_based') %>%
  #   mutate(n_iter = i)
  
  # while (is.na(total_auc_knowledge_based[3,]$.estimate)) {
  #   
  #   total_auc_knowledge_based <- z %>%
  #     group_by(clinical) %>%
  #     sample_frac(0.8) %>%
  #     ungroup() %>%
  #     group_by(reliability_score) %>%
  #     roc_auc(clinical, cnvscore_score_knowledge_based) %>%
  #     mutate(tool = 'cnvscore_score_knowledge_based') %>%
  #     mutate(n_iter = i)
  #   
  # }
  
  
  total_auc <- total_auc %>%
    bind_rows(
  x %>%
    group_by(clinical) %>%
    sample_frac(0.8) %>%
    ungroup() %>%
    group_by(reliability_score) %>%
    roc_auc(clinical, cnvscore_score) %>%
    mutate(tool = 'cnvscore_score') %>%
    mutate(n_iter = i)
    ) %>%
    bind_rows(
      x %>%
        group_by(clinical) %>%
        sample_frac(0.8) %>%
        ungroup() %>%
        group_by(reliability_score) %>%
        roc_auc(clinical, structure_score) %>%
        mutate(tool = 'structure_code') %>%
        mutate(n_iter = i)
    )
    # bind_rows(total_auc_both) %>%
    # bind_rows(total_auc_knowledge_based)
    
  }
  
  # sanity_check <- total_auc %>% 
  #   filter(tool == 'cnvscore_score_both') %>%
  #   filter(reliability_score %in% c(1,3)) %>%
  #   group_by(reliability_score) %>%
  #   summarise(min_estimate = min(.estimate)) %>%
  #   pull(min_estimate) %>%
  #   sum()
  # 
  # sanity_check2 <- total_precision_sd %>%
  #   filter(tool == 'result_cnvscore_knowledge_based') %>%
  #   filter(tile_general %in% c(3)) %>%
  #     pull(precision) %>%
  #   sum()
  
  # CHECK THIS 
 # if (sanity_check2 == 0) {
 #   
 #   total_precision_sd <- total_precision_sd %>%
 #     mutate(precision = ifelse(tool == 'result_cnvscore_knowledge_based' & tile_general == 3,
 #                                sample(c(0.01, 0.02)), precision))
 #   
 # }

  # if (sanity_check == 2) {
  #   
  #   tag_sign_auc <- total_auc %>%
  #     filter(!tool %in% 'cnvscore_score_both') %>%
  #     select(-n_iter) %>%
  #     group_by(tool) %>%
  #     wilcox_test(.estimate ~ reliability_score) %>%
  #     group_by(tool) %>%
  #     slice(2) %>%
  #     bind_rows(tibble(tool = 'result_cnvscore_both', p = 1)) %>%
  #     mutate(tag = paste(tool, paste0('(p.value = ', p, ')'))) %>%
  #     pull(tag) %>%
  #     paste(collapse = ' \n') %>%
  #     paste0('\n', .)
  #   
  #   tag_sign_precision <-  total_precision_sd %>%
  #     filter(!tool %in% 'result_cnvscore_both') %>%
  #     select(-n_iter) %>%
  #     group_by(tool) %>%
  #     # ggplot(aes(as.factor(tile_general), precision)) +
  #     # geom_boxplot() + facet_wrap(vars(tool))
  #     wilcox_test(precision ~ tile_general) %>%
  #     group_by(tool) %>%
  #     slice(2) %>%
  #     bind_rows(tibble(tool = 'result_cnvscore_both', p = 1)) %>%
  #     mutate(tag = paste(tool, paste0('(p.value = ', p, ')'))) %>%
  #     pull(tag) %>%
  #     paste(collapse = ' \n') %>%
  #     paste0('\n', .)
  #   
  #   
  #   tag_sign_npv <-  total_npv_sd %>%
  #     filter(tile_general %in% c(1,3)) %>%
  #     filter(!tool %in% 'result_cnv') %>%
  #     select(-n_iter) %>%
  #     group_by(tool) %>%
  #     wilcox_test(npv ~ tile_general) %>%
  #     group_by(tool) %>%
  #     slice(1) %>%
  #     # bind_rows(tibble(tool = 'result_cnvscore_both', p = 1)) %>%
  #     mutate(tag = paste(tool, paste0('(p.value = ', p, ')'))) %>%
  #     pull(tag) %>%
  #     paste(collapse = ' \n') %>%
  #     paste0('\n', .)
  # } else {
  
  
  tag_sign_auc <- total_auc %>%
    select(-n_iter) %>%
    group_by(tool) %>%
    wilcox_test(.estimate ~ reliability_score) %>%
    group_by(tool) %>%
    slice(2) %>%
    mutate(tag = paste(tool, paste0('(p.value = ', p, ')'))) %>%
    pull(tag) %>%
    paste(collapse = ' \n') %>%
    paste0('\n', .)

  tag_sign_precision <-  total_precision_sd %>%
    select(-n_iter) %>%
    group_by(tool) %>%
    wilcox_test(precision ~ tile_general) %>%
    group_by(tool) %>%
    slice(2) %>%
    # select(tool, p.adj) %>%
    mutate(tag = paste(tool, paste0('(p.value = ', p, ')'))) %>%
    pull(tag) %>%
    paste(collapse = ' \n') %>%
    paste0('\n', .)

  
 tag_sign_npv <-  total_npv_sd %>%
    filter(tile_general %in% c(1,3)) %>%
    select(-n_iter) %>%
    group_by(tool) %>%
    wilcox_test(npv ~ tile_general) %>%
    group_by(tool) %>%
    slice(1) %>%
    # select(tool, p.adj) %>%
    mutate(tag = paste(tool, paste0('(p.value = ', p, ')'))) %>%
    pull(tag) %>%
    paste(collapse = ' \n') %>%
    paste0('\n', .)
 
  # }
 
 auc_plot <- total_auc %>%
   group_by(tool, reliability_score) %>%
   summarise(
     mean = mean(.estimate),
     sd = sd(.estimate),
     n = n(),
     se = sd / sqrt(n)
   ) %>%
   ggplot(aes(as.factor(reliability_score), mean)) +
   geom_smooth(aes(group = tool, fill = tool, color = tool)) +
   geom_point(aes(fill = tool), shape = 21, color = 'black') +
   scale_y_continuous(label = percent) +
   coord_cartesian(ylim = c(limit_y,1)) +
   geom_pointrange(aes(ymin = mean - se, ymax = mean + se)) +
   theme_minimal() +
   labs(title = paste(tag, 'dataset (AUROC) - ', tag2),
        subtitle = paste0(tmp_sub, tag_sign_auc),
        x = 'Reliability score (percentile(residual))', y = 'AUROC') +
   theme(axis.text.x = element_text(angle = 45, hjust=1),
         plot.subtitle = element_text(size = 10))
 
 
 precision_sd <- total_precision_sd %>%
   group_by(tool, tile_general) %>%
   summarise(
     mean = mean(precision),
     sd = sd(precision),
     n = n(),
     se = sd / sqrt(n)
   ) %>%
   ggplot(aes(as.factor(tile_general), mean)) +
   geom_smooth(aes(group = tool, fill = tool, color = tool)) +
   geom_point(aes(fill = tool), shape = 21, color = 'black') +
   scale_y_continuous(label = percent) +
   coord_cartesian(ylim = c(limit_y,1)) +
   geom_pointrange(aes(ymin = mean - se, ymax = mean + se)) +
   theme_minimal() +
   labs(title = paste(tag, 'dataset (Precision) - ', tag2),
        subtitle = paste0(tmp_sub, tag_sign_precision),
        x = 'Reliability score (percentile(residual))', y = 'Precision (%)') +
   theme(axis.text.x = element_text(angle = 45, hjust=1),
         plot.subtitle = element_text(size = 10))
 
  
  
  npv_sd <- total_npv_sd %>%
    group_by(tool, tile_general) %>%
    summarise(
      mean = mean(npv),
      sd = sd(npv),
      n = n(),
      se = sd / sqrt(n)
    ) %>%
    ggplot(aes(as.factor(tile_general), mean)) +
    geom_smooth(aes(group = tool, fill = tool, color = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(limit_y,1)) +
    geom_pointrange(aes(ymin = mean - se, ymax = mean + se)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (NPV) - ', tag2),
         subtitle = paste0(tmp_sub, tag_sign_npv),
         x = 'Reliability score (percentile(residual))', y = 'NPV (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          plot.subtitle = element_text(size = 10))

  precision_sd +
    npv_sd +
    auc_plot +
    plot_layout(nrow = 1)

}


# ------------------------------------------------------------------------------
# PLOT PRECISION AND NPV IMPROVED
# ------------------------------------------------------------------------------


plot_rates2 <- function(x, tag, tag2, 
                        input_patho = 0.5, 
                        input_benign = 0.5, 
                        fixed_interval = TRUE,
                        limit_y = 0.75,
                        n_split = 20, 
                        tmp_sub) {
  
  
  # x <- almost_decipher
  # tag <- 'DECIPHER'
  # tag2 <- 'Model (unbiased features)'
  # fixed_interval <- 'fixed-interval'
  # n_split <- 3
  # what_is_pathogenic <- 0.5
  # what_is_benign <- 0.5
  # limit_y <- 0.6
  # tmp_sub <- 'test'
  

  
  what_is_pathogenic <- input_patho
  what_is_benign <- input_benign

  conditional_interval <- function(fixed_interval, a) {
    
    if (fixed_interval == 'fixed-interval') {
      
      # cut_interval(a, n = 10)
      cut_interval(a, n = n_split)
      
      
    } else {
      
      ntile(a, n = n_split)
      
    }

  }
  
  


  precision_sd <- x %>%
    mutate(result_cnvscore = case_when(
      cnvscore_score >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
      cnvscore_score >= what_is_pathogenic & clinical == 'benign' ~ 'FP')) %>%
    mutate(result_structure = case_when(
      structure_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
      structure_score >= 0.5 & clinical == 'benign' ~ 'FP')) %>%
    select(reliability_score, contains('result')) %>%
    pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TP', 'FP')) %>%
    mutate(tile_general = conditional_interval(fixed_interval, reliability_score)) %>%
    count(tile_general, tool, result) %>%
    pivot_wider(id_cols = c(tile_general, tool),
                names_from = 'result', values_from = 'n') %>%
    group_by(tile_general, tool) %>%
    mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
    mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
    mutate(precision = TP / (TP + FP)) %>%
    ungroup() %>%
    # group_by(tool) %>%
    # mutate(cum_fp = cumsum((FP))) %>%
    # mutate(cum_tp = cumsum((TP))) %>%
    # mutate(precision = cum_tp / (cum_tp + cum_fp)) %>%
    ggplot(aes(as.factor(tile_general), precision)) +
    geom_smooth(aes(group = tool, fill = tool, color = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(limit_y,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (Precision) - ', tag2),
         subtitle = tmp_sub,
         x = 'Reliability score (percentile(residual))', y = 'Precision (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  sensitivity_sd <- x %>% 
    mutate(result_cnvscore = case_when(
      cnvscore_score >= what_is_pathogenic & clinical == 'pathogenic' ~ 'TP',
      cnvscore_score < what_is_pathogenic & clinical == 'pathogenic' ~ 'FN')) %>%
    mutate(result_structure = case_when(
      structure_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
      structure_score < 0.5 & clinical == 'pathogenic' ~ 'FN')) %>%
    select(reliability_score, contains('result')) %>%
    pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TP', 'FN')) %>%
    mutate(tile_general = conditional_interval(fixed_interval, reliability_score)) %>%
    count(tile_general, tool, result) %>%
    pivot_wider(id_cols = c(tile_general, tool), 
                names_from = 'result', values_from = 'n') %>%
    group_by(tile_general, tool) %>%
    mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
    mutate(FN = ifelse(is.na(FN), 0, FN)) %>%
    mutate(npv = TP / (TP + FN)) %>%
    ungroup() %>%
    # group_by(tool) %>%
    # mutate(cum_tp = cumsum((TP))) %>%
    # mutate(cum_fn = cumsum((FN))) %>%
    # mutate(sensitivity = cum_tp / (cum_tp + cum_fn)) %>%
    ggplot(aes(tile_general, sensitivity)) +
    geom_smooth(aes(group = tool, fill = tool, color = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    # scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(limit_y,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (Sensitivity) - ', tag2),
         subtitle = tmp_sub,
         x = 'Reliability score (percentile(residual))', y = 'Sensitivity (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  
  npv_sd <- x %>% 
    mutate(result_cnvscore = case_when(
      cnvscore_score < what_is_benign & clinical == 'benign' ~ 'TN',
      cnvscore_score < what_is_benign & clinical == 'pathogenic' ~ 'FN')) %>%
    mutate(result_structure = case_when(
      structure_score < 0.5 & clinical == 'benign' ~ 'TN',
      structure_score < 0.5 & clinical == 'pathogenic' ~ 'FN')) %>%
    select(reliability_score, contains('result')) %>%
    pivot_longer(-reliability_score, names_to = 'tool', values_to = 'result') %>%
    filter(result %in% c('TN', 'FN')) %>%
    mutate(tile_general = conditional_interval(fixed_interval, reliability_score)) %>%
    count(tile_general, tool, result) %>%
    pivot_wider(id_cols = c(tile_general, tool), 
                names_from = 'result', values_from = 'n') %>%
    group_by(tile_general, tool) %>%
    mutate(FN = ifelse(is.na(FN), 0, FN)) %>%
    mutate(TN = ifelse(is.na(TN), 0, TN)) %>%
    mutate(npv = TN / (TN + FN)) %>%
    ungroup() %>%
    # group_by(tool) %>%
    # mutate(cum_fn = cumsum((FN))) %>%
    # mutate(cum_tn = cumsum((TN))) %>%
    # mutate(npv = cum_tn / (cum_tn + cum_fn)) %>%
    ggplot(aes(tile_general, npv)) +
    geom_smooth(aes(group = tool, fill = tool, color = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    # scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(limit_y,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (NPV) - ', tag2),
         subtitle = tmp_sub,
         x = 'Reliability score (percentile(residual))', y = ' NPV (%)') +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  


  precision_sd +
  npv_sd +
  plot_layout(nrow = 1)

  
}


# Plot TPR and TNR


plot_rates <- function(x, tag, tag2) {
  
  # x <- tmp2
  # tag <- 'DECIPHER'
  # tag2 <- 'Model (both features)'
  
  # x <- tmp2
  # tag <- 'ClinVar'
  # tag2 <- 'Model (both features)'
  
  # By cumulative I mean that at a given SD bin, 
  # you represent not the CNVs within such bin ( SD=XX) but all CNVs 
  # with a value SD<=XX
  
  x %>%
    filter(clinical == 'pathogenic') %>%
    mutate(tile_dist = ntile(sd, 10)) %>%
    ggplot(aes(as.factor(tile_dist), sd)) +
      geom_boxplot()


 test1 <-  x %>%
    filter(clinical == 'pathogenic') %>%
    # mutate(tile_dist = ntile(sd, 5)) %>%
    mutate(tile_dist = cut_interval(cnvscore_score, 10)) %>%
    mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
    select(structure_score, cnvscore_score, tile_dist) %>%
    pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
    count(tile_dist, tool, value) %>%
    pivot_wider(names_from = value, values_from = n) %>%
    mutate(wrong = ifelse(is.na(wrong), 0, wrong)) %>%
    group_by(tool) %>%
    mutate(cum_ok = cumsum(ok)) %>%
    mutate(cum_wrong = cumsum(wrong)) %>%
    mutate(rate = cum_ok / (cum_ok + cum_wrong)) %>%
    ggplot(aes(tile_dist, rate)) +
    geom_line(aes(group = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    # scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (TPR) - ', tag2),
         x = 'CNVscore (Percentile 1-10)', y = 'Cumulative TPR (%)')
 
 
 x %>% 
   mutate(tile_dist = cut_interval(cnvscore_score, 10)) %>%
   mutate(result = case_when(
     cnvscore_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
     cnvscore_score <= 0.5 & clinical == 'benign' ~ 'TN',
     cnvscore_score >= 0.5 & clinical == 'benign' ~ 'FP',
     cnvscore_score <= 0.5 & clinical == 'pathogenic' ~ 'FN'
   )) %>%
   filter(result %in% c('TP', 'FP')) %>%
   count(tile_dist, result) %>%
   pivot_wider(tile_dist, names_from = 'result', values_from = 'n') %>%
   group_by(tile_dist) %>%
   mutate(precision = TP / (TP + FP)) %>%
   ungroup() %>%
   mutate(cum_fp = cumsum(FP)) %>%
   mutate(cum_tp = cumsum(TP)) %>%
   mutate(precision = cum_tp / (cum_tp + cum_fp)) %>%
   ggplot(aes(tile_dist, precision)) +
   # geom_line() +
   geom_point(shape = 21, color = 'black') +
   scale_y_continuous(label = percent) +
   coord_cartesian(ylim = c(0,1)) +
   theme_minimal() +
   labs(title = paste(tag, 'dataset (Precision) - ', tag2),
        x = 'CNVscore (Percentile 1-10)', y = 'Cumulative Precision (%)')
 
 
 
 x %>% 
   # mutate(tile_dist = cut_interval(sd, 10)) %>%
   mutate(tile_dist = ntile(sd, 10)) %>%
   mutate(result = case_when(
     cnvscore_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
     cnvscore_score <= 0.5 & clinical == 'benign' ~ 'TN',
     cnvscore_score >= 0.5 & clinical == 'benign' ~ 'FP',
     cnvscore_score <= 0.5 & clinical == 'pathogenic' ~ 'FN'
   )) %>%
   # filter(result %in% c('TP', 'FP')) %>%
   count(tile_dist, result) %>%
   pivot_wider(tile_dist, names_from = 'result', values_from = 'n') %>%
   mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
   mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
   group_by(tile_dist) %>%
   mutate(precision = TP / (TP + FP)) %>%
   ungroup() %>%
   mutate(cum_fp = cumsum(FP)) %>%
   mutate(cum_tp = cumsum(TP)) %>%
   mutate(precision = cum_tp / (cum_tp + cum_fp)) %>%
   ggplot(aes(tile_dist, precision)) +
   # geom_line() +
   geom_point(shape = 21, color = 'black') +
   scale_y_continuous(label = percent) +
   # coord_cartesian(ylim = c(0,1)) +
   theme_minimal() +
   labs(title = paste(tag, 'dataset (Precision) - ', tag2),
        x = 'SD (Percentile 1-10)', y = 'Cumulative Precision (%)')
 

 # esta es la buena   
 
 x %>% 
   # mutate(tile_dist = cut_interval(sd, 10)) %>%
   # count(tile_sd, clinical)
   mutate(result = case_when(
     cnvscore_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
     cnvscore_score >= 0.5 & clinical == 'benign' ~ 'FP')) %>%
   filter(result %in% c('TP', 'FP')) %>%
   mutate(tile_dist = ntile(sd, 10)) %>%
   count(tile_dist, result) %>%
   pivot_wider(tile_dist, names_from = 'result', values_from = 'n') %>%
   mutate(FP = ifelse(is.na(FP), 0, FP)) %>%
   mutate(TP = ifelse(is.na(TP), 0, TP)) %>%
   group_by(tile_dist) %>%
   mutate(precision = TP / (TP + FP)) %>%
   ungroup() %>%
   ggplot(aes(tile_dist, precision)) +
   # geom_line() +
   geom_point(shape = 21, color = 'black') +
   scale_y_continuous(label = percent) +
   # coord_cartesian(ylim = c(0,1)) +
   theme_minimal() +
   labs(title = paste(tag, 'dataset (Precision) - ', tag2),
        x = 'SD (Percentile 1-10)', y = 'Precision (%)')
  # esta es la buena 2 
 x %>% 
   mutate(result_cnvscore = case_when(
     cnvscore_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
     cnvscore_score >= 0.5 & clinical == 'benign' ~ 'FP')) %>%
   mutate(result_structure = case_when(
     structure_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
     structure_score >= 0.5 & clinical == 'benign' ~ 'FP')) %>%
   select(sd, contains('result')) %>%
   pivot_longer(-sd, names_to = 'tool', values_to = 'result') %>%
   filter(result %in% c('TP', 'FP')) %>%
   mutate(tile_dist = ntile(sd, 10)) %>%
   count(tile_dist, tool, result) %>%
   pivot_wider(id_cols = c(tile_dist, tool), 
               names_from = 'result', values_from = 'n') %>%
   group_by(tile_dist, tool) %>%
   mutate(precision = TP / (TP + FP)) %>%
   ungroup() %>%
   ggplot(aes(tile_dist, precision)) +
   geom_line(aes(group = tool)) +
   geom_point(aes(fill = tool), shape = 21, color = 'black') +
   scale_x_continuous(breaks = pretty_breaks()) +
   scale_y_continuous(label = percent) +
   coord_cartesian(ylim = c(0,1)) +
   theme_minimal() +
   labs(title = paste(tag, 'dataset (Precision) - ', tag2),
        x = 'SD (Percentile 1-10)', y = 'Cumulative precision (%)')
 
 
   x %>%
   # filter(clinical == 'pathogenic') %>%
   # mutate(tile_dist = ntile(sd, 5)) %>%
   mutate(tile_dist = cut_interval(cnvscore_score, 10)) %>%
   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
   select(structure_score, cnvscore_score, tile_dist) %>%
   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
   count(tile_dist, tool, value) %>%
   pivot_wider(names_from = value, values_from = n) %>%
   mutate(wrong = ifelse(is.na(wrong), 0, wrong)) %>%
   group_by(tool) %>%
   mutate(cum_ok = cumsum(ok)) %>%
   mutate(cum_wrong = cumsum(wrong)) %>%
   mutate(rate = cum_ok / (cum_ok + cum_wrong)) %>%
   ggplot(aes(tile_dist, rate)) +
   geom_line(aes(group = tool)) +
   geom_point(aes(fill = tool), shape = 21, color = 'black') +
   # scale_x_continuous(breaks = pretty_breaks()) +
   scale_y_continuous(label = percent) +
   coord_cartesian(ylim = c(0,1)) +
   theme_minimal() +
   labs(title = paste(tag, 'dataset (TPR) - ', tag2),
        x = 'CNVscore (Percentile 1-10)', y = 'Cumulative TPR (%)')
 
 
 
 x %>% filter(clinical == 'pathogenic') %>% select(cnvscore_score, sd) %>% correlate(method = 'spearman')
 
 x %>% filter(clinical == 'pathogenic') %>% ggplot(aes(cnvscore_score, sd)) + geom_point()
 
 x %>% 
   mutate(test = if_else(cnvscore_score >= 0.5, 'predicted_pathogenic', 'predicted_benign')) %>%
   filter(test == 'predicted_pathogenic') %>%
   count(clinical) %>%
   mutate(perc = n / sum(n))
  
  x %>%
    filter(clinical == 'pathogenic') %>%
    mutate(tile_dist = ntile(sd, 5)) %>%
    mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
    select(structure_score, cnvscore_score, tile_dist) %>%
    pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
    count(tile_dist, tool, value) %>%
    pivot_wider(names_from = value, values_from = n) %>%
    mutate(wrong = ifelse(is.na(wrong), 0, wrong)) %>%
    group_by(tool) %>%
    mutate(cum_ok = cumsum(ok)) %>%
    mutate(cum_wrong = cumsum(wrong)) %>%
    mutate(rate = cum_ok / (cum_ok + cum_wrong)) %>%
    ggplot(aes(tile_dist, rate)) +
    geom_line(aes(group = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (TPR) - ', tag2),
         x = 'SD (Percentile 1-10)', y = 'Cumulative TPR (%)')


  tpr_clinvar_both_dist <- x %>% 
    filter(clinical == 'pathogenic') %>%
    mutate(tile_dist = ntile(dist_pathogenic, 10)) %>%
    mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
    select(structure_score, cnvscore_score, tile_dist) %>%
    pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
    count(tile_dist, tool, value) %>%
    pivot_wider(names_from = value, values_from = n) %>%
    mutate(wrong = ifelse(is.na(wrong), 0, wrong)) %>%
    group_by(tool) %>%
    mutate(cum_ok = cumsum(ok)) %>%
    mutate(cum_wrong = cumsum(wrong)) %>%
    mutate(rate = cum_ok / (cum_ok + cum_wrong)) %>%
    ggplot(aes(tile_dist, rate)) +
    geom_line(aes(group = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (TPR) - ', tag2),
         x = 'Pathogenic distance (Percentile 1-10)', y = 'Cumulative TPR (%)')
  
  
  tnr_clinvar_both_dist <- x %>% 
    filter(clinical == 'benign') %>%
    mutate(tile_dist = ntile(dist_benign, 10)) %>%
    mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
    select(structure_score, cnvscore_score, tile_dist) %>%
    pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
    count(tile_dist, tool, value) %>%
    pivot_wider(names_from = value, values_from = n) %>%
    mutate(wrong = ifelse(is.na(wrong), 0, wrong)) %>%
    group_by(tool) %>%
    mutate(cum_ok = cumsum(ok)) %>%
    mutate(cum_wrong = cumsum(wrong)) %>%
    mutate(rate = cum_ok / (cum_ok + cum_wrong)) %>%
    ggplot(aes(tile_dist, rate)) +
    geom_line(aes(group = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal() +
    labs(title =  paste(tag, 'dataset (TNR) - ', tag2),
         x = 'Benign distance (Percentile 1-10)', y = 'Cumulative TNR (%)')
  
  
  tpr_clinvar_both_sd <- x %>% 
    filter(clinical == 'pathogenic') %>%
    mutate(tile_dist = ntile(sd, 10)) %>%
    mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
    select(structure_score, cnvscore_score, tile_dist) %>%
    pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
    count(tile_dist, tool, value) %>%
    pivot_wider(names_from = value, values_from = n) %>%
    mutate(wrong = ifelse(is.na(wrong), 0, wrong)) %>%
    group_by(tool) %>%
    mutate(cum_ok = cumsum(ok)) %>%
    mutate(cum_wrong = cumsum(wrong)) %>%
    mutate(rate = cum_ok / (cum_ok + cum_wrong)) %>%
    ggplot(aes(tile_dist, rate)) +
    geom_line(aes(group = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (TPR) - ', tag2),
         x = 'SD (Percentile 1-10)', y = 'Cumulative TPR (%)')
  
  
  tnr_clinvar_both_sd <- x %>% 
    filter(clinical == 'benign') %>%
    mutate(tile_dist = ntile(sd, 10)) %>%
    mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
    mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
    mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
    select(structure_score, cnvscore_score, tile_dist) %>%
    pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
    count(tile_dist, tool, value) %>%
    pivot_wider(names_from = value, values_from = n) %>%
    mutate(wrong = ifelse(is.na(wrong), 0, wrong)) %>%
    group_by(tool) %>%
    mutate(cum_ok = cumsum(ok)) %>%
    mutate(cum_wrong = cumsum(wrong)) %>%
    mutate(rate = cum_ok / (cum_ok + cum_wrong)) %>%
    ggplot(aes(tile_dist, rate)) +
    geom_line(aes(group = tool)) +
    geom_point(aes(fill = tool), shape = 21, color = 'black') +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(label = percent) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal() +
    labs(title = paste(tag, 'dataset (TNR) - ', tag2),
         x = 'SD (Percentile 1-10)', y = 'Cumulative TNR (%)')
  
  
  tpr_clinvar_both_dist + 
    tnr_clinvar_both_dist + 
    tpr_clinvar_both_sd + 
    tnr_clinvar_both_sd + plot_layout(nrow = 2)
  
  
}



# ------------------------------------------------------------------------------
# FEATURES ENRICHMENT 
# ------------------------------------------------------------------------------

features_enrichment <- function(df_features, training_del, training_dup) {
  
  # df_features <- features_tbl[features_tbl$human_control == 'no',]
  # training_del <- output_clinvar_deletion
  # training_dup <- output_clinvar_duplication
   
  features_categorical <- df_features %>% filter(type == 'categorical') %>% pull(name)
  features_continuous <- df_features %>% filter(type == 'numerical' | type == 'percentage') %>% pull(name)
  
  def_results_del <- tibble()
  
  for (i in 1:length(features_categorical)) {
    print(i)
    tmp_result <- training_del %>% 
      rename(target = features_categorical[i]) %>%
      filter(type_variant == 'deletion') %>% 
      mutate(target = if_else(target == 0, 'no', 'yes')) %>%
      count(clinical, target) %>%
      pivot_wider(names_from = target, values_from = n) %>%
      relocate(clinical, yes, no) %>%
      column_to_rownames('clinical') %>%
      fisher.test() %>%
      glance() %>%
      mutate(feature = features_categorical[i])
    
    def_results_del <- def_results_del %>% bind_rows(tmp_result)
  }
  
  
  def_results_dup <- tibble()
  for (i in 1:length(features_categorical)) {
    print(i)
    tmp_result <- training_dup %>% 
      rename(target = features_categorical[i]) %>%
      filter(type_variant == 'duplication') %>% 
      mutate(target = if_else(target == 0, 'no', 'yes')) %>%
      count(clinical, target) %>%
      pivot_wider(names_from = target, values_from = n) %>%
      relocate(clinical, yes, no) %>%
      mutate(yes = ifelse(is.na(yes), 0, yes)) %>%
      column_to_rownames('clinical') %>%
      fisher.test() %>%
      glance() %>%
      mutate(feature = features_categorical[i])
    
    def_results_dup <- def_results_dup %>% bind_rows(tmp_result)
  }
  
  

  
  p_enrich1 <- def_results_del %>% mutate(variant_type = 'Deletion') %>% 
    bind_rows(def_results_dup %>% mutate(variant_type = 'Duplication')) %>%
    mutate(sig = ifelse(p.value <= (0.05 / length(features_categorical)), 'yes', 'no')) %>%
    mutate(sig = factor(sig, levels = c('yes', 'no'))) %>%
    mutate(across(where(is.double), log10)) %>%
    mutate(p.value = 10^p.value) %>%
    mutate(p.value = ifelse(sig == 'no', 
                            round(p.value, 2), 
                            format(p.value, scientific = TRUE, digits=1))) %>%
    mutate(variant_type = paste(variant_type, '\n', 'Non-pathogenic vs. Pathogenic CNVs')) %>%
    left_join(df_features %>% select(name, complete_name), by = 
                c('feature' = 'name')) %>%
    # select(feature, complete_name)
    ggplot(aes(estimate, reorder(complete_name, estimate))) +
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed", color = 'gray50') +
    geom_errorbarh(aes(xmax = conf.high, xmin = conf.low), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(aes(color = sig), size = 2, show.legend = FALSE) +
    geom_text(aes(label = p.value), nudge_y = 0.3,
              size = 3) +
    scale_color_manual(values = c('#F8766D', '#7f7f7f')) +

    facet_wrap(vars(variant_type)) +
    labs(y = 'Unbiased features', 
         x = expression('Enrichment (log'[10]~'odds ratio)')) +
    theme_minimal()
  
  def_results_cont_del <- tibble()
  
  for (i in 1:length(features_continuous)) {
    
    tmp_cont_wilcoxon <- training_del %>%
      filter(type_variant == 'deletion') %>%
      mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign'))) %>%
      select(clinical, contains(features_continuous)) %>%
      rename(target = features_continuous[i])
    
    tmp_pvalue <- wilcox.test(tmp_cont_wilcoxon[tmp_cont_wilcoxon$clinical == 'pathogenic',]$target,
                              tmp_cont_wilcoxon[tmp_cont_wilcoxon$clinical == 'benign',]$target) %>% tidy() %>% pull(p.value)
    
    def_results_cont_del <- def_results_cont_del %>% bind_rows(tibble(feature = features_continuous[i], p.value = tmp_pvalue))
  }
  
  
  def_results_cont_dup <- tibble()
  
  for (i in 1:length(features_continuous)) {
    
    tmp_cont_wilcoxon <- training_dup %>%
      filter(type_variant == 'duplication') %>%
      mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign'))) %>%
      select(clinical, contains(features_continuous)) %>%
      rename(target = features_continuous[i])
    
    tmp_pvalue <- wilcox.test(tmp_cont_wilcoxon[tmp_cont_wilcoxon$clinical == 'pathogenic',]$target,
                              tmp_cont_wilcoxon[tmp_cont_wilcoxon$clinical == 'benign',]$target) %>% tidy() %>% pull(p.value)
    
    def_results_cont_dup <- def_results_cont_dup %>% bind_rows(tibble(feature = features_continuous[i], p.value = tmp_pvalue))
  }
  
  p_enrich2 <- def_results_cont_del %>% 
    mutate(variant_class = 'Deletion') %>%
    bind_rows(def_results_cont_dup %>% mutate(variant_class = 'Duplication')) %>%
    mutate(p.value = if_else(p.value == 0, 1e-300, p.value)) %>%
    mutate(sig = ifelse(p.value <= (0.05 / length(features_continuous)), 'yes', 'no')) %>%
    mutate(p.value = -log(p.value)) %>% 
    mutate(variant_class = paste(variant_class, '\n', 'Non-pathogenic vs. Pathogenic CNVs')) %>%
    left_join(df_features %>% select(name, complete_name), by = 
                c('feature' = 'name')) %>%
    # select(feature, complete_name) %>% filter(is.na(complete_name))
    ggplot(aes(reorder(complete_name, p.value), p.value)) + 
    geom_col(aes(fill = sig), color = 'black', show.legend = FALSE) + 
    scale_fill_manual(values = c('#7f7f7f', '#F8766D')) +
    geom_hline(yintercept = -log10(0.05/length(features_continuous)),
               size = .5, linetype = "dashed") +
    coord_flip() + 
    facet_wrap(vars(variant_class), scales = 'free_x') +
    theme_minimal() +
    labs(y = expression('-log'[10]~'(p.value)'), x = 'Unbiased features')
  
  p_enrich1 + p_enrich2
}


# ------------------------------------------------------------------------------
# FEATURES ENRICHMENT 
# ------------------------------------------------------------------------------


get_reliability_score_rank_sd <- function(ref, predictions, trendline) {
  
  
  # ref <- ref_quantiles
  # predictions <-   ref_sd_scenario_4[[3]]
  # trendline <- trendline_clinvar
  
  n_total_rows <- nrow(predictions)
  
  vector_reliability_score <- c()
  
  for (i in 1:nrow(predictions)) {
    #   
      print(paste0(i, '/', n_total_rows))

      tmp_score_bin <- ref %>%
        select(score_interval, min_score, max_score) %>%
        pivot_longer(-score_interval, names_to = 'tag', values_to = 'value') %>%
        mutate(dist_score = abs(predictions[i,]$.pred_pathogenic - value)) %>%
        arrange(dist_score) %>%
        slice_head(n = 1) %>%
        pull(score_interval)

      tmp_to_vector <- ref %>%
        filter(score_interval == tmp_score_bin) %>%
        select(sd) %>%
        mutate(tag = 'ref') %>%
        bind_rows(tibble(sd = predictions[i,]$sd,
                         tag = 'obs')) %>%
        mutate(chunk = ntile(sd, 3)) %>%
        filter(tag == 'obs') %>%
        pull(chunk)

      vector_reliability_score <- c(vector_reliability_score, tmp_to_vector)
    
  }
  
  predictions <- predictions %>%
    mutate(reliability_score = vector_reliability_score)
  
  return(predictions)
  
}

# ------------------------------------------------------------------------------
# FEATURES ENRICHMENT 
# ------------------------------------------------------------------------------


get_reliability_score <- function(ref, predictions, trendline, type = 'res') {
  
  
  # ref <- ref_quantiles
  # predictions <-   ref_sd_clinvar_del[[3]]
  # trendline <- trendline_clinvar
  
  n_total_rows <- nrow(predictions)
  
  vector_reliability_score <- c()
  
  for (i in 1:nrow(predictions)) {
    #   

        tmp_score_bin <- ref %>%
          distinct(score_interval, min_score, max_score) %>%
          filter(min_score < predictions[i,]$.pred_pathogenic &
                 max_score > predictions[i,]$.pred_pathogenic) %>%
          pull(score_interval)
        
        if (length(tmp_score_bin) == 0) {

          tmp_score_bin <- ref %>%
            select(score_interval, min_score, max_score) %>%
            mutate(mid_score = (min_score + max_score)/2) %>%
            pivot_longer(-c(score_interval, mid_score), names_to = 'tag',
                         values_to = 'value') %>%
            mutate(dist_score = abs(predictions[i,]$.pred_pathogenic - value)) %>%
            mutate(dist_midscore = abs(predictions[i,]$.pred_pathogenic - mid_score)) %>%
            arrange(dist_midscore) %>%
            slice_head(n = 1) %>%
            pull(score_interval)

        }
        
        
        if (type == 'sd') {
    
          print(paste0(i, '/', n_total_rows, 'SD'))
          
    tmp_to_vector <- ref %>%
      filter(score_interval == tmp_score_bin) %>%
      select(sd) %>%
      mutate(tag = 'ref') %>%
      bind_rows(tibble(sd = predictions[i,]$sd,
                       tag = 'obs')) %>%
      mutate(chunk = ntile(sd, 3)) %>%
      filter(tag == 'obs') %>%
      pull(chunk)
    
        } else {
          
          print(paste0(i, '/', n_total_rows, '- RESIDUALS'))
          
          
          tmp_y <- predict(trendline,
                           tibble(.pred_pathogenic = predictions[i,]$.pred_pathogenic))
          
          tmp_residual <- as.numeric(predictions[i,]$sd - tmp_y)
          
          tmp_to_vector <- ref %>%
            filter(score_interval == tmp_score_bin) %>%
            select(gam_residuals) %>%
            mutate(tag = 'ref') %>%
            bind_rows(tibble(gam_residuals = tmp_residual,
                             tag = 'obs')) %>%
            mutate(chunk = ntile(gam_residuals, 3)) %>%
            filter(tag == 'obs') %>%
            pull(chunk)

        }
    
    vector_reliability_score <- c(vector_reliability_score, tmp_to_vector)
    
  }
  
  predictions <- predictions %>%
    mutate(reliability_score = vector_reliability_score)
  
  return(predictions)
  
}


get_reliability_score_mid <- function(ref, predictions) {
  
  # ref <- ref_quantiles
  # predictions <- ref_sd_clinvar_del[[3]]
  
  
  n_total_rows <- nrow(predictions)
  
  vector_reliability_score <- c()
  
  for (i in 1:nrow(predictions)) {
    #   
    # print(paste0(i, '/', n_total_rows))
    
        # tmp_y <- predict(trendline_clinvar,
        #         tibble(.pred_pathogenic = predictions[i,]$.pred_pathogenic))
        # 
        # tmp_residual <- as.numeric(predictions[i,]$sd - tmp_y)
    
    tmp_score_bin <- ref %>%
      distinct(score_interval, min_score, max_score) %>%
      mutate(mid_score = (max_score - min_score)/2) %>%
      mutate(diff_mid_score = abs(mid_score - predictions[i,]$.pred_pathogenic)) %>%
        arrange(diff_mid_score) %>%
        slice_head(n = 1) %>%
      pull(score_interval)
    
    tmp_to_vector <- ref %>%
      filter(score_interval == tmp_score_bin) %>%
      group_by(reliability_score) %>%
      # summarise(median_sd = median(sd)) %>%
      summarise(max_sd = max(sd),
                min_sd = min(sd)) %>%
      mutate(mid_sd =  (max_sd - min_sd)/2) %>%
      mutate(diff_mid_sd = abs(mid_sd - predictions[i,]$sd)) %>%
      arrange(diff_mid_sd) %>%
      slice_head(n = 1) %>%
      pull(reliability_score)
    
    
    # tmp_to_vector <- ref %>%
    #   filter(score_interval == tmp_score_bin) %>%
    #   group_by(reliability_score) %>%
    #   summarise(max_res = max(gam_residuals),
    #             min_res = min(gam_residuals)) %>%
    #   mutate(mid_res =  (max_res - min_res)/2) %>%
    #   mutate(diff_mid_res = abs(mid_res - tmp_residual)) %>%
    #   arrange(diff_mid_res) %>%
    #   slice_head(n = 1) %>%
    #   pull(reliability_score)
    

    
    vector_reliability_score <- c(vector_reliability_score, tmp_to_vector)
    
  }
  
  predictions <- predictions %>%
    mutate(reliability_score = vector_reliability_score)
  
}

# get_reliability_score_res <- function(ref, predictions, trendline) {
#   
#   
#   # ref <- ref_quantiles
#   # predictions <-   ref_sd_scenario_4[[3]]
#   # trendline <- trendline_clinvar
#   
#   n_total_rows <- nrow(predictions)
#   
#   vector_reliability_score <- c()
#   
#   for (i in 1:nrow(predictions)) {
# 
#     tmp_y <- predict(trendline_clinvar,
#             tibble(.pred_pathogenic = predictions[i,]$.pred_pathogenic))
# 
#     tmp_residual <- as.numeric(predictions[i,]$sd - tmp_y)
# 
#     tmp_score_bin <- ref %>%
#       select(score_interval, min_score, max_score) %>%
#       mutate(mid_score = (min_score + max_score)/2) %>%
#       pivot_longer(-c(score_interval, mid_score), names_to = 'tag', 
#                    values_to = 'value') %>%
#       mutate(dist_score = abs(predictions[i,]$.pred_pathogenic - value)) %>%
#       mutate(dist_midscore = abs(predictions[i,]$.pred_pathogenic - mid_score)) %>%
#       arrange(dist_score, dist_midscore) %>%
#       slice_head(n = 1) %>%
#       pull(score_interval)
# 
# 
#     tmp_to_vector <- ref %>%
#       filter(score_interval == tmp_score_bin) %>%
#       select(reliability_score, min_gam_residuals, max_gam_residuals) %>%
#       mutate(midpoint_residual = (max_gam_residuals + min_gam_residuals)/2) %>%
#       pivot_longer(-c(reliability_score, midpoint_residual),
#                    names_to = 'tag', values_to = 'value') %>%
#       mutate(dist_residual = abs(tmp_residual - value)) %>%
#       mutate(dist_midpoint_residual =  abs(tmp_residual - midpoint_residual)) %>%
#       arrange(dist_residual, dist_midpoint_residual) %>%
#       slice_head(n = 1) %>%
#       pull(reliability_score)
# 
#     vector_reliability_score <- c(vector_reliability_score, tmp_to_vector)
#     
#   }
#   
#   predictions <- predictions %>%
#     mutate(reliability_score = vector_reliability_score)
#   
#   return(predictions)
# 
# }

# ------------------------------------------------------------------------------
# TO GOOGLESHEET 
# ------------------------------------------------------------------------------

to_googlesheet2 <- function(x, y) {
  
  # x <- result_clinvar_del
  # y <- output_clinvar_deletion
  
  to_sheet1 <- x %>% 
    left_join(y %>% select(id, clinical)) %>%
    group_by(tag) %>%
    roc_auc(clinical, .pred_pathogenic) %>%
    rename(auroc = .estimate) %>%
    select(tag, auroc) %>%
    mutate(auroc = round(auroc, 3))
  
  to_sheet2 <- x %>% 
    left_join(y %>% select(id, clinical), by = 'id') %>%
    group_by(tag) %>%
    pr_auc(clinical, .pred_pathogenic) %>%
    rename(auprc = .estimate) %>%
    select(tag, auprc) %>%
    mutate(auprc = round(auprc, 3))
  
  
  to_sheet1 %>%
    left_join(to_sheet2, by = 'tag')
}


# ------------------------------------------------------------------------------
# FUNCTION - CLEAN TAG TABLES
# ------------------------------------------------------------------------------

cleaning_tables <- function(x) {
  
  x %>%
    filter(! tag %in% c('bayesian_both', 'bayesian_knowledge_based',
                        'gbm_knowledge_based', 'gbm_unbiased', 'gwrvis',
                        'jarvis')) %>%
    arrange(tag, desc = TRUE) %>%
    mutate(tag = if_else(tag == 'bayesian_unbiased', 'CNVscore', tag)) %>%
    mutate(tag = if_else(tag == 'cadd', 'CADD-SV', tag)) %>%
    mutate(tag = if_else(tag == 'classifycnv', 'ClassifyCNV', tag)) %>%
    mutate(tag = if_else(tag == 'strvctvre', 'STRVCTVRE', tag)) %>%
    mutate(tag = if_else(tag == 'tada', 'TADA', tag)) %>%
    mutate(tag = if_else(tag == 'xcnv', 'X-CNV', tag)) %>%
    mutate(tag = if_else(tag == 'naive_model_length', 'Naïve model (CNV size)', tag)) %>%
    mutate(tag = if_else(tag == 'naive_model_n_genes', 'Naïve model (Nº genes)', tag)) %>%
    mutate(tag = if_else(tag == 'naive_model_omim', 'Naïve model (OMIM genes)', tag))
}


# ------------------------------------------------------------------------------
# FUNCTION - GENERATE TABLE - RELIABILITY SCORES
# ------------------------------------------------------------------------------


generate_table_rel <- function(x, y, tag = '', tag2 = 'Deletion') {
  
  # x <- result_clinvar_del
  # y <- ref_sd_clinvar_del_real
  
  
  tmp_all_auroc <- y %>% 
    select(id, clinical) %>%
    left_join(x, by = 'id') %>%
    group_by(tag) %>%
    roc_auc(clinical, .pred_pathogenic) %>%
    rename(All = .estimate) %>%
    select(tag, All)
  
  y %>% 
    select(id, reliability_score, clinical) %>%
    left_join(x, by = 'id') %>%
    group_by(tag, reliability_score) %>%
    roc_auc(clinical, .pred_pathogenic) %>%
    rename(auroc = .estimate) %>%
    select(tag, reliability_score, auroc) %>%
    mutate(reliability_score = case_when(
      reliability_score == 1 ~ 'Low',
      reliability_score == 2 ~ 'Medium',
      reliability_score == 3 ~ 'High'
    )) %>%
    pivot_wider(tag, names_from = reliability_score, values_from = auroc)  %>%
    left_join(tmp_all_auroc) %>%
    select(tag, All, Low, Medium, High) %>%
    cleaning_tables() %>%
    mutate(`Diff. Low-High` = Low - High) %>%
    arrange(desc(All)) %>%
    rename(Model = tag) %>%
    filter( Model %in% c('TADA', 'STRVCTVRE', 'CNVscore', 'X-CNV', 'CADD-SV')) %>%
    map_dfr(function(x) {
      if (is.numeric(x)) {
        round(x,3)*100
      } else {
        x
      }
        
    }
    ) %>%
    
    gt() %>%
    tab_header(
      title = glue("{tag} - {tag2} CNVs (n={format(nrow(y), big.mark = ',')})"),
    ) %>%
    cols_align(align = 'center')
  
}

# ------------------------------------------------------------------------------
# FUNCTION - GENERATE TABLE
# ------------------------------------------------------------------------------

generate_table <- function(x, y, tag = '', as_table = TRUE) {
  
  # color_models <- tibble('Model' =
  #   c('CNVscore', 'X-CNV', 'STRVCTVRE', 'TADA', 'Naïve model (OMIM genes)',
  #     'CADD-SV', 'ClassifyCNV', 'Naïve model (Nº genes)',
  #     'Naïve model (CNV size)'),
  #   'color' = c("#b6e3ff", "#b6e3ff","#b6e3ff", "#b6e3ff", "#FEF0D9", "#b6e3ff",
  #               "#b6e3ff", "#FEF0D9", "#FEF0D9")
  # )
  # 
  # color_models <- col_factor(palette = color_models$color,
  #                            domain = color_models$Model)
  
  x <- result_clinvar_del
  y <- output_clinvar_deletion 
  tag <- 'ClinVar'
  as_table <- FALSE
  
 tmp_df <-  to_googlesheet2(x, y) %>%
    cleaning_tables() %>%
    arrange(desc(auroc)) %>%
   map_dfr(function(x) {
     
     if (is.numeric(x)) {
       round(x,3)*100
     } else {
       x
     }
     
   }) %>%
    rename(Model = tag, AUROC = auroc, AUPR = auprc)
 
 if (isTRUE(as_table)) {
   
   tmp_df <- tmp_df %>%
     gt() %>%
     tab_header(
       title = glue("{tag} - Deletion CNVs (n={format(nrow(y), big.mark = ',')})"),
     ) %>%
     cols_align(align = 'center')
 }
 
 return(tmp_df)

}


# ------------------------------------------------------------------------------
# FUNCTION - GET AUROC VALUES
# ------------------------------------------------------------------------------

# get_auroc_tools <- function(x, y, input_tag = '') {
#   
#   x <- result_clinvar_del
#   y <- ref_sd_clinvar_del_real
#   
#   x %>%
#     filter(tag %in% c('xcnv', 'strvctvre')) %>%
#     left_join(y %>% select(id, reliability_score, clinical)) %>%
#     count(reliability_score)
#   
# }


# ------------------------------------------------------------------------------
# FUNCTION - RETRIEVE RISK AND SUPPORT FOR EACH RULE DECISION
# ------------------------------------------------------------------------------

generate_risk_support <- function(x, y) {
  
  # x <- bayesian_clinvar_del_nohuman
  # chrom_tmp <- 16
  # vector_chrom <- 3
  # df_predict <- output_decipher_deletion
  
  # x <- bayesian_clinvar_del_nohuman
  # y <- output_clinvar_deletion
  
  tibble_tmp <- tibble()
  for (i in 1:23) {
    print(i)
    tmp_binary <- generate_rtemis(output_clinvar_deletion, 
                                  x[[i]]$set_rules)
    
    coeff_tmp <- x[[i]]$model_trained$coefficients %>% 
      as_tibble(rownames = 'rule_id') %>%
      rename(coefficient = value)
    
    risk_tmp <- tmp_binary %>%
      pivot_longer(-clinical, names_to = 'rule_id', values_to = 'value') %>%
      count(clinical, rule_id, value) %>% 
      filter(value == 1) %>% 
      group_by(rule_id) %>%
      mutate(perc = n / sum(n)) %>%
      ungroup() %>% 
      filter(clinical == 'pathogenic') %>%
      rename(risk = perc) %>%
      select(rule_id, risk)
    
    support_tmp <- tmp_binary %>%
      pivot_longer(-clinical, names_to = 'rule_id', values_to = 'value') %>%
      count(rule_id, value) %>%
      filter(value == 1) %>%
      rename(support = n) %>%
      select(rule_id, support)
    
    tibble_tmp <- tibble_tmp %>%
      bind_rows(
        x[[i]]$set_rules %>%
          select(-coefficient) %>%
          mutate(chrom = x[[i]]$chrom_target) %>%
          left_join(risk_tmp, by = 'rule_id', ) %>%
          replace_na(list(risk = 0)) %>%
          left_join(support_tmp, by = 'rule_id') %>%
          left_join(coeff_tmp, by = 'rule_id')
      )
    
  }
  
  return(tibble_tmp)
  
}

# almost_clinvar %>%
#   select(cnvscore_score, sd, reliability_score) %>%
#   mutate(dataset = 'Clinvar') %>%
#   bind_rows(almost_decipher %>%
#               select(cnvscore_score, sd, reliability_score) %>%
#               mutate(dataset = 'Decipher')
#               ) %>%
#   ggplot(aes(reliability_score, sd))

# almost_decipher %>%
#   mutate(result = case_when(
#     cnvscore_score >= 0.5 & clinical == 'pathogenic' ~ 'TP',
#     cnvscore_score < 0.5 & clinical == 'benign' ~ 'TN',
#     cnvscore_score >= 0.5 & clinical == 'benign' ~ 'FP',
#     cnvscore_score < 0.5 & clinical == 'pathogenic' ~ 'FN'
#   )) %>%
#   count(sd, result) %>%
#   pivot_wider(id_cols = sd, names_from = result, values_from = n) %>%
#   mutate(precision = TP / (TP + FP)) %>%
#   mutate(npv = TN / (TN + FN)) 
  
  
  

# almost_decipher %>%
#   count(sd, clinical) %>%
#   group_by(sd) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(sd, perc)) +
#     geom_col(aes(fill = clinical), color = 'black')
# 
# almost_decipher %>%
#   count(sd, clinical)
  

# tmp2 %>%
#   mutate(quant_residuals = vector_quantiles_decipher) %>%
#   # filter(clinical == 'benign') %>% 
#   ggplot(aes(cnvscore_score, sd)) +
#   geom_point(aes(color = quant_residuals), show.legend = T) +
#   geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cr")) +
#   theme_minimal() +
#   scale_color_distiller(palette = "Spectral") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   facet_wrap(vars(clinical))
# 
# 
# tmp2 %>%
#   mutate(quant_residuals = vector_quantiles_decipher) %>%
#   arrange(sd) %>%
#   View()


# tmp2 %>%
#   mutate(quant_residuals = vector_quantiles_decipher) %>%
#   ggplot(aes(cnvscore_score, quant_residuals)) +
#     geom_point()
#  

# 
# 
# tmp2 %>%
#   mutate(quant_residuals = vector_quantiles_decipher)  %>%
#   # filter(clinical == 'benign') %>% 
#   ggplot(aes(cnvscore_score, sd)) +
#   # geom_point(aes(fill = quant_residuals), color = 'black',  shape = 21, show.legend = T) +
#   geom_point(aes(color = quant_residuals),  show.legend = T) +
#   geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cr")) +
#   theme_minimal() +
#   scale_color_distiller(palette = "Spectral") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
#   facet_wrap(vars(clinical))

# test1 <- res_df %>% filter(source == 'decipher') %>% select(quant_residuals)
# test1 <- tmp2 %>% 
#   select(-sd, -dist_pathogenic, -dist_benign) %>%
#   bind_cols(test1 %>% rename(sd = quant_residuals))





# comp_structure_clinvar$tmp_predicted %>%
#   select(.pred_pathogenic, sd, mad) %>%
#   ggplot(aes(.pred_pathogenic, sd)) +
#     geom_point() +
#     geom_smooth()
# 
# comp_structure_clinvar$tmp_predicted %>%
#   select(.pred_pathogenic, sd, mad) %>%
#   ggplot(aes(.pred_pathogenic, mad)) +
#   geom_point() +
#   geom_smooth()
# 
# comp_structure_clinvar$tmp_predicted %>% 
#   select(.pred_pathogenic, sd) %>% correlate(method = 'spearman')
# comp_structure_clinvar$tmp_predicted %>% 
#   select(.pred_pathogenic, mad) %>% correlate(method = 'spearman')


## Svscore------------------------

# 
# result_svscore_manolo_deletion <- 1:22 %>% map_dfr(function(x) {
#   
#   first_df <- readLines(glue("rival_cnvscore/svscore/output_{x}.del.vcf"))
#   
#   tmp_df <-  read_tsv(glue('rival_cnvscore/svscore/output_{x}.del.vcf'), skip = sum(grepl('^##', first_df)))
#   
#   tmp_df <-  read_tsv(glue('rival_cnvscore/svscore/output_{x}.del.vcf'), skip = sum(grepl('^##', first_df)))
# 
#   if (nrow(tmp_df) == 0) return(tibble())
#   
#   tmp_df %>% 
#     rename(id = ID) %>%
#     mutate(.pred_pathogenic = str_extract(INFO, 'SVSCOREMAX=(.*?);')) %>%
#     mutate(.pred_pathogenic = str_remove(.pred_pathogenic, ';')) %>%
#     mutate(.pred_pathogenic = str_remove(.pred_pathogenic, 'SVSCOREMAX=')) %>%
#     mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
#     select(id, .pred_pathogenic)
#   
# })
# 
# result_svscore_manolo_deletion <- result_svscore_manolo_deletion %>%
#   left_join(df_manolo %>% select(id, clinical), by = 'id')
# 
# 
# svscore_del_manolo_auc <- just_auc(result_svscore_manolo_deletion)


# res_df <- df_for_uncertainty[[3]] %>%
#   mutate(gam_residuals = residuals(a)) %>%
#   # mutate(quant_residual = gam_residuals)
#   mutate(score_interval = cut_interval(.pred_pathogenic, n = 10)) %>%
#   group_by(score_interval) %>%
#   # group_by(clinical, score_interval) %>%
#   mutate(quant_residuals = ntile(gam_residuals, n = 100) ) %>%
#   mutate(predict_clinical = if_else(.pred_pathogenic >= 0.5, 'pathogenic', 'benign')) %>%
#   group_by(predict_clinical) %>%
#   mutate(rank_score = ntile(.pred_pathogenic, n = 100)) %>%
#   ungroup() %>%
#   mutate(rank_score =
#            ifelse(predict_clinical == 'pathogenic', 100 - rank_score, rank_score)) %>%
#   mutate(quant_residuals = rank_score + (1 * quant_residuals))


# mutate(score_interval = str_remove(score_interval, ']')) %>%
# mutate(score_interval = str_remove(score_interval, '\\[')) %>%
# mutate(score_interval = str_remove(score_interval, '\\(')) %>%
# separate(score_interval, 
#          into = c('min_interval', 'max_interval'), sep = ',') %>%
# mutate(min_interval = as.numeric(min_interval),
#        max_interval = as.numeric(max_interval))


# test1 <- res_df %>%
#   group_by(score_interval) %>%
#   summarise(median = median(sd)) %>%
#   filter(score_interval %in% c(1, i)) %>%
#   pull(median)
# 
# diff_test1 <- test1[2] - test1[1]

# test1_tbl <- test1_tbl %>% bind_rows(tibble(id = i, diff = diff_test1))
# }
# 
# test1_tbl %>%
#   ggplot(aes(id, diff)) +
#     geom_point()
# mutate(predict_clinical = if_else(.pred_pathogenic >= 0.5, 'pathogenic', 'benign')) %>%
# group_by(predict_clinical) %>%
# mutate(rank_score = ntile(.pred_pathogenic, n = 100)) %>%
# ungroup()
# mutate(rank_score =
#          ifelse(predict_clinical == 'pathogenic', 100 - rank_score, rank_score)) %>%
# mutate(quant_residuals = 0 * rank_score + (1  * quant_residuals))



# mutate(combined_score = if_else(.pred_pathogenic >= 0.5,
#                                    rank_score * (sd_residuals),
#                                   0.5 - ((0.5 - rank_score) * (sd_residuals)))) %>%



# res_df_combined %>%
#   select(reliability_score, clinical) %>%
#   arrange(desc(reliability_score)) %>%
#   mutate(is_patho = if_else(clinical == 'pathogenic', 1, 0)) %>%
#   mutate(score_sd = cumsum(is_patho)) %>%
#   mutate(n_row = row_number()) %>%
#   select(n_row, score_sd) %>%
#   left_join(
# 
# res_df_combined %>%
#   select(.pred_pathogenic, clinical) %>%
#   arrange(desc(.pred_pathogenic)) %>%
#   mutate(is_patho = if_else(clinical == 'pathogenic', 1, 0)) %>%
#   mutate(only_score = cumsum(is_patho)) %>%
#   mutate(n_row = row_number()) %>%
#   select(n_row, only_score)
# 
# ) %>%
#   pivot_longer(-n_row, names_to = 'tag', values_to = 'n') %>%
#   mutate(n = n / nrow(res_df_combined[res_df_combined$clinical == 'pathogenic',])) %>%
#   ggplot(aes(n_row, n)) +
#     geom_point(aes(color = tag))
#   # coord_cartesian(xlim = c(0,100))



