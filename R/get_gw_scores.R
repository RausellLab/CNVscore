
##
# From TSV to BED format:
# gunzip whole_genome_SNVs.tsv.gz && cat whole_genome_SNVs.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' | awk -v n=1 '{print $1, $2-n, $3, $7}' OFS='\t' | tail -n +3 > whole_genome_SNVs.bed
# Find intersects: 
# parallel --jobs 22 'bedtools intersect -sorted -wa -wb -a /data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{}.bed  -b /data/non-coding/CADD_v.1.6/whole_genome_SNVs.bed > /data-cbl/frequena_data/cnvscore/cadd_feature/output_chrom_{}.bed' ::: {1..22}

cadd_max <- 1:22 %>% map_dfr(function(x) {
  
  
  tmp_length <- coord_chrom_hg19 %>% filter(chrom == x) %>% pull(length)
  
  tmp_chunks <- tibble('chrom' = x, 'start' = seq(1, tmp_length, by = 100), 'end' = seq(1, tmp_length, by = 100) + 99)
  # We do this because it's out of the genomic boundaries
  tmp_chunks <- tmp_chunks %>% head(-1)
  
  tmp_chunks <- tmp_chunks %>% 
    mutate(start = start - 1) 
  
  return(tmp_chunks)
  })



1:22 %>% map(function(x) {
  
    cadd_max %>%
    filter(chrom == x) %>%
    write_tsv(glue('/data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{x}.bed'), col_names = FALSE)
  
})



tic()
cadd_max <- 1:22 %>% map_dfr(function(x) {
  
  print(x)
  tmp_df <- read_tsv(glue('cadd_feature/output_chrom_{x}.bed'), col_names = FALSE)
  tmp_df2 <- tmp_df %>%
    select(-X4, -X5, -X6) %>%
    group_by(X1, X2, X3) %>%
    summarise(max_cadd = max(X7)) %>%
    ungroup()
  
  return(tmp_df2)
  
})
toc()

cadd_max <- cadd_max %>% rename(chrom = X1, start = X2, end = X3)

write_tsv(cadd_max, 'cadd_feature/result_cadd.tsv')


cadd_max_sure <- cadd_max %>% mutate(start = start + 1)


# GERP
# gunzip /data/non-coding/CADD_v1.3/annotations/gerp/gerp_scores.tsv.gz
# cat gerp_scores.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' | awk -v n=1 '{print $1, $2-n, $3, $5}' OFS='\t' > gerp_scores.bed
# Find intersects: 
# parallel --jobs 22 'bedtools intersect -sorted -wa -wb -a /data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{}.bed  -b /data/non-coding/gerp/gerp_scores.bed > /data-cbl/frequena_data/cnvscore/gerp_feature/output_chrom_{}.bed' ::: {1..22}

gerp_raw <- 1:22 %>% map_dfr(function(x) {
  
  print(x)
  tmp_df <- read_tsv(glue('gerp_feature/output_chrom_{x}.bed'), col_names = FALSE)
  tmp_df2 <- tmp_df %>%
    select(-X4, -X5, -X6) %>%
    group_by(X1, X2, X3) %>%
    summarise(max_gerp = max(X7)) %>%
    ungroup()
  
  return(tmp_df2)
  
})

gerp_max_sure <- gerp_raw %>% 
  rename(chrom = X1, start = X2, end = X3) %>% 
  mutate(chrom = as.character(chrom)) %>%
  mutate(start = start + 1)


##
# From TSV to BED format:
# gunzip whole_genome_SNVs.tsv.gz && cat whole_genome_SNVs.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' | awk -v n=1 '{print $1, $2-n, $3, $7}' OFS='\t' | tail -n +3 > whole_genome_SNVs.bed
# Find intersects: 
# parallel --jobs 22 'bedtools intersect -sorted -wa -wb -a /data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{}.bed  -b /data/non-coding/CADD_v.1.6/whole_genome_SNVs.bed > /data-cbl/frequena_data/cnvscore/cadd_feature/output_chrom_{}.bed' ::: {1..22}




## remot-gw


remot_cnvscore 

