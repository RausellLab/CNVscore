


# ------------------------------------------------------------------------------
# LOAD LIBRARIES
# ------------------------------------------------------------------------------

library(rules)
library(valr)
library(recipes)
library(future)
library(grid)
library(tidyverse)
library(tictoc)
library(ontologyIndex)
library(corrr)
library(furrr)
library(tune)
library(gbm)
library(parsnip)
library(yardstick)
library(glue)
library(patchwork)
library(rsample)
library(xrftest)
library(xrf)
library(rules)
library(dials)
library(workflows)
library(rstanarm)
library(bayestestR)
library(prettydoc)
library(applicable)
library(reticulate)
library(ggpubr)
library(overlapping)
library(corrplot)
library(broom)
library(ranger)
library(rtemis)
library(enrichR)
library(chromPlot)
library(NoiseFiltersR)
library(rtracklayer)
library(liftOver)
library(readxl)
library(ggridges)
library(googlesheets4)
library(gghighlight)
library(see)
library(FactoMineR)
library(factoextra)
library(mgcv)
library(rstatix)
library(MetBrewer)

set.seed(123)


rename <- dplyr::rename
slice <- dplyr::slice
select <- dplyr::select


google_calc_results <- gs4_create(paste0("cnvscore_results_", Sys.Date()), sheets = 'index')


## CREATE FOLDER

# current_date <- '02_04_22'




# save(pnull,
#      expression_features,
#      recomb,
#      lads,
#      gene_density_tbl,
#      genes_promoter,
#      pubmed_df,
#      string_db,
#      para_genes,
#      ensembl_reg,
#      hotspot,
#      hars,
#      eds,
#      hu_map,
#      problematic_regions,
#      region_gaps,
#      ucne,
#      crispr_score_df,
#      haplo_triplo_genes,
#      hu_map,
#      prot_complex,
#      loeuf_score,
#      cadd_max,
#      gerp_max,
#      remot_cnvscore,
#      file = 'local_data_features.RData')

# system('gunzip -c ../cnvxplorer/local_data.RData.gz > local_data.RData')
load('../cnvxplorer/local_data.RData.gz')
# load('bancco_app/local_data_features.RData.gz')
load('local_data_features.RData')

source('R/cnvscore_functions.R')

coord_chrom_hg19 <- read_tsv('https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',
                             col_names = c('chrom', 'length'), col_types = 
                               list(chrom = col_character(),length = col_double())) %>%
  filter(nchar(chrom) < 6) %>% 
  filter(!str_detect(chrom, 'chrM')) %>%
  mutate(chrom = str_remove(chrom, 'chr'))


coord_cytobands <- chromPlot::hg_cytoBandIdeo %>% 
  mutate(Start = Start + 1) %>% 
  rename(chrom = Chrom, start = Start, end = End)


# ------------------------------------------------------------------------------
# FEATURES FORMULE
# ------------------------------------------------------------------------------

features_tbl <- read_tsv('features.tsv')


human_control <- reformulate(termlabels = features_tbl %>% 
                               filter(human_control == 'yes') %>% pull(name), response = 'clinical')

human_no_control <- reformulate(termlabels = features_tbl %>% 
                                  filter(human_control == 'no') %>% pull(name), response = 'clinical')

formule_total <- reformulate(termlabels = features_tbl %>% pull(name), response = 'clinical')

vector_yes_human_control <- features_tbl %>% filter(human_control == 'yes') %>% pull(name)
vector_no_human_control <- features_tbl %>% filter(human_control == 'no') %>% pull(name)
vector_total_human_control <- features_tbl %>% pull(name)

# ------------------------------------------------------------------------------
# DOWNLOAD CNV files
# ------------------------------------------------------------------------------

# download.file('https://www.cell.com/cms/10.1016/j.neuron.2017.06.010/attachment/e1fbd380-d3f4-4c6e-be59-0e25162697f2/mmc2.xlsx',
#               destfile = 'tourette.xlsx')
# 
# download.file('https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/csv/supporting_variants_for_nstd152.csv.gz',
#               destfile = 'supporting_variants_for_nstd152.csv.gz')
# 
# download.file('https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-021-02382-3/MediaObjects/13059_2021_2382_MOESM2_ESM.xlsx',
#               destfile = 'tibet.xlsx')
# 
# download.file('https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/csv/supporting_variants_for_nstd162.csv.gz', 
#               destfile = 'supporting_variants_for_nstd162.csv.gz')
# 
# download.file('https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.bed.gz',
#               'gnomad_v2.1_sv.controls_only.sites.bed.gz')
# 
# download.file('https://decipher.sanger.ac.uk/files/downloads/population_cnv_grch37.txt.gz', destfile = 'population_cnv_grch37.txt.gz')
# 
# download.file('http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt', destfile = 'GRCh37_hg19_variants_2020-02-25.txt')
# 
# 
# download.file('https://github.com/DecodeGenetics/LRS_SV_sets/blob/master/ont_sv_high_confidence_SVs.sorted.vcf.gz?raw=true', destfile = 'ont_sv_high_confidence_SVs.sorted.vcf.gz')
# 
# download.file('https://raw.githubusercontent.com/DecodeGenetics/LRS_SV_sets/master/ont_sv_high_confidence_tandemdup.csv', destfile = 'ont_sv_high_confidence_tandemdup.csv')
# 
# # RUN LINUX
# 'wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2021-10.txt.gz ; gunzip variant_summary_2021-10.txt.gz'




## Icelandic study

# 1 "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"                       
# 2 "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">"             
# 3 "##INFO=<ID=TRRBEGIN,Number=1,Type=Integer,Description=\"Begin position of the Tandem Repeat Region the structural variant is in
# 4 "##INFO=<ID=TRREND,Number=1,Type=Integer,Description=\"End position of the Tandem Repeat Region the structural variant is in
# 5 "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" 

beyter_et_al <- read_tsv('data/ont_sv_high_confidence_SVs.sorted.vcf.gz', skip = 32) %>%
  select(`#CHROM`, POS, INFO) %>%
  separate(INFO, into = c('end', 'variant_class', 'length', 'trr_begin', 'trr_end'), sep = ';') %>%
  mutate(variant_class = str_remove(variant_class, 'SVTYPE=')) %>%
  mutate(end = str_remove(end, 'END=')) %>%
  mutate(trr_begin = str_remove(trr_begin, 'TRRBEGIN=')) %>%
  mutate(trr_end = str_remove(trr_end, 'TRREND=')) %>%
  na_if('.') %>%
  mutate(across(starts_with('trr'), as.numeric)) %>%
  filter(variant_class == 'DEL') %>%
  mutate(variant_class = 'deletion') %>%
  rename(chrom = `#CHROM`, start = POS) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  mutate(end = as.numeric(end)) %>%
  mutate(length_cnv = end - start + 1) %>%
  filter(length_cnv >= 50) %>%
  select(-starts_with('trr'), -length) %>%
  bind_rows(read_tsv('data/ont_sv_high_confidence_tandemdup.csv') %>% 
  select(ref_chrom, ref_begin, ref_end) %>%
  distinct() %>%
  rename(chrom = ref_chrom, start = ref_begin, end = ref_end) %>%
  mutate(length_cnv = end - start + 1) %>%
  mutate(variant_class = 'duplication')
  ) %>%
  mutate(id_liftover = row_number()) %>%
  mutate(chrom = str_remove(chrom, 'chr'))


from_hg38_tohg19 = import.chain('/data-cbl/liftover/hg38ToHg19.over.chain')


granges_total_df <- beyter_et_al %>% 
  select(id_liftover, chrom, start, end) %>% 
  as.data.frame() %>% 
  GRanges()

seqlevelsStyle(granges_total_df) = "UCSC"  # necessary

total_liftover = liftOver(granges_total_df, from_hg38_tohg19)

beyter_after_liftover <- total_liftover %>% 
  as_tibble() %>%
  count(id_liftover) %>%
  filter(n == 1) %>%
  pull(id_liftover)

beyter_et_al <- beyter_et_al %>%
  filter(id_liftover %in% beyter_after_liftover) %>%
  left_join(total_liftover %>% as_tibble() %>% select(id_liftover, width)) %>%
  mutate(is_same = if_else(length_cnv == width, 'yes', 'no')) %>%
  filter(is_same == 'yes') %>%
  mutate(source = 'beyter_et_al') %>%
  mutate(clinical = 'benign') %>%
  mutate(id = 'is_gonna_be_removed') %>%
  select(id, chrom, start, end, variant_class, length_cnv, source, clinical)
  

submission_summary <- read_tsv('data/submission_summary_2021-10.txt', skip = 15)
submission_summary_only <- read_tsv('data/variant_summary_2021-10.txt') %>% 
  mutate(Start = as.numeric(Start)) %>%
  mutate(Stop = as.numeric(Stop)) %>%
  mutate(length_cnv = Stop - Start + 1) %>% 
  filter(length_cnv >= 50) %>% 
  filter(Type %in% c('copy number gain', 'copy number loss', 'Deletion', 'Duplication')) %>% 
  # filter(ClinicalSignificance %in% 'Pathogenic') %>%
  # filter(ClinicalSignificance %in% c('Benign', 'Likely benign', 'Uncertain significance', 'Likely pathogenic', 'Pathogenic')) %>
  filter(Assembly %in% c('GRCh37')) %>% 
  filter(OriginSimple == 'germline') %>% 
  # filter(str_detect(OriginSimple, 'germline') | str_detect(OriginSimple, 'not provided') | str_detect(OriginSimple, 'unknown')) %>%
  filter(str_detect(ReviewStatus, 'submitter') | str_detect(ReviewStatus, 'reviewed by expert panel')) %>% 
  select(VariationID, Chromosome, Start, Stop, Type, ClinicalSignificance, LastEvaluated, ReviewStatus, NumberSubmitters, length_cnv) %>% 
  rename(chrom = Chromosome, start = Start, end = Stop, type_variant = Type, 
         clinical = ClinicalSignificance, score = ReviewStatus) %>% 
  # mutate(id_tmp = paste0(row_number(), '_clinvar')) %>%
  # distinct(chrom, start, end, .keep_all = TRUE) %>%
  filter(chrom != 'Y') %>% 
  mutate(clinical = tolower(clinical)) %>% 
  rename(variant_class = type_variant) %>%
  mutate(variant_class = tolower(variant_class)) %>%
  mutate(variant_class = if_else(str_detect(variant_class, 'copy number gain'), 'duplication', variant_class)) %>%
  mutate(variant_class = if_else(str_detect(variant_class, 'copy number loss'), 'deletion', variant_class)) %>%
  mutate(source = 'clinvar') 

submission_nigreisy <- submission_summary %>% 
  select(Submitter, `#VariationID` ) %>%
  rename(VariationID = `#VariationID`) %>%
  left_join(submission_summary_only %>% 
              select(chrom, start, end, VariationID, variant_class), by = 'VariationID') %>%
  na.omit() %>%
  select(chrom, start, end, Submitter, variant_class) %>%
  rename(submitter = Submitter) %>%
  # there are 216 variants with at least >1 submitter
  # submission_nigreisy %>% count(chrom, start, end, variant_class) %>% filter(n > 1)
  group_by(chrom, start, end, variant_class) %>%
  sample_n(1) %>%
  ungroup() %>%
  distinct()


clinvar_cnvs_hg37 <- read_tsv('data/variant_summary_2021-10.txt') %>% 
  # filter(!str_detect(LastEvaluated, '2021')) %>%
  mutate(Start = as.numeric(Start)) %>%
  mutate(Stop = as.numeric(Stop)) %>%
  mutate(length_cnv = Stop - Start + 1) %>% 
  filter(length_cnv >= 50) %>% 
  filter(Type %in% c('copy number gain', 'copy number loss', 'Deletion', 'Duplication')) %>% 
  # count(ClinicalSignificance) %>% arrange(desc(n))
  # filter(ClinicalSignificance %in% 'Pathogenic') %>%
  # filter(ClinicalSignificance %in% c('Benign', 'Likely benign', 'Uncertain significance', 'Likely pathogenic', 'Pathogenic')) %>
  filter(Assembly %in% c('GRCh37')) %>% 
  filter(OriginSimple == 'germline') %>% 
  # filter(str_detect(OriginSimple, 'germline') | str_detect(OriginSimple, 'not provided') | str_detect(OriginSimple, 'unknown')) %>%
  filter(str_detect(ReviewStatus, 'submitter') | str_detect(ReviewStatus, 'reviewed by expert panel')) %>% 
  select(Chromosome, Start, Stop, Type, ClinicalSignificance, LastEvaluated, ReviewStatus, NumberSubmitters, length_cnv) %>% 
  rename(chrom = Chromosome, start = Start, end = Stop, type_variant = Type, 
         clinical = ClinicalSignificance, score = ReviewStatus) %>% 
  # mutate(id_tmp = paste0(row_number(), '_clinvar')) %>%
  # distinct(chrom, start, end, .keep_all = TRUE) %>%
  filter(chrom != 'Y') %>% 
  mutate(clinical = tolower(clinical)) %>% 
  rename(variant_class = type_variant) %>%
  mutate(variant_class = tolower(variant_class)) %>%
  mutate(variant_class = if_else(str_detect(variant_class, 'copy number gain'), 'duplication', variant_class)) %>%
  mutate(variant_class = if_else(str_detect(variant_class, 'copy number loss'), 'deletion', variant_class)) %>%
  mutate(source = 'clinvar') 


pre_clinvar <- clinvar_cnvs_hg37 %>% count(chrom, start, end, variant_class) %>% filter(n > 1)


problematic_clinvar <- tibble()
for (i in 1:nrow(pre_clinvar)) {
  
  print(i)
  
  c <- clinvar_cnvs_hg37 %>% filter(chrom == pre_clinvar[i,]$chrom,
                               start == pre_clinvar[i,]$start,
                               end == pre_clinvar[i,]$end,
                               variant_class == pre_clinvar[i,]$variant_class)
  
  d <- c %>% filter(clinical == 'pathogenic')
  
  if (nrow(d) == 0) next
  if (nrow(c) == nrow(d)) next
  
  problematic_clinvar <- problematic_clinvar %>% bind_rows(tibble(chrom = pre_clinvar[i,]$chrom,
                                                                  start = pre_clinvar[i,]$start,
                                                                  end = pre_clinvar[i,]$end,
                                                                  variant_class = pre_clinvar[i,]$variant_class))

}

clinvar_cnvs_hg37 <- clinvar_cnvs_hg37 %>% anti_join(problematic_clinvar, by = c('chrom', 'start', 'end'))

clinvar_stacked_plot <- clinvar_cnvs_hg37

clinvar_cnvs_hg37 <- clinvar_cnvs_hg37 %>% filter(clinical == 'pathogenic')

# DECIPHER

tmp_decipher <- read_tsv('/data-cbl/decipher_data/decipher-cnvs-grch37-2020-12-06.txt', skip = 1) %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  select(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  filter(!contribution %in% c('None')) %>%
  # filter(pathogenicity %in% c('Pathogenic',
  #                             'Likely pathogenic',
  #                             'Unknown'
  # ) & (! contribution %in% c('None'))) %>%
  filter(inheritance %in% 'De novo') %>%
  filter(genotype == 'Heterozygous') %>%
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>%
  filter(variant_class %in% c('Deletion', 'Duplication')) %>%
  mutate(variant_class = tolower(variant_class)) %>%
  rename(clinical = pathogenicity)


pre_decipher <- tmp_decipher %>% count(chrom, start, end, variant_class) %>% filter(n > 1)


problematic_decipher <- tibble()

for (i in 1:nrow(pre_decipher)) {
  
  print(i)
  
  c <-  tmp_decipher %>% filter(chrom == pre_decipher[i,]$chrom,
                               start == pre_decipher[i,]$start,
                               end == pre_decipher[i,]$end,
                               variant_class == pre_decipher[i,]$variant_class)
  
  # d <- c %>% filter(str_detect(clinical, 'Pathogenic|pathogenic'))
  d <- c %>% filter(str_detect(clinical, 'Pathogenic|pathogenic|Unknown'))
  
  if (nrow(d) == 0) next
  if (nrow(c) == nrow(d)) next
  
  problematic_decipher <- problematic_decipher %>% bind_rows(tibble(chrom = pre_decipher[i,]$chrom,
                                                                  start = pre_decipher[i,]$start,
                                                                  end = pre_decipher[i,]$end,
                                                                  variant_class = pre_decipher[i,]$variant_class))
}


tmp_decipher <- tmp_decipher %>% anti_join(problematic_decipher, by = c('chrom', 'start', 'end'))

decipher_stacked_plot <- tmp_decipher

tmp_decipher <- tmp_decipher %>% filter(clinical %in% c('Pathogenic', 'Unknown', 'Likely pathogenic'))
# tmp_decipher <- tmp_decipher %>% filter(clinical %in% c('Pathogenic', 'Likely pathogenic'))
tmp_decipher <- tmp_decipher  %>% mutate(clinical = 'pathogenic')

tmp_decipher_only_patho <- read_tsv('/data-cbl/decipher_data/decipher-cnvs-grch37-2020-12-06.txt', skip = 1) %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  select(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>% 
  filter(pathogenicity %in% c('Pathogenic') & (! contribution %in% c('None'))) %>%
  filter(inheritance %in% 'De novo') %>%
  filter(genotype == 'Heterozygous') %>%
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>%
  filter(variant_class %in% c('Deletion', 'Duplication')) %>%
  mutate(variant_class = tolower(variant_class)) %>%
  rename(clinical = pathogenicity) %>%
  mutate(clinical = 'pathogenic')


# decipher Control https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd183/
tmp_decipher_control <- read_tsv('data/population_cnv_grch37.txt.gz', 
                                 col_names = TRUE, col_types = list(chr = col_character())) %>%
  rename(id = `#population_cnv_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>%
  mutate(chrom = as.character(chrom)) %>%
  mutate(chrom = if_else(chrom == 23, 'X', chrom)) %>%
  mutate(source = 'decipher_control') %>%
  filter(frequency <= 0.01) %>%
  filter(deletion_observations == 0 | duplication_observations == 0) %>%
  filter(!(deletion_observations == 0 & duplication_observations == 0)) %>%
  mutate(variant_class = if_else(deletion_observations > 0, 'deletion', 'duplication')) %>%
  mutate(clinical = 'benign')

# Chaisson et al.

chaisson_et_al <- read_csv('data/supporting_variants_for_nstd152.csv.gz') %>%
  rename_with(~ tolower(str_replace_all(., ' ', '_')))

chaisson_et_al <- chaisson_et_al %>%
  filter(assembly_name == 'GRCh37.p13') %>%
  mutate(length_cnv = end - start + 1) %>%
  filter(length_cnv >= 50) %>%
  filter(remap_score == 1) %>%
  rename(chrom = chromosome, variant_class = variant_call_type ) %>%
  select(chrom, start, end, variant_class) %>%
  mutate(chrom = str_remove(chrom, '\\|.*$')) %>%
  filter(chrom != 'Unplaced') %>%
  filter(variant_class %in% c('deletion', 'duplication')) %>%
  na.omit() %>%
  mutate(source = 'chaisson_et_al',
         clinical = 'benign') %>%
  mutate(id_tmp = paste(row_number(), 'chaisson_et_al', sep = '_'))


# Audano

audano_et_al <- read_csv('data/supporting_variants_for_nstd162.csv.gz') %>%
  rename_with(~ tolower(str_replace_all(., ' ', '_')))


audano_et_al <- audano_et_al %>% 
  filter(assembly_name == 'GRCh37.p13') %>% 
  filter(remap_score == 1) %>%
  mutate(length_cnv = end - start + 1) %>% 
  filter(length_cnv >= 50) %>% 
  rename(chrom = chromosome, variant_class = variant_call_type ) %>%
  select(chrom, start, end, variant_class) %>%
  mutate(chrom = str_remove(chrom, '\\|.*$')) %>%
  filter(chrom != 'Unplaced') %>%
  filter(variant_class %in% c('deletion', 'duplication')) %>%
  na.omit() %>%
  mutate(source = 'audano_et_al',
         clinical = 'benign') %>%
  distinct() %>% 
  mutate(id_tmp = paste(row_number(), 'audano_et_al', sep = '_'))


# gnomAD

tmp_gnomad <-  read_tsv('data/gnomad_v2.1_sv.controls_only.sites.bed.gz', col_types = cols(.default = "c")) %>%
  filter(SVTYPE %in% c('DEL', 'DUP')) %>%
  filter(FILTER == 'PASS') %>%
  rename(id = name) %>%
  mutate(id  = str_remove(id, 'gnomAD-SV_v2.1_')) %>%
  mutate(source = 'gnomad_v2.1') %>%
  rename(chrom = `#chrom`) %>%
  mutate(start = as.double(start), end = as.double(end), AF = as.double(AF)) %>%
  mutate(start = start + 1) %>%
  filter(AF <= 0.01) %>%
  select(id, chrom, start, end, svtype, -AF) %>%
  rename(variant_class = svtype) %>%
  mutate(variant_class = if_else(str_detect(variant_class, 'DEL'), 'Deletion', 'Duplication')) %>%
  mutate(source = 'gnomad_v2.1') %>%
  mutate(variant_class = tolower(variant_class)) %>%
  mutate(clinical = 'benign')

# dgv check gold standard
tmp_dgv <- read_tsv('GRCh37_hg19_variants_2020-02-25.txt',
                    col_types = list(chr = col_character())) %>%
  filter(varianttype == 'CNV') %>%
  filter(observedgains == 0 | observedlosses == 0) %>%
  filter(!(observedlosses == 0 & observedgains == 0)) %>%
  mutate(observedtotal = observedgains + observedlosses) %>%
  mutate(frequency = observedtotal / (samplesize * 2)) %>%
  filter(frequency <= 0.01) %>%
  filter(variantsubtype %in% c('deletion', 'duplication', 'loss', 'gain')) %>%
  rename(id = variantaccession, chrom = chr, variant_type = variantsubtype) %>%
  filter(!str_detect(reference, 'gnomAD')) %>%
  mutate(source = 'dgv',
         chrom = as.character(chrom))  %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  mutate(variant_class = if_else(str_detect(variant_type, 'loss'), 'deletion', 'duplication')) %>%
  mutate(clinical = 'benign') %>%
  select(id, chrom, start, end, source, variant_class, clinical) 

input_check_cnv <- tmp_decipher %>%
  bind_rows(tmp_gnomad,
            beyter_et_al,
            tmp_decipher_control,
            tmp_dgv,
            chaisson_et_al,
            audano_et_al,
            clinvar_cnvs_hg37)  %>%
  mutate(length_cnv = end - start  + 1) %>% 
  filter(length_cnv > 49) %>%
  mutate(id_tmp = row_number()) %>%
  # bind_rows(clinvar_cnvs_hg37 %>% select(chrom, source, variant_class)
  # mutate(pathogenicity = if_else(is.na(pathogenicity), 'benign', pathogenicity)) %>%
  # rename(clinical = pathogenicity) %>%
  filter(chrom != 'Y') %>%
  select(id, chrom, start, end, variant_class, clinical, source, length_cnv, id_tmp, LastEvaluated, NumberSubmitters)


preproc_flow <- tibble()
preproc_flow <- preproc_flow %>% bind_rows(fill_flow(input_check_cnv, 'start'))


# 0. Overlap with problematic regions------------------------------------------------------------------------------

remove_ids_overlap <- input_check_cnv %>% 
  bed_coverage(problematic_regions) %>% 
  select(id_tmp, .frac) %>% filter(.frac >= 0.3) %>% pull(id_tmp)

input_check_cnv <- input_check_cnv %>% filter(!id_tmp %in% remove_ids_overlap)

preproc_flow <- preproc_flow %>% bind_rows(fill_flow(input_check_cnv, 'problematic_regions'))

# 1. LiftOver step------------------------------------------------------------------------------

split_ids_decipher <- report_split_cnvs(input_check_cnv)

input_check_cnv <- input_check_cnv %>% filter(id_tmp %in% split_ids_decipher)

preproc_flow <- preproc_flow %>% bind_rows(fill_flow(input_check_cnv, 'liftover'))


input_check_cnv_del <- input_check_cnv %>% filter(variant_class == 'deletion')
input_check_cnv_dup <- input_check_cnv %>% filter(variant_class == 'duplication')


# 2. Remove identical CNVs------------------------------------------------------------------------------

input_check_cnv_del <- input_check_cnv_del %>%
  distinct(chrom, start, end, clinical, .keep_all = TRUE)

input_check_cnv_dup <- input_check_cnv_dup %>%
  distinct(chrom, start, end, clinical, .keep_all = TRUE)

preproc_flow <- preproc_flow %>% bind_rows(fill_flow(bind_rows(input_check_cnv_del, input_check_cnv_dup), 'identical_cnvs'))

# 3. Reciprocal overlap------------------------------------------------------------------------------

input_check_cnv_del_benign <- input_check_cnv_del %>% filter(clinical == 'benign')
input_check_cnv_dup_benign <- input_check_cnv_dup %>% filter(clinical == 'benign')
input_check_cnv_del_pathogenic_clinvar <- input_check_cnv_del %>% filter(clinical == 'pathogenic' & source == 'clinvar')
input_check_cnv_dup_pathogenic_clinvar <- input_check_cnv_dup %>% filter(clinical == 'pathogenic' & source == 'clinvar')
input_check_cnv_del_pathogenic_decipher <- input_check_cnv_del %>% filter(clinical == 'pathogenic' & source == 'decipher')
input_check_cnv_dup_pathogenic_decipher <- input_check_cnv_dup %>% filter(clinical == 'pathogenic' & source == 'decipher')


remove_input_check_cnv_del_benign <- reciprocal_overlap(input_check_cnv_del_benign)
remove_input_check_cnv_dup_benign <- reciprocal_overlap(input_check_cnv_dup_benign)
remove_input_check_cnv_del_pathogenic_clinvar <- reciprocal_overlap(input_check_cnv_del_pathogenic_clinvar)
remove_input_check_cnv_dup_pathogenic_clinvar <- reciprocal_overlap(input_check_cnv_dup_pathogenic_clinvar)
remove_input_check_cnv_del_pathogenic_decipher <- reciprocal_overlap(input_check_cnv_del_pathogenic_decipher)
remove_input_check_cnv_dup_pathogenic_decipher <- reciprocal_overlap(input_check_cnv_dup_pathogenic_decipher)


input_check_cnv_del_benign <- input_check_cnv_del_benign %>% filter(! id_tmp %in% remove_input_check_cnv_del_benign) 
input_check_cnv_dup_benign <- input_check_cnv_dup_benign %>% filter(! id_tmp %in% remove_input_check_cnv_dup_benign) 
input_check_cnv_del_pathogenic_clinvar <- input_check_cnv_del_pathogenic_clinvar %>% filter(! id_tmp %in% remove_input_check_cnv_del_pathogenic_clinvar) 
input_check_cnv_dup_pathogenic_clinvar <- input_check_cnv_dup_pathogenic_clinvar %>% filter(! id_tmp %in% remove_input_check_cnv_dup_pathogenic_clinvar) 
input_check_cnv_del_pathogenic_decipher <- input_check_cnv_del_pathogenic_decipher %>% filter(! id_tmp %in% remove_input_check_cnv_del_pathogenic_decipher) 
input_check_cnv_dup_pathogenic_decipher <- input_check_cnv_dup_pathogenic_decipher %>% filter(! id_tmp %in% remove_input_check_cnv_dup_pathogenic_decipher) 


preproc_flow <- preproc_flow %>% 
  bind_rows(fill_flow(bind_rows(input_check_cnv_del_benign, 
                                input_check_cnv_del_pathogenic_clinvar,
                                input_check_cnv_dup_pathogenic_clinvar,
                                input_check_cnv_del_pathogenic_decipher,
                                input_check_cnv_dup_pathogenic_decipher,
                                input_check_cnv_dup_benign), 'reciprocal_overlap'))


# Check pathogenic - benign overlap ------------------------------------------------------------------------------

second_overlap_patho_del_clinvar <- overlap_benign_pathogenic(input_check_cnv_del_pathogenic_clinvar, input_check_cnv_del_benign)
second_overlap_patho_del_decipher <- overlap_benign_pathogenic(input_check_cnv_del_pathogenic_decipher, input_check_cnv_del_benign)
second_overlap_patho_dup_clinvar <- overlap_benign_pathogenic(input_check_cnv_dup_pathogenic_clinvar, input_check_cnv_dup_benign)
second_overlap_patho_dup_decipher <- overlap_benign_pathogenic(input_check_cnv_dup_pathogenic_decipher, input_check_cnv_dup_benign)

input_check_cnv_del_pathogenic_clinvar <- input_check_cnv_del_pathogenic_clinvar %>% 
  filter(!id_tmp %in% second_overlap_patho_del_clinvar$id_tmp_patho)

input_check_cnv_del_pathogenic_decipher <- input_check_cnv_del_pathogenic_decipher %>% 
  filter(!id_tmp %in% second_overlap_patho_del_decipher$id_tmp_patho)

input_check_cnv_dup_pathogenic_clinvar <- input_check_cnv_dup_pathogenic_clinvar %>% 
  filter(!id_tmp %in% second_overlap_patho_dup_clinvar$id_tmp_patho)

input_check_cnv_dup_pathogenic_decipher <- input_check_cnv_dup_pathogenic_decipher %>% 
  filter(!id_tmp %in% second_overlap_patho_dup_decipher$id_tmp_patho)

input_check_cnv_del_benign <- input_check_cnv_del_benign %>% 
  filter(!id_tmp %in% c(second_overlap_patho_del_clinvar$id_tmp_benign, second_overlap_patho_del_decipher$id_tmp_benign ))

input_check_cnv_dup_benign <- input_check_cnv_dup_benign %>% 
  filter(!id_tmp %in% c(second_overlap_patho_dup_clinvar$id_tmp_benign, second_overlap_patho_dup_decipher$id_tmp_benign ))

preproc_flow <- preproc_flow %>% 
  bind_rows(fill_flow(bind_rows(input_check_cnv_del_benign, 
                                input_check_cnv_del_pathogenic_clinvar,
                                input_check_cnv_dup_pathogenic_clinvar,
                                input_check_cnv_del_pathogenic_decipher,
                                input_check_cnv_dup_pathogenic_decipher,
                                input_check_cnv_dup_benign), 'conflict_pathogenic_benign'))

                  
# Split clinvar (before - after Jan 2021)------------------------------------------------------------------------------

input_check_cnv_del_pathogenic_clinvar_after_20 <- input_check_cnv_del_pathogenic_clinvar %>% 
  filter(str_detect(LastEvaluated, '2021') & NumberSubmitters == 1) %>% select(-LastEvaluated, -NumberSubmitters)

input_check_cnv_del_pathogenic_clinvar <- input_check_cnv_del_pathogenic_clinvar %>% 
  filter(!str_detect(LastEvaluated, '2021')) %>% select(-LastEvaluated, -NumberSubmitters)

preproc_flow <- preproc_flow %>% 
  bind_rows(fill_flow(bind_rows(input_check_cnv_del_benign, 
                                input_check_cnv_del_pathogenic_clinvar,
                                input_check_cnv_dup_pathogenic_clinvar,
                                input_check_cnv_del_pathogenic_decipher,
                                input_check_cnv_dup_pathogenic_decipher,
                                input_check_cnv_dup_benign), 'remove_before_2021'))


# Figure 2 - length_distribution------------------------------------------------------------------------------


density_length_df <- bind_rows(input_check_cnv_del_benign, 
                               input_check_cnv_del_pathogenic_clinvar,
                               input_check_cnv_dup_pathogenic_clinvar,
                               input_check_cnv_del_pathogenic_decipher,
                               input_check_cnv_dup_pathogenic_decipher,
                               input_check_cnv_dup_benign) %>%
  select(source, length_cnv, clinical, variant_class)

density_length_df <- density_length_df %>%
  mutate(source = case_when(
    source == 'gnomad_v2.1' ~ 'gnomAD-SV',
    source == 'decipher' ~ 'DECIPHER',
    source == 'decipher_control' ~ 'DECIPHER Control Set',
    source == 'dgv' ~ 'DGV',
    source == 'clinvar' ~ 'ClinVar',
    source == 'beyter_et_al' ~ 'Beyter et al., 2021',
    source == 'chaisson_et_al' ~ 'Chaisson et al., 2019',
    source == 'audano_et_al' ~ 'Audano et al., 2019',
    
  )) %>%
  # mutate(tag3 = case_when(
  #   str_detect(source, 'ClinVar') ~ 'Pathogenic CNVs',
  #   str_detect(source, 'DECIPHER Control Set') ~ 'Non-pathogenic CNVs',
  #   str_detect(source, 'DECIPHER') ~ 'Pathogenic CNVs',
  #   str_detect(source, 'DGV') ~ 'Non-pathogenic CNVs',
  #   str_detect(source, 'gnomAD-SV') ~ 'Non-pathogenic CNVs',
  #   str_detect(source, 'Beyter et al., 2021') ~ 'Non-pathogenic CNVs',
  #   str_detect(source, 'Chaisson et al., 2019') ~ 'Non-pathogenic CNVs',
  #   str_detect(source, 'Audano et al., 2019') ~ 'Non-pathogenic CNVs'
  # )) %>%
mutate(tag5 = case_when(
  str_detect(source, 'ClinVar') ~ 2,
  str_detect(source, 'DECIPHER Control Set') ~ 3,
  str_detect(source, 'DECIPHER') ~ 1,
  str_detect(source, 'DGV') ~ 4,
  str_detect(source, 'gnomAD-SV') ~ 5,
  str_detect(source, 'Beyter et al., 2021') ~ 6,
  str_detect(source, 'Chaisson et al., 2019') ~ 7,
  str_detect(source, 'Audano et al., 2019') ~ 8
))


median_length <- density_length_df  %>%
  group_by(source) %>%
  summarise(median_length_cnv = median(length_cnv))

total_n_source <- density_length_df %>% count(source)

figure_length_distribution_before_matching <- density_length_df %>% 
  left_join(median_length, by = 'source') %>%
  left_join(total_n_source, by = 'source') %>%
  mutate(source = paste0(source, ' (', n , ')')) %>%
  mutate(variant_class = ifelse(variant_class == 'deletion', 'Deletion', 'Duplication')) %>%
  ggplot(aes(x = log10(length_cnv), y = reorder(source, -tag5))) + 
  stat_density_ridges(aes(fill = clinical), alpha = 0.8, quantiles = 2, quantile_lines = TRUE, show.legend = FALSE) +
  labs(x = 'Log10(CNV length)', y = 'Dataset') +

  scale_fill_manual(values = c('#4b6319', '#dc3545')) +
  facet_wrap(vars(variant_class)) +
  theme_minimal() +
  theme( axis.text=element_text(size=12,face="bold"),
    axis.title=element_text(size=14,face="bold"))



figure_length_distribution_before_matching_variant_class <- density_length_df %>% 
  left_join(median_length, by = 'source') %>%
  left_join(total_n_source, by = 'source') %>%
  mutate(source = paste0(source, ' (', n , ')')) %>%
  ggplot(aes(x = log10(length_cnv), y = reorder(source, median_length_cnv))) + 
  # geom_density_ridges() + 
  stat_density_ridges(aes(fill = source), alpha = 0.4, quantiles = 2, quantile_lines = TRUE, show.legend = FALSE) +
  # facet_wrap(vars(variant_class)) +
  labs(title = 'Length distribution across sources and CNV types', x = 'Log10(CNV length)', y = 'Source') +
  facet_wrap(vars(variant_class)) +
  theme_minimal()



# Figure 3 - Stacked barplot------------------------------------------------------------------------------

annotation_stacked <- bind_rows(input_check_cnv_del_benign, 
                                input_check_cnv_dup_benign,
          input_check_cnv_del_pathogenic_clinvar,
          input_check_cnv_dup_pathogenic_clinvar,
          input_check_cnv_del_pathogenic_decipher,
          input_check_cnv_dup_pathogenic_decipher) %>% 
  left_join(decipher_stacked_plot %>% rename(clinical2 = clinical), 
            by = c("id", "chrom", "start", "end", "variant_class", "source")) %>%
  mutate(clinical = case_when(
    clinical == 'benign' ~ 'benign',
    !is.na(clinical2) ~ clinical2,
    TRUE ~ 'Pathogenic'
  )) %>%
  select(-id) %>%
  check_cnv_v2(factor_clinical = FALSE)

tmp_only_omim_outside <- annotation_stacked %>%
  filter(omim == 0) %>%
  select(id, chrom, start, end) %>%
  bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
  filter(.overlap > 0) %>%
  pull(id.x) %>%
  unique()


df_genes_inside <- annotation_stacked %>%
  filter(omim == 1) %>%
  bed_intersect(hgcn_genes) %>%
  select(id.x, gene.y) %>%
  rename(id = id.x, gene_inside = gene.y)

df_genes_outside <- annotation_stacked %>%
  filter(omim == 1) %>%
  bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
  select(id.x, gene.y) %>%
  rename(id = id.x, gene_outside = gene.y)

unique_id_genes <- unique(df_genes_inside$id)

result_ids <- c()

for (i in 1:length(unique_id_genes)) {
  
  how_many_genes_outside_no_inside <- df_genes_outside %>% 
    filter(id == unique_id_genes[i] ) %>%
    filter(!gene_outside %in% df_genes_inside[df_genes_inside$id == unique_id_genes[i],]$gene_inside) %>%
    nrow()
  
  if (how_many_genes_outside_no_inside > 0)  result_ids <- c(unique_id_genes[i],result_ids)
  
  
}



tmp_inside_outside_omim <- annotation_stacked %>%
  filter(omim == 1) %>%
  select(id, chrom, start, end) %>%
  bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
  filter(.overlap > 0) %>%
  pull(id.x) %>%
  unique()


remove_small_clinicals <- annotation_stacked %>% 
  mutate(tag2 = paste(source, '-', clinical)) %>% 
  count(tag2) %>% 
  arrange(desc(n)) %>%
  filter(n < 22) %>%
  pull(tag2)

show_n <- annotation_stacked %>% 
  mutate(tag2 = paste(source, '-', clinical)) %>% 
  count(tag2)


figure_stacked_barplot <- annotation_stacked %>% 
  mutate(target_omim_only_outside = if_else(id %in% tmp_only_omim_outside, 'yes', 'no')) %>%
  mutate(tag = case_when(
    id %in% result_ids ~ 'OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV',
    omim == 0 & target_omim_only_outside == 'yes' ~ 'Enhancers targeted of OMIM genes outside the CNV',
    n_genes > 0 & omim == 0 ~ 'Other protein-coding genes targeted',
    n_genes == 0 ~ 'No protein-coding genes targeted',
    omim > 0 ~ 'OMIM genes targeted'
  )) %>% 
  # mutate(source = if_else(source == 'clinvar', 'ClinVar', source)) %>%
  # mutate(source = if_else(source == 'decipher', 'DECIPHER', source)) %>%
  mutate(tag2 = paste(source, '-', clinical)) %>%  
  filter(!tag2 %in% remove_small_clinicals) %>%
  # left_join(show_n, by = 'tag2') %>%
  # mutate(tag2 = glue('{tag2} ({n})')) %>%
  count(tag, tag2, source) %>% 
  group_by(tag2, source) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(tag = factor(tag, c('No protein-coding genes targeted', 'Other protein-coding genes targeted', 'Enhancers targeted of OMIM genes outside the CNV','OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV', 'OMIM genes targeted'))) %>%
  mutate(tag3 = case_when(
    str_detect(source, 'decipher_control') ~ 'General population',
    str_detect(source, 'decipher') ~ 'DECIPHER',
    str_detect(source,'clinvar') ~ 'ClinVar',
    TRUE ~ 'General population'
  )) %>%
  mutate(tag2 = str_remove(tag2, ' - benign')) %>%
  mutate(tag2 = str_to_title(tag2)) %>%   
  # filter((tag3 == 'ClinVar' & str_detect(tag2, 'Pathogenic')) | 
  #          (tag3 == 'DECIPHER' & str_detect(tag2, 'Pathogenic |Unknown')) |
  #          (tag3 == 'General population')) %>%
  mutate(tag4 = if_else(tag3 == 'General population', 'Benign CNVs', 'Pathogenic CNVs')) %>%
  mutate(tag4 = factor(tag4, levels = c('Pathogenic CNVs', 'Benign CNVs'))) %>% 
  ungroup() %>%
  # mutate(tag2 = factor(tag2)) %>% 
  mutate(tag2 = str_replace(tag2, 'Decipher', 'DECIPHER')) %>%
  mutate(tag2 = str_replace(tag2, 'DECIPHER_control', 'DECIPHER Control Set')) %>%
  mutate(tag2 = str_replace(tag2, 'Dgv', 'DGV')) %>%
  mutate(tag2 = str_replace(tag2, 'Gnomad_v2.1', 'gnomAD-SV')) %>%
  mutate(tag2 = str_replace(tag2, 'Beyter_et_al', 'Beyter et al., 2021')) %>%
  mutate(tag2 = str_replace(tag2, 'Audano_et_al', 'Audano et al., 2019')) %>%
  mutate(tag2 = str_replace(tag2, 'Chaisson_et_al', 'Chaisson et al., 2019')) %>% 
  mutate(tag2 = factor(tag2, levels = c('Clinvar - Pathogenic', 'DECIPHER - Pathogenic', 'DECIPHER - Likely Pathogenic', 'DECIPHER - Unknown',
                                        'Audano et al., 2019', 'Beyter et al., 2021', 'Chaisson et al., 2019', 'DECIPHER Control Set',
                                        'DGV', 'gnomAD-SV'))) %>%
  # mutate(tag2 = fct_inseq(tag2, tag5)) %>%
  # mutate(tag2 = str_remove(tag2, '\\([^)]*\\)')) %>%
  ggplot(aes(tag2, perc, group = tag)) +
  geom_col(aes(fill = tag), color = 'black') +
  scale_y_continuous(label = percent) +
  # geom_label(aes(label = paste0(100*round(perc, 2), '%' , ' (', n, ')')), size = 3, position = position_stack(vjust = 0.5)) +
  # scale_fill_manual(values = c(hue_pal()(4)[2], '#Ecd26d', '#Eca9a9', hue_pal()(1), '#DD2E44')) +
  facet_wrap(vars(tag4), scale = 'free_x') +
  theme_minimal() +
  labs(fill = 'Category', y = 'Percentage', x = '') +
  scale_fill_met_d('Archambault') +
  theme_minimal() +
  theme(axis.text=element_text(size= 8,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = -45, hjust = 0))



figure_stacked_per_type_variant <- annotation_stacked %>% 
  mutate(target_omim_only_outside = if_else(id %in% tmp_only_omim_outside, 'yes', 'no')) %>%
  mutate(tag = case_when(
    id %in% result_ids ~ 'OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV',
    omim == 0 & target_omim_only_outside == 'yes' ~ 'Enhancers targeted of OMIM genes outside the CNV',
    n_genes > 0 & omim == 0 ~ 'Other protein-coding genes targeted',
    n_genes == 0 ~ 'No protein-coding genes targeted',
    omim > 0 ~ 'OMIM genes targeted'
  )) %>% 
  # mutate(source = if_else(source == 'clinvar', 'ClinVar', source)) %>%
  # mutate(source = if_else(source == 'decipher', 'DECIPHER', source)) %>%
  mutate(tag2 = paste(source, '-', clinical)) %>%  
  filter(!tag2 %in% remove_small_clinicals) %>%
  # left_join(show_n, by = 'tag2') %>%
  # mutate(tag2 = glue('{tag2} ({n})')) %>%
  count(tag, tag2, source, type_variant) %>% 
  group_by(tag2, source, type_variant) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(tag = factor(tag, c('No protein-coding genes targeted', 'Other protein-coding genes targeted', 'Enhancers targeted of OMIM genes outside the CNV','OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV', 'OMIM genes targeted'))) %>%
  mutate(tag3 = case_when(
    str_detect(source, 'decipher_control') ~ 'General population',
    str_detect(source, 'decipher') ~ 'DECIPHER',
    str_detect(source,'clinvar') ~ 'ClinVar',
    TRUE ~ 'General population'
  )) %>%
  mutate(tag2 = str_remove(tag2, ' - benign')) %>%
  mutate(tag2 = str_to_title(tag2)) %>%   
  # filter((tag3 == 'ClinVar' & str_detect(tag2, 'Pathogenic')) | 
  #          (tag3 == 'DECIPHER' & str_detect(tag2, 'Pathogenic |Unknown')) |
  #          (tag3 == 'General population')) %>%
  mutate(tag4 = if_else(tag3 == 'General population', 'Benign CNVs', 'Pathogenic CNVs')) %>%
  mutate(tag4 = factor(tag4, levels = c('Pathogenic CNVs', 'Benign CNVs'))) %>% 
  ungroup() %>%
  # mutate(tag2 = factor(tag2)) %>% 
  mutate(tag2 = str_replace(tag2, 'Decipher', 'DECIPHER')) %>%
  mutate(tag2 = str_replace(tag2, 'DECIPHER_control', 'DECIPHER Control Set')) %>%
  mutate(tag2 = str_replace(tag2, 'Dgv', 'DGV')) %>%
  mutate(tag2 = str_replace(tag2, 'Gnomad_v2.1', 'gnomAD-SV')) %>%
  mutate(tag2 = str_replace(tag2, 'Beyter_et_al', 'Beyter et al., 2021')) %>%
  mutate(tag2 = str_replace(tag2, 'Audano_et_al', 'Audano et al., 2019')) %>%
  mutate(tag2 = str_replace(tag2, 'Chaisson_et_al', 'Chaisson et al., 2019')) %>% 
  mutate(tag2 = factor(tag2, levels = c('Clinvar - Pathogenic', 'DECIPHER - Pathogenic', 'DECIPHER - Likely Pathogenic', 'DECIPHER - Unknown',
                                        'Audano et al., 2019', 'Beyter et al., 2021', 'Chaisson et al., 2019', 'DECIPHER Control Set',
                                        'DGV', 'gnomAD-SV'))) %>%
  # mutate(tag2 = fct_inseq(tag2, tag5)) %>%
  # mutate(tag2 = str_remove(tag2, '\\([^)]*\\)')) %>%
  mutate(type_variant = str_to_title(type_variant)) %>%
  mutate(tag4 = paste0(type_variant, ' - ', tag4)) %>%
  ggplot(aes(tag2, perc, group = tag)) +
  geom_col(aes(fill = tag), color = 'black') +
  scale_y_continuous(label = percent) +
  facet_wrap(vars(tag4), scale = 'free_x') +
  theme_minimal() +
  labs(fill = 'Category', y = 'Percentage', x = '') +
  scale_fill_met_d('Archambault') +
  facet_wrap(vars(tag4), scales = 'free') +
  theme_minimal() +
  theme(axis.text=element_text(size= 8,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_text(angle = -45, hjust = 0))


  

# Matching by length------------------------------------------------------------------------------

# Deletion


clinvar_match_deletion <- matching_length(bin_length = 100, input_check_cnv_del_pathogenic_clinvar %>% 
                                            bind_rows(input_check_cnv_del_benign %>% select(-LastEvaluated, -NumberSubmitters))
)

remove_from_clinvar_to_decipher_and_ind_del <- clinvar_match_deletion %>% filter(clinical == 'benign') %>% pull(id_tmp)

decipher_match_deletion <- matching_length(bin_length = 100, input_check_cnv_del_pathogenic_decipher %>% 
                                            bind_rows(input_check_cnv_del_benign %>% select(-LastEvaluated, -NumberSubmitters) %>%
                                                        filter(!id_tmp %in% remove_from_clinvar_to_decipher_and_ind_del))
)

clinvar_20_match <-  matching_length(bin_length = 100, input_check_cnv_del_pathogenic_clinvar_after_20 %>% 
                                       bind_rows(input_check_cnv_del_benign %>% select(-LastEvaluated, -NumberSubmitters) %>%
                                                   filter(!id_tmp %in% remove_from_clinvar_to_decipher_and_ind_del))
)

# Duplications

clinvar_match_duplication <- matching_length(bin_length = 100, input_check_cnv_dup_pathogenic_clinvar %>% 
                                            bind_rows(input_check_cnv_dup_benign %>% select(-LastEvaluated, -NumberSubmitters))
)

remove_from_clinvar_to_decipher_and_ind_dup <- clinvar_match_duplication %>% filter(clinical == 'benign') %>% pull(id_tmp)

decipher_match_duplication <- matching_length(bin_length = 100, input_check_cnv_dup_pathogenic_decipher %>% 
                                             bind_rows(input_check_cnv_dup_benign %>% select(-LastEvaluated, -NumberSubmitters) %>%
                                                         filter(!id_tmp %in% remove_from_clinvar_to_decipher_and_ind_dup))
)


preproc_flow <- preproc_flow %>% 
  bind_rows(fill_flow(bind_rows(clinvar_match_deletion, 
                                clinvar_match_duplication,
                                decipher_match_deletion,
                                decipher_match_duplication), 'match_by_length'))




# Figure 3 after matching------------------------------------------------------------------------------


annotation_stacked2 <- bind_rows(clinvar_match_deletion, 
                                decipher_match_deletion) %>% 
  left_join(decipher_stacked_plot %>% rename(clinical2 = clinical), 
            by = c("id", "chrom", "start", "end", "variant_class", "source")) %>%
  mutate(clinical = case_when(
    clinical == 'benign' ~ 'benign',
    !is.na(clinical2) ~ clinical2,
    TRUE ~ 'Pathogenic'
  )) %>% 
  check_cnv_v2(factor_clinical = FALSE)

tmp_only_omim_outside <- annotation_stacked2 %>%
  filter(omim == 0) %>%
  select(id, chrom, start, end) %>%
  bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
  filter(.overlap > 0) %>%
  pull(id.x) %>%
  unique()


df_genes_inside <- annotation_stacked2 %>%
  filter(omim == 1) %>%
  bed_intersect(hgcn_genes) %>%
  select(id.x, gene.y) %>%
  rename(id = id.x, gene_inside = gene.y)

df_genes_outside <- annotation_stacked2 %>%
  filter(omim == 1) %>%
  bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
  select(id.x, gene.y) %>%
  rename(id = id.x, gene_outside = gene.y)

unique_id_genes <- unique(df_genes_inside$id)

result_ids <- c()

for (i in 1:length(unique_id_genes)) {
  
  how_many_genes_outside_no_inside <- df_genes_outside %>% 
    filter(id == unique_id_genes[i] ) %>%
    filter(!gene_outside %in% df_genes_inside[df_genes_inside$id == unique_id_genes[i],]$gene_inside) %>%
    nrow()
  
  if (how_many_genes_outside_no_inside > 0)  result_ids <- c(unique_id_genes[i],result_ids)
  
  
}



tmp_inside_outside_omim <- annotation_stacked2 %>%
  filter(omim == 1) %>%
  select(id, chrom, start, end) %>%
  bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
  filter(.overlap > 0) %>%
  pull(id.x) %>%
  unique()


remove_small_clinicals <- annotation_stacked2 %>% 
  mutate(tag2 = paste(source, '-', clinical)) %>% 
  count(tag2) %>% 
  arrange(desc(n)) %>%
  filter(n < 22) %>%
  pull(tag2)

show_n <- annotation_stacked2 %>% 
  mutate(tag2 = paste(source, '-', clinical)) %>% 
  count(tag2)



figure_stacked_barplot2 <- annotation_stacked2 %>% 
  mutate(target_omim_only_outside = if_else(id %in% tmp_only_omim_outside, 'yes', 'no')) %>%
  mutate(tag = case_when(
    id %in% result_ids ~ 'OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV',
    omim == 0 & target_omim_only_outside == 'yes' ~ 'Enhancers targeted of OMIM genes outside the CNV',
    n_genes > 0 & omim == 0 ~ 'Other protein-coding genes targeted',
    n_genes == 0 ~ 'No protein-coding genes targeted',
    omim > 0 ~ 'OMIM genes targeted'
  )) %>% 
  # mutate(source = if_else(source == 'clinvar', 'ClinVar', source)) %>%
  # mutate(source = if_else(source == 'decipher', 'DECIPHER', source)) %>%
  mutate(tag2 = paste(source, '-', clinical)) %>%
  filter(!tag2 %in% remove_small_clinicals) %>%
  left_join(show_n, by = 'tag2') %>%
  mutate(tag2 = glue('{tag2} ({n})')) %>%
  count(tag, tag2, source) %>%
  group_by(tag2, source) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(tag = factor(tag, c('No protein-coding genes targeted', 'Other protein-coding genes targeted', 'Enhancers targeted of OMIM genes outside the CNV','OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV', 'OMIM genes targeted'))) %>%
  mutate(tag3 = case_when(
    str_detect(source, 'decipher_control') ~ 'General population',
    str_detect(source, 'decipher') ~ 'DECIPHER',
    str_detect(source,'clinvar') ~ 'ClinVar',
    TRUE ~ 'General population'
  )) %>%
  mutate(tag2 = str_remove(tag2, ' - benign')) %>%
  mutate(tag2 = str_to_title(tag2)) %>%
  filter((tag3 == 'ClinVar' & str_detect(tag2, 'Pathogenic')) | 
           (tag3 == 'DECIPHER' & str_detect(tag2, 'Pathogenic |Unknown')) |
           (tag3 == 'General population')) %>%
  mutate(tag4 = if_else(tag3 == 'General population', 'Benign CNVs', 'Pathogenic CNVs')) %>%
  mutate(tag4 = factor(tag4, levels = c('Pathogenic CNVs', 'Benign CNVs'))) %>%
  mutate(tag5 = case_when(
    str_detect(tag2, 'Clinvar - Pathogenic') ~ 1,
    # str_detect(tag2, 'Clinvar - Likely Pathogenic') ~ 3,
    str_detect(tag2, 'Decipher - Pathogenic') ~ 2,
    str_detect(tag2, 'Decipher - Likely Pathogenic') ~ 3,
    str_detect(tag2, 'Decipher - Unknown') ~ 4,
    str_detect(tag2, 'Decipher_control') ~ 5,
    str_detect(tag2, 'Dgv') ~ 6,
    str_detect(tag2, 'Gnomad_v2.1') ~ 10,
    str_detect(tag2, 'Chaisson') ~ 9,
    str_detect(tag2, 'Audano') ~ 8,
    str_detect(tag2, 'Beyter') ~ 7
  )) %>%
  ungroup() %>%
  mutate(tag2 = factor(tag2)) %>%
  mutate(tag2 = str_replace(tag2, 'Decipher', 'DECIPHER')) %>%
  mutate(tag2 = str_replace(tag2, 'Decipher_control', 'DECIPHER Control Set')) %>%
  mutate(tag2 = str_replace(tag2, 'Dgv', 'DGV')) %>%
  mutate(tag2 = str_replace(tag2, 'Gnomad_v2.1', 'gnomAD-SV')) %>%
  mutate(tag2 = str_replace(tag2, 'Beyter_et_al', 'Beyter et al., 2021')) %>%
  mutate(tag2 = str_replace(tag2, 'Audano_et_al', 'Audano et al., 2019')) %>%
  mutate(tag2 = str_replace(tag2, 'Chaisson_et_al', 'Chaisson et al., 2019')) %>%
  mutate(tag2 = fct_reorder(tag2, tag5)) %>%
  ggplot(aes(tag2, perc, group = tag)) +
  geom_col(aes(fill = tag), color = 'black') +
  scale_y_continuous(label = percent) +
  geom_label(aes(label = paste0(100*round(perc, 2), '%' , ' (', n, ')')), size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c(hue_pal()(4)[2], '#Ecd26d', '#Eca9a9', hue_pal()(1), '#DD2E44')) +
  facet_wrap(vars(tag4), scale = 'free_x') +
  theme_minimal() +
  labs(fill = 'Category', y = 'Percentage') +
  theme(axis.text.x=element_text(angle = -45, hjust = 0))




# N. Length distribution after matching------------------------------------------------------------------------------


density_length_df_after <- bind_rows(clinvar_match_deletion, 
                               decipher_match_deletion,
                               clinvar_match_duplication,
                               decipher_match_duplication) %>%
  select(source, length_cnv, clinical)


median_length <- density_length_df_after  %>%
  group_by(source) %>%
  summarise(median_length_cnv = median(length_cnv))

total_n_source <- density_length_df_after %>% count(source)


figure_length_distribution_after_matching <- density_length_df_after %>% 
  left_join(median_length, by = 'source') %>%
  left_join(total_n_source, by = 'source') %>%
  mutate(source = paste0(source, ' (', n , ')')) %>%
  mutate(clinical = fct_rev(clinical)) %>%
  ggplot(aes(x = log10(length_cnv), y = reorder(source, median_length_cnv))) + 
  # geom_density_ridges() + 
  stat_density_ridges(aes(fill = clinical), alpha = 0.4, quantiles = 2, quantile_lines = TRUE, show.legend = TRUE) +
  # facet_wrap(vars(variant_class)) +
  labs(title = 'Length distribution across sources and CNV types', x = 'Log10(CNV length)', y = 'Source') +
  # scale_fill_branded() +
  theme_minimal() +
  theme(legend.position="top")


figure_length_distribution_after_matching2 <- density_length_df_after %>% 
  mutate(clinical = str_to_title(clinical)) %>%
  ggplot(aes(log10(length_cnv), clinical)) + 
  # geom_density(aes(fill = clinical)) +
  stat_density_ridges(aes(fill = clinical), alpha = 0.4, quantiles = 2, quantile_lines = TRUE, show.legend = FALSE) +
  scale_fill_manual(values = c('#4b6319', '#dc3545')) +
  labs(x = 'Log10(CNV length)', y = '') +
  theme_minimal() +
  theme(legend.position="right")





# N. Annotation------------------------------------------------------------------------------


output_clinvar_deletion <- check_cnv_v2(clinvar_match_deletion %>% select(-id))
output_clinvar_duplication <- check_cnv_v2(clinvar_match_duplication  %>% select(-id))
output_decipher_deletion <- check_cnv_v2(decipher_match_deletion  %>% select(-id))
output_decipher_duplication <- check_cnv_v2(decipher_match_duplication  %>% select(-id))
output_clinvar_20 <- check_cnv_v2(clinvar_20_match  %>% select(-id))


# output_clinvar_omim_deletion <- check_cnv_v2(check_omim_matched_del)
# N. Nigreisy writing------------------------------------------------------------------------------


# df_match_deletion_only_pathogenic <- match_patho_benign(output_clinvar_deletion, max(output_clinvar_deletion$length_cnv))
# df_match_duplication_only_pathogenic <- match_patho_benign(output_clinvar_duplication, max(output_clinvar_duplication$length_cnv))
# 
# df_match_deletion_decipher_only_pathogenic <- match_patho_benign(output_decipher_deletion, max(output_decipher_deletion$length_cnv))
# df_match_duplication_decipher_only_pathogenic <- match_patho_benign(output_decipher_duplication, max(output_decipher_duplication$length_cnv))
# 
# df_match_clinvar_20_match <- match_patho_benign(output_clinvar_20, max(output_clinvar_20$length_cnv))
# 
# output_clinvar_20 %>%
#   select(-c(id.x, id.y)) %>% 
#   left_join(submission_nigreisy, by = c('chrom', 'start', 'end', 'type_variant' = 'variant_class')) %>%
#   write_tsv('nigreisy_data/23_01_22/clinvar_21_deletion.tsv')
#   
#   
# 
# output_clinvar_deletion %>% 
#   select(-c(id.x, id.y)) %>% 
#   left_join(submission_nigreisy, by = c('chrom', 'start', 'end', 'type_variant' = 'variant_class')) %>%
#   write_tsv('nigreisy_data/23_01_22/clinvar_deletion.tsv')
# 
# output_clinvar_duplication %>% 
#   select(-c(id.x, id.y, LastEvaluated.x, LastEvaluated.y,
#             NumberSubmitters.x, NumberSubmitters.y)) %>%
#   left_join(submission_nigreisy, by = c('chrom', 'start', 'end', 'type_variant' = 'variant_class')) %>%
#   write_tsv('nigreisy_data/23_01_22/clinvar_duplication.tsv')
# 
# output_decipher_deletion %>% 
#   select(-c(id.x, id.y, LastEvaluated.x, LastEvaluated.y,
#             NumberSubmitters.x, NumberSubmitters.y)) %>%
#   write_tsv('nigreisy_data/23_01_22/decipher_deletion.tsv')
# 
# output_decipher_duplication %>% 
#   select(-c(id.x, id.y, LastEvaluated.x, LastEvaluated.y,
#             NumberSubmitters.x, NumberSubmitters.y)) %>%
#   write_tsv('nigreisy_data/23_01_22/decipher_duplication.tsv')
# 
# 
# df_match_deletion_only_pathogenic %>% write_tsv('nigreisy_data/23_01_22/match_patho_benign_clinvar_deletion.tsv')
# df_match_duplication_only_pathogenic %>% write_tsv('nigreisy_data/23_01_22/match_patho_benign_clinvar_duplication.tsv')
# 
# df_match_deletion_decipher_only_pathogenic %>% write_tsv('nigreisy_data/23_01_22/match_patho_benign_decipher_deletion.tsv')
# df_match_duplication_decipher_only_pathogenic %>% write_tsv('nigreisy_data/23_01_22/match_patho_benign_decipher_duplication.tsv')
# 
# df_match_clinvar_20_match %>% write_tsv('nigreisy_data/23_01_22/match_patho_benign_clinvar_21_deletion.tsv')


# Figure 1------------------------------------------------------------------------------

preproc_flow <- bind_rows(preproc_flow[1,], preproc_flow)


preproc_flow[1,]$n_del_patho_clinvar <- problematic_clinvar %>% 
  left_join(pre_clinvar) %>% 
  filter(variant_class == 'deletion') %>% 
  summarise(total_n = sum(n)) %>% 
  pull() %>%
  sum(preproc_flow[1,]$n_del_patho_clinvar)

preproc_flow[1,]$n_dup_patho_clinvar <- problematic_clinvar %>% 
  left_join(pre_clinvar) %>% 
  filter(variant_class == 'deletion') %>% 
  summarise(total_n = sum(n)) %>% 
  pull() %>%
  sum(preproc_flow[1,]$n_dup_patho_clinvar)

preproc_flow[1,]$n_del_patho_decipher <- problematic_decipher %>% 
  left_join(pre_decipher) %>% 
  filter(variant_class == 'deletion') %>% 
  summarise(total_n = sum(n)) %>% 
  pull() %>%
  sum(preproc_flow[1,]$n_del_patho_decipher)

preproc_flow[1,]$n_dup_patho_decipher <- problematic_decipher %>% 
  left_join(pre_decipher) %>% 
  filter(variant_class == 'deletion') %>% 
  summarise(total_n = sum(n)) %>% 
  pull() %>%
  sum(preproc_flow[1,]$n_dup_patho_decipher)

preproc_flow[2,]$tag <- 'conflict_clinical_assessment'
  

preproc_flow_plot <- preproc_flow %>% slice(1:7)

n_removed <- preproc_flow_plot$n_del_patho_clinvar[1] - preproc_flow_plot$n_del_patho_clinvar[nrow(preproc_flow_plot)]
  
p_flow1 <- tibble(tag = c('CNVs removed', 'CNVs OK'),
       n = c(n_removed, preproc_flow_plot$n_del_patho_clinvar[nrow(preproc_flow_plot)] )) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(x = '', y = perc)) +
    geom_col(aes(fill = tag), color = 'black') +
  geom_label(aes(group = tag, label = paste0(100 * round(perc, 3), '%', ' (', n, ')')), size = 4, position = position_stack(vjust = .5)) +
  labs(x = 'Pathogenic CNVs', fill = 'Legend', y = 'Percentage of CNVs', 
       # title = paste0('Pathogenic deletion CNVs', ' (n = ', n_removed, ')'))
       title = paste0('Pathogenic deletion CNVs', ' (n = ', preproc_flow_plot$n_del_patho_clinvar[1], ')')) +
  scale_y_continuous(label = percent) +
  scale_fill_manual(values = rev(hue_pal()(2))) +
  theme_lucid()


p_flow2 <- preproc_flow_plot %>% 
  mutate(n_del_patho = n_del_patho_clinvar) %>% 
  select(tag, n_del_patho) %>%
  mutate(diff_del_patho = c(0,diff(n_del_patho))) %>%
  mutate(diff_del_patho = abs(diff_del_patho)) %>%
  mutate(perc_diff_del_patho = diff_del_patho / sum(diff_del_patho)) %>%
  select(-n_del_patho) %>%

  filter(tag != 'start') %>%
  #merging 4 rows into 2 categories
  mutate(tag2 = case_when(
    tag == 'conflict_clinical_assessment' ~ 'conflict_clinical_assessment',
    tag == 'problematic_regions' ~ 'problematic_regions',
    tag == 'liftover' ~ 'liftover',
    tag == 'identical_cnvs' ~ 'reciprocal_overlap',
    tag == 'reciprocal_overlap' ~ 'reciprocal_overlap',
    tag == 'conflict_pathogenic_benign' ~ 'conflict_clinical_assessment'
  )) %>%
  group_by(tag2) %>%
  mutate(perc_diff_del_patho2 = sum(perc_diff_del_patho)) %>%
  mutate(diff_del_patho2 = sum(diff_del_patho)) %>%
  ungroup() %>%
  filter(! tag %in% c('conflict_pathogenic_benign', 'identical_cnvs')) %>%
  # mutate(id = row_number()) %>%
  # mutate(id = paste0('Step ', id, ':')) %>%
  # mutate(id2 = row_number()) %>%
  # mutate(tag = paste(id, tag)) %>%
  mutate(tag = fct_reorder(tag, -diff_del_patho)) %>%
  ggplot(aes('', perc_diff_del_patho2)) +
  geom_col(aes(fill = tag), color = 'black') +
  scale_y_continuous(label = percent) +
  geom_label(aes(group = tag, label = paste0(tag, ' - ', 100 * round(perc_diff_del_patho2, 3), '%', ' (', diff_del_patho2, ')')), size = 4, position = position_stack(vjust = .5)) +
  theme_lucid() +
  labs(x = 'Pathogenic CNVs', y = 'Percentage of CNVs removed', fill = 'Legend',
       title = paste0('Pathogenic deletion CNVs removed', ' (n = ', n_removed, ')'))

figure_qc_cnv_removed <- p_flow1 + p_flow2



# + 100 ??
(nrow(output_clinvar_deletion) + nrow(output_decipher_deletion) + nrow(output_clinvar_20)) / 2

# Features enrichment---------------------------------------------

# 
p_enrich1_nobias <- features_enrichment(features_tbl[features_tbl$human_control == 'no',],
                                        output_clinvar_deletion,
                                        output_clinvar_duplication)






# Counting across chromosomes---------------------------------------------

coverage_df <- coord_chrom_hg19 %>% 
  rename(end = length) %>% 
  mutate(start = 1) %>%
  bed_coverage(output_clinvar_deletion %>% filter(clinical == 'pathogenic')) %>%
  rename(patho_del = .frac) %>%
  select(chrom, patho_del) %>%
  left_join(
    
    coord_chrom_hg19 %>% 
      rename(end = length) %>% 
      mutate(start = 1) %>%
      bed_coverage(output_clinvar_deletion %>% filter(clinical == 'benign')) %>%
      rename(benign_del = .frac) %>%
      select(chrom, benign_del)
  ) %>%
  left_join(
    
    coord_chrom_hg19 %>% 
      rename(end = length) %>% 
      mutate(start = 1) %>%
      bed_coverage(output_clinvar_duplication %>% filter(clinical == 'pathogenic')) %>%
      rename(patho_dup = .frac) %>%
      select(chrom, patho_dup)
  ) %>% left_join(
    
    
    coord_chrom_hg19 %>% 
      rename(end = length) %>% 
      mutate(start = 1) %>%
      bed_coverage(output_clinvar_duplication %>% filter(clinical == 'benign')) %>%
      rename(benign_dup = .frac) %>%
      select(chrom, benign_dup)
  ) %>%
  pivot_longer(-chrom, names_to = 'category', values_to = 'coverage') %>%
  separate(category, into = c('clinical', 'type_variant')) %>%
  mutate(chrom = factor(chrom, levels = c(1:22, 'X'))) %>%
  mutate(clinical = if_else(clinical == 'patho', 'pathogenic', clinical)) %>%
  mutate(type_variant = if_else(type_variant == 'del', 'deletion', 'duplication'))


counting_df <- output_clinvar_deletion %>%
  count(chrom, clinical, type_variant) %>%
  bind_rows(output_clinvar_duplication %>% count(chrom, clinical, type_variant)) %>%
  right_join(coverage_df, by = c('chrom', 'clinical', 'type_variant')) %>%
  mutate(chrom = factor(chrom, levels = c(1:22, 'X'))) %>%
  replace_na(list(n = 0, coverage = 0)) %>%
  group_by(type_variant) %>%
  mutate(n_perc = n / sum(n)) %>%
  mutate(coverage = 100 * coverage)

counting_df %>% filter(type_variant == 'deletion') %>% sheet_write(google_calc_results, sheet = "counting_training_deletions")
counting_df %>% filter(type_variant == 'duplication') %>% sheet_write(google_calc_results, sheet = "counting_training_duplications")

figure_coverage <- coverage_df %>%
  ggplot(aes(chrom, coverage)) +
    geom_col(aes(fill = clinical), color = 'black', position = 'dodge') +
    scale_y_continuous(label = percent) +
  scale_fill_manual(values = rev(hue_pal()(2)))+
  facet_wrap(vars(type_variant)) +
  theme_lucid()





# Figure - Correlation between models---------------------------------------------

for (i in c(seq(2,22), 'X')) {
  
  print(i)
  train_split <- output_clinvar_deletion %>% 
    filter(type_variant == 'deletion') %>% 
    filter(!chrom %in% c(i, 1))
  
  test_split <- output_clinvar_deletion %>% filter(type_variant == 'deletion') %>%
    filter(chrom == 1)
  
  prev_tmp_model <- rand_forest() %>%
    set_mode('classification') %>%
    set_engine('ranger') %>%
    fit(human_no_control, data = train_split)
  
  tmp_result <- predict(prev_tmp_model, test_split, type = 'prob') %>%
    select(.pred_pathogenic)
  
  colnames(tmp_result) <- paste0('Chromosome ', i)
  
  if (i == 2) {
    
    result_corr <- tmp_result
  } else {
    
    result_corr <- bind_cols(result_corr, tmp_result)
  }
  
}

figure_correlation <- result_corr %>%
correlate(method = 'pearson') %>% 
  mutate(term = factor(term, levels = as.character(paste0('Chromosome ', c(seq(1,22), 'X'))))) %>%
rplot(print_cor = TRUE, legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))



# ------------------------------------------------------------------------------
# 16 models (6 naive models (logistic) + 4 GBM +  6 Bayesian)
# ------------------------------------------------------------------------------

# Load models

bayesian_clinvar_del_nohuman <- readRDS('models/bayesian_clinvar_del_nohuman.RData')
bayesian_clinvar_dup_nohuman <- readRDS('models/bayesian_clinvar_dup_nohuman.RData')
bayesian_clinvar_del_human <- readRDS('models/bayesian_clinvar_del_human.RData')
bayesian_clinvar_dup_human <- readRDS('models/bayesian_clinvar_dup_human.RData')
bayesian_clinvar_del_both <- readRDS('models/bayesian_clinvar_del_both.RData')
bayesian_clinvar_dup_both <- readRDS('models/bayesian_clinvar_dup_both.RData')


# ------------------------------------------------------------------------------
# NAIVE MODELS
# ------------------------------------------------------------------------------

plan('multiprocess', workers = 5)

# Deletion CNVs------------------------------------------------------------------




# bayesian_clinvar_del_omim <- rtemis_step1(input_tbl = output_clinvar_omim_deletion,
#                                           tag_variant = 'deletion',
#                                           vector_features = vector_no_human_control,
#                                           tag_features <- 'unbiased (only CNVs mapping OMIM genes)',
#                                           input_prior = 'hs',
#                                           nc = 23,
#                                           input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 


tic()

logistic_clinvar_del_length <- chrom_aware(input_tbl = output_clinvar_deletion,
                                    tag_variant = 'deletion',
                                    tag_formule = 'clinical ~ length',
                                    model_name = 'logistic',
                                    formule_model = clinical ~ length_cnv,
                                    list_hyper = NULL)

logistic_clinvar_del_n_genes <- chrom_aware(input_tbl = output_clinvar_deletion,
                                           tag_variant = 'deletion',
                                           tag_formule = 'clinical ~ n_genes',
                                           model_name = 'logistic',
                                           formule_model = clinical ~ n_genes,
                                           list_hyper = NULL)

logistic_clinvar_del_omim <- chrom_aware(input_tbl = output_clinvar_deletion,
                                            tag_variant = 'deletion',
                                            tag_formule = 'clinical ~ omim',
                                            model_name = 'logistic',
                                            formule_model = clinical ~ omim,
                                            list_hyper = NULL)


# Duplication CNVs------------------------------------------------------------------------------


logistic_clinvar_dup_length <- chrom_aware(input_tbl = output_clinvar_duplication,
                                           tag_variant = 'duplication',
                                           tag_formule = 'clinical ~ length',
                                           model_name = 'logistic',
                                           formule_model = clinical ~ length_cnv,
                                           list_hyper = NULL)

logistic_clinvar_dup_n_genes <- chrom_aware(input_tbl = output_clinvar_duplication,
                                            tag_variant = 'duplication',
                                            tag_formule = 'clinical ~ n_genes',
                                            model_name = 'logistic',
                                            formule_model = clinical ~ n_genes,
                                            list_hyper = NULL)

logistic_clinvar_dup_omim <- chrom_aware(input_tbl = output_clinvar_duplication,
                                         tag_variant = 'duplication',
                                         tag_formule = 'clinical ~ omim',
                                         model_name = 'logistic',
                                         formule_model = clinical ~ omim,
                                         list_hyper = NULL)


# ------------------------------------------------------------------------------
# GRADIENT BOOSTING - ClinVar
# ------------------------------------------------------------------------------


# Deletion CNVs------------------------------------------------------------------------------


gbm_clinvar_del_nohuman <- chrom_aware(input_tbl = output_clinvar_deletion,
                                       tag_variant = 'deletion',
                                       tag_formule = 'unbiased_approach',
                                       model_name = 'gbm',
                                       formule_model = human_no_control,
                                       list_hyper = NULL)

gbm_clinvar_del_human <- chrom_aware(input_tbl = output_clinvar_deletion,
                                     tag_variant = 'deletion',
                                     tag_formule = 'knowledge-based',
                                     model_name = 'gbm',
                                     formule_model = human_control,
                                     list_hyper = NULL)


# Duplication CNVs------------------------------------------------------------------------------


gbm_clinvar_dup_nohuman <- chrom_aware(input_tbl = output_clinvar_duplication,
                                       tag_variant = 'duplication',
                                       tag_formule = 'unbiased_approach',
                                       model_name = 'gbm',
                                       formule_model = human_no_control,
                                       list_hyper = NULL)

gbm_clinvar_dup_human <- chrom_aware(input_tbl = output_clinvar_deletion,
                                     tag_variant = 'duplication',
                                     tag_formule = 'knowledge-based',
                                     model_name = 'gbm',
                                     formule_model = human_control,
                                     list_hyper = NULL)



# ------------------------------------------------------------------------------
# BAYESIAN MODEL - ClinVar
# ------------------------------------------------------------------------------


# bayesian_clinvar_del_both <- rtemis_step1(input_tbl = output_clinvar_deletion,
#                                           tag_variant = 'deletion',
#                                           vector_features = vector_total_human_control,
#                                           tag_features <- 'both',
#                                           input_prior = 'hs',
#                                           nc = 23,
#                                           input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 
# bayesian_clinvar_dup_both <- rtemis_step1(input_tbl = output_clinvar_duplication,
#                                           tag_variant = 'duplication',
#                                           vector_features = vector_total_human_control,
#                                           tag_features <- 'both',
#                                           input_prior = 'hs',
#                                           nc = 23,
#                                           input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 


# Deletion CNVs------------------------------------------------------------------------------
# bayesian_clinvar_del_nohuman <- rtemis_step1(input_tbl = output_clinvar_deletion,
#                                              tag_variant = 'deletion',
#                                              vector_features = vector_no_human_control,
#                                              tag_features <- 'no_human',
#                                              input_prior = 'hs',
#                                              nc = 23,
#                                              input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 
# 
# bayesian_decipher_del_nohuman <- rtemis_step1(input_tbl = output_decipher_deletion,
#                                              tag_variant = 'deletion',
#                                              vector_features = vector_no_human_control,
#                                              tag_features <- 'no_human',
#                                              input_prior = 'hs',
#                                              nc = 23,
#                                              input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 
# bayesian_decipher_dup_nohuman <- rtemis_step1(input_tbl = output_decipher_duplication,
#                                               tag_variant = 'duplication',
#                                               vector_features = vector_no_human_control,
#                                               tag_features <- 'no_human',
#                                               input_prior = 'hs',
#                                               nc = 23,
#                                               input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 
# 

# 
# 
# 
# bayesian_clinvar_del_human <- rtemis_step1(input_tbl = output_clinvar_deletion,
#                                            tag_variant = 'deletion',
#                                            vector_features = vector_yes_human_control,
#                                            tag_features <- 'no_human',
#                                            input_prior = 'hs',
#                                            nc = 23,
#                                            input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 
# Duplication CNVs------------------------------------------------------------------------------
# 
# 
# bayesian_clinvar_dup_nohuman <- rtemis_step1(input_tbl = output_clinvar_duplication,
#                                              tag_variant = 'duplication',
#                                              vector_features = vector_no_human_control,
#                                              tag_features <- 'no_human',
#                                              input_prior = 'hs',
#                                              nc = 23,
#                                              input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
# 
# bayesian_clinvar_dup_human <- rtemis_step1(input_tbl = output_clinvar_duplication,
#                                            tag_variant = 'duplication',
#                                            vector_features = vector_yes_human_control,
#                                            tag_features <- 'no_human',
#                                            input_prior = 'hs',
#                                            nc = 23,
#                                            input_hyper = tibble(trees = 500, depth = 2, min_n = 1))






# GENERAL STEP-------------------

tmp_decipher_del_setting_general <- check_cnv_v2(bind_rows(input_check_cnv_del_pathogenic_decipher,
                                                           input_check_cnv_del_benign %>%
                                                            filter(!id_tmp %in% remove_from_clinvar_to_decipher_and_ind_del)) %>% select(-id))

tmp_clinvar_del_setting_general <- check_cnv_v2(bind_rows(input_check_cnv_del_pathogenic_clinvar,
                                                           input_check_cnv_del_benign) %>% select(-id))
# Before 4 -> Now 1
# Before 2 -> Now 2
# Before 1 -> Now 3
# Before 3 -> Now 4
## Important! Check bin_length = 1,000

# Scenario #1 patho (>0 OMIM) - non_patho (>0 OMIM)
# Scenario #2 patho (0 OMIM) - non_patho (0 OMIM)
# Scenario #3 patho (0 OMIM) - non_patho (>0 OMIM)
# Scenario #4 patho (0 genes) - non_patho (0 genes)

# DATASET DECIPHER SETTING I-----------
# Scenario #1 patho (>0 OMIM) - non_patho (>0 OMIM)


decipher_setting_1 <- tmp_decipher_del_setting_general %>%
  filter(omim == 1) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e3)

clinvar_setting_1 <- tmp_clinvar_del_setting_general %>%
  filter(omim == 1) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e2)


# DATASET DECIPHER SETTING II-----------
# Scenario #2 patho (0 OMIM) - non_patho (0 OMIM)

decipher_setting_2 <- tmp_decipher_del_setting_general %>%
  filter(omim == 0) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e3)

clinvar_setting_2 <- tmp_clinvar_del_setting_general %>%
  filter(omim == 0) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e2)

# DATASET DECIPHER SETTING III-----------
# Scenario #3 patho (0 OMIM) - non_patho (>0 OMIM)

decipher_setting_3 <- tmp_decipher_del_setting_general %>%
  filter(clinical == 'pathogenic' & omim == 0 | clinical == 'benign' & omim == 1) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e3)

clinvar_setting_3 <- tmp_clinvar_del_setting_general %>%
  filter(clinical == 'pathogenic' & omim == 0 | clinical == 'benign' & omim == 1) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e2)

# DATASET DECIPHER SETTING IV-----------
# Scenario #4 patho (0 genes) - non_patho (0 genes)


decipher_setting_4 <- tmp_decipher_del_setting_general %>%
  filter(n_genes == 0) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e3)

clinvar_setting_4 <- tmp_clinvar_del_setting_general %>%
  filter(n_genes == 0) %>%
  mutate(id_tmp = row_number()) %>%
  matching_length(bin_length = 1e2)



# Benchmark analysis-------------------

# test55 <-  predict_chrom_aware_rtemis(bayesian_decipher_del_nohuman, output_clinvar_deletion, 'deletion', 'unbiased approach')
# test56 <-  predict_chrom_aware_rtemis(bayesian_decipher_dup_nohuman, output_clinvar_duplication, 'duplication', 'unbiased approach')

# ClinVar ------------------------------------------------------------------------------

result_clinvar_del <- get_results(output_clinvar_deletion, 'DEL', 'ClinVar')

result_clinvar_dup <- get_results(output_clinvar_duplication, 'DUP', 'Clinvar')

# ClinVar independent dataset ------------------------------------------------------------------------------

result_clinvar_ind_del <- get_results(output_clinvar_20, 'DEL', 'ClinVar ind.(> Jan 21)')

# DECIPHER TOTAL ------------------------------------------------------------------------------

result_decipher_del <- get_results(output_decipher_deletion, 'DEL', 'DECIPHER total')

result_decipher_dup <- get_results(output_decipher_duplication, 'DUP', 'DECIPHER total')

# DECIPHER special I ------------------------------------------------------------------------------

result_decipher_setting_1 <- get_results(decipher_setting_1, 'DEL', 'DECIPHER - patho (>0 OMIM) - non_patho (>0 OMIM)')
result_clinvar_setting_1 <- get_results(clinvar_setting_1, 'DEL', 'DECIPHER - patho (>0 OMIM) - non_patho (>0 OMIM)')

# DECIPHER special II ---Reqeuan---------------------------------------------------------------------------

result_decipher_setting_2 <- get_results(decipher_setting_2, 'DEL', 'DECIPHER - patho (0 OMIM) - non_patho (0 OMIM)')
result_clinvar_setting_2 <- get_results(clinvar_setting_2, 'DEL', 'DECIPHER - patho (0 OMIM) - non_patho (0 OMIM)')

# DECIPHER noncoding ------------------------------------------------------------------------------

result_decipher_setting_3 <- get_results(decipher_setting_3, 'DEL', 'DECIPHER - patho (0 OMIM) - non_patho (>0 OMIM)')
result_clinvar_setting_3 <- get_results(clinvar_setting_3, 'DEL', 'DECIPHER - patho (0 OMIM) - non_patho (>0 OMIM)')

# DECIPHER special IV ------------------------------------------------------------------------------

result_decipher_setting_4 <- get_results(decipher_setting_4, 'DEL', 'DECIPHER - patho (0 genes) - non_patho (0 genes)')
result_clinvar_setting_4 <- get_results(clinvar_setting_4, 'DEL', 'DECIPHER - patho (0 genes) - non_patho (0 genes)')



# Before 4 -> Now 1
# Before 2 -> Now 2
# Before 1 -> Now 3
# Before 3 -> Now 4

toc()

# ------------------------------------------------------------------------------
# BANCCO results
# ------------------------------------------------------------------------------

# Description of the number of CNVs across athogeninc / benign, separately for deletions / duplications
# In case we have it: distribution lengths for pathogeninc / benign, separately for deletions / duplications
#  AUROC and PR values on the global sets for CNVscore and the 3 naive models (in case we have it)
#  % of CNVs across uncertainty bins, split by deletions and duplications
# AUROC split per uncertainty bin

# bancco_df <- read_tsv('bancco_results/file_results_BANCCO_2022_05_31.tsv') %>%
#   mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign')))

bancco_decipher <- read_tsv('bancco_results/file_results_BANCCO_2022_06_04.tsv') %>%
  mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign')))





bancco_del_decipher_rel <- get_reliability_score_mid(ref_quantiles_decipher_del, bancco_decipher %>% filter(str_detect(tag, 'deletion - decipher - bayesian - unbiased')))
bancco_dup_decipher_rel <- get_reliability_score_mid(ref_quantiles_decipher_dup, bancco_decipher %>% filter(str_detect(tag, 'duplication - decipher - bayesian - unbiased')))

bancco_del_rel <- get_reliability_score_mid(ref_quantiles, bancco_df %>% filter(str_detect(tag, 'deletion - bayesian - unbiased')))
bancco_dup_rel <- get_reliability_score_mid(ref_quantiles_dup, bancco_df %>% filter(str_detect(tag, 'duplication - bayesian - unbiased')))


bancco_decipher2 <- bancco_decipher %>%
  filter(str_detect(tag, 'unbiased')) %>%
  mutate(tag2 = ifelse(str_detect(tag, 'deletion'), 'Deletion', 'Duplication')) %>%
  mutate(tag = ifelse(str_detect(tag, 'decipher'), 'CNVscore trained on DECIPHER',
                      'CNVscore trained on ClinVar'))

n_bancco_deletions_pathogenic <- bancco_decipher2 %>% filter(tag2 == 'Deletion', str_detect(tag, 'ClinVar')) %>% count(clinical) %>% filter(clinical == 'pathogenic') %>% pull(n) %>% format(big.mark = ',')
n_bancco_deletions_benign <- bancco_decipher2 %>% filter(tag2 == 'Deletion', str_detect(tag, 'ClinVar')) %>% count(clinical) %>% filter(clinical == 'benign') %>% pull(n) %>% format(big.mark = ',')
n_bancco_duplications_pathogenic <- bancco_decipher2 %>% filter(tag2 == 'Duplication', str_detect(tag, 'ClinVar')) %>% count(clinical) %>% filter(clinical == 'pathogenic') %>% pull(n) %>% format(big.mark = ',')
n_bancco_duplications_benign <- bancco_decipher2 %>% filter(tag2 == 'Duplication', str_detect(tag, 'ClinVar')) %>% count(clinical) %>% filter(clinical == 'benign') %>% pull(n) %>% format(big.mark = ',')


p1_bancco <- bancco_decipher2 %>%
  filter(tag2 == 'Deletion') %>%
  group_by(tag, tag2) %>%
  roc_curve(clinical, .pred_pathogenic) %>%
  ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  theme_roc() +
  theme(legend.title = element_blank(), legend.position = 'top')
        # plot.title = element_text(size = 10)) +
  # ggtitle(glue('Bancco dataset ({n_bancco_deletions_pathogenic} pathogenic and {n_bancco_deletions_benign} non-pathogenic deletions)'))


p11_bancco <- bancco_decipher2 %>%
  filter(tag2 == 'Deletion') %>%
  group_by(tag, tag2) %>%
  pr_curve(clinical, .pred_pathogenic) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  # ggtitle(paste0(str_remove(x[[8]], pattern = "(\\d{1}).*"), ' (merged:', x[[12]], ')', '')) +
  theme_pr() +
  theme(legend.title = element_blank(), legend.position = 'top')
  # facet_wrap(vars(tag2), nrow = 2)

p111_bancco <- bancco_del_decipher_rel %>% 
  count(reliability_score) %>% 
  mutate(perc = n / sum(n)) %>%
  mutate(tag = 'CNVscore trained on DECIPHER') %>%
  bind_rows(
    
    bancco_del_rel %>% 
      count(reliability_score) %>% 
      mutate(perc = n / sum(n)) %>%
      mutate(tag = 'CNVscore trained on ClinVar')
  ) %>%
  mutate(reliability_score = case_when(
    reliability_score == 1 ~ 'Lowly uncertain',
    reliability_score == 2 ~ 'Moderately uncertain',
    reliability_score == 3 ~ 'Highly uncertain'
  )) %>%
  mutate(reliability_score = factor(reliability_score, levels = c('Highly uncertain', 'Moderately uncertain', 'Lowly uncertain'))) %>%
  ggplot(aes(tag, perc)) +
  geom_col(aes(fill = reliability_score), color = 'black') +
  scale_y_continuous(label = percent) +
  # scale_fill_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Greens', direction = -1) +
  # scale_fill_manual(values = c('#F8766D', '#619CFF', '#00BA38')) +
  theme_minimal() +
  labs(fill = 'Uncertainty level', y = '% of CNVs across uncertainty CNVscore levels', x = '') +
  theme(axis.text.x = element_text(size = 9))


p2_bancco <- bancco_decipher2 %>%
  filter(tag2 == 'Duplication') %>%
  group_by(tag, tag2) %>%
  roc_curve(clinical, .pred_pathogenic) %>%
  ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  theme_roc() +
  theme(legend.title = element_blank(), legend.position = 'top')
  # ggtitle(glue('Bancco dataset ({n_bancco_duplications_pathogenic} pathogenic and {n_bancco_duplications_benign} non-pathogenic duplications)'))


p22_bancco <- bancco_decipher2 %>%
  filter(tag2 == 'Duplication') %>%
  group_by(tag, tag2) %>%
  pr_curve(clinical, .pred_pathogenic) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  # ggtitle(paste0(str_remove(x[[8]], pattern = "(\\d{1}).*"), ' (merged:', x[[12]], ')', '')) +
  theme_pr() +
  theme(legend.title = element_blank(), legend.position = 'top')
# facet_wrap(vars(tag2), nrow = 2)

p222_bancco <- bancco_dup_decipher_rel %>% 
  count(reliability_score) %>% 
  mutate(perc = n / sum(n)) %>%
  mutate(tag = 'CNVscore trained on DECIPHER') %>%
  bind_rows(
    
    bancco_del_rel %>% 
      count(reliability_score) %>% 
      mutate(perc = n / sum(n)) %>%
      mutate(tag = 'CNVscore trained on ClinVar')
  ) %>%
  mutate(reliability_score = case_when(
    reliability_score == 1 ~ 'Lowly uncertain',
    reliability_score == 2 ~ 'Moderately uncertain',
    reliability_score == 3 ~ 'Highly uncertain'
  )) %>%
  mutate(reliability_score = factor(reliability_score, levels = c('Highly uncertain', 'Moderately uncertain', 'Lowly uncertain'))) %>%
  ggplot(aes(tag, perc)) +
  geom_col(aes(fill = reliability_score), color = 'black') +
  scale_y_continuous(label = percent) +
  # scale_fill_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Greens', direction = -1) +
  # scale_fill_manual(values = c('#F8766D', '#619CFF', '#00BA38')) +
  theme_minimal() +
  labs(fill = 'Uncertainty level', y = '% of CNVs across uncertainty CNVscore levels', x = '') +
  theme(axis.text.x = element_text(size = 9))


p111_bancco + p1_bancco + p11_bancco 
p222_bancco + p2_bancco + p22_bancco 

result_bancco_auroc <- bancco_decipher2 %>%
  group_by(tag, tag2) %>%
  roc_auc(clinical, .pred_pathogenic) %>%
  rename(auroc_total = .estimate) %>%
  select(tag2, auroc_total)

result_bancco_prauc <- bancco_decipher2 %>%
  group_by(tag, tag2) %>%
  pr_auc(clinical, .pred_pathogenic) %>%
  rename(aupr_total = .estimate) %>%
  select(tag2, aupr_total)


# Send to googlesheet ------------------------------------------------------------------------------


bind_rows(
output_clinvar_deletion %>% count(clinical) %>% mutate(source = 'clinvar_del'),
output_clinvar_duplication %>% count(clinical) %>% mutate(source = 'clinvar_dup'),
output_clinvar_20 %>% count(clinical) %>% mutate(source = 'clinvar_ind_del'),
output_decipher_deletion %>% count(clinical) %>% mutate(source = 'decipher_del'),
output_decipher_duplication %>% count(clinical) %>% mutate(source = 'decipher_dup'),
decipher_setting_1 %>% count(clinical) %>% mutate(source = 'decipher_setting_1'),
decipher_setting_2 %>% count(clinical) %>% mutate(source = 'decipher_setting_2'),
decipher_setting_3 %>% count(clinical) %>% mutate(source = 'decipher_setting_3'),
decipher_setting_4 %>% count(clinical) %>% mutate(source = 'decipher_setting_4'),
clinvar_setting_1 %>% count(clinical) %>% mutate(source = 'clinvar_setting_1'),
clinvar_setting_2 %>% count(clinical) %>% mutate(source = 'clinvar_setting_2'),
clinvar_setting_3 %>% count(clinical) %>% mutate(source = 'clinvar_setting_3'),
clinvar_setting_4 %>% count(clinical) %>% mutate(source = 'clinvar_setting_4'),

) %>%
  pivot_wider(id_cols = source, names_from = clinical, values_from = n) %>%
  sheet_write(google_calc_results, sheet = "Supp. Table 5")


to_googlesheet2(result_clinvar_del, output_clinvar_deletion) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 6")
result_clinvar_del_realscore %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 7")
to_googlesheet(result_clinvar_dup) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 8")
to_googlesheet(result_clinvar_ind_del) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 9")
to_googlesheet2(result_decipher_del, output_decipher_deletion) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 10")
result_decipher_del_realscore %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 11")
to_googlesheet(result_decipher_dup) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 12")
to_googlesheet(result_decipher_setting_1) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 13")
to_googlesheet(result_decipher_setting_2) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 14")
to_googlesheet(result_decipher_setting_3) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 15")
to_googlesheet2(result_decipher_setting_4, decipher_setting_4) %>% 
  sheet_write(google_calc_results, sheet = "Supp. Table 16")
to_googlesheet2(result_clinvar_setting_1, clinvar_setting_1) %>% gt()
  sheet_write(google_calc_results, sheet = "Supp. Table 17")
to_googlesheet2(result_clinvar_setting_2, clinvar_setting_2) %>%  gt()
  sheet_write(google_calc_results, sheet = "Supp. Table 18") 
to_googlesheet2(result_clinvar_setting_3, clinvar_setting_3) %>% gt()
  sheet_write(google_calc_results, sheet = "Supp. Table 19")
to_googlesheet2(result_clinvar_setting_4, clinvar_setting_4) %>% gt()
  sheet_write(google_calc_results, sheet = "Supp. Table 20")

ref_sd_clinvar_del_real %>% count(reliability_score) %>%
  mutate(Dataset = 'ClinVar - Deletion CNVs') %>%
  bind_rows(ref_sd_decipher_del_real %>% count(reliability_score) %>%
              mutate(Dataset = 'DECIPHER - Deletion CNVs')) %>%
  rename(`Uncertainty level` = reliability_score) %>%
  select(Dataset, `Uncertainty level`, n) %>%
  sheet_write(google_calc_results, sheet = "Supp. Table 21")

ref_sd_clinvar_del_real %>% count(reliability_score) %>%
  mutate(Dataset = 'ClinVar - Deletion CNVs') %>%
  bind_rows(ref_sd_decipher_del_real %>% count(reliability_score) %>%
              mutate(Dataset = 'DECIPHER - Deletion CNVs')) %>%
  rename(`Uncertainty level` = reliability_score) %>%
  select(Dataset, `Uncertainty level`, n) %>%
  gt()



# Table formatting ------------------------------------------------------------------------------
# install.packages( "gt", repos = "https://mran.revolutionanalytics.com/snapshot/2020-09-22/")

library(gt)

# Table 2
left_join(
generate_table(result_clinvar_del, output_clinvar_deletion, tag = 'ClinVar', as_table = FALSE) %>%
  rename(`AUROC (ClinVar)` = AUROC,
         `AUPR (ClinVar)` = AUPR),
generate_table(result_clinvar_ind_del, output_clinvar_20, tag = 'ClinVar (> Jan 2020)', as_table = FALSE) %>%
  rename(`AUROC (ClinVar ind.)` = AUROC,
         `AUPR (ClinVar ind.)` = AUPR),
by = 'Model'
) %>%
  left_join(generate_table(result_decipher_del, output_decipher_deletion, tag = 'DECIPHER', as_table = FALSE) %>%
              rename(`AUROC (DECIPHER)` = AUROC,
                     `AUPR (DECIPHER)` = AUPR), by = 'Model') %>%
  mutate(Model = ifelse(Model == 'TADA', 'TADA (*)', Model)) %>%
  mutate(Model = ifelse(Model == 'X-CNV', 'X-CNV (*)', Model)) %>%
  write_excel_csv('paper.xlsx')

generate_table(result_clinvar_dup, output_clinvar_duplication, tag = 'ClinVar', as_table = FALSE)
generate_table(result_decipher_dup, output_decipher_duplication, tag = 'DECIPHER', as_table = FALSE)


to_googlesheet(result_clinvar_dup) %>%
  rename(`AUROC (ClinVar)` = roc_auc,
         `AUPR (ClinVar)` = pr_auc) %>%
  left_join(to_googlesheet(result_decipher_dup) %>%
              rename(`AUROC (DECIPHER)` = roc_auc,
                     `AUPR (DECIPHER)` = pr_auc)) %>%
  map_dfr(function(x) {
    if (is.numeric(x)) {
      round(x,3)*100
    } else {
      x
    }
    
  }) %>%
  rename(Model = tag) %>%
  write_excel_csv('paper2.xlsx')

# Table 3
generate_table_rel(result_clinvar_del, res_df, tag = 'ClinVar')
generate_table_rel(result_decipher_del, ref_sd_decipher_del_real, tag = 'DECIPHER')

generate_table_rel(result_clinvar_dup, res_df_dup, tag = 'ClinVar', tag2 = 'Duplication')
generate_table_rel(result_decipher_dup, ref_sd_decipher_dup_real, tag = 'DECIPHER', tag2 = 'Duplication')



# Table 4

generate_table(result_decipher_setting_1, decipher_setting_1, tag = 'Scenario #1', as_table = FALSE) %>%
  select(Model, AUROC) %>%
  rename(`AUROC (Scenario #1)` = AUROC) %>%
  left_join(
generate_table(result_decipher_setting_2, decipher_setting_2, tag = 'Scenario #2', as_table = FALSE) %>%
  select(Model, AUROC) %>%
  rename(`AUROC (Scenario #2)` = AUROC),
by = 'Model') %>%
  left_join(
generate_table(result_decipher_setting_3, decipher_setting_3, tag = 'Scenario #3', as_table = FALSE) %>%
  select(Model, AUROC) %>%
  rename(`AUROC (Scenario #3)` = AUROC),
by = 'Model') %>%
  left_join(
generate_table(result_decipher_setting_4, decipher_setting_4, tag = 'Scenario #4', as_table = FALSE) %>%
  select(Model, AUROC) %>%
  rename(`AUROC (Scenario #4)` = AUROC),
by = 'Model'
) %>%
  write_excel_csv('paper3.xlsx')




generate_table(result_clinvar_setting_1, clinvar_setting_1, tag = 'Scenario #1', as_table = FALSE) %>%
  select(Model, AUROC) %>%
  rename(`AUROC (Scenario #1)` = AUROC) %>%
  left_join(
    generate_table(result_clinvar_setting_2, clinvar_setting_2, tag = 'Scenario #2', as_table = FALSE) %>%
      select(Model, AUROC) %>%
      rename(`AUROC (Scenario #2)` = AUROC),
    by = 'Model') %>%
  left_join(
    generate_table(result_clinvar_setting_3, clinvar_setting_3, tag = 'Scenario #3', as_table = FALSE) %>%
      select(Model, AUROC) %>%
      rename(`AUROC (Scenario #3)` = AUROC),
    by = 'Model') %>%
  left_join(
    generate_table(result_clinvar_setting_4, clinvar_setting_4, tag = 'Scenario #4', as_table = FALSE) %>%
      select(Model, AUROC) %>%
      rename(`AUROC (Scenario #4)` = AUROC),
    by = 'Model'
  ) %>%
  write_excel_csv('paper4.xlsx')





 

# AUROC- Perc. Reliab. 1 ------------------------------------------------------------------------------

ref_sd_clinvar_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_deletion, 'deletion', 'unbiased approach')
ref_sd_clinvarind_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_20, 'deletion', 'unbiased approach')
ref_sd_decipher_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_decipher_deletion, 'deletion', 'unbiased approach')
ref_sd_scenario_1 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, decipher_setting_1, 'deletion', 'unbiased approach')
ref_sd_scenario_2 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, decipher_setting_2, 'deletion', 'unbiased approach')
ref_sd_scenario_3 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, decipher_setting_3, 'deletion', 'unbiased approach')
ref_sd_scenario_4 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, decipher_setting_4, 'deletion', 'unbiased approach')
ref_sd_scenario_1_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, clinvar_setting_1, 'deletion', 'unbiased approach')
ref_sd_scenario_2_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, clinvar_setting_2, 'deletion', 'unbiased approach')
ref_sd_scenario_3_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, clinvar_setting_3, 'deletion', 'unbiased approach')
ref_sd_scenario_4_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, clinvar_setting_4, 'deletion', 'unbiased approach')

ref_sd_clinvar_dup <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_nohuman, output_clinvar_duplication, 'duplication', 'unbiased approach')
ref_sd_decipher_dup <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_nohuman, output_decipher_duplication, 'duplication', 'unbiased approach')

# check 
ref_sd_clinvar_del_real <- get_reliability_score_mid(ref_quantiles, ref_sd_clinvar_del[[3]])
ref_sd_clinvarind_del_real <- get_reliability_score_mid(ref_quantiles, ref_sd_clinvarind_del[[3]])
ref_sd_decipher_del_real <- get_reliability_score_mid(ref_quantiles, ref_sd_decipher_del[[3]])
ref_sd_scenario_1_real_decipher <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_1[[3]])
ref_sd_scenario_2_real_decipher <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_2[[3]])
ref_sd_scenario_3_real_decipher <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_3[[3]]) 
ref_sd_scenario_4_real_decipher <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_4[[3]]) 
ref_sd_scenario_1_real_clinvar <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_1_clinvar[[3]])
ref_sd_scenario_2_real_clinvar <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_2_clinvar[[3]])
ref_sd_scenario_3_real_clinvar <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_3_clinvar[[3]]) 
ref_sd_scenario_4_real_clinvar <- get_reliability_score_mid(ref_quantiles, ref_sd_scenario_4_clinvar[[3]]) 

# ref_sd_clinvar_del_real <- get_reliability_score_mid(ref_quantiles, ref_sd_clinvar_del[[3]])
ref_sd_decipher_dup_real <- get_reliability_score_mid(ref_quantiles_dup, ref_sd_decipher_dup[[3]])



ref_sd_auc <- bind_rows(
generate_table(result_clinvar_del, 
               output_clinvar_deletion, tag = 'ClinVar', FALSE) %>%
  mutate(tag2 = 'ClinVar'),
generate_table(result_clinvar_ind_del, 
               output_clinvar_20, tag = 'ClinVar (> Jan 2020)', FALSE) %>%
  mutate(tag2 = 'ClinVar (> Jan 2020)'),
generate_table(result_decipher_del, 
               output_decipher_deletion, tag = 'DECIPHER', FALSE) %>%
  mutate(tag2 = 'DECIPHER'),
generate_table(result_decipher_setting_1, 
               decipher_setting_1, tag = 'Scenario #1', FALSE) %>%
  mutate(tag2 = 'Scenario #1 - DECIPHER'),
generate_table(result_decipher_setting_2, 
               decipher_setting_2, tag = 'Scenario #2', FALSE) %>%
  mutate(tag2 = 'Scenario #2 - DECIPHER'),
generate_table(result_decipher_setting_3, 
               decipher_setting_3, tag = 'Scenario #3', FALSE) %>%
  mutate(tag2 = 'Scenario #3 - DECIPHER'),
generate_table(result_decipher_setting_4, 
               decipher_setting_4, tag = 'Scenario #4', FALSE) %>%
  mutate(tag2 = 'Scenario #4 - DECIPHER'),
generate_table(result_clinvar_setting_1, 
               clinvar_setting_1, tag = 'Scenario #1', FALSE) %>%
  mutate(tag2 = 'Scenario #1 - ClinVar'),
generate_table(result_clinvar_setting_2, 
               clinvar_setting_2, tag = 'Scenario #2', FALSE) %>%
  mutate(tag2 = 'Scenario #2 - ClinVar'),
generate_table(result_clinvar_setting_3, 
               clinvar_setting_3, tag = 'Scenario #3', FALSE) %>%
  mutate(tag2 = 'Scenario #3 - ClinVar'),
generate_table(result_clinvar_setting_4, 
               clinvar_setting_4, tag = 'Scenario #4', FALSE) %>%
  mutate(tag2 = 'Scenario #4 - ClinVar')
) %>%
  filter(Model %in% c('X-CNV', 'STRVCTVRE', 'CNVscore', 'TADA')) %>%
  select(-AUPR)

  


test100 <- bind_rows(ref_sd_clinvar_del_real %>% mutate(tag2 = 'ClinVar'), 
                     ref_sd_clinvarind_del_real  %>% mutate(tag2 = 'ClinVar (> Jan 2020)'), 
                     ref_sd_decipher_del_real %>% mutate(tag2 = 'DECIPHER'),
                     ref_sd_scenario_1_real_decipher %>% mutate(tag2 = 'Scenario #1 - DECIPHER'),
                     ref_sd_scenario_2_real_decipher %>% mutate(tag2 = 'Scenario #2 - DECIPHER'),
                     ref_sd_scenario_3_real_decipher %>% mutate(tag2 = 'Scenario #3 - DECIPHER'),
                     ref_sd_scenario_4_real_decipher %>% mutate(tag2 = 'Scenario #4 - DECIPHER'),
                     ref_sd_scenario_1_real_clinvar %>% mutate(tag2 = 'Scenario #1 - ClinVar'),
                     ref_sd_scenario_2_real_clinvar %>% mutate(tag2 = 'Scenario #2 - ClinVar'), 
                     ref_sd_scenario_3_real_clinvar %>% mutate(tag2 = 'Scenario #3 - ClinVar'), 
                     ref_sd_scenario_4_real_clinvar %>% mutate(tag2 = 'Scenario #4 - ClinVar') 

          ) %>%
  count(tag2, reliability_score) %>%
  group_by(tag2) %>%
  mutate(perc = n / sum(n)) %>%
  left_join(ref_sd_auc, by = 'tag2') %>%
  mutate(source = if_else(str_detect(tag2, 'ClinVar'), 'ClinVar', 'DECIPHER'))


plot_reli1 <- test100 %>%
  group_by(tag2) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(reliability_score = fct_rev(as.factor(reliability_score))) %>%
  ggplot(aes(tag2, perc)) +
    geom_col(aes(fill = as.factor(reliability_score)),
             color = 'black') +
  scale_y_continuous(label = percent) +
  labs(x = 'Dataset', y = 'Percentage', fill = 'Uncertainty level') +
  theme_minimal()


library(ggrepel)

plot_reli3 <- test100 %>%
  mutate(AUROC = AUROC / 100) %>%
  filter(!str_detect(tag2, 'Scenario #4')) %>%
  mutate(tag2 = str_remove(tag2, ' - ClinVar| - DECIPHER')) %>%
  mutate(tag2 = str_remove(tag2, 'Scenario ')) %>%
  filter(reliability_score == 3) %>%
  ggplot(aes(perc, AUROC)) +
  geom_smooth(aes(group = Model), method = 'lm',  level = 0.90) +
  geom_text_repel(aes(label = tag2),
             size = 4
             ) +
  geom_point(aes(fill = source), color = 'black', shape = 21, size = 3) +
  scale_x_continuous(label = percent) +
  coord_cartesian(ylim  = c(0.45,1)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  labs(x = 'Percentage of CNVs associated with high uncertainty') +
  facet_wrap(vars(Model))


plot_reli1 + plot_reli2

plot_reli3



# Write models ------------------------------------------------------------------------------

# saveRDS(bayesian_clinvar_del_nohuman, 'api/bayesian_clinvar_del_nohuman.RData')
# saveRDS(bayesian_clinvar_dup_nohuman, 'api/bayesian_clinvar_dup_nohuman.RData')
# saveRDS(bayesian_clinvar_del_human, 'api/bayesian_clinvar_del_human.RData')
# saveRDS(bayesian_clinvar_dup_human, 'api/bayesian_clinvar_dup_human.RData')
# saveRDS(bayesian_clinvar_del_both, 'api/bayesian_clinvar_del_both.RData')
# saveRDS(bayesian_clinvar_dup_both, 'bayesian_clinvar_dup_both.RData')
# saveRDS(bayesian_clinvar_del_omim, 'bayesian_clinvar_del_omim.RData')

# saveRDS(bayesian_decipher_del_nohuman, 'bayesian_decipher_del_nohuman.RData')
# saveRDS(bayesian_decipher_dup_nohuman, 'bayesian_decipher_dup_nohuman.RData')


# TSNE ------------------------------------------------------------------------------


selected_features <- c(features_tbl[features_tbl$human_control == 'no',]$name, 
                       'tag', 'reliability_score', 'clinical')

a <- output_clinvar_deletion %>% 
  mutate(tag = 'training_data') %>%
  bind_rows(output_decipher_deletion %>%
              mutate(tag = 'decipher_dataset') %>%
              left_join(realscore_decipher_del %>% 
                          select(id, reliability_score), by = 'id')) %>%
  select(all_of(selected_features)) %>%
  distinct() %>%
  mutate(id_tsne = row_number())
  

tsne_a <- Rtsne(a %>% select(-c(tag, reliability_score, clinical)) %>% scale())

tSNE_df <- tsne_a$Y %>% 
  as_tibble() %>%
  rename(tSNE1="V1",
         tSNE2="V2") %>%
  mutate(id_tsne = row_number())

tSNE_df %>%
  left_join(a %>% select(id_tsne, tag, clinical, reliability_score), by="id_tsne") %>%
  # mutate(reliability_score = ifelse(is.na(reliability_score), 4, reliability_score)) %>%
  filter(reliability_score %in% c(1,2,3, NA)) %>%
  mutate(reliability_score = ifelse(is.na(reliability_score), 'ClinVar', paste('DECIPHER - ','Reliab_score:', reliability_score))) %>%
  mutate(reliability_score = fct_rev(reliability_score)) %>%
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = as.factor(reliability_score)), color = 'black', 
             shape = 21, size = 2) +
  # facet_wrap(vars(tag)) +
  scale_fill_manual(values = c(hue_pal()(3), '#b3b1b1')) +
  # stat_ellipse(aes(color = clinical)) +
  theme_minimal()

# Two case studies ------------------------------------------------------------------------------

realscore_decipher_del


# Figure N. Comparison STRVCTVRE ~ unbiased approach (sd) ------------------------------------------------------------------------------

from_hg38_to_hg19 = import.chain('/data-cbl/liftover/hg38ToHg19.over.chain')

comp_structure_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_deletion, 'deletion', 'unbiased approach')
comp_structure_both_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_del_both, output_clinvar_deletion, 'deletion', 'both approach')

comp_structure_clinvarind <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_20, 'deletion', 'unbiased approach')

comp_structure_decipher <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_decipher_deletion, 'deletion', 'unbiased approach')
comp_structure_both_decipher <- predict_chrom_aware_rtemis(bayesian_clinvar_del_both, output_decipher_deletion, 'deletion', 'both approach')

comp_structure_clinvar1 <- comp_structure_clinvar[[3]] %>%
  left_join(output_clinvar_deletion %>% select(chrom, start, end, id)) %>%
  rename(cnvscore_score = .pred_pathogenic) %>%
  select(chrom, start, end, cnvscore_score, sd, clinical, id)

comp_structure_both_clinvar1 <- comp_structure_both_clinvar[[3]] %>%
  left_join(output_clinvar_deletion %>% select(chrom, start, end, id)) %>%
  rename(cnvscore_score = .pred_pathogenic) %>%
  select(chrom, start, end, cnvscore_score, sd, clinical, id)

comp_structure_clinvarind1 <- comp_structure_clinvarind[[3]] %>%
  left_join(output_clinvar_20 %>% select(chrom, start, end, id)) %>%
  rename(cnvscore_score = .pred_pathogenic) %>%
  select(chrom, start, end, cnvscore_score, sd, clinical, id)

comp_structure_decipher1 <- comp_structure_decipher[[3]] %>%
  left_join(output_decipher_deletion %>% select(chrom, start, end, id)) %>%
  rename(cnvscore_score = .pred_pathogenic) %>%
  # mutate(sd = mad) %>%
  select(chrom, start, end, cnvscore_score, sd, clinical, id, mad)

comp_structure_both_decipher1 <- comp_structure_both_decipher[[3]] %>%
  left_join(output_decipher_deletion %>% select(chrom, start, end, id)) %>%
  rename(cnvscore_score = .pred_pathogenic) %>%
  select(chrom, start, end, cnvscore_score, sd, clinical, id)


# ClinVar dataset


tmp1 <- read_tsv('structure_tmp/output_clinvar.bed', 
              col_names = c('chrom', 'start', 'end')) %>%
              mutate(start = start + 1) %>%
              mutate(chrom = str_remove(chrom, 'chr')) %>%
              rename(.pred_pathogenic = X5) %>%
              mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
              mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
              rename(structure_score = .pred_pathogenic) %>%
              select(chrom, start, end, structure_score) %>%
              GRanges()

seqlevelsStyle(tmp1) = "UCSC" 

tmp2 <- liftOver(tmp1, from_hg38_to_hg19) %>%
  as_tibble() %>%
  rename(chrom = seqnames) %>%
  select(chrom, start, end, structure_score) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  left_join(comp_structure_clinvar1) %>%
  na.omit() %>%
  select(-c('chrom', 'start', 'end'))
  

tmp_corr_score <- tmp2 %>% 
  select(structure_score, cnvscore_score) %>% 
  correlate(method = 'spearman') %>%
  slice(1) %>%
  .[[3]] %>%
  round(2)


structure_corr1 <- tmp2 %>% 
  select(structure_score, cnvscore_score) %>%
  ggplot(aes(structure_score, cnvscore_score)) +
  theme_minimal() +
    geom_point() +
    geom_smooth() +
  ggtitle('ClinVar dataset',
          subtitle = glue('Spearman correlation: {tmp_corr_score}')) +
  labs(x = 'STRVCTVRE score', y = 'Unbiased CNVscore')


structure_auc1 <- 1:10 %>% map_df(function(x) {
  
  
  tmp2 %>% 
    mutate(tile_sd = ntile(sd, 10)) %>%
    filter(tile_sd == x) %>%
    pivot_longer(-c(sd, tile_sd, clinical), names_to = 'tool', values_to = 'score') %>%
    group_by(tool) %>%
    roc_auc(clinical, score) %>%
    mutate(tile_sd = x)
  
}) %>%
  ggplot(aes(tile_sd, .estimate)) +
  geom_line(aes(group = tool)) +
  geom_point(aes(fill = tool), shape = 21, color = 'black') +
  theme_minimal() +
  labs(title = 'ClinVar dataset') +
  tmp2 %>%
    mutate(tile_sd = ntile(sd, 10)) %>%
    count(clinical, tile_sd) %>%
    group_by(tile_sd) %>%
    mutate(perc = n / sum(n)) %>%
    ggplot(aes(tile_sd, perc)) +
      geom_col(aes(fill = clinical), color = 'black') +
      scale_y_continuous(label = percent) +
      theme_minimal() +
    labs(x = 'SD (Percentile 1-10)', y = 'ROC-AUC')




# ClinVar independent dataset
tmp1 <- read_tsv('structure_tmp/output_clinvar_ind.bed', 
                 col_names = c('chrom', 'start', 'end')) %>%
  mutate(start = start + 1) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  rename(.pred_pathogenic = X5) %>%
  mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
  mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
  rename(structure_score = .pred_pathogenic) %>%
  select(chrom, start, end, structure_score) %>%
  GRanges()

seqlevelsStyle(tmp1) = "UCSC" 

tmp2 <- liftOver(tmp1, from_hg38_to_hg19) %>%
  as_tibble() %>%
  rename(chrom = seqnames) %>%
  select(chrom, start, end, structure_score) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  left_join(comp_structure_clinvarind1) %>%
  na.omit() %>%
  select(-c('chrom', 'start', 'end'))


tmp_corr_score <- tmp2 %>% 
  select(structure_score, cnvscore_score) %>% 
  correlate(method = 'spearman') %>%
  slice(1) %>%
  .[[3]] %>%
  round(2)


structure_corr2 <- tmp2 %>% 
  select(structure_score, cnvscore_score) %>%
  ggplot(aes(structure_score, cnvscore_score)) +
  theme_minimal() +
  geom_point() +
  geom_smooth() +
  ggtitle('ClinVar (> Jan 2020) dataset',
          subtitle = glue('Spearman correlation: {tmp_corr_score}')) +
  labs(x = 'STRVCTVRE score', y = 'Unbiased CNVscore')


structure_auc2 <- 1:10 %>% map_df(function(x) {
  
  
  tmp2 %>% 
    mutate(tile_sd = ntile(sd, 10)) %>%
    filter(tile_sd == x) %>%
    pivot_longer(-c(sd, tile_sd, clinical), names_to = 'tool', values_to = 'score') %>%
    group_by(tool) %>%
    roc_auc(clinical, score) %>%
    mutate(tile_sd = x)
  
}) %>%
  ggplot(aes(tile_sd, .estimate)) +
  geom_line(aes(group = tool)) +
  geom_point(aes(fill = tool), shape = 21, color = 'black') +
  theme_minimal() +
  labs(title = 'ClinVar independent dataset (> Jan 2020)') +
  tmp2 %>%
  mutate(tile_sd = ntile(sd, 10)) %>%
  count(clinical, tile_sd) %>%
  group_by(tile_sd) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(tile_sd, perc)) +
  geom_col(aes(fill = clinical), color = 'black') +
  scale_y_continuous(label = percent) +
  theme_minimal() +
  labs(x = 'SD (Percentile 1-10)', y = 'ROC-AUC')

# DECIPHER dataset
tmp1 <- read_tsv('structure_tmp/output_decipher.bed', 
                 col_names = c('chrom', 'start', 'end')) %>%
  mutate(start = start + 1) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  rename(.pred_pathogenic = X5) %>%
  mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
  mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
  rename(structure_score = .pred_pathogenic) %>%
  select(chrom, start, end, structure_score) %>%
  GRanges()

seqlevelsStyle(tmp1) = "UCSC" 

tmp2 <- liftOver(tmp1, from_hg38_to_hg19) %>%
  as_tibble() %>%
  rename(chrom = seqnames) %>%
  select(chrom, start, end, structure_score) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  left_join(comp_structure_decipher1) %>%
  na.omit() %>%
  select(-c('chrom', 'start', 'end'))

tmp_corr_score <- tmp2 %>% 
  select(structure_score, cnvscore_score) %>% 
  correlate(method = 'spearman') %>%
  slice(1) %>%
  .[[3]] %>%
  round(2)


structure_corr3 <- tmp2 %>% 
  select(structure_score, cnvscore_score) %>%
  ggplot(aes(structure_score, cnvscore_score)) +
  theme_minimal() +
  geom_point() +
  geom_smooth() +
  ggtitle('DECIPHER dataset',
          subtitle = glue('Spearman correlation: {tmp_corr_score}')) +
  labs(x = 'STRVCTVRE score', y = 'Unbiased CNVscore')


structure_auc3 <- 1:10 %>% map_df(function(x) {
  
  
  tmp2 %>% 
    mutate(tile_sd = ntile(sd, 10)) %>%
    filter(tile_sd == x) %>%
    pivot_longer(-c(id, sd, tile_sd, clinical), names_to = 'tool', values_to = 'score') %>%
    group_by(tool) %>%
    roc_auc(clinical, score) %>%
    mutate(tile_sd = x)
  
}) %>%
  ggplot(aes(tile_sd, .estimate)) +
  geom_line(aes(group = tool)) +
  geom_point(aes(fill = tool), shape = 21, color = 'black') +
  theme_minimal() +
  labs(title = 'DECIPHER dataset')

structure_decipher_sd_barplot <- tmp2 %>%
  mutate(tile_sd = ntile(sd, 10)) %>%
  count(clinical, tile_sd) %>%
  group_by(tile_sd) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(tile_sd, perc)) +
  geom_col(aes(fill = clinical), color = 'black') +
  scale_y_continuous(label = percent) +
  theme_minimal() +
  labs(x = 'SD (Percentile 1-10)', y = 'ROC-AUC')


# ggsave(glue("figures/{current_date}/auc_strvctvre_unbiased_cnvscore.png"),
#        structure_auc1 + structure_auc2 + structure_auc3 + plot_layout(nrow = 3)
#        , width = 17, height = 9.6, dpi = 300, units = "in", device='png')



# Figure N. Comparison STRVCTVRE ~ unbiased + total (distance) ------------------------------------------------------------------------------

## ClinVar

# Unbiased model

distances_tbl <- get_distances(
  training_tbl = output_clinvar_deletion,
  input_tbl = output_clinvar_deletion,
  formule_input = human_no_control, 
  threshold = 0.95
)


comp_structure_clinvar2 <- comp_structure_clinvar1 %>%
  left_join(distances_tbl %>% select(id, dist_pathogenic, dist_benign), by = 'id')


tmp1 <- read_tsv('structure_tmp/output_clinvar.bed', 
                 col_names = c('chrom', 'start', 'end')) %>%
  mutate(start = start + 1) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  rename(.pred_pathogenic = X5) %>%
  mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
  mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
  rename(structure_score = .pred_pathogenic) %>%
  select(chrom, start, end, structure_score) %>%
  GRanges()

seqlevelsStyle(tmp1) = "UCSC" 

plot_rates_clinvar_unbiased <- liftOver(tmp1, from_hg38_to_hg19) %>%
  as_tibble() %>%
  rename(chrom = seqnames) %>%
  select(chrom, start, end, structure_score) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  left_join(comp_structure_clinvar2) %>%
  na.omit() %>%
  select(-c('chrom', 'start', 'end'))


# Total model 

# distances_tbl_both <- get_distances(
#   training_tbl = output_clinvar_deletion,
#   input_tbl = output_clinvar_deletion,
#   formule_input = formule_total, 
#   threshold = 0.95
# )
# 
# 
# comp_structure_both_clinvar <- comp_structure_both_clinvar1 %>%
#   left_join(distances_tbl_both %>% select(id, dist_pathogenic, dist_benign), by = 'id')
# 
# 
# tmp1 <- read_tsv('structure_tmp/output_clinvar.bed', 
#                  col_names = c('chrom', 'start', 'end')) %>%
#   mutate(start = start + 1) %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   rename(.pred_pathogenic = X5) %>%
#   mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
#   mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
#   rename(structure_score = .pred_pathogenic) %>%
#   select(chrom, start, end, structure_score) %>%
#   GRanges()
# 
# seqlevelsStyle(tmp1) = "UCSC" 
# 
# tmp2 <- liftOver(tmp1, from_hg38_to_hg19) %>%
#   as_tibble() %>%
#   rename(chrom = seqnames) %>%
#   select(chrom, start, end, structure_score) %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   left_join(comp_structure_both_clinvar) %>%
#   na.omit() %>%
#   select(-c('chrom', 'start', 'end'))


## DECIPHER

# Unbiased model

distances_tbl <- get_distances(
  training_tbl = output_clinvar_deletion,
  input_tbl = output_decipher_deletion,
  formule_input = human_no_control, 
  threshold = 0.95
)


comp_structure_decipher2 <- comp_structure_decipher1 %>%
  left_join(distances_tbl %>% select(id, dist_pathogenic, dist_benign), by = 'id')


tmp1 <- read_tsv('structure_tmp/output_decipher.bed', 
                 col_names = c('chrom', 'start', 'end')) %>%
  mutate(start = start + 1) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  rename(.pred_pathogenic = X5) %>%
  mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
  mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
  rename(structure_score = .pred_pathogenic) %>%
  select(chrom, start, end, structure_score) %>%
  GRanges()

seqlevelsStyle(tmp1) = "UCSC" 

plot_rates_decipher_unbiased <- liftOver(tmp1, from_hg38_to_hg19) %>%
  as_tibble() %>%
  rename(chrom = seqnames) %>%
  select(chrom, start, end, structure_score) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  left_join(comp_structure_decipher2) %>%
  na.omit() %>%
  select(-c('chrom', 'start', 'end'))


## ClinVar INDEPENDENT

distances_tbl <- get_distances(
  training_tbl = output_clinvar_deletion,
  input_tbl = output_clinvar_20,
  formule_input = human_no_control, 
  threshold = 0.95
)


comp_structure_clinvarind <- comp_structure_clinvarind1 %>%
  left_join(distances_tbl %>% select(id, dist_pathogenic, dist_benign), by = 'id')


tmp1 <- read_tsv('structure_tmp/output_clinvar_ind.bed', 
                 col_names = c('chrom', 'start', 'end')) %>%
  mutate(start = start + 1) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  rename(.pred_pathogenic = X5) %>%
  mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
  mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
  rename(structure_score = .pred_pathogenic) %>%
  select(chrom, start, end, structure_score) %>%
  GRanges()

seqlevelsStyle(tmp1) = "UCSC" 

plot_rates_clinvarind_unbiased <- liftOver(tmp1, from_hg38_to_hg19) %>%
  as_tibble() %>%
  rename(chrom = seqnames) %>%
  select(chrom, start, end, structure_score) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  left_join(comp_structure_clinvarind) %>%
  na.omit() %>%
  select(-c('chrom', 'start', 'end'))


# supp. figure 3 . Comparison gbm and Bayesian approach ------------------------------------------------------------------------------

comp_gbm_bay1 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_deletion, 'deletion', 'unbiased approach')

comp_gbm_bay2 <- predict_chrom_aware(gbm_clinvar_del_nohuman, output_clinvar_deletion)

comp_gbm_bay3 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_decipher_deletion, 'deletion', 'unbiased approach')

comp_gbm_bay4 <- predict_chrom_aware(gbm_clinvar_del_nohuman, output_decipher_deletion)


p_comp_gbm_bay_clinvar <- comp_gbm_bay1$tmp_roc_curve %>%
  bind_rows(comp_gbm_bay2$tmp_roc_curve) %>%
  mutate(tag = str_replace(tag, 'deletion - bayesian - unbiased approach - ', 'CNVscore - ')) %>%
  mutate(tag = str_replace(tag, 'deletion - gbm - unbiased_approach', 'Gradient boosting model')) %>%
  ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  theme_roc() +
  theme(legend.text=element_text(size=12), legend.title = element_text(size = 15),
        legend.position="top") +
  labs(color = 'Model') +
  # theme(legend.text=element_text(size=12), legend.title = element_text(size = 15)) +
  comp_gbm_bay1$tmp_pr_curve %>%
  bind_rows(comp_gbm_bay2$tmp_pr_curve) %>%
  mutate(tag = str_replace(tag, 'deletion - bayesian - unbiased approach - ', 'CNVscore - ')) %>%
  mutate(tag = str_replace(tag, 'deletion - gbm - unbiased_approach', 'Gradient boosting model')) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  labs(color = 'Model') +
  theme_pr() +
  theme(legend.text=element_text(size=12), legend.title = element_text(size = 15),
        legend.position="top")

p_comp_gbm_bay_decipher <- comp_gbm_bay3$tmp_roc_curve %>%
  bind_rows(comp_gbm_bay4$tmp_roc_curve) %>%
  mutate(tag = str_replace(tag, 'deletion - bayesian - unbiased approach - ', 'CNVscore - ')) %>%
  mutate(tag = str_replace(tag, 'deletion - gbm - unbiased_approach', 'Gradient boosting model')) %>%
ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  # ggtitle('Comparison gradient boosting and Bayesian approach - DECIPHER dataset') +
  theme_roc() +
  labs(color = 'Model') +
  theme(legend.text=element_text(size=12), legend.title = element_text(size = 15),
        legend.position="top") +
  comp_gbm_bay3$tmp_pr_curve %>%
  bind_rows(comp_gbm_bay4$tmp_pr_curve) %>%
  mutate(tag = str_replace(tag, 'deletion - bayesian - unbiased approach - ', 'CNVscore - ')) %>%
  mutate(tag = str_replace(tag, 'deletion - gbm - unbiased_approach', 'Gradient boosting model')) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = tag, color = tag),  show.legend = TRUE, size = 1.5) +
  labs(color = 'Model') +
  theme_pr() +
  theme(legend.text=element_text(size=12), legend.title = element_text(size = 15),
        legend.position="top")
  


p2_comp_gbm_bay_clinvar <- tibble('GBM model' = comp_gbm_bay1$tmp_predicted$.pred_pathogenic, 
       'Bayesian model' = comp_gbm_bay2$tmp_predicted$.pred_pathogenic) %>%
  ggplot(aes(`GBM model`, `Bayesian model`)) +
  geom_bin2d(bins = 80) +
  theme_minimal() +
  ggtitle('ClinVar dataset') +
  tibble('GBM model' = comp_gbm_bay3$tmp_predicted$.pred_pathogenic, 
         'Bayesian model' = comp_gbm_bay4$tmp_predicted$.pred_pathogenic) %>%
  ggplot(aes(`GBM model`, `Bayesian model`)) +
  geom_bin2d(bins = 80) +
  ggtitle('DECIPHER dataset') +
  theme_minimal() +
  theme(legend.position="top")
  

p_comp_gbm_bay_clinvar
p_comp_gbm_bay_decipher

# Figure N. PCA ------------------------------------------------------------------------------


training_pca <- PCA(output_clinvar_deletion[, names(output_clinvar_deletion) %in% 
                                              features_tbl[features_tbl$human_control == 'no',]$name], 
                    graph = FALSE)


p1_pca_training <- fviz_pca_ind(training_pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = output_clinvar_deletion$clinical, # color by groups
             palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)


p2_pca_training <- fviz_eig(training_pca, addlabels = TRUE, ylim = c(0, 50))





# Figure N. GW-PCA ------------------------------------------------------------------------------


chrom_gw <- tibble()
for (i in c(1:22, 'X')) {
  
  chrom_selected <- as.character(i)
  length_chrom_selected <- coord_chrom_hg19 %>% filter(chrom == chrom_selected) %>% pull(length)
  
  # 2e4
  # 1e3
  
  chrom_tmp <- tibble(chrom = chrom_selected, 
                      start = seq(1, length_chrom_selected, 1e3)) %>% 
    mutate(end = start + median(output_clinvar_deletion$length_cnv)  - 1) %>%
    slice(-n()) %>%
    mutate(
      clinical = sample(c('benign', 'pathogenic'), size = nrow(.), replace = TRUE),
      variant_class = 'deletion')
  
  chrom_gw <- chrom_gw %>% bind_rows(chrom_tmp)
}

chrom_gw <- chrom_gw %>% 
  mutate(id_tmp = row_number()) %>%
  mutate(length_cnv = end - start + 1) %>%
  mutate(source = 'chrom_gw')

chrom_gw[1:50,] %>%
  mutate(mid_point = (end + start)/2) %>%
  ggplot(aes(mid_point, id_tmp)) +
  geom_pointrange(aes(xmin = start, xmax = end))


remove_whole_ids <- chrom_gw %>% 
  bed_coverage(problematic_regions) %>% 
  select(id_tmp, .frac) %>% filter(.frac >= 0.3) %>% pull(id_tmp)

chrom_gw <- chrom_gw %>% filter(!id_tmp %in% remove_whole_ids) 

plan('multisession', workers = 10)

annotated_gw_df <- c(1:22, 'X') %>% map_dfr(function(x) {
  
  annotated_gw_df <- check_cnv_v2(chrom_gw %>% filter(chrom == x))
  print(x)
  
  return(annotated_gw_df)
})


predicted_gw_df <- c(1:22, 'X') %>% map_dfr(function(x) {
  
  predicted_gw_df <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, 
                                                annotated_gw_df %>% filter(chrom == x), 
                                                'deletion', 'unbiased approach')[[3]]
  print(x)
  
  return(predicted_gw_df)
})


distances_tbl <- get_distances(
  training_tbl = output_clinvar_deletion,
  input_tbl = annotated_gw_df,
  formule_input = human_no_control, 
  threshold = 0.95
)


target_gw_df <- distances_tbl %>% 
  mutate(ood = if_else(dist_pathogenic >= 90 & dist_benign >= 90 , 'yes', 'no')) %>%
  select(id, ood) %>%
  right_join(annotated_gw_df, by = 'id')


p0_ood <- predicted_gw_df %>% 
  select(id, chrom, sd) %>%
  left_join(target_gw_df %>% select(id, ood), by = 'id') %>%
  mutate(chrom = factor(chrom, levels = c(1:22, 'X'))) %>%
  ggplot(aes(chrom, sd)) + 
  geom_boxplot(aes(fill = ood)) +
  scale_fill_discrete(labels = c('In-of-distribution (IOD)', 'Out-of-distribution (OOD)')) +
  theme_minimal() +
  theme(legend.position="top", plot.title = element_text(size = 11)) +
  labs(x = 'Chromosome', y = 'Standard deviations (sd)', 
       title = 'SD across chromosomes and region type',
       fill = 'Legend')


p1_ood <- target_gw_df %>%
  filter(ood == 'yes') %>%
  select(chrom, start, end) %>%
  bed_merge() %>%
  mutate(length_chunk = end - start + 1) %>%
  group_by(chrom) %>%
  summarise(total_length = sum(length_chunk)) %>%
  left_join(coord_chrom_hg19, by = 'chrom') %>%
  summarise(total_ood_across_chrom = sum(total_length),
            total_across_chrom = sum(length)) %>%
  pivot_longer(everything(), names_to = 'tag', values_to = 'value') %>%
  mutate(perc = value / sum(value)) %>%
  ggplot(aes('', perc)) +
  geom_col(aes(fill = tag), color = 'black') +
  geom_label(aes(label = paste0(round(100*perc, 2), '%'), group = tag), position = position_stack(vjust = 0.5)) +
  scale_fill_discrete(labels = c('In-of-distribution (IOD)', 'Out-of-distribution (OOD)')) +
  scale_y_continuous(label = percent) +
  theme_minimal() +
  theme(legend.position="top", plot.title = element_text(size = 11)) +
  labs(x = 'Human genome', y = 'Coverage', fill = 'Legend',
       title = 'Coverage of OOD regions across the human genome')



p2_ood <- target_gw_df %>%
  filter(ood == 'yes') %>%
  select(chrom, start, end) %>%
  bed_merge() %>%
  mutate(length_chunk = end - start + 1) %>%
  group_by(chrom) %>%
  summarise(total_length = sum(length_chunk)) %>%
  left_join(coord_chrom_hg19, by = 'chrom') %>%
  mutate(coverage = total_length / length) %>%
  mutate(chrom = factor(chrom, levels = c(1:22, 'X'))) %>%
  ggplot(aes(chrom, coverage)) +
  geom_col(color = 'black', fill = 'steelblue') +
  scale_y_continuous(label = percent) +
  theme_minimal() +
  theme(plot.title = element_text(size = 11)) +
  labs(x = 'Chromosome', y = 'Coverage', title = 'Coverage OOD regions across chromosomes')



p3_ood <- target_gw_df %>%
  select(ood, n_genes) %>%
  rename(tag = ood) %>%
  mutate(tag = if_else(tag == 'no', 'IOD', 'OOD')) %>%
  mutate(tag = paste('Genome-wide -', tag)) %>%
  bind_rows(output_clinvar_deletion %>% 
              rename(tag = clinical) %>%
              select(tag, n_genes) %>%
              mutate(tag = paste('Training -', tag))
  ) %>%
  mutate(tag = factor(tag, 
                      levels = c('Training - benign', 'Training - pathogenic', 'Genome-wide - IOD', 'Genome-wide - OOD'))) %>%
  ggplot(aes(tag, n_genes)) +
  geom_boxplot(aes(fill = tag), show.legend = FALSE) +
  theme_lucid() +
  scale_fill_material() +
  scale_y_continuous(breaks= pretty_breaks()) +
  labs(x = 'Category', y = 'Nº genes') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


p4_ood <- target_gw_df %>%
  select(ood, max_cadd) %>%
  rename(tag = ood) %>%
  mutate(tag = if_else(tag == 'no', 'IOD', 'OOD')) %>%
  mutate(tag = paste('Genome-wide -', tag)) %>%
  bind_rows(output_clinvar_deletion %>% 
              rename(tag = clinical) %>%
              select(tag, max_cadd) %>%
              mutate(tag = paste('Training -', tag))
  ) %>%
  mutate(tag = factor(tag, 
                      levels = c('Training - benign', 'Training - pathogenic', 'Genome-wide - IOD', 'Genome-wide - OOD'))) %>%
  ggplot(aes(tag, max_cadd)) +
  geom_boxplot(aes(fill = tag), show.legend = FALSE) +
  theme_lucid() +
  scale_fill_material() +
  scale_y_continuous(breaks= pretty_breaks()) +
  labs(x = 'Category', y = 'Maximum CADD score') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))



p5_ood <- target_gw_df %>%
  select(ood, omim) %>%
  rename(tag = ood) %>%
  mutate(tag = if_else(tag == 'no', 'IOD', 'OOD')) %>%
  mutate(tag = paste('Genome-wide -', tag)) %>%
  bind_rows(output_clinvar_deletion %>% 
              rename(tag = clinical) %>%
              select(tag, omim) %>%
              mutate(tag = paste('Training -', tag))
  ) %>%
  mutate(tag = factor(tag, 
                      levels = c('Training - benign', 'Training - pathogenic', 'Genome-wide - IOD', 'Genome-wide - OOD'))) %>%
  count(tag, omim) %>%
  group_by(tag) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(omim = factor(omim)) %>%
  ggplot(aes(tag, perc)) +
  geom_col(aes(fill = omim), color = 'black', show.legend = TRUE) +
  theme_lucid() +
  scale_fill_material_d() +
  geom_label(aes(label = paste0(round(100*perc,2), '%', '(', n , ')'), 
                 group = omim), size = 4, position = position_stack(vjust = .5)) +
  scale_y_continuous(label = percent) +
  labs(x = 'Category', y = 'omim-associated genes?') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


p6_ood <- target_gw_df %>%
  filter(ood == 'yes') %>%
  select(chrom, start, end) %>%
  bed_merge() %>%
  mutate(length_chunk = end - start + 1) %>%
  group_by(chrom) %>%
  summarise(total_length = sum(length_chunk)) %>%
  left_join(coord_chrom_hg19, by = 'chrom') %>%
  mutate(coverage = total_length / length) %>%
  select(chrom, coverage) %>%
  left_join(hgcn_genes %>% 
              filter(pLI < 90) %>%
              count(chrom) %>%
              left_join(coord_chrom_hg19) %>%
              mutate(gene_density = n / length)) %>%
  ggplot(aes(gene_density, coverage)) +
  geom_label(aes(label = chrom)) +
  # geom_point() +
  theme_lucid() +
  labs(y = 'Coverage out-of-distribution (OOD) regions', x = 'Gene density (Nº genes (pLI < 90) / chromosome length)')

analysis_ood <- p1_ood + p2_ood + p0_ood
analysis_ood2 <- p3_ood + p4_ood + p5_ood
analysis_ood3 <- p6_ood






# Figure N. Uncertainty------------------------------------------------------------------------------


df_for_uncertainty <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_deletion, 'deletion', 'unbiased approach')

sd_df <- df_for_uncertainty[[3]] %>%
  mutate(prediction = if_else(.pred_pathogenic >= 0.5, 'pathogenic', 'benign')) %>%
  mutate(is_correct = if_else(clinical == prediction, 'yes', 'no')) %>%
  mutate(dist_one_zero = if_else(.pred_pathogenic >= 0.5, 1 - .pred_pathogenic, abs(0 - .pred_pathogenic ))) %>%
  mutate(tile_sd = ntile(sd, 10)) %>% 
  mutate(rank_score = dense_rank(dist_one_zero)) %>%
  mutate(rank_sd = dense_rank(sd)) %>%
  mutate(evaluation_priority = rank_score + 20 * rank_sd) %>%
  select(-c('rank_score', 'rank_sd')) %>%
  mutate(tile_evalprior = ntile(evaluation_priority, 10))

spearman_correlation <- sd_df %>% select(dist_one_zero, sd) %>% correlate(method = 'spearman') %>% head() %>% pull(sd) %>% .[1]

p1_uncertainty <- df_for_uncertainty[[3]] %>%
  ggplot(aes(.pred_pathogenic, sd)) +
  geom_point() +
  labs(title = glue('Spearman correlation = {round(spearman_correlation, 2)}')) +
  theme_minimal()

df_for_uncertainty[[3]] %>%
  
  ggplot(aes(.pred_pathogenic, sd)) +
  geom_point() +
  # labs(title = glue('Spearman correlation = {round(spearman_correlation, 2)}')) +
  facet_wrap(vars(clinical)) +
  theme_minimal()

df_for_uncertainty[[3]] %>%
  mutate(tile_sd = ntile())
  ggplot(aes(sd, .pred_pathogenic)) +
  geom_boxplot() +
  # labs(title = glue('Spearman correlation = {round(spearman_correlation, 2)}')) +
  facet_wrap(vars(clinical)) +
  theme_minimal()


p2_uncertainty <- sd_df %>%
  count(tile_sd, clinical) %>%
  group_by(tile_sd) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(tile_sd, perc)) +
  geom_col(aes(fill = clinical), color = 'black') +
  scale_y_continuous(label = percent) +
  theme_minimal()

p3_uncertainty <- sd_df %>%
  ggplot(aes(sd)) +
  geom_density(aes(fill = clinical), alpha = 0.4) +
  theme_minimal()

p4_uncertainty <- sd_df %>%
  count(tile_sd, clinical) %>%
  pivot_wider(id_cols = tile_sd, names_from = clinical, values_from = n) %>%
  mutate(Pathogenic_CNVs = cumsum(pathogenic / sum(pathogenic))) %>%
  mutate(Benign_CNVs = cumsum(benign / sum(benign))) %>%
  select(-c(pathogenic, benign)) %>%
  pivot_longer(-tile_sd, names_to = 'tag', values_to = 'value') %>%
  mutate(tile_sd = as.factor(tile_sd)) %>%
  ggplot(aes(tile_sd, value)) +
  geom_line(aes(group = tag, color = tag), size = 2) +
  geom_point(size = 2) +
  scale_y_continuous(label = percent) +
  scale_color_manual(values = rev(hue_pal()(2)))+
  theme_minimal() +
  labs(x = 'Standard deviation (Percentile 1-10)', y = 'Percentage of CNVs')


p5_uncertainty <- sd_df %>%
  mutate(tile_sd = ntile(sd, 10)) %>%
  count(tile_sd, is_correct) %>%
  # count(clinical, tile_sd, is_correct) %>%
  # group_by(clinical, tile_sd) %>%
  group_by(tile_sd) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(is_correct = factor(is_correct, levels = c('yes', 'no'))) %>%
  ggplot(aes(tile_sd, perc)) +
  geom_col(aes(fill = is_correct), color = 'black') +
  scale_y_continuous(label = percent) +
  geom_label(aes(label = n), position = position_stack(vjust = .5)) +
  scale_fill_manual(values = rev(hue_pal()(2)))+
  ggtitle('% correctly predicted across standard deviation (sd)') +
  theme_minimal() +
  theme(plot.title = element_text(size = 11))



p6_uncertainty <- sd_df %>%
  mutate(dist_one_zero = if_else(.pred_pathogenic >= 0.5, 1 - .pred_pathogenic, abs(0 - .pred_pathogenic ))) %>%
  mutate(dist_one_zero = ntile(dist_one_zero, 10)) %>%
  count(dist_one_zero, is_correct) %>%
  group_by(dist_one_zero) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(is_correct = factor(is_correct, levels = c('yes', 'no'))) %>%
  ggplot(aes(dist_one_zero, perc)) +
  geom_col(aes(fill = is_correct), color = 'black') +
  scale_y_continuous(label = percent) +
  geom_label(aes(label = n), position = position_stack(vjust = .5)) +
  scale_fill_manual(values = rev(hue_pal()(2)))+
  ggtitle('% correctly predicted across pathogenic prob. [1 - score || score - 1]') +
  theme_minimal() +
  theme(plot.title = element_text(size = 11))





p7_uncertainty <- sd_df %>%
  ggplot(aes(sd)) +
  geom_histogram(aes(fill = is_correct), alpha = 0.8, bins = 10, color = 'black') +
  theme_minimal()

p8_uncertainty <- sd_df %>%
  ggplot(aes(sd)) +
  geom_histogram(aes(fill = is_correct), alpha = 0.8, bins = 10, color = 'black') +
  facet_wrap(vars(clinical)) +
  theme_minimal()

tmp_id_decipher_clinical <- read_tsv('/data-cbl/decipher_data/decipher-cnvs-grch37-2020-12-06.txt', skip = 1)

tmp_id_decipher_clinical <- tmp_id_decipher_clinical %>% 
  rename(id = `# patient_id`, clinical2 = pathogenicity) %>%
  select(id, clinical2, contribution) %>%
  filter(clinical2 %in% c('Pathogenic', 'Unknown', 'Likely pathogenic')) %>%
  mutate(id = as.character(id))

p9_uncertainty <- output_decipher_deletion %>% select(id.x, id) %>% 
  left_join(tmp_id_decipher_clinical, by = c('id.x' = 'id')) %>%
  distinct() %>%
  right_join(sd_df, by = 'id') %>%
  mutate(clinical2 = if_else(is.na(clinical2), 'Benign', clinical2)) %>%
  mutate(clinical2 = factor(clinical2, levels = c('Benign', 'Pathogenic', 'Likely pathogenic', 'Unknown'))) %>%
  ggplot(aes(sd, y = clinical2)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = clinical2), alpha = 0.6, show.legend = FALSE, size = 1.25) +
  scale_fill_viridis_d() +
  xlab('Standard deviation (sd)') +
  ylab('Clinical assessment') +
  labs(title = 'SD across across clinical interpretations in DECIPHER') +
  theme_ridges()



p10_uncertainty <- output_decipher_deletion %>% select(id.x, id) %>% 
  left_join(tmp_id_decipher_clinical, by = c('id.x' = 'id')) %>%
  distinct() %>%
  right_join(sd_df, by = 'id') %>%
  filter(!is.na(contribution)) %>%
  mutate(contribution = factor(contribution, levels = c('Full', 'Uncertain', 'Partial', 'Unknown'))) %>%
  ggplot(aes(sd, y = contribution)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, aes(fill = contribution), alpha = 0.6, show.legend = FALSE, size = 1.25) +
  scale_fill_viridis_d() +
  xlab('Standard deviation (sd)') +
  ylab('Clinical contribution') +
  labs(title = 'SD across across clinical contribution in DECIPHER') +
  theme_ridges()


p11_uncertainty <- sd_df %>%
  ggplot(aes(sd)) +
  geom_density(aes(fill = is_correct), alpha = 0.4) +
  facet_wrap(vars(clinical)) +
  theme_minimal()


figure_p1_uncertainty <- (p1_uncertainty + p3_uncertainty + p11_uncertainty) / p8_uncertainty
  
figure_p2_uncertainty <- (p5_uncertainty + p6_uncertainty) / (p9_uncertainty + p10_uncertainty)



# Figure N. applicability  -----------------------------------------------------------------------------


# ClinVar
distances_tbl <- get_distances(
  training_tbl = output_clinvar_deletion,
  input_tbl = output_clinvar_deletion,
  formule_input = human_no_control, 
  threshold = 0.95
)


df_ood <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, 
                                     output_clinvar_deletion , 'deletion', 'human_nocontrol')[[3]]

p1_app_clinvar <- df_ood %>% 
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  mutate(sd = percent_rank(sd)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  ggplot(aes(sd, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(vars(tag)) +
  labs(x = 'Standard deviation (sd)', y = 'Distance') +
  theme_minimal()


p1_1_app <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  mutate(tag = paste(tag, '-', clinical, 'CNVs', '(true label)')) %>%
  group_by(tag) %>%
  mutate(sd = percent_rank(sd)) %>%
  ggplot(aes(sd, value)) +
  geom_point(aes(fill = clinical), shape = 21, color =  'black') +
  geom_smooth(aes(color = clinical), show.legend = F) +
  facet_wrap(vars(tag)) +
  labs(title = 'ClinVar dataset', x = 'Standard deviation (sd)', y = 'Distance', fill = 'True label') +
  theme_minimal()

p1_2_app <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  mutate(tag = paste(tag, '-', clinical, 'CNVs', '(true label)')) %>%
  group_by(tag) %>%
  mutate(sd = ntile(sd, 10)) %>%
  ggplot(aes(as.factor(sd), value)) +
  geom_boxplot(aes(fill = clinical), shape = 21, color =  'black') +
  facet_wrap(vars(tag)) +
  labs(title = 'ClinVar dataset', x = 'Standard deviation (sd)', y = 'Distance', fill = 'True label') +
  theme_minimal()


p1_3_app <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  mutate(tag = paste(tag, '-', clinical, 'CNVs', '(true label)')) %>%
  group_by(tag) %>%
  ggplot(aes(sd, value)) +
  geom_point(aes(fill = clinical), shape = 21, color =  'black') +
  geom_smooth(aes(color = clinical), show.legend = F) +
  coord_cartesian(ylim = c(0,1)) +
  facet_wrap(vars(tag)) +
  labs(title = 'ClinVar dataset', x = 'Standard deviation (sd)', y = 'Distance', fill = 'True label') +
  theme_minimal()


  

p2_app_clinvar <- df_ood %>% 
  left_join(distances_tbl, by = 'id') %>%
  mutate(tag = case_when(
    dist_pathogenic >= 95 & dist_benign >= 95 ~ 'ood',
    dist_pathogenic <= 50 & dist_benign <= 50 ~ 'class_overlap',
    TRUE ~ 'ok'
  )) %>% 
  mutate(tag = factor(tag, levels = c('ok', 'class_overlap', 'ood'))) %>%
  ggplot(aes(clinical, sd)) +
  geom_boxplot(aes(fill = tag)) +
  theme_minimal()


p3_app_clinvar <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(ood = if_else(dist_pathogenic >= 95 & dist_benign >= 95 , 'yes', 'no')) %>%
  mutate(prediction = if_else(.pred_pathogenic >= 0.5, 'pathogenic', 'benign')) %>%
  mutate(is_correct = if_else(clinical == prediction, 'yes', 'no')) %>%
  count(ood, is_correct, clinical) %>%
  group_by(ood, clinical) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(ood, perc)) +
  geom_col(aes(fill = is_correct), color = 'black') +
  scale_y_continuous(label = percent) +
  facet_wrap(vars(clinical)) +
  labs(x = 'OOD', y = 'Percentage') +
  theme_minimal()


# df_ood %>%
#   left_join(distances_tbl, by = 'id') %>%
#   mutate(tile_sd = ntile(sd, 2)) %>%
#   mutate(ood = if_else(dist_pathogenic >= 95 & dist_benign >= 95 , 'yes', 'no')) %>%
#   mutate(prediction = if_else(.pred_pathogenic >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(is_correct = if_else(clinical == prediction, 'yes', 'no')) %>%
#   mutate(sd_ood = paste(ood, tile_sd)) %>%
#   count(sd_ood, is_correct, clinical) %>%
#   group_by(sd_ood, clinical) %>%
#   mutate(perc = n / sum(n)) %>%
#   mutate(sd_ood = factor(sd_ood, levels = c('no 1', 'yes 1', 'no 2', 'yes 2', 'no 3', 'yes 3'))) %>%
#   ggplot(aes(sd_ood, perc)) +
#   geom_col(aes(fill = is_correct), color = 'black') +
#   geom_label(aes(label = paste0(round(100*perc, 2), '%', '(', n, ')'), group = is_correct), position = position_stack(vjust = 0.5)) +
#   scale_y_continuous(label = percent) +
#   facet_wrap(vars(clinical)) +
#   labs(x = 'OOD', y = 'Percentage') +
#   theme_minimal()




# DECIPHER
distances_tbl <- get_distances(
  training_tbl = output_clinvar_deletion,
  input_tbl = output_decipher_deletion,
  formule_input = human_no_control, 
  threshold = 0.95
)


df_ood <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, 
                                     output_decipher_deletion , 'deletion', 'human_nocontrol')[[3]]

p1_app_decipher <- df_ood %>% 
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  mutate(sd = percent_rank(sd)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  ggplot(aes(sd, value)) +
  geom_point() +
  geom_smooth() +
  facet_grid(vars(tag)) +
  labs(x = 'Standard deviation (sd)', y = 'Distance') +
  theme_minimal()



p1_1_app <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  mutate(tag = paste(tag, '-', clinical, 'CNVs', '(true label)')) %>%
  group_by(tag) %>%
  mutate(sd = percent_rank(sd)) %>%
  ggplot(aes(sd, value)) +
  geom_point(aes(fill = clinical), shape = 21, color =  'black') +
  geom_smooth(aes(color = clinical), show.legend = F) +
  facet_wrap(vars(tag)) +
  labs(title = 'DECIPHER dataset', x = 'Standard deviation (sd)', y = 'Distance', fill = 'True label') +
  theme_minimal()

p1_2_app <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  mutate(tag = paste(tag, '-', clinical, 'CNVs', '(true label)')) %>%
  group_by(tag) %>%
  mutate(sd = ntile(sd, 10)) %>%
  ggplot(aes(as.factor(sd), value)) +
  geom_boxplot(aes(fill = clinical), shape = 21, color =  'black') +
  facet_wrap(vars(tag)) +
  labs(title = 'DECIPHER dataset', x = 'Standard deviation (sd)', y = 'Distance', fill = 'True label') +
  theme_minimal()



df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  mutate(tag = paste(tag, '-', clinical, 'CNVs', '(true label)')) %>%
  group_by(tag) %>%
  mutate(sd = cut_interval(sd, 10)) %>%
  ggplot(aes(as.factor(sd), value)) +
  geom_boxplot(aes(fill = clinical), shape = 21, color =  'black') +
  facet_wrap(vars(tag), scales = 'free') +
  labs(title = 'DECIPHER dataset', x = 'Standard deviation (sd)', y = 'Distance', fill = 'True label') +
  theme_minimal()


p1_3_app <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(dist_pathogenic = percent_rank(dist_pathogenic)) %>%
  mutate(dist_benign = percent_rank(dist_benign)) %>%
  select(sd, clinical, dist_benign, dist_pathogenic) %>%
  pivot_longer(-c(sd, clinical), names_to = 'tag', values_to = 'value') %>%
  mutate(tag = paste(tag, '-', clinical, 'CNVs', '(true label)')) %>%
  group_by(tag) %>%
  ggplot(aes(sd, value)) +
  geom_point(aes(fill = clinical), shape = 21, color =  'black') +
  geom_smooth(aes(color = clinical), show.legend = F) +
  coord_cartesian(ylim = c(0,1)) +
  facet_wrap(vars(tag)) +
  labs(title = 'DECIPHER dataset', x = 'Standard deviation (sd)', y = 'Distance', fill = 'True label') +
  theme_minimal()



p2_app_decipher <- df_ood %>% 
  left_join(distances_tbl, by = 'id') %>%
  mutate(tag = case_when(
    dist_pathogenic >= 95 & dist_benign >= 95 ~ 'ood',
    dist_pathogenic <= 50 & dist_benign <= 50 ~ 'class_overlap',
    TRUE ~ 'ok'
  )) %>% 
  mutate(tag = factor(tag, levels = c('ok', 'class_overlap', 'ood'))) %>%
  ggplot(aes(clinical, sd)) +
  geom_boxplot(aes(fill = tag)) +
  theme_minimal()


p3_app_decipher <- df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(ood = if_else(dist_pathogenic >= 95 & dist_benign >= 95 , 'yes', 'no')) %>%
  mutate(prediction = if_else(.pred_pathogenic >= 0.5, 'pathogenic', 'benign')) %>%
  mutate(is_correct = if_else(clinical == prediction, 'yes', 'no')) %>%
  count(ood, is_correct, clinical) %>%
  group_by(ood, clinical) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(ood, perc)) +
  geom_col(aes(fill = is_correct), color = 'black') +
  scale_y_continuous(label = percent) +
  facet_wrap(vars(clinical)) +
  labs(x = 'OOD', y = 'Percentage') +
  theme_minimal()


df_ood %>%
  left_join(distances_tbl, by = 'id') %>%
  mutate(tile_sd = ntile(sd, 2)) %>%
  mutate(ood = if_else(dist_pathogenic >= 95 & dist_benign >= 95 , 'yes', 'no')) %>%
  mutate(prediction = if_else(.pred_pathogenic >= 0.5, 'pathogenic', 'benign')) %>%
  mutate(is_correct = if_else(clinical == prediction, 'yes', 'no')) %>%
  mutate(sd_ood = paste(ood, tile_sd)) %>%
  count(sd_ood, is_correct, clinical) %>%
  group_by(sd_ood, clinical) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(sd_ood = factor(sd_ood, levels = c('no 1', 'yes 1', 'no 2', 'yes 2', 'no 3', 'yes 3'))) %>%
  ggplot(aes(sd_ood, perc)) +
  geom_col(aes(fill = is_correct), color = 'black') +
  geom_label(aes(label = paste0(round(100*perc, 2), '%', '(', n, ')'), group = is_correct), position = position_stack(vjust = 0.5)) +
  scale_y_continuous(label = percent) +
  facet_wrap(vars(clinical)) +
  labs(x = 'OOD', y = 'Percentage') +
  theme_minimal()


# sd_df %>%
#   left_join(distances_tbl %>% select(id, dist_total), by = 'id') %>%
#   mutate(tile_dist_total = ntile(dist_total, 10)) %>%
#   ggplot(aes(factor(tile_dist_total), sd)) +
#     geom_boxplot()
# 
# sd_df %>%
#   left_join(distances_tbl %>% select(id, dist_total), by = 'id') %>%
#   mutate(tile_dist_total = ntile(dist_total, 10)) %>%
#   ggplot(aes(factor(tile_dist_total), dist_one_zero)) +
#   geom_boxplot()


# sd_df %>%
#   left_join(distances_tbl %>% select(id, dist_total), by = 'id') %>%
#   mutate(tile_dist_total = ntile(dist_total, 10)) %>%
#   filter(tile_dist_total %in% c(9,10)) %>%
#   group_by(tile_sd) %>%
#   count(tile_sd, is_correct) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(factor(tile_sd), perc)) +
#   geom_col(aes(fill = is_correct))


# sd_df %>%
#   left_join(distances_tbl %>% select(id, dist_total, dist_pathogenic, dist_benign), by = 'id') %>%
#   mutate(tile_dist_total = ntile(dist_total, 10)) %>%
#   filter(tile_dist_total %in% c(9,10)) %>%
#   ggplot(aes(dist_pathogenic, dist_benign)) +
#     geom_point()

# sd_df %>%
#   left_join(distances_tbl %>% select(id, dist_pathogenic, dist_benign), by = 'id') %>%
#   mutate(ad = case_when(
#     dist_pathogenic >= 50 & dist_benign >= 50 ~ 'ood',
#     dist_pathogenic <= 50 & dist_benign <= 50 ~ 'class_overlap',
#     TRUE ~ 'ok'
#   )) %>%
#   filter(ad == 'ood') %>%
#   mutate(tile_dist_one_zero = ntile(dist_one_zero, 10)) %>%
#   mutate(tile_sd = ntile(sd, 10)) %>%
#   # pivot_longer(cols = starts_with('tile'), names_to = 'measure', values_to = 'value') %>%
#   count(tile_dist_one_zero, is_correct) %>%
#   group_by(tile_dist_one_zero) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(factor(tile_dist_one_zero), perc)) +
#     geom_col(aes(fill = is_correct))
#   



# sd_df %>%
#   # filter(clinical == 'pathogenic') %>%
#   left_join(distances_tbl %>% select(id, dist_total), by = 'id') %>%
#   mutate(tile_dist_total = ntile(dist_total, 10)) %>%
#   # group_by(tile_dist_total) %>%
#   count(tile_dist_total, is_correct, clinical) %>%
#   group_by(tile_dist_total, clinical) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(factor(tile_dist_total), perc)) +
#   geom_col(aes(fill = is_correct)) +
#   facet_wrap(vars(clinical))


# distances_tbl %>%
#   left_join(sd_df, by = 'id') %>%
#   select_if(is.double) %>%
#   correlate()


# distances_tbl %>%
#   left_join(sd_df, by = 'id') %>%
#   ggplot(aes(dist_pathogenic, dist_benign)) +
#     geom_point(aes(fill = tile_evalprior), shape = 21, color = 'black') +
#     scale_fill_viridis_c() +
#     theme_minimal()

# p1_dist_evalprior <- distances_tbl %>%
#   left_join(sd_df, by = 'id') %>%
#   filter(clinical == 'benign') %>%
#   mutate(tile_evalprior = ntile(evaluation_priority, 10)) %>%
#   mutate(tile_dist_benign = ntile(dist_benign, 100)) %>%
#   mutate(tile_dist_pathogenic = ntile(dist_pathogenic, 100)) %>%
#   select(tile_evalprior, tile_dist_pathogenic, tile_dist_benign) %>%
#   pivot_longer(-tile_evalprior, names_to = 'clinical', values_to = 'value') %>%
#   ggplot(aes(as.factor(tile_evalprior), value)) +
#     geom_boxplot(aes(fill = clinical)) +
#     scale_fill_manual(values = rev(hue_pal()(2))) +
#     ggtitle('Benign CNVs') +
#     theme_minimal()
# 
# p2_dist_evalprior <- distances_tbl %>%
#   left_join(sd_df, by = 'id') %>%
#   filter(clinical == 'pathogenic') %>%
#   mutate(tile_evalprior = ntile(evaluation_priority, 10)) %>%
#   mutate(tile_dist_benign = ntile(dist_benign, 100)) %>%
#   mutate(tile_dist_pathogenic = ntile(dist_pathogenic, 100)) %>%
#   select(tile_evalprior, tile_dist_pathogenic, tile_dist_benign) %>%
#   pivot_longer(-tile_evalprior, names_to = 'clinical', values_to = 'value') %>%
#   ggplot(aes(as.factor(tile_evalprior), value)) +
#   geom_boxplot(aes(fill = clinical)) +
#   scale_fill_manual(values = rev(hue_pal()(2))) +
#   ggtitle('Pathogenic CNVs') +
#   theme_minimal()
# 
# p1_dist_evalprior + p2_dist_evalprior
# 
# origin_tbl <- get_distances(
#   training_tbl = output_clinvar_deletion,
#   input_tbl = output_decipher_deletion,
#   formule_input = human_no_control, 
#   threshold = 0.95,
#   only_origin = TRUE
# )
# 
# 
# pca_output <- output_clinvar_deletion %>%
#   select(any_of(all.vars(human_no_control)[-1])) %>%
#   prcomp(scale = TRUE, center = TRUE)
# 
# 
# pca_output$x %>%
#   as_tibble() %>%
#   select(PC1, PC2) %>%
#   bind_cols(output_clinvar_deletion %>% select(clinical)) %>%
#   ggplot(aes(PC1, PC2)) +
#   geom_point(aes(fill = clinical), shape = 21, color = 'black') +
#   theme_minimal()
# 
# origin_tbl %>%
#   filter(clinical == 'pathogenic') %>%
#   ggplot(aes(pc1_origin, pc2_origin)) +
#   geom_point(aes(fill = clinical), shape = 21, color = 'black') +
#   theme_minimal()
#   
# 
# p1_appli_training <- pca_output$x %>%
#   as_tibble() %>%
#   select(PC1, PC2) %>%
#   bind_cols(output_clinvar_deletion %>% select(clinical)) %>%
#   ggplot(aes(PC1, PC2)) +
#   geom_point(aes(fill = clinical), shape = 21, color = 'black') +
#   theme_minimal()
# 
# 
# 
# output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   select(clinical, pc1_pathogenic, pc2_pathogenic) %>%
#   filter(pc2_pathogenic < 25) %>%
#   # bind_cols(pca_score %>% select(distance_pctl) %>% rename(applicability = distance_pctl)) %>%
#   # mutate(applicability = if_else(applicability >= 95, 'No >= 95 Percentile', 'Yes')) %>%
#   ggplot(aes(pc1_pathogenic, pc2_pathogenic)) +
#   geom_point(aes(fill = clinical), shape = 21, color = 'black') +
#   labs(fill = 'Applicability?') +
#   theme_minimal()
# 
# 
# output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   pivot_longer(starts_with('dist')) %>%
#   ggplot(aes(value)) +
#   geom_density(aes(fill = clinical), alpha = 0.4) +
#   facet_wrap(vars(name)) +
#   theme_minimal()
# 
# 
# p1_distance <- output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   mutate(pctl_dist = ntile(dist_pathogenic, 10)) %>%
#   count(clinical, pctl_dist) %>%
#   group_by(pctl_dist) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(pctl_dist, perc)) +
#   geom_col(aes(fill = clinical), color = 'black') +
#   theme_minimal()
# 
# p2_distance <- output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   mutate(pctl_dist = ntile(dist_benign, 10)) %>%
#   count(clinical, pctl_dist) %>%
#   group_by(pctl_dist) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(pctl_dist, perc)) +
#   geom_col(aes(fill = clinical), color = 'black') +
#   theme_minimal()
# 
# p1_distance + p2_distance
# 
# 
# 
# output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   ggplot(aes(dist_pathogenic, dist_benign)) +
#     geom_vline(xintercept = 50) +
#     geom_hline(yintercept = 50) +
#     geom_point(aes(fill = clinical), color = 'black', shape = 21) +
#   theme_minimal()
# 
# output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   mutate(tile_dist_pathogenic = ntile(dist_pathogenic, 100)) %>%
#   mutate(tile_dist_benign = ntile(dist_benign, 100)) %>%
#   mutate(tag = case_when(
#     tile_dist_pathogenic >= 50 & tile_dist_benign >= 50 ~ 'ood',
#     tile_dist_pathogenic >= 50 & tile_dist_benign <= 50 ~ 'benign',
#     tile_dist_pathogenic <= 50 & tile_dist_benign >= 50 ~ 'pathogenic',
#     tile_dist_pathogenic <= 50 & tile_dist_benign <= 50 ~ 'class_overlap'
#   )) %>%
#   count(tag, clinical) %>%
#   group_by(tag) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(tag, perc )) +
#     geom_col(aes(fill = clinical), color = 'black') +
#     geom_label(aes(label = paste0(round(100*perc,2), '%', '(', n , ')'), group = clinical), size = 5, position = position_stack(vjust = .5)) +
#     theme_minimal()
# 
# output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   left_join(sd_df) %>%
#   mutate(tile_dist_pathogenic = ntile(dist_pathogenic, 100)) %>%
#   mutate(tile_dist_benign = ntile(dist_benign, 100)) %>%
#   mutate(tag = case_when(
#     tile_dist_pathogenic >= 50 & tile_dist_benign >= 50 ~ 'ood',
#     tile_dist_pathogenic >= 50 & tile_dist_benign <= 50 ~ 'benign',
#     tile_dist_pathogenic <= 50 & tile_dist_benign >= 50 ~ 'pathogenic',
#     tile_dist_pathogenic <= 50 & tile_dist_benign <= 50 ~ 'class_overlap'
#   )) %>%
#   ggplot(aes(x = mad, y = tag)) +
#   stat_density_ridges(quantile_lines = TRUE, aes(fill = tag), quantiles = 2, alpha = 0.6, show.legend = FALSE, size = 1.25) +
#   scale_fill_viridis_d() +
#   facet_wrap(vars(clinical)) +
#   theme_minimal()
# 
# # me gusta
# output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   left_join(sd_df) %>%
#   mutate(tile_dist_pathogenic = ntile(dist_pathogenic, 100)) %>%
#   mutate(tile_dist_benign = ntile(dist_benign, 100)) %>%
#   mutate(tag = case_when(
#     tile_dist_pathogenic >= 50 & tile_dist_benign >= 50 ~ 'ood',
#     tile_dist_pathogenic >= 50 & tile_dist_benign <= 50 ~ 'benign',
#     tile_dist_pathogenic <= 50 & tile_dist_benign >= 50 ~ 'pathogenic',
#     tile_dist_pathogenic <= 50 & tile_dist_benign <= 50 ~ 'class_overlap'
#   )) %>%
#   ggplot(aes(x = tag, y = sd)) +
#   geom_boxplot(aes(fill = tag)) +
#   scale_fill_viridis_d() +
#   facet_wrap(vars(clinical)) +
#   theme_minimal()
#   
# output_decipher_deletion %>% 
#   select(id, clinical) %>% 
#   left_join(distances_tbl) %>%
#   left_join(sd_df) %>%
#   mutate(tile_dist_pathogenic = ntile(dist_pathogenic, 100)) %>%
#   mutate(tile_dist_benign = ntile(dist_benign, 100)) %>%
#   mutate(tag = case_when(
#     tile_dist_pathogenic >= 50 & tile_dist_benign >= 50 ~ 'ood',
#     tile_dist_pathogenic >= 50 & tile_dist_benign <= 50 ~ 'benign',
#     tile_dist_pathogenic <= 50 & tile_dist_benign >= 50 ~ 'pathogenic',
#     tile_dist_pathogenic <= 50 & tile_dist_benign <= 50 ~ 'class_overlap',
#     TRUE ~ 'na'
#   )) %>%
#   # group_by(clinical) %>%
#   mutate(rank_score = dense_rank(.pred_pathogenic)) %>%
#   mutate(rank_sd = 0.5 * desc(dense_rank(mad))) %>%
#   # ungroup() %>%
#   mutate(ep = if_else(.pred_pathogenic >= 0.5, rank_score + rank_sd, rank_score - rank_sd)) %>%
#   select(clinical, ep, rank_score, tag) %>%
#   pivot_longer(cols = c(ep, rank_score), names_to = 'score_type', values_to = 'value') %>%
#   group_by(score_type, tag) %>%
#   roc_auc(clinical, value) %>%
#   ggplot(aes(tag, .estimate)) +
#     geom_line(aes(group = score_type)) +
#     geom_point(aes(fill = score_type), shape = 21, color = 'black', size = 5) +
#   geom_label(aes(label = round(.estimate, 2), fill = score_type), size = 5) +
#   # geom_label(aes(label = round(.estimate, 2))) +
#   theme_minimal()

# MAYBE DELETE
# Figure N. Interpretability (maybe delete) ------------------------------------------------------------------------------

which_model_chosen <- bayesian_clinvar_del_nohuman[[1]]
which_database_chosen <- output_clinvar_deletion



enriched_rules <- rules_rtemis(which_model_chosen$set_rules, 
                               which_model_chosen$model_trained, 
                               which_database_chosen)



p1_rule <- enriched_rules %>%
  mutate(support = support / nrow(which_database_chosen)) %>%
  mutate(tag = if_else(risk >= 0.5, 'red', 'blue')) %>%
  mutate(high = sample(c('yes', 'no'), size = nrow(.), replace = TRUE, prob = c(0.02, 0.99))) %>%
  mutate(high = if_else(high == 'yes' & support > 0.25 & (risk >= 0.60 | risk <= 0.40), 'yes', 'no')) %>%
  ggplot(aes(risk, support)) +
  geom_point(aes(fill = tag), shape = 21, size = 2) +
  gghighlight(high == 'yes', label_key = description) +
  scale_y_continuous(label = percent) +
  scale_x_continuous(label = percent) +
  geom_vline(aes(xintercept = 0.5), linetype = 'dashed') +
  scale_x_continuous(labels=c("100%","75%","50%","75%","100%")) +
  scale_fill_manual(values = c('green', 'red')) +
  theme_minimal() +
  labs(x = 'Risk (% of training pathogenic CNVs following the rule)', y = 'Support (% of training CNVs following the rule)',
       title = glue('Set of rules (n = {nrow(enriched_rules)}) extracted from the model (chrom1) - Training dataset (Deletions - ClinVar)'))



topn_risk <- enriched_rules %>% filter(support >= 1000) %>% arrange(desc(risk)) %>%
  head(20) %>%
  pull(rule)

topn_norisk <- enriched_rules %>% filter(support >= 1000) %>% arrange(risk) %>%
  head(20) %>%
  pull(rule)

terms_selected <- c(topn_risk, topn_norisk)


order_by_risk <- enriched_rules %>% 
  arrange(desc(risk)) %>% 
  mutate(id_row = row_number()) %>%
  select(rule, id_row)

# enriched_rules %>%
#   filter(rule %in% terms_selected) %>%
#   left_join(order_by_risk, by = 'rule') %>%
#   mutate(risk_inverse = 1 - risk) %>%
#   pivot_longer(cols = starts_with('risk'), values_to = 'value', names_to = 'tag') %>%
#   ggplot(aes(reorder(description, id_row), value)) +
#   geom_col(aes(fill = tag), color = 'black') +
#   coord_flip() +
#   # scale_fill_manual(values = rev(hue_pal()(2)))+
#   theme_minimal()



# enriched_rules %>%
#   filter(support > 100) %>%
#   # filter(rule %in% terms_selected) %>%
#   mutate(tag = if_else(risk > 0.5, 'benign', 'pathogenic')) %>%
#   ggplot(aes(risk, reorder(description, risk))) +
#   geom_col(aes(fill = tag), color = 'black') +
#   scale_fill_manual(values = rev(hue_pal()(2)))+
#   theme_minimal(base_size = 10)

enriched_rules %>%
  ggplot(aes(estimate, risk)) +
  geom_point(shape = 21, color = 'black', aes(fill = support), size = 4) +
  scale_fill_viridis_c() +
  geom_hline(yintercept = 0.5, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(label = percent) +
  scale_x_continuous(limits = c(-0.5,0.5)) +
  labs(x = 'Coefficient (Bayesian model)', y = 'Risk (% pathogenic CNVs)') +
  theme_minimal()


# ggsave(glue("figures/{current_date}/p1_rule.png"),
#        p1_rule, width = 17, height = 9.6, dpi = 300, units = "in", device='png')


# Uncertainty II part (maybe delete) ------------------------------------------------------------------------------


result_toleration_roc_curve <- tibble()
result_toleration_pr_curve <- tibble()


number_chunks <- 10

for (i in 1:number_chunks) {

  ### ROC curves

  tmp_roc_auc <- sd_df %>%
    mutate(tile_split_sd = ntile(sd, number_chunks)) %>%
    # group_by(clinical) %>%
    # mutate(tile_split_sd = ntile(sd, 10)) %>%
    # ungroup() %>%
    filter(tile_split_sd %in% c(1:i)) %>% 
    mutate(clinical = as.factor(clinical)) %>%
    roc_auc(clinical, .pred_pathogenic) %>%
    mutate(tolerated_sd = i) %>%
    pull(.estimate) %>%
    round(3)

  tmp_roc_curve <- sd_df %>%
    # group_by(clinical) %>%
    mutate(tile_split_sd = ntile(sd, number_chunks)) %>%
    # ungroup() %>%
    filter(tile_split_sd %in% c(1:i)) %>%
    mutate(clinical = as.factor(clinical)) %>%
    roc_curve(clinical, .pred_pathogenic) %>%
    mutate(tolerated_sd = glue("{i} - AUC: {tmp_roc_auc}"))

  result_toleration_roc_curve <- bind_rows(result_toleration_roc_curve, tmp_roc_curve)


  ### PR curves

 how_many_benign <-  sd_df %>%
    select(clinical, sd) %>%
    mutate(tile_split_sd = ntile(sd, number_chunks)) %>%
    count(clinical, tile_split_sd) %>%
    pivot_wider(names_from = clinical, values_from = n) %>%
    mutate(cum_pathogenic = cumsum(pathogenic),
           cum_benign = cumsum(benign)) %>%
    mutate(max_n_clinical = case_when(
      cum_pathogenic > cum_benign ~ 'pathogenic',
      TRUE ~ 'benign'
    )) %>%
    select(tile_split_sd, cum_pathogenic, cum_benign, max_n_clinical) %>%
   filter(tile_split_sd %in% c(1:i)) %>%
   slice_tail(n = 1) %>%
   pull(cum_pathogenic)


  tmp_tmp_pr <- sd_df %>%
    mutate(tile_split_sd = ntile(sd, number_chunks)) %>%
    filter(tile_split_sd %in% c(1:i)) %>%
    filter(clinical == 'benign') %>%
    slice_sample(n = how_many_benign)

  tmp_tmp_pr <- sd_df %>%
    mutate(tile_split_sd = ntile(sd, number_chunks)) %>%
    filter(tile_split_sd %in% c(1:i)) %>%
    filter(clinical == 'pathogenic') %>%
    bind_rows(tmp_tmp_pr)

  tmp_pr_auc <- tmp_tmp_pr %>%
    mutate(clinical = as.factor(clinical)) %>%
    pr_auc(clinical, .pred_pathogenic) %>%
    mutate(tolerated_sd = i) %>%
    pull(.estimate) %>%
    round(3)

  tmp_pr_curve <- tmp_tmp_pr %>%
    mutate(clinical = as.factor(clinical)) %>%
    pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tolerated_sd = glue("{i} - AUC: {tmp_pr_auc}"))

  result_toleration_pr_curve <- bind_rows(result_toleration_pr_curve, tmp_pr_curve)

}




p1_toleration_roc_curve <- result_toleration_roc_curve %>%
  mutate(tolerated_sd = as.factor(tolerated_sd)) %>%
  mutate(tolerated_sd = fct_inorder(tolerated_sd)) %>%
  ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(group = tolerated_sd, color = tolerated_sd),  show.legend = TRUE) +
  scale_color_viridis_d() +
  theme_roc()

p2_toleration_pr_curve <- result_toleration_pr_curve %>%
  mutate(tolerated_sd = as.factor(tolerated_sd)) %>%
  mutate(tolerated_sd = fct_inorder(tolerated_sd)) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = tolerated_sd, color = tolerated_sd),  show.legend = TRUE) +
  scale_color_viridis_d() +
  theme_pr()


p1_toleration_roc_curve + p2_toleration_pr_curve


result_toleration_evalprior_roc_curve <- tibble()
result_toleration_evalprior_pr_curve <- tibble()



for (i in 1:10) {
  
  ### ROC curves
  tmp_roc_auc <- sd_df %>%
    mutate(tile_dist_one_zero = ntile(dist_one_zero, 10)) %>%
    mutate(tile_sd = ntile(sd, 10)) %>%
    mutate(new_score = tile_dist_one_zero + tile_sd) %>%
    mutate(new_score = ntile(new_score, 10)) %>%
    # filter(new_score %in% c(1:i)) %>% count(clinical, prediction)
    filter(tile_evalprior %in% c(1:i)) %>% count(clinical, prediction)
    mutate(clinical = as.factor(clinical)) %>% 
    roc_auc(clinical, .pred_pathogenic) %>%
    mutate(tolerated_sd = i) %>%
    pull(.estimate) %>%
    round(3)
  
  tmp_roc_curve <- sd_df %>%
    filter(tile_evalprior %in% c(1:i)) %>%
    mutate(clinical = as.factor(clinical)) %>%
    roc_curve(clinical, .pred_pathogenic) %>%
    mutate(evaluation_score = glue("{i} - AUC: {tmp_roc_auc}"))
  
  result_toleration_evalprior_roc_curve <- bind_rows(result_toleration_evalprior_roc_curve, 
                                                     tmp_roc_curve)
  
  
  ### PR curves
  
  how_many_benign <-  sd_df %>%
    count(clinical, tile_evalprior) %>%
    pivot_wider(names_from = clinical, values_from = n) %>%
    mutate(cum_pathogenic = cumsum(pathogenic),
           cum_benign = cumsum(benign)) %>%
    mutate(max_n_clinical = case_when(
      cum_pathogenic > cum_benign ~ 'pathogenic',
      TRUE ~ 'benign'
    )) %>%
    select(tile_evalprior, cum_pathogenic, cum_benign, max_n_clinical) %>%
    filter(tile_evalprior %in% c(1:i)) %>%
    slice_tail(n = 1) %>%
    pull(cum_pathogenic)

 
  tmp_tmp_pr <- sd_df %>%
    filter(tile_evalprior %in% c(1:i)) %>%
    filter(clinical == 'benign') %>%
    slice_sample(n = how_many_benign)
  # 
  tmp_tmp_pr <- sd_df %>%
    filter(tile_evalprior %in% c(1:i)) %>%
    filter(clinical == 'pathogenic') %>%
    bind_rows(tmp_tmp_pr)
  # 
  tmp_pr_auc <- tmp_tmp_pr %>%
    mutate(clinical = as.factor(clinical)) %>%
    pr_auc(clinical, .pred_pathogenic) %>%
    mutate(tolerated_sd = i) %>%
    pull(.estimate) %>%
    round(3)
  # 
  tmp_pr_curve <- tmp_tmp_pr %>%
    mutate(clinical = as.factor(clinical)) %>%
    pr_curve(clinical, .pred_pathogenic) %>%
    mutate(evaluation_score = glue("{i} - AUC: {tmp_pr_auc}"))
  
  result_toleration_evalprior_pr_curve <- bind_rows(result_toleration_evalprior_pr_curve, tmp_pr_curve)
  

  
}

p1_toleration_eval_prior_roc_curve <- result_toleration_evalprior_roc_curve %>%
  mutate(evaluation_score = as.factor(evaluation_score)) %>%
  mutate(evaluation_score = fct_inorder(evaluation_score)) %>%
  ggplot(aes(1 - specificity, sensitivity)) +
  geom_path(aes(group = evaluation_score, color = evaluation_score),  show.legend = TRUE) +
  scale_color_viridis_d() +
  theme_roc() +
  ggtitle('ROC curve - Evaluation score (Beta = 20)')

p2_toleration_eval_prior_pr_curve <- result_toleration_evalprior_pr_curve %>%
  mutate(evaluation_score = as.factor(evaluation_score)) %>%
  mutate(evaluation_score = fct_inorder(evaluation_score)) %>%
  ggplot(aes(recall, precision)) +
  geom_path(aes(group = evaluation_score, color = evaluation_score),  show.legend = TRUE) +
  scale_color_viridis_d() +
  theme_pr() +
  ggtitle('PR curve - Evaluation score (Beta = 20)')

p1_toleration_eval_prior_roc_curve + p2_toleration_eval_prior_pr_curve
# (p1_toleration_eval_prior_roc_curve + p2_toleration_eval_prior_pr_curve) / (p1_toleration_roc_curve + p2_toleration_pr_curve)


### Test


a <- expand.grid(y = seq(0,30, 0.1))
result_expan <- tibble()

for (i in 1:nrow(a)) {

tmp_auc <- sd_df %>%
  mutate(rank_score = dense_rank(.pred_pathogenic)) %>%
  mutate(rank_sd = dense_rank(sd)) %>%
mutate(evaluation_priority = if_else(.pred_pathogenic >= 0.5, rank_score - a$y[i] * rank_sd,
                                     rank_score + a$y[i] * rank_sd)) %>%
  roc_auc(clinical, evaluation_priority) %>%
    pull(.estimate) %>%
    round(2)

result_expan <- result_expan %>% bind_rows(tibble(y = a$y[i], result_auc = tmp_auc))

}


result_expan %>% ggplot(aes(y, result_auc)) + geom_point()


sd_df %>%
  ggplot(aes(.pred_pathogenic)) +
  geom_density(aes(fill = clinical), alpha = 0.4)


sd_df %>%
  mutate(rank_score = dense_rank(.pred_pathogenic)) %>%
  mutate(rank_sd = dense_rank(sd)) %>%
  mutate(evaluation_priority = ifelse(.pred_pathogenic >= 0.5,  0.5 * rank_score - 0.5 *rank_sd,
                                       rank_score)) %>%
  roc_auc(clinical, evaluation_priority)


sd_df %>%
  mutate(tile_dist_zero = ntile(dist_one_zero, 100)) %>%
  mutate(tile_sd = ntile(sd, 100)) %>%
  mutate(final_score = tile_dist_zero + tile_sd) %>%
  mutate(tile_final_score = ntile(final_score, 10)) %>%
  roc_auc(clinical, tile_sd)


try_betas <- function(x) {
  
  
  p1_tmp_patho <- sd_df %>% 
    filter(clinical == 'pathogenic') %>%
    mutate(rank_score = desc(dense_rank(.pred_pathogenic))) %>%
    mutate(rank_sd = dense_rank(sd)) %>%
    mutate(evaluation_priority = rank_score + x * rank_sd) %>% 
    mutate(tile_evalprior = ntile(evaluation_priority, 5)) %>% 
    count(is_correct, tile_evalprior) %>%
    group_by(tile_evalprior) %>%
    mutate(perc = n / sum(n)) %>%
    ggplot(aes(tile_evalprior, perc)) +
    geom_col(aes(fill = is_correct), color = 'black') +
    geom_label(aes(label = paste0(round(100*perc,2), '%', '(', n , ')'), group = is_correct), size = 3, position = position_stack(vjust = .5)) +
    scale_y_continuous(label = percent) +
    ggtitle(glue('Beta =  {x}')) +
    theme_minimal()
  
  
  p2_tmp_benign <- sd_df %>% 
    filter(clinical == 'benign') %>%
    mutate(rank_score = dense_rank(dist_one_zero)) %>%
    mutate(rank_sd = dense_rank(sd)) %>%
    mutate(evaluation_priority = rank_score + x * rank_sd) %>%
    mutate(tile_evalprior = ntile(evaluation_priority, 5)) %>%
    count(clinical, is_correct, tile_evalprior) %>%
    group_by(clinical, tile_evalprior) %>%
    mutate(perc = n / sum(n)) %>%
    ggplot(aes(tile_evalprior, perc)) +
    geom_col(aes(fill = is_correct), color = 'black') +
    # facet_wrap(vars(clinical)) +
    geom_label(aes(label = paste0(round(100*perc,2), '%', '(', n , ')'), group = is_correct), size = 3, position = position_stack(vjust = .5)) +
    scale_y_continuous(label = percent) +
    ggtitle(glue('Beta =  {x}')) +
    theme_minimal()
  
  
  p1_tmp_patho + p2_tmp_benign
  
  
}

p0_beta <- try_betas(0)
p1_beta <- try_betas(1)
p2_beta <- try_betas(10)
p3_beta <- try_betas(20)
p4_beta <- try_betas(30)
p5_beta <- try_betas(50)
p6_beta <- try_betas(100)


p0_beta + p1_beta + p2_beta + p3_beta + p4_beta + p5_beta  + plot_layout(nrow = 6, byrow = TRUE)







# Tsne (maybe delete) ------------------------------------------------------------------------------
# library(Rtsne)
# 
# tmp_tsne2 <- output_clinvar_deletion %>% 
#   select(-c('start', 'end', 'length_cnv', 'id')) %>% 
#   select_if(function(col) is.double(col) | all(col == .$clinical)) %>% 
#   distinct() %>%
#   mutate(tag = 'training')
# 
# tmp_tsne <- output_decipher_deletion %>% 
#   left_join(sd_df %>% select(evaluation_priority, tile_evalprior, sd, id), by = c('id')) %>%
#   select(-c('start', 'end', 'length_cnv', 'id')) %>% 
#   select_if(function(col) is.numeric(col) | all(col == .$clinical)) %>% 
#   distinct() %>%
#   mutate(tag = 'test')
# 
# 
# tmp_tsne_def <- tmp_tsne %>% bind_rows(tmp_tsne2)
# 
# 
# # está mal faltan 2 observaciones
# tsne_result <- Rtsne(as.matrix(tmp_tsne_def %>% 
#                                  select(!contains('NumberSubmitter')) %>%
#                                  select(-c('clinical', 'tag', 'sd', 'tile_evalprior', 'evaluation_priority')) %>% distinct()),
#                     perplexity = 40)
# 
# 
# tsne_result_plot <-  tibble(x = tsne_result$Y[,1], 
#                             y = tsne_result$Y[,2],
#                             tag = as.factor(tmp_tsne_def$tag),
#                             sd = tmp_tsne_def$sd,
#                             evaluation_priority = tmp_tsne_def$evaluation_priority,
#                             tile_evalprior = tmp_tsne_def$tile_evalprior,
#                             clinical = as.factor(tmp_tsne_def$clinical)) %>%
#   mutate(tile_sd = ntile(sd, 5)) %>%
#   mutate(tile_evalprior = ntile(tile_evalprior, 5))
# 
# 
# # centroid_pathogenic <- tsne_result_plot %>% 
# #   filter(clinical == 'pathogenic') %>%
# #   summarise(x = mean(x), y = mean(y)) %>% 
# #   mutate(clinical = 'pathogenic')
# # 
# # centroid_benign <- tsne_result_plot %>% 
# #   filter(clinical == 'benign') %>%
# #   summarise(x = mean(x), y = mean(y)) %>% 
# #   mutate(clinical = 'benign')
# 
# 
# # centroid_patho_kmeans <- kmeans(tsne_result_plot %>% filter(clinical == 'pathogenic') %>% select(x, y), 1, iter.max = 10, nstart = 1)
# # centroid_benign_kmeans <- kmeans(tsne_result_plot %>% filter(clinical == 'benign') %>% select(x, y), 1, iter.max = 10, nstart = 1)
# 
# centroid_kmeans <- kmeans(tsne_result_plot %>% 
#                             filter(tag == 'training') %>% 
#                             select(x, y), 2, iter.max = 10, nstart = 1)$centers %>% 
#   as_tibble() %>% mutate(clinical = c('benign', 'pathogenic'))
# 
# tsne_result_plot %>% 
#   filter(tag == 'test') %>%
#   mutate(dist_patho = sqrt((x - centroid_kmeans$x[2])^2 + (y - centroid_kmeans$y[2])^2)) %>%
#   mutate(dist_benign = sqrt((x - centroid_kmeans$x[1])^2 + (y - centroid_kmeans$y[1])^2)) %>%
#   mutate(closest_dist = if_else(dist_patho > dist_benign, dist_benign, dist_patho)) %>%
#   select(closest_dist, evaluation_priority) %>%
#   correlate()
# 
# tsne_result_plot %>% 
#   filter(tag == 'test') %>%
#   mutate(dist_patho = sqrt((x - centroid_kmeans$x[2])^2 + (y - centroid_kmeans$y[2])^2)) %>%
#   mutate(dist_benign = sqrt((x - centroid_kmeans$x[1])^2 + (y - centroid_kmeans$y[1])^2)) %>%
#   mutate(closest_dist = if_else(dist_patho > dist_benign, dist_benign, dist_patho)) %>%
#   select(closest_dist, evaluation_priority) %>%
#   ggplot(aes(closest_dist, evaluation_priority)) +
#     geom_point() +
#     geom_smooth()
# 
# tsne_result_plot %>% 
#   mutate(dist_patho = sqrt((x - centroid_kmeans$x[2])^2 + (y - centroid_kmeans$y[2])^2)) %>%
#   mutate(dist_benign = sqrt((x - centroid_kmeans$x[1])^2 + (y - centroid_kmeans$y[1])^2)) %>%
#   ggplot(aes(dist_patho, dist_benign)) +
#     geom_point(aes(color = clinical)) +
#   # geom_hline(yintercept = 0) +
#   # geom_vline(xintercept = 0) +
#   theme_minimal()
# 
# 
# tsne_result_plot %>% 
#   mutate(dist_patho = sqrt((x - centroid_kmeans$x[2])^2 + (y - centroid_kmeans$y[2])^2)) %>%
#   mutate(dist_benign = sqrt((x - centroid_kmeans$x[1])^2 + (y - centroid_kmeans$y[1])^2)) %>%
#   ggplot(aes(dist_patho, dist_benign)) +
#   geom_point(aes(fill = tile_evalprior), shape = 21) +
#   scale_fill_viridis_c()
# 
# 
# tsne_result_plot %>% 
#   mutate(dist_patho = sqrt((x - centroid_kmeans$x[2])^2 + (y - centroid_kmeans$y[2])^2)) %>%
#   mutate(dist_benign = sqrt((x - centroid_kmeans$x[1])^2 + (y - centroid_kmeans$y[1])^2)) %>%
#   mutate(dist_patho = ntile(dist_patho, 3)) %>%
#   mutate(dist_benign = ntile(dist_benign, 3)) %>%
#   mutate(dist_patho = case_when(
#     dist_patho == 1 ~ 'close_pathogenic',
#     dist_patho == 2 ~ 'mid_pathogenic',
#     dist_patho == 3 ~ 'far_pathogenic'
#   )) %>%
#   mutate(dist_benign = case_when(
#     dist_benign == 1 ~ 'close_benign',
#     dist_benign == 2 ~ 'mid_benign',
#     dist_benign == 3 ~ 'far_benign'
#   )) %>%
#   mutate(summary_dist = paste0(dist_patho, '-', dist_benign)) %>%
#   count(tile_evalprior, summary_dist) %>% 
#   group_by(tile_evalprior) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(tile_evalprior, perc)) +
#     geom_col(aes(fill = summary_dist), color = 'black') +
#   geom_label(aes(label = paste0(round(100*perc,2), '%', '(', n , ')'), group = summary_dist), size = 3, position = position_stack(vjust = .5)) +
#   theme_minimal()
#   
# 
#   
#  tsne_result_plot %>%
#   filter(tag == 'training') %>%
#   ggplot(aes(x, y)) +
#     geom_point(aes(fill = clinical), shape = 21, color = 'black', size = 3) +
#   # scale_x_continuous(limits = c(-60, 60)) +
#   # scale_y_continuous(limits = c(-60, 60)) +
#   stat_ellipse(aes(color = clinical)) +
#   # geom_point(data = centroid_total, aes(x, y, fill = clinical), shape = 21, size = 8) +
#   geom_point(data = centroid_test, aes(x, y, fill = clinical), shape = 21, size = 8) +
#   theme_minimal()
# 
# 
# # p2_eval_priority <- tsne_result_plot %>%
# #   filter(tag == 'test') %>%
# #   ggplot(aes(x, y)) +
# #   geom_point(aes(fill = tile_sd), shape = 21, color = 'black') +
# #   scale_fill_distiller(palette = "Spectral") + 
# #   # scale_y_continuous(limits = c(-60, 60)) +
# #   # scale_x_continuous(limits = c(-60, 60)) +
# #   theme_minimal()
# 
# p3_eval_priority <- tsne_result_plot %>%
#   filter(tag == 'test') %>%
#   ggplot(aes(x, y)) +
#   geom_point(aes(fill = tile_evalprior), shape = 21, color = 'black') +
#   scale_fill_distiller(palette = "Spectral") + 
#   # scale_y_continuous(limits = c(-60, 60)) +
#   # scale_x_continuous(limits = c(-60, 60)) +
#   theme_minimal()
# 
# p1_eval_priority + p3_eval_priority








# Non-coding (maybe delete) ------------------------------------------------------------------------------


p1_mode_reg <- get_mode_noncoding(decipher_match_deletion, 
                                  bayesian_clinvar_del_nohuman, input_tag = 'Deletion - ClinVar')

# ------------------------------------------------------------------------------
# GW - MAIN ANALYSIS
# ------------------------------------------------------------------------------



result_gw_bayesian_del_human <- predict_chrom_aware_rtemis(bayesian_total_del_human,
                                                           output_gw_df,
                                                           'deletion', 'human_control', only_table = TRUE)

result_gw_bayesian_del_nohuman <- predict_chrom_aware_rtemis(bayesian_total_del_nohuman, 
                                                             output_gw_df, 
                                                             'deletion', 'human_no_control', only_table = TRUE)



n_genes_df <- output_gw_df %>% select(id, chrom, start, end) %>%
  bed_intersect(hgcn_genes) %>%
  count(chrom, start.x, end.x) %>%
  rename(start = start.x, end = end.x, n_genes = n) %>%
  select(chrom, start, end, n_genes)


main_gw <-  result_gw_bayesian_del_nohuman %>%
  rename(pred_nohuman = .pred_pathogenic) %>%
  left_join(output_gw_df %>% select(start, end, id, chrom), by = c('chrom', 'id')) %>%
  left_join(result_gw_bayesian_del_human %>% 
              rename(pred_human = .pred_pathogenic) %>% select(chrom, id, sd, pred_human), by = c('id', 'chrom')) %>%
  left_join(n_genes_df, by = c('chrom', 'start', 'end')) %>%
  mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes)) %>%
  pivot_longer(cols = starts_with('pred'), names_to = 'category', values_to = 'prediction') %>%
  mutate(mid_point = (start + end) / 2) %>%
  mutate(chrom = factor(chrom, levels = as.character(c(seq(1,22), 'X')))) %>%
  rename('sd_nohuman' = sd.x, 'sd_human' = sd.y) %>%
  filter(sd_human < 0.1 & sd_nohuman < 0.1) # 50% FILTERED OUT


# Exploratory analysis ------------------------------------------------------------------------------

p1_gw <- main_gw %>%
  select(chrom, id, category, prediction) %>%
  pivot_wider(id_cols = c('id', 'chrom'), names_from = 'category', values_from = 'prediction') %>%
  # filter(chrom == '22') %>%
  mutate(category = case_when(
    pred_nohuman > 0.70 & pred_human < 0.30 ~ 'pyu',
    pred_human > 0.70 & pred_nohuman < 0.30 ~ 'nf',
    pred_human > 0.70 & pred_nohuman > 0.70 ~ 'ok_up',
    pred_human < 0.30 & pred_nohuman < 0.30 ~ 'ok_down',
    TRUE ~ 'not_useful'
  )) %>%
  ggplot(aes(pred_human, pred_nohuman)) +
  geom_point(aes(fill = category ), shape = 21, color = 'black', size = 3) +
  geom_hline(yintercept = c(0.3, 0.7), linetype="dashed") +
  geom_vline(xintercept = c(0.3, 0.7), linetype="dashed") +
  theme_minimal() +
  ggtitle('Example - chrom 22')

p2_gw <- main_gw %>%
  select(chrom, id, category, prediction) %>%
  pivot_wider(id_cols = c('id', 'chrom'), names_from = 'category', values_from = 'prediction') %>%
  # filter(chrom == '22') %>%
  mutate(category = case_when(
    pred_nohuman > 0.70 & pred_human < 0.30 ~ 'pyu',
    pred_human > 0.70 & pred_nohuman < 0.30 ~ 'nf',
    pred_human > 0.70 & pred_nohuman > 0.70 ~ 'ok_up',
    pred_human < 0.30 & pred_nohuman < 0.30 ~ 'ok_down',
    TRUE ~ 'not_useful'
  )) %>%
  count(category) %>%
  mutate(perc = round(n / sum(n), 3)*100) %>%
  ggplot(aes(reorder(category,-perc), perc)) +
  geom_col(aes(fill = category), color = 'black') +
  geom_label(aes(label = paste0(perc, '%')), position = position_stack(vjust = .5)) +
  theme_minimal()


pubmed_gw <- pubmed_df %>% 
  # filter(chrom == '22') %>%
  mutate(mid_point = (start + end) / 2) %>%
  mutate(hits_del = scale(hits_del))

p3_gw <- main_gw %>%
  pivot_wider(id_cols = -category, names_from = category, values_from = prediction) %>%
  group_by(chrom) %>%
  mutate(pred_human = scale(pred_human)) %>%
  mutate(pred_nohuman = scale(pred_nohuman)) %>%
  mutate(n_genes = scale(n_genes)) %>%
  pivot_longer(cols = c('n_genes', 'pred_nohuman', 'pred_human'), names_to = 'category', values_to = 'prediction') %>%
  filter(category != 'n_genes') %>%
  ggplot(aes(mid_point, prediction)) +
  geom_smooth(aes(group = category, color = category), method = 'gam') +
  # geom_line(data = pubmed_gw, aes(mid_point, hits_del), size = 2, alpha = 0.6) +
  facet_wrap(vars(chrom), scales = 'free') +
  theme_minimal()

(p1_gw + p2_gw) / p3_gw


# ------------------------------------------------------------------------------
# GW - PFAM 
# ------------------------------------------------------------------------------

download.file('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ucscGenePfam.txt.gz',
              destfile = 'ucscGenePfam.txt.gz')

pfam_hg19 <- read_tsv('ucscGenePfam.txt.gz', col_names = FALSE) %>%
  rename(chrom = X2, start = X3, end = X4, domain = X5) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  select(chrom, start, end, domain)

file.remove('ucscGenePfam.txt.gz')

main2_gw <- main_gw %>% 
  pivot_wider(id_cols = c('id', 'chrom', 'start', 'end'), names_from = 'category', values_from = 'prediction') %>%
  mutate(category = case_when(
    pred_nohuman > 0.70 & pred_human < 0.30 ~ 'pyu',
    pred_human > 0.70 & pred_nohuman < 0.30 ~ 'nf',
    pred_human > 0.70 & pred_nohuman > 0.70 ~ 'ok_up',
    pred_human < 0.30 & pred_nohuman < 0.30 ~ 'ok_down',
    TRUE ~ 'not_useful')) %>%
  mutate(category2 = if_else(category %in% c('ok_up', 'pyu'), 'pathogenic_nohuman', 'benign_nohuman'))

run_odds_pfam <- function(x, tag) {
  
  vector_domains <- pfam_hg19 %>% select(domain) %>% distinct() %>% pull()
  
  # x <-  main_gw %>%
  #   filter(category == 'pyu')
  
  x_genes <- x %>%
    bed_intersect(hgcn_genes %>% select(chrom, start, end, gene) %>%
                    mutate(length_gene = end - start + 1)) %>%
    mutate(perc_overlap = .overlap / length_gene.y) %>%
    filter(perc_overlap > 0.5) %>%
    select(gene.y, chrom, start.y, end.y) %>%
    rename(gene = gene.y, start = start.y, end = end.y) %>%
    distinct()
  
  for (i in 1:length(vector_domains)) {
    
    print(i)
    y_genes <- pfam_hg19 %>% 
      filter(domain == vector_domains[i]) %>%
      bed_intersect(hgcn_genes) %>%
      pull(gene.y) %>%
      unique()
    
    
    
    first_value <- sum(x_genes$gene %in% y_genes)
    second_value <- length(y_genes) - first_value
    third_value <- nrow(x_genes) - first_value
    fourth_value <- hgcn_genes %>% filter(! gene %in% x_genes$gene) %>% filter(! gene %in% y_genes) %>% nrow()
    
    
    result_tmp <- fisher.test( matrix(c(first_value,second_value,third_value,fourth_value), nrow = 2, byrow = TRUE)) %>% 
      glance() %>%
      select(-method, -alternative) %>%
      mutate(domain = vector_domains[i]) %>%
      mutate(tag = tag)
    
    if (i == 1) {
      result <- result_tmp
    } else {
      result <- result %>% bind_rows(result_tmp)
    }
    
  }
  
  return(result)
}


pfam_pathogenic  <-  run_odds_pfam(main2_gw %>% filter(category2 == 'pathogenic_nohuman'), tag = 'pathogenic_nohuman')
pfam_benign  <-  run_odds_pfam(main2_gw %>% filter(category2 != 'pathogenic_nohuman'), tag = 'benign_nohuman')

fdr_pvalues_pathogenic <- p.adjust(pfam_pathogenic$p.value, method = 'BH', n = nrow(pfam_pathogenic))
fdr_pvalues_benign <- p.adjust(pfam_benign$p.value, method = 'BH', n = nrow(pfam_benign))

pfam_pathogenic <- pfam_pathogenic %>% mutate(fdr_pvalues = fdr_pvalues_pathogenic)
pfam_benign <- pfam_benign %>% mutate(fdr_pvalues = fdr_pvalues_benign)

pfam_pathogenic %>% 
  bind_rows(pfam_benign) %>%
  filter(fdr_pvalues < 0.05) %>%
  mutate(across(where(is.double), log10)) %>% 
  ggplot(aes(estimate, reorder(domain, estimate))) +
  # geom_point(aes(fill = sig), shape = 21, color = 'black') +
  geom_pointrange(aes(xmin=conf.low, xmax=conf.high), size = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_minimal(base_size = 10) +
  labs(title = 'Enrichment of Pfam domain families (BH correction - p.value < 0.05)', y = 'Pfam domains',
       x = 'log10(odds ratio)', color = 'Significant after Bonferroni correction') +
  facet_wrap(vars(tag), scales = 'free')


# ------------------------------------------------------------------------------
# GW - MERGED INTERVALS 
# ------------------------------------------------------------------------------

bed_merge( main2_gw %>% filter(category == 'pyu')) %>% mutate(length_cnv = end - start + 1) %>% summarise(total = sum(length_cnv)/ 1e6)

# ------------------------------------------------------------------------------
# GW - GENE CATEGORIES 
# ------------------------------------------------------------------------------


# genes_pyu_entrez <- main_gw %>%
#   filter(category == 'pyu') %>%
#   bed_intersect(hgcn_genes %>% select(chrom, start, end, entrez_id) %>%
#                   mutate(length_gene = end - start + 1)) %>%
#   mutate(perc_overlap = .overlap / length_gene.y) %>%
#   filter(perc_overlap > 0.9) %>%
#   pull(entrez_id.y) %>%
#   unique()

# genes_ok_down <- main2_gw %>%
#   filter(category == 'ok_down') %>%
#   bed_intersect(hgcn_genes %>% select(chrom, start, end, gene) %>%
#                   mutate(length_gene = end - start + 1)) %>%
#   mutate(perc_overlap = .overlap / length_gene.y) %>%
#   filter(perc_overlap > 0.9) %>%
#   pull(gene.y) %>%
#   unique()
# 
# genes_ok_up <- main2_gw %>%
#   filter(category == 'ok_up') %>%
#   bed_intersect(hgcn_genes %>% select(chrom, start, end, gene) %>%
#                   mutate(length_gene = end - start + 1)) %>%
#   mutate(perc_overlap = .overlap / length_gene.y) %>%
#   filter(perc_overlap > 0.9) %>%
#   pull(gene.y) %>%
#   unique()
# 
# genes_nf <- main2_gw %>%
#   filter(category == 'nf') %>%
#   bed_intersect(hgcn_genes %>% select(chrom, start, end, gene) %>%
#                   mutate(length_gene = end - start + 1)) %>%
#   mutate(perc_overlap = .overlap / length_gene.y) %>%
#   filter(perc_overlap > 0.9) %>%
#   pull(gene.y) %>%
#   unique()

run_total_odds <- function(x, tag) {
  
  tmp_df <- bind_rows(
    run_odds(x, samochadenovo, tag = 'denovo', tag2 = tag),
    run_odds(x, dl_no_disease, tag = 'dl_163', tag2 = tag),
    run_odds(x, proteome_placenta$gene, tag = 'placenta', tag2 = tag),
    run_odds(x, male_infertility$gene, tag = 'male_infertility', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(clinvar == 'Yes') %>% pull(gene), tag = 'clinvar', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(gwas == 'Yes') %>% pull(gene), tag = 'GWAS', tag2 = tag),
    run_odds(x, para_genes$gene, tag = 'paral', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(disease == 'Yes') %>% pull(gene), tag = 'disease', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(omim == 'Yes') %>% pull(gene), tag = 'omim', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(fusil == 'CL (Cellular lethal)') %>% pull(gene), tag ='cellular_lethal', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(fusil == 'DL (Developmental lethal)') %>% pull(gene), tag = 'devel_lethal', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(fusil == 'SV (Subviable)') %>% pull(gene), tag = 'subviable', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(fusil == 'VP (Viable with phenotype)') %>% pull(gene), tag = 'viable_with_pheno', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(fusil == 'VN (Viable with no phenotype)') %>% pull(gene), tag = 'viable_no_pheno', tag2 = tag),
    run_odds(x, hgcn_genes %>% filter(hi >= 90) %>% pull(gene), tag = 'high_HI', tag2 = tag)
  )
  
  fdr_pvalues_tmp <- p.adjust(tmp_df$p.value, method = 'BH', n = nrow(tmp_df))
  tmp_df <- tmp_df %>% mutate(fdr_pvalues = fdr_pvalues_tmp)
  return(tmp_df)
}


genes_pathogenic_nohuman <- main2_gw %>%
  filter(category2 == 'pathogenic_nohuman') %>%
  bed_intersect(hgcn_genes %>% select(chrom, start, end, gene) %>%
                  mutate(length_gene = end - start + 1)) %>%
  mutate(perc_overlap = .overlap / length_gene.y) %>%
  filter(perc_overlap > 0.5) %>%
  pull(gene.y) %>%
  unique()


genes_benign_nohuman <- main2_gw %>%
  filter(category2 == 'benign_nohuman') %>%
  bed_intersect(hgcn_genes %>% select(chrom, start, end, gene) %>%
                  mutate(length_gene = end - start + 1)) %>%
  mutate(perc_overlap = .overlap / length_gene.y) %>%
  filter(perc_overlap > 0.5) %>%
  pull(gene.y) %>%
  unique()


bind_rows(
  run_total_odds(genes_pathogenic_nohuman, 'pathogenic_nohuman'),
  run_total_odds(genes_benign_nohuman, 'benign_nohuman')
) %>%
  mutate(pvalue_lower_than_0.05 = if_else(fdr_pvalues < 0.05, 'yes', 'no')) %>%
  # filter(fdr_pvalues < 0.05) %>%
  mutate(across(where(is.double), log10)) %>% 
  ggplot(aes(estimate, tag)) +
  geom_pointrange(aes(estimate, xmin = conf.low, xmax = conf.high,
                      color = pvalue_lower_than_0.05)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_minimal() +
  theme_minimal(base_size = 18) +
  facet_grid(cols = vars(tag2), scales = 'free') +
  labs(title = 'Enrichment of gene categories (BH correction - p.value < 0.05)', y = 'Gene categories',
       x = 'log10(odds ratio)', color = 'Significant ')

# ------------------------------------------------------------------------------
# GW - ENRICHMENT ANALYSES
# ------------------------------------------------------------------------------




enrich_categories <- function(x) {
  
  # x <- genes_pathogenic_nohuman
  list_databases <- c('GO_Biological_Process_2018', 'Reactome_2016')
  
  result_enrichr <- enrichr(x, list_databases)
  result_enrichr <- result_enrichr %>% map_dfr(~ .x)
  result_enrichr %>%
    as_tibble() %>% 
    filter(Adjusted.P.value < 0.05) %>%
    mutate(Adjusted.P.value = -log10(Adjusted.P.value)) %>%
    arrange(desc(Adjusted.P.value)) %>%
    slice_head(n = 10) %>%
    ggplot(aes(reorder(Term, Adjusted.P.value), Adjusted.P.value)) +
    geom_col(aes(fill = Adjusted.P.value), color = 'black') +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(y = '-log10(p.value.adjusted)', x = 'GO-BP & Reactome') + 
    coord_flip()
}

p1_enrich <- enrich_categories(genes_pathogenic_nohuman) + labs(title = 'Pathogenic regions (nohuman model)')
p2_enrich <- enrich_categories(genes_benign_nohuman) + labs(title = 'Benign regions (nohuman model)')

p1_enrich + p2_enrich + plot_layout(nrow = 2)

# ------------------------------------------------------------------------------
# GW - PUBMED HITS
# ------------------------------------------------------------------------------
library(ggpubr)



pubmed_main2_gw <- main2_gw %>%
  bed_intersect(pubmed_df, suffix = c('', '.y')) %>%
  rename(hits_del = hits_del.y, hits_dup = hits_dup.y) %>%
  filter(0.2 < (.overlap / 1e5)) %>%
  select(id, category2, hits_del, hits_dup) %>%
  distinct()



p1_pubmed <- pubmed_main2_gw %>%
  ggplot(aes(category2, hits_del)) +
  geom_boxplot(aes(fill = category2)) +
  theme_minimal() +
  stat_compare_means()

p2_pubmed <- pubmed_main2_gw %>%
  ggplot(aes(hits_del)) +
  geom_density(aes(fill = category2), alpha = 0.6)+
  theme_minimal()

p3_pubmed <- pubmed_main2_gw %>%
  ggplot(aes(category2, hits_dup)) +
  geom_boxplot(aes(fill = category2)) +
  theme_minimal() +
  stat_compare_means()

p4_pubmed <- pubmed_main2_gw %>%
  ggplot(aes(hits_dup)) +
  geom_density(aes(fill = category2), alpha = 0.6)+
  theme_minimal()

(p1_pubmed + p2_pubmed) / (p3_pubmed + p4_pubmed)

# ------------------------------------------------------------------------------
# GW - DECIPHER - UNCERTAIN CNVs
# ------------------------------------------------------------------------------
library(regioneR)

tmp_decipher_uncertain <- read_tsv('/data-cbl/decipher_data/decipher-cnvs-grch37-2020-12-06.txt', skip = 1) %>%
  mutate(length = end - start + 1) %>%
  filter(length >= 50) %>%
  select(-length) %>%
  mutate(source = 'decipher') %>%
  rename(id = `# patient_id`, chrom = chr) %>%
  mutate(id = as.character(id)) %>% 
  filter(pathogenicity %in% c('Uncertain') & (! contribution %in% c('None'))) %>%
  filter(inheritance %in% 'De novo') %>%
  # filter(inheritance %in% c('De novo', 'Unknown')) %>%
  filter(genotype == 'Heterozygous') %>%
  mutate(phenotypes = str_replace_all(phenotypes, '\\|', '<br>')) %>%
  filter(variant_class %in% c('Deletion', 'Duplication')) %>%
  mutate(variant_class = tolower(variant_class)) %>%
  rename(clinical = pathogenicity) %>%
  mutate(clinical = 'uncertain') %>%
  select(chrom, start, end)


perm_pathogenic <- overlapPermTest(A = main2_gw %>% 
                                     filter(category2 == 'pathogenic_nohuman') %>% 
                                     select(chrom, start, end) %>%
                                     as.data.frame() %>% 
                                     mutate(chrom = paste0('chr', chrom)) %>%
                                     toGRanges(), 
                                   B = tmp_decipher_uncertain %>% 
                                     mutate(chrom = paste0('chr', chrom)) %>% 
                                     as.data.frame() %>% 
                                     toGRanges(),
                                   ntimes = 1e3)


perm_benign <- overlapPermTest(A = main2_gw %>% 
                                 filter(category2 == 'benign_nohuman') %>% 
                                 select(chrom, start, end) %>%
                                 as.data.frame() %>% 
                                 mutate(chrom = paste0('chr', chrom)) %>%
                                 toGRanges(), 
                               B = tmp_decipher_uncertain %>% 
                                 mutate(chrom = paste0('chr', chrom)) %>% 
                                 as.data.frame() %>% 
                                 toGRanges(),
                               ntimes = 1e3)

plot(perm_pathogenic)
plot(perm_benign)



# ------------------------------------------------------------------------------
# GW - CNV SYNDROMES
# ------------------------------------------------------------------------------

# cnv_syndromes_clingen_hg37 <- read_tsv('https://ftp.clinicalgenome.org/ClinGen%20recurrent%20CNV%20.bed%20file%20V1.1-hg19.bed',
#                                        skip = 1, col_names = FALSE)
# 
# cnv_syndromes_clingen_hg37 <- cnv_syndromes_clingen_hg37 %>%
#   # select(X1, X2, X3, X4) %>%
#   rename(chrom = X1,
#          start = X2,
#          end = X3,
#          syndrome_name = X4
#   ) %>%
#   mutate(start = 1 + start)  %>% # .bed format 0-based
#   mutate(source = 'clingen') %>%
#   select(chrom, start, end) %>%
#   mutate(id = row_number())
# 
# cnv_syndromes_clingen_hg37 %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   mutate(length_cnv = end - start + 1) %>%
#   bed_intersect(tmp_decipher) %>%
#   mutate(coverage = (.overlap / length_cnv.x)*100) %>% View()
#   filter(coverage > 90) %>% 
#   select(id.x) %>%
#   distinct() %>%
#   pull()
# 
# 
# 
# 
# perm_syndrome_pathogenic <- overlapPermTest(A = main2_gw %>% 
#                                      filter(category2 == 'pathogenic_nohuman') %>% 
#                                      select(chrom, start, end) %>%
#                                      as.data.frame() %>% 
#                                      mutate(chrom = paste0('chr', chrom)) %>%
#                                      toGRanges(), 
#                                    B = cnv_syndromes_clingen_hg37 %>% 
#                                      as.data.frame() %>% 
#                                      toGRanges(),
#                                    ntimes = 1e3)
# 
# 
# perm_syndrome_benign <- overlapPermTest(A = main2_gw %>% 
#                                  filter(category2 == 'benign_nohuman') %>% 
#                                  select(chrom, start, end) %>%
#                                  as.data.frame() %>% 
#                                  mutate(chrom = paste0('chr', chrom)) %>%
#                                  toGRanges(), 
#                                B = cnv_syndromes_clingen_hg37 %>% 
#                                  as.data.frame() %>% 
#                                  toGRanges(),
#                                ntimes = 1e3)
# 
# plot(perm_syndrome_pathogenic)
# plot(perm_syndrome_benign)
# 
# 



# ------------------------------------------------------------------------------
# GW - ALLELE FREQUENCY TEST
# ------------------------------------------------------------------------------


download.file('https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.bed.gz',
              'gnomad_v2.1_sv.controls_only.sites.bed.gz')


gnomad_af_test <- read_tsv('gnomad_v2.1_sv.controls_only.sites.bed.gz', col_types = cols(.default = "c")) %>%
  filter(SVTYPE %in% c('DEL', 'DUP')) %>%
  filter(FILTER == 'PASS') %>%
  rename(id = name) %>%
  mutate(id  = str_remove(id, 'gnomAD-SV_v2.1_')) %>%
  mutate(source = 'gnomad_v2.1') %>%
  rename(chrom = `#chrom`) %>%
  mutate(start = as.double(start), end = as.double(end), AF = as.double(AF)) %>%
  mutate(start = start + 1) %>%
  select(chrom, start, end, id, AF) %>%
  mutate(AF = ifelse(AF > 0.5, 1 - AF, AF))

# file.remove('gnomad_v2.1_sv.controls_only.sites.bed.gz')


main2_gw %>%
  bed_intersect(gnomad_af_test) %>%
  mutate(length_cnv = end.y - start.y + 1) %>%
  filter(.overlap / length_cnv > 0.9) %>%
  group_by(id.x) %>%
  summarise(max_af = min(AF.y)) %>%
  rename(id = id.x) %>%
  right_join(main2_gw) %>%
  mutate(max_af = log10(max_af)) %>%
  ggplot(aes(category2, max_af)) +
  geom_boxplot(aes(fill = category2)) +
  # scale_y_log10() +
  theme_minimal()


main2_gw %>%
  bed_intersect(gnomad_af_test) %>%
  mutate(length_cnv = end.y - start.y + 1) %>%
  # filter(.overlap / length_cnv > 0.9) %>%
  group_by(id.x) %>%
  # CAREFUL HERE
  summarise(max_af = min(AF.y)) %>%
  rename(id = id.x) %>%
  right_join(main2_gw) %>%
  # mutate(max_af = log10(max_af)) %>%
  ggplot(aes(category2, max_af)) +
  geom_boxplot(aes(fill = category2)) +
  # scale_x_log10() +
  theme_minimal()


main2_gw %>%
  bed_intersect(gnomad_af_test) %>%
  mutate(length_cnv = end.y - start.y + 1) %>%
  # filter(.overlap / length_cnv > 0.9) %>%
  group_by(id.x) %>%
  summarise(max_af = max(AF.y)) %>%
  mutate(tag_max_af = ntile(max_af, 10)) %>%
  rename(id = id.x) %>%
  right_join(main2_gw) %>%
  count(category2, tag_max_af) %>%
  group_by(tag_max_af) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(tag_max_af, perc)) +
    geom_col(aes(fill = category2))

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# GENERATE DATASETS FOR WHOLE GENOME ANALYSIS
# ------------------------------------------------------------------------------

# url_end <- 'https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/hmg/26/3/10.1093_hmg_ddw405/1/ddw405_Supp.zip?Expires=1616637355&Signature=b4T154U17FAN4GzkNjL-hW1RuUxrStcuMqaW1fo-zQ6lfmcfxnPx57nvV4vCcc-3TK6btuK~nQxCTzKeeqQQHdffpmzI52BvJGkD6R8yRyxJA0P4bbXJQqwcgodbCbygCIa0XoU30BMqO6H~vWuExK-kG5~cZ6QqCK-BJ1aGmZ4lvkp4CCHvVJTKcd2XndLwOJr-fsYcBeArVaEKvJ0Oejv2hkxessSS6-oeDZqaQBxl5i~J-STl4klfeObILVWz3codCG1qK8cC~DFUbmDyFvm~pDMRwjiWleQYzR8P0ZiTq~v6HAwiJEBerz9fdOlzF8BJkeIdN51MgrKEpRcEBA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA'
# 
# system('mkdir end_genes')
# download.file(url_end, destfile = 'end_genes/list_gene.zip')
# system('unzip end_genes/list_gene.zip')
# end_genes <- read_xlsx('ddw405-suppl_data/SupplTablesS1toS5SpataroetalHMG.xlsm', sheet = 2) 
# system('rm -r ddw405-suppl_data')
# system('rm -r end_genes')
# 
# end_genes <- end_genes %>% 
#   filter(Group == 'END') %>%
#   pull(`Gene Name`)

download.file('https://www.medrxiv.org/content/medrxiv/early/2021/05/04/2021.05.01.21256465/DC1/embed/media-1.xlsx?download=true',
              destfile = 'delete.xlsx')

male_infertility <- read_xlsx('delete.xlsx')

male_infertility <- male_infertility %>% rename(score = `Conclusion in 2020`, gene = `HGNC gene name`) %>% select(gene, score) %>%
  filter(score %in% c('Definitive', 'Strong', 'Moderate')) %>% select(gene)

file.remove('delete.xlsx' )


proteome_placenta <- read_tsv('https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Tissue_Protein_Expression_from_Human_Proteome_Map',
                              col_names = FALSE)

proteome_placenta <- proteome_placenta %>%
  filter(str_detect(X1, 'placenta')) %>%
  select(-X2) %>%
  pivot_longer(-X1, names_to = 'delete', values_to = 'gene') %>%
  na.omit() %>% 
  select(gene) %>%
  mutate(gene = str_remove(gene, ',.*$') ) %>%
  distinct()

# https://www.nature.com/articles/s41467-020-14284-2
download.file('https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-14284-2/MediaObjects/41467_2020_14284_MOESM7_ESM.xlsx',
              destfile = '41467_2020_14284_MOESM7_ESM.xlsx')
dl_no_disease <- read_xlsx('41467_2020_14284_MOESM7_ESM.xlsx', skip = 2) %>%
  pull(`Gene symbol`)


denovo_genes <- denovo %>% 
  filter(PrimaryPhenotype == 'developmentalDisorder') %>%
  # filter(PrimaryPhenotype != 'control') %>% 
  select(Gene) %>% 
  distinct() %>% 
  pull()

# High pLI

high_pli_genes <- hgcn_genes %>% filter(pLI > 90) %>% pull(gene)
# high_pli_genes <- hgcn_genes %>% filter(ccr > 0) %>% pull(gene)


# PFAM domains



download.file('https://s3.us-east-2.amazonaws.com/pathoscore-data/samocha/samochadenovo.xlsx', 'samochadenovo.xlsx')

samochadenovo <- read_xlsx('samochadenovo.xlsx', sheet = 2, col_types = 'text') %>%
  filter(dataset == 'id_ddd') %>%
  rename(gene = VEP_gene, start = pos) %>%
  # mutate(end = start) %>%
  # select(chrom, start, end, gene) %>%
  select(gene) %>%
  distinct() %>%
  pull()


file.remove('samochadenovo.xlsx')


samocha_new <- read_tsv('https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2832-5/MediaObjects/41586_2020_2832_MOESM3_ESM.txt')
samocha_new <- samocha_new %>% select(symbol) %>% distinct()

# ------------------------------------------------------------------------------
# MANHATTAN PLOT WHOLE GENOME
# ------------------------------------------------------------------------------


plot_data <- main_gw %>%  
  mutate(chrom = factor(chrom, levels = as.character(seq(1,22)))) %>%
  group_by(chrom) %>%
  summarise(chr_len=as.numeric(max(mid_point))) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(main_gw, ., by=c("chrom"="chrom")) %>%
  # arrange(chrom, mid_point) %>%
  mutate( BPcum=as.numeric(mid_point + tot))

axisdf <- plot_data %>% 
  mutate(chrom = factor(chrom, levels = as.character(seq(1,22)))) %>%
  group_by(chrom) %>% 
  summarize(center=(max(BPcum) + min(BPcum)) / 2)

plot_data %>%
  filter(category == 'pyu') %>%
  mutate(chrom = factor(chrom, levels = as.character(seq(1,22)))) %>%
  ggplot(aes(x=BPcum, y=pred_patho_control)) + 
  geom_point( aes(color= category), alpha=0.8, size=1.3, show.legend = FALSE) + 
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center) + 
  scale_y_continuous(expand = c(0, 0)) +
  # geom_smooth(formula = 'y ~ x', method = 'loess') +
  theme_minimal()

plot_data %>%
  # filter(!category %in% c('not_useful', 'ok_down')) %>%
  filter(category == 'pyu') %>%
  mutate(chrom = factor(chrom, levels = as.character(seq(1,22)))) %>%
  ggplot(aes(x=mid_point, y=pred_patho_control)) + 
  geom_point( aes(fill= category), color = 'black', alpha=0.8, size=1.3, show.legend = FALSE, shape = 21 ) + 
  # scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  # scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center) + 
  scale_y_continuous(expand = c(0, 0)) +
  geom_smooth(aes(group = category, color = category), formula = 'y ~ x', method = 'loess') +
  theme_minimal() +
  facet_wrap(vars(chrom), scales = 'free')

plot_data <- plot_data %>%
  bed_intersect(hgcn_genes) %>%
  count(chrom, start.x, end.x) %>%
  rename(start = start.x, end = end.x, n_genes = n) %>%
  right_join(plot_data) %>%
  mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes)) %>%
  group_by(chrom) %>%
  mutate(n_genes = scale(n_genes)) %>%
  mutate(pred_patho_control = scale(pred_patho_control)) %>%
  mutate(pred_patho_no_control = scale(pred_patho_no_control))

plot_data <- plot_data %>%
  bed_intersect(pubmed_df) %>%
  group_by(chrom, mid_point.x) %>%
  summarise(n_hits_del = sum(hits_del.y)) %>%
  rename(mid_point = mid_point.x) %>%
  select(chrom, mid_point, n_hits_del) %>%
  group_by(chrom) %>%
  mutate(n_hits_del = scale(n_hits_del)) %>%
  # rename(n_hits_del = `n_hits_del[,1]`)
  right_join(plot_data) %>%
  mutate(n_hits_del = ifelse(is.na(n_hits_del), 0, n_hits_del))




plot_data %>%
  # filter(category == 'pyu') %>%
  select(chrom, mid_point, category, pred_patho_control, pred_patho_no_control, 
         n_genes
         # n_hits_del
  ) %>%
  pivot_longer(cols = -c('chrom', 'mid_point', 'category'), names_to = 'type', values_to = 'value') %>%
  mutate(chrom = factor(chrom, levels = as.character(seq(1,22)))) %>%
  ggplot(aes(x=mid_point, y=value)) + 
  # geom_point( aes(fill= category), color = 'black', alpha=0.8, size=1.3, show.legend = FALSE, shape = 21 ) + 
  # scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  # scale_x_continuous(label = axisdf$chrom, breaks= axisdf$center) + 
  scale_y_continuous(expand = c(0, 0)) +
  geom_smooth(aes(group = type, color = type), formula = 'y ~ x', method = 'loess') +
  theme_minimal() +
  facet_wrap(vars(chrom), scales = 'free')


test_count_per_cytoband <- coord_cytobands %>%
  rowwise() %>%
  mutate(n_cnvs = valr::bed_intersect(tibble('chrom' = chrom, 'start' = start, 'end' = end), input_check_cnv_deletion) %>% nrow()) %>%
  ungroup() 

test_count_per_cytoband %>%
  arrange(chrom, start) %>%
  # ggplot(aes(reorder(Name, -n_cnvs), n_cnvs)) +
  ggplot(aes(Name, n_cnvs)) +
  geom_col(aes(fill = chrom), color = 'black', show.legend = FALSE) +
  facet_wrap(~ chrom, scales = 'free') +
  theme_minimal() +
  ggtitle(paste('Distribution of CNVs across genome (counts by chromosome) -')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))





# ------------------------------------------------------------------------------
# CALCULATE RESIDUALS DISTANCE QUANTILE PACO
# ------------------------------------------------------------------------------


# ClinVar
# plots_clinvar_unbiased
# write_tsv(plot_rates_decipher_unbiased, file = 'plot_rates_decipher_with_unknown.tsv')

# cut_interval, ntile(residuals) - 
# Best combinations - ntile100) - ntile(100) - fixed=FALSE, n = 20
# not bad ntile(500), ntile(10), FALSE
# good: ntile(500) ntile(10) FALSE - 20
# good: ntile(600) ntile(10) FALSE - 10
# the best so far: ntile(700) - ntile(10) - FALSE - 20
# the best so far (unknown): ntile(700) - ntile(10) - FALSE - 10
# the best so far (unknown): ntile(100) - ntile(5) - FALSE - 10--5-30
# the best so far NPV (unknown): ntile(10) - ntile(10) - FALSE - 10--5-30
# the best so far!!! (unknown): ntile(10) - ntile(5) - FALSE - 10
# the best so far!!!!!! (unknown): ntile(5) - ntile(10) - FALSE - 5-10-20-30
# 5-5 DOESN'T WORK WELL
# 10-10 doesn't work either
# the best so far!!!!!! (unknown): ntile(10) - ntile(3) - TRUE - must be 3
# the best so far!!!!!! (unknown): ntile(5) - ntile(3) - TRUE - must be 3 - better Prec - Plain NPV
# (same) sensitivity is higher for cnvscore than structure


# test1_tbl <- tibble()
# for (i in seq(5, 100)) {


# tmp_both_clinvar_prescore <- predict_chrom_aware_rtemis(bayesian_clinvar_del_both, output_clinvar_deletion, 'deletion', 'both approach')
# tmp_both_decipher_prescore <- predict_chrom_aware_rtemis(bayesian_clinvar_del_both, output_decipher_deletion, 'deletion', 'both approach')
# tmp_both_clinvarind_prescore <- predict_chrom_aware_rtemis(bayesian_clinvar_del_both, output_clinvar_20, 'deletion', 'both approach')
# 
# tmp_human_clinvar_prescore <- predict_chrom_aware_rtemis(bayesian_clinvar_del_human, 
#                                                          output_clinvar_deletion, 'deletion', 'human')
# tmp_human_decipher_prescore <- predict_chrom_aware_rtemis(bayesian_clinvar_del_human, 
#                                                           output_decipher_deletion, 'deletion', 'human')
# tmp_human_clinvarind_prescore <- predict_chrom_aware_rtemis(bayesian_clinvar_del_human, 
#                                                             output_clinvar_20, 'deletion', 'human')
# 
# tmp_unbiased_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_deletion, 'deletion', 'unbiased')
# tmp_unbiased_decipher <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_decipher_deletion, 'deletion', 'unbiased')
# tmp_unbiased_clinvarind <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_20, 'deletion', 'unbiased')
# 
# tmp_unbiased_decipher_challenge1 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, decipher_setting_1, 'deletion', 'unbiased')
# tmp_unbiased_decipher_challenge2 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, decipher_setting_2, 'deletion', 'unbiased')
# tmp_unbiased_decipher_challenge3 <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, decipher_setting_3, 'deletion', 'unbiased')

tmp_dup_clinvar <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_nohuman, output_clinvar_duplication, 'deletion', 'unbiased')


split_score <- 10
split_residuals <- 3





# Deletion CNVs----------

ref_sd_clinvar_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_deletion, 'deletion', 'unbiased approach')
ref_sd_clinvar_dup <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_nohuman, output_clinvar_duplication, 'duplication', 'unbiased approach')


trendline_clinvar <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
                         data= ref_sd_clinvar_del[[3]])


res_df <- ref_sd_clinvar_del[[3]] %>%
  mutate(gam_residuals = residuals(trendline_clinvar)) %>%
  mutate(score_interval = ntile(.pred_pathogenic, n = 10)) %>%
  group_by(score_interval) %>%
  mutate(reliability_score = ntile(gam_residuals, 3)) %>%
  ungroup()


score_intervals_df <- res_df %>% 
  group_by(score_interval) %>% 
  summarise(max_score = max(.pred_pathogenic), 
            min_score = min(.pred_pathogenic)) %>%
  arrange(score_interval)

ref_quantiles <- res_df %>%
  left_join(score_intervals_df, by = 'score_interval') %>%
  mutate(score_interval = as.character(score_interval)) %>%
  select(score_interval, min_score, max_score, reliability_score, sd, gam_residuals)


split_score_dup <- 3

# Duplication CNVs----------

trendline_clinvar_dup <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
                         data= ref_sd_clinvar_dup[[3]])

res_df_dup <- ref_sd_clinvar_dup[[3]] %>%
  mutate(gam_residuals = residuals(trendline_clinvar_dup)) %>%
  mutate(score_interval = ntile(.pred_pathogenic, n = 3)) %>%
  group_by(score_interval) %>%
  mutate(reliability_score = ntile(gam_residuals, 3)) %>%
  ungroup()


score_intervals_df_dup <- res_df_dup %>% 
  group_by(score_interval) %>% 
  summarise(max_score = max(.pred_pathogenic), 
            min_score = min(.pred_pathogenic)) %>%
  arrange(score_interval)


ref_quantiles_dup <- res_df_dup %>%
  left_join(score_intervals_df_dup, by = 'score_interval') %>%
  mutate(score_interval = as.character(score_interval)) %>%
  select(score_interval, min_score, max_score, reliability_score, sd, gam_residuals)


#----------

# 
ref_quantiles %>%
  ggplot(aes(score_interval, sd)) +
    geom_boxplot(aes(fill = as.factor(reliability_score)))

ref_quantiles %>%
  ggplot(aes(score_interval, sd)) +
  geom_point(aes(fill = as.factor(reliability_score)))

ref_quantiles %>% select(.pred_pathogenic, reliability_score) %>% correlate()

res_df %>%
  mutate(reliability_score = factor(reliability_score)) %>%
  mutate(reliability_score = case_when(
    reliability_score == 1 ~ 'Lowly uncertain',
    reliability_score == 2 ~ 'Moderately uncertain',
    reliability_score == 3 ~ 'Highly uncertain'
  )) %>%
  mutate(reliability_score = factor(reliability_score, levels = c('Highly uncertain','Moderately uncertain','Lowly uncertain'))) %>%
ggplot(aes(.pred_pathogenic, sd)) +
geom_point(aes(fill = reliability_score), shape = 21, alpha = 0.6) +
  scale_fill_brewer(palette = 'Set1') +
  theme_minimal() +
  labs(x = 'CNVscore', y = 'Standard deviation (sd)', fill = 'Uncertainty level') +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15))


res_df_dup %>%
  mutate(reliability_score = factor(reliability_score)) %>%
  mutate(reliability_score = case_when(
    reliability_score == 1 ~ 'Lowly uncertain',
    reliability_score == 2 ~ 'Moderately uncertain',
    reliability_score == 3 ~ 'Highly uncertain'
  )) %>%
  mutate(reliability_score = factor(reliability_score, levels = c('Highly uncertain','Moderately uncertain','Lowly uncertain'))) %>%
  ggplot(aes(.pred_pathogenic, sd)) +
  geom_point(aes(fill = reliability_score), shape = 21, alpha = 0.6) +
  scale_fill_brewer(palette = 'Set1') +
  theme_minimal() +
  labs(x = 'CNVscore', y = 'Standard deviation (sd)', fill = 'Uncertainty level') +
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15))
  
   
# theme_minimal()

p2_rate <- res_df %>%
  ggplot(aes(.pred_pathogenic, reliability_score)) +
  geom_point(aes(color = as.factor(score_interval))) +
  theme_minimal()

res_df %>%
  ggplot(aes(reliability_score, .pred_pathogenic)) +
  geom_point(aes(color = as.factor(reliability_score))) +
  theme_minimal()

res_df %>%
  ggplot(aes(.pred_pathogenic)) +
  geom_density(aes(fill = as.factor(reliability_score))) +
  theme_minimal()

res_df_dup %>%
  ggplot(aes(.pred_pathogenic)) +
  geom_density(aes(fill = as.factor(reliability_score))) +
  theme_minimal()


p3_rate <- res_df %>%
  ggplot(aes(.pred_pathogenic, sd)) +
  geom_point(aes(color = reliability_score),  show.legend = T) +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cr")) +
  theme_minimal() +
  scale_color_distiller(palette = "Spectral") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

almost_clinvar <- plot_rates_clinvar_unbiased %>%
  left_join(res_df %>% select(gam_residuals, id, reliability_score), by = 'id') %>%
  na.omit()


# Reference knowledge-based

# ref_knowledge_based_quantiles <- res_df_knowledge_based %>%
#   group_by(score_interval, reliability_score) %>%
#   summarise(max_sd = max(sd),
#             min_sd = min(sd)
#   ) %>%
#   ungroup() %>%
#   left_join(res_df_knowledge_based %>% 
#               group_by(score_interval) %>% 
#               summarise(max_score = max(cnvscore_score_knowledge_based), 
#                         min_score = min(cnvscore_score_knowledge_based)), by = 'score_interval') %>%
#   mutate(score_interval = as.character(score_interval))

# Reference general - PACO VERSION II






# ref_quantiles <- res_df %>%
#   group_by(score_interval, reliability_score) %>%
#   summarise(max_gam_residuals = max(gam_residuals),
#             min_gam_residuals = min(gam_residuals)
#   ) %>%
#   ungroup() %>%
#   left_join(res_df %>% 
#               group_by(score_interval) %>% 
#               summarise(max_score = max(.pred_pathogenic), 
#                         min_score = min(.pred_pathogenic)), by = 'score_interval') %>%
#   mutate(score_interval = as.character(score_interval)) %>%
#   select(score_interval, reliability_score, max_score, min_score, max_gam_residuals , min_gam_residuals)
# mutate(mid_score = (max_score + min_score) / 2) %>%
# mutate(mid_sd = (max_sd + min_sd) / 2)

# Reference general - USING ANTONIO VERSION

# ref_quantiles <- res_df %>%
  # group_by(score_interval, reliability_score) %>%
  # summarise(max_gam_residuals = max(gam_residuals),
  #           min_gam_residuals = min(gam_residuals)
  # ) %>%
  # ungroup() %>%
  # left_join(res_df %>% 
  #             group_by(score_interval) %>% 
  #             summarise(max_score = max(.pred_pathogenic), 
  #                       min_score = min(.pred_pathogenic)), by = 'score_interval') %>%
  # mutate(score_interval = as.character(score_interval)) %>%
  # select(score_interval, min_score, max_score, sd)
  # mutate(mid_score = (max_score + min_score) / 2) %>%
  # mutate(mid_sd = (max_sd + min_sd) / 2)



  

# Reference total model

# ref_both_quantiles <- res_df_both %>%
#   group_by(score_interval, reliability_score) %>%
#   summarise(max_sd = max(sd),
#             min_sd = min(sd)
#   ) %>%
#   ungroup() %>%
#   left_join(res_df_both %>% 
#               group_by(score_interval) %>% 
#               summarise(max_score = max(cnvscore_score_both), 
#                         min_score = min(cnvscore_score_both)), by = 'score_interval') %>%
#   mutate(score_interval = as.character(score_interval))

# Reference combined score (sd + score)

# ref_combined_quantiles <- res_df_combined %>%
#   group_by(score_interval, reliability_score) %>%
#   summarise(max_sd = max(sd),
#             min_sd = min(sd)
#   ) %>%
#   ungroup() %>%
#   left_join(res_df_combined %>% 
#               group_by(score_interval) %>% 
#               summarise(max_score = max(cnvscore_score_combined), 
#                         min_score = min(cnvscore_score_combined)), by = 'score_interval') %>%
#   mutate(score_interval = as.character(score_interval))


# Result combined 
# vector_quantiles_combined_decipher <- c()
# 
# for (i in 1:nrow(plot_rates_decipher_unbiased)) {
#   
#   tmp_score_interval <- ref_combined_quantiles %>%
#     mutate(mid_score = (max_score + min_score) / 2) %>%
#     mutate(dist_mid_score = 
#              abs(plot_rates_decipher_unbiased[i,]$cnvscore_score - mid_score)) %>%
#     arrange(dist_mid_score) %>%
#     slice(1) %>%
#     pull(score_interval) %>% 
#     as.character()
#   
#   tmp_to_vector <- ref_combined_quantiles %>%
#     filter(score_interval == tmp_score_interval) %>%
#     mutate(mid_sd = (max_sd + min_sd) / 2) %>%
#     mutate(dist_mid_sd = abs(plot_rates_decipher_unbiased[i,]$sd - mid_sd)) %>%
#     arrange(dist_mid_sd) %>%
#     slice(1) %>%
#     pull(reliability_score)
#   
#   vector_quantiles_combined_decipher <- c(vector_quantiles_combined_decipher, tmp_to_vector)
#   
# }
# 
# almost_combined_decipher <- plot_rates_decipher_unbiased %>%
#   mutate(reliability_score = vector_quantiles_combined_decipher)

# almost_combined_decipher %>% roc_auc(clinical, cnvscore_score)
# almost_combined_decipher %>% roc_auc(clinical, reliability_score)


# DECIPHER - Knowledge-based model

# 
# vector_quantiles_knowledge_based_decipher <- c()
# 
# for (i in 1:nrow(tmp_human_decipher_prescore[[3]])) {
#   
#   tmp_score_interval <- ref_both_quantiles %>%
#     mutate(mid_score = (max_score + min_score) / 2) %>%
#     mutate(dist_mid_score = 
#              abs(tmp_human_decipher_prescore[[3]][i,]$.pred_pathogenic - mid_score)) %>%
#     arrange(dist_mid_score) %>%
#     slice(1) %>%
#     pull(score_interval) %>% 
#     as.character()
#   
#   tmp_to_vector <- ref_both_quantiles %>%
#     filter(score_interval == tmp_score_interval) %>%
#     mutate(mid_sd = (max_sd + min_sd) / 2) %>%
#     mutate(dist_mid_sd = abs(tmp_human_decipher_prescore[[3]][i,]$sd - mid_sd)) %>%
#     arrange(dist_mid_sd) %>%
#     slice(1) %>%
#     pull(reliability_score)
#   
#   vector_quantiles_knowledge_based_decipher <- c(vector_quantiles_knowledge_based_decipher, 
#                                                  tmp_to_vector)
#   
# }
# 
# almost_knowledge_based_decipher <- tmp_human_decipher_prescore[[3]] %>%
#   mutate(reliability_score = vector_quantiles_both_decipher)


# DECIPHER - Both model
# 
# vector_quantiles_both_decipher <- c()
# 
# for (i in 1:nrow(tmp_both_decipher_prescore[[3]])) {
#   
#   tmp_score_interval <- ref_both_quantiles %>%
#     mutate(mid_score = (max_score + min_score) / 2) %>%
#     mutate(dist_mid_score = 
#              abs(tmp_both_decipher_prescore[[3]][i,]$.pred_pathogenic - mid_score)) %>%
#     arrange(dist_mid_score) %>%
#     slice(1) %>%
#     pull(score_interval) %>% 
#     as.character()
#   
#   tmp_to_vector <- ref_both_quantiles %>%
#     filter(score_interval == tmp_score_interval) %>%
#     mutate(mid_sd = (max_sd + min_sd) / 2) %>%
#     mutate(dist_mid_sd = abs(tmp_both_decipher_prescore[[3]][i,]$sd - mid_sd)) %>%
#     arrange(dist_mid_sd) %>%
#     slice(1) %>%
#     pull(reliability_score)
#   
#   vector_quantiles_both_decipher <- c(vector_quantiles_both_decipher, tmp_to_vector)
#   
# }
# 
# almost_both_decipher <- tmp_both_decipher_prescore[[3]] %>%
#   mutate(reliability_score = vector_quantiles_both_decipher)

# DECIPHER - Unbiased model
# 
# vector_quantiles_decipher <- c()
# 
# for (i in 1:nrow(plot_rates_decipher_unbiased)) {
#   
#   tmp_score_interval <- ref_quantiles %>%
#     mutate(mid_score = (max_score + min_score) / 2) %>%
#     mutate(dist_mid_score = 
#              abs(plot_rates_decipher_unbiased[i,]$cnvscore_score - mid_score)) %>%
#     arrange(dist_mid_score) %>%
#     slice(1) %>%
#     pull(score_interval) %>% 
#     as.character()
#   
#   tmp_to_vector <- ref_quantiles %>%
#     filter(score_interval == tmp_score_interval) %>%
#     mutate(mid_sd = (max_sd + min_sd) / 2) %>%
#     mutate(dist_mid_sd = abs(plot_rates_decipher_unbiased[i,]$sd - mid_sd)) %>%
#     arrange(dist_mid_sd) %>%
#     slice(1) %>%
#     pull(reliability_score)
#   
#   vector_quantiles_decipher <- c(vector_quantiles_decipher, tmp_to_vector)
#   
# }
# 
# almost_decipher <- plot_rates_decipher_unbiased %>%
#   mutate(reliability_score = vector_quantiles_decipher)

# DECIPHER - Scenario #1
# 
# vector_quantiles_decipher_challenge1 <- c()
# 
# for (i in 1:nrow(tmp_unbiased_decipher_challenge1[[3]])) {
#   
#   tmp_score_interval <- ref_quantiles %>%
#     mutate(mid_score = (max_score + min_score) / 2) %>%
#     mutate(dist_mid_score = 
#              abs(tmp_unbiased_decipher_challenge1[[3]][i,]$.pred_pathogenic - mid_score)) %>%
#     arrange(dist_mid_score) %>%
#     slice(1) %>%
#     pull(score_interval) %>% 
#     as.character()
#   
#   tmp_to_vector <- ref_quantiles %>%
#     filter(score_interval == tmp_score_interval) %>%
#     mutate(mid_sd = (max_sd + min_sd) / 2) %>%
#     mutate(dist_mid_sd = abs(tmp_unbiased_decipher_challenge1[[3]][i,]$sd - mid_sd)) %>%
#     arrange(dist_mid_sd) %>%
#     slice(1) %>%
#     pull(reliability_score)
#   
#   vector_quantiles_decipher_challenge1 <- c(vector_quantiles_decipher_challenge1, tmp_to_vector)
#   
# }


# DECIPHER - Scenario #2

vector_quantiles_decipher_challenge2 <- c()

for (i in 1:nrow(tmp_unbiased_decipher_challenge2[[3]])) {
  
  tmp_score_interval <- ref_quantiles %>%
    mutate(mid_score = (max_score + min_score) / 2) %>%
    mutate(dist_mid_score = 
             abs(tmp_unbiased_decipher_challenge2[[3]][i,]$.pred_pathogenic - mid_score)) %>%
    arrange(dist_mid_score) %>%
    slice(1) %>%
    pull(score_interval) %>% 
    as.character()
  
  tmp_to_vector <- ref_quantiles %>%
    filter(score_interval == tmp_score_interval) %>%
    mutate(mid_sd = (max_sd + min_sd) / 2) %>%
    mutate(dist_mid_sd = abs(tmp_unbiased_decipher_challenge2[[3]][i,]$sd - mid_sd)) %>%
    arrange(dist_mid_sd) %>%
    slice(1) %>%
    pull(reliability_score)
  
  vector_quantiles_decipher_challenge2 <- c(vector_quantiles_decipher_challenge2, tmp_to_vector)
  
}

# DECIPHER - Scenario #3

vector_quantiles_decipher_challenge3 <- c()

for (i in 1:nrow(tmp_unbiased_decipher_challenge3[[3]])) {
  
  tmp_score_interval <- ref_quantiles %>%
    mutate(mid_score = (max_score + min_score) / 2) %>%
    mutate(dist_mid_score = 
             abs(tmp_unbiased_decipher_challenge3[[3]][i,]$.pred_pathogenic - mid_score)) %>%
    arrange(dist_mid_score) %>%
    slice(1) %>%
    pull(score_interval) %>% 
    as.character()
  
  tmp_to_vector <- ref_quantiles %>%
    filter(score_interval == tmp_score_interval) %>%
    mutate(mid_sd = (max_sd + min_sd) / 2) %>%
    mutate(dist_mid_sd = abs(tmp_unbiased_decipher_challenge2[[3]][i,]$sd - mid_sd)) %>%
    arrange(dist_mid_sd) %>%
    slice(1) %>%
    pull(reliability_score)
  
  vector_quantiles_decipher_challenge3 <- c(vector_quantiles_decipher_challenge3, tmp_to_vector)
  
}



# ClinVar independent dataset - Unbiased model

vector_quantiles_clinvarind <- c()

for (i in 1:nrow(plot_rates_clinvarind_unbiased)) {
  
  tmp_score_interval <- ref_quantiles %>%
    mutate(mid_score = (max_score + min_score) / 2) %>%
    mutate(dist_mid_score = 
             abs(plot_rates_clinvarind_unbiased[i,]$cnvscore_score - mid_score)) %>%
    arrange(dist_mid_score) %>%
    slice(1) %>%
    pull(score_interval) %>% 
    as.character()
  
  tmp_to_vector <- ref_quantiles %>%
    filter(score_interval == tmp_score_interval) %>%
    mutate(mid_sd = (max_sd + min_sd) / 2) %>%
    mutate(dist_mid_sd = abs(plot_rates_clinvarind_unbiased[i,]$sd - mid_sd)) %>%
    arrange(dist_mid_sd) %>%
    slice(1) %>%
    pull(reliability_score)
  
  vector_quantiles_clinvarind <- c(vector_quantiles_clinvarind, tmp_to_vector)
  
}

almost_clinvarind <- plot_rates_clinvarind_unbiased %>%
  mutate(reliability_score = vector_quantiles_clinvarind)



# ClinVar independent dataset - Knowledge-based model

vector_quantiles_clinvarind_knowledge_based <- c()

for (i in 1:nrow(tmp_human_clinvarind_prescore[[3]])) {
  
  tmp_score_interval <- ref_both_quantiles %>%
    mutate(mid_score = (max_score + min_score) / 2) %>%
    mutate(dist_mid_score = 
             abs(tmp_human_clinvarind_prescore[[3]][i,]$.pred_pathogenic - mid_score)) %>%
    arrange(dist_mid_score) %>%
    slice(1) %>%
    pull(score_interval) %>% 
    as.character()
  
  tmp_to_vector <- ref_both_quantiles %>%
    filter(score_interval == tmp_score_interval) %>%
    mutate(mid_sd = (max_sd + min_sd) / 2) %>%
    mutate(dist_mid_sd = abs(tmp_human_clinvarind_prescore[[3]][i,]$sd - mid_sd)) %>%
    arrange(dist_mid_sd) %>%
    slice(1) %>%
    pull(reliability_score)
  
  vector_quantiles_clinvarind_knowledge_based <- c(vector_quantiles_clinvarind_knowledge_based, tmp_to_vector)
  
}

almost_clinvarind_knowledge_based <- tmp_human_clinvarind_prescore[[3]] %>%
  mutate(reliability_score = vector_quantiles_clinvarind_knowledge_based)



# ClinVar independent dataset - Both model

vector_quantiles_clinvarind_both <- c()

for (i in 1:nrow(tmp_both_clinvarind_prescore[[3]])) {
  
  tmp_score_interval <- ref_both_quantiles %>%
    mutate(mid_score = (max_score + min_score) / 2) %>%
    mutate(dist_mid_score = 
             abs(tmp_both_clinvarind_prescore[[3]][i,]$.pred_pathogenic - mid_score)) %>%
    arrange(dist_mid_score) %>%
    slice(1) %>%
    pull(score_interval) %>% 
    as.character()
  
  tmp_to_vector <- ref_both_quantiles %>%
    filter(score_interval == tmp_score_interval) %>%
    mutate(mid_sd = (max_sd + min_sd) / 2) %>%
    mutate(dist_mid_sd = abs(tmp_both_clinvarind_prescore[[3]][i,]$sd - mid_sd)) %>%
    arrange(dist_mid_sd) %>%
    slice(1) %>%
    pull(reliability_score)
  
  vector_quantiles_clinvarind_both <- c(vector_quantiles_clinvarind_both, tmp_to_vector)
  
}

almost_clinvarind_both <- tmp_both_clinvarind_prescore[[3]] %>%
  mutate(reliability_score = vector_quantiles_clinvarind_both)


# Special plot IRENE

plot1_reli_challenges <- almost_decipher %>% 
  # filter(clinical == 'pathogenic') %>%
  select(clinical, reliability_score) %>%
  mutate(dataset = 'DECIPHER') %>%
  bind_rows(
    almost_clinvar %>%
      select(clinical, reliability_score) %>%
      mutate(dataset = 'ClinVar')
  ) %>%
  bind_rows(
    almost_clinvarind %>%
      select(clinical, reliability_score) %>%
      mutate(dataset = 'ClinVar >Jan 2021')
  ) %>%
  
  bind_rows(
    decipher_setting_1 %>% 
      mutate(reliability_score = vector_quantiles_decipher_challenge1) %>%
      mutate(dataset = 'Scenario #1') %>%
      select(clinical, reliability_score, dataset)

  ) %>%
  bind_rows(
    decipher_setting_2 %>% 
      mutate(reliability_score = vector_quantiles_decipher_challenge2) %>%
      mutate(dataset = 'Scenario #2') %>%
      select(clinical, reliability_score, dataset)
  ) %>%
  bind_rows(
    decipher_setting_3 %>% 
      mutate(reliability_score = vector_quantiles_decipher_challenge3) %>%
      mutate(dataset = 'Scenario #3') %>%
      select(clinical, reliability_score, dataset)
  ) %>%
  count(clinical, reliability_score, dataset) %>%
  group_by(clinical, dataset) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(reliability_score = factor(reliability_score)) %>%
  filter(!is.na(clinical)) %>%
  ggplot(aes(dataset, perc, group = clinical)) +
    scale_y_continuous(label = percent) +
    geom_col(aes(fill = reliability_score), color = 'black') +
    geom_label(aes(label = paste0(100*round(perc, 2), '%' , ' (', n, ')')), size = 3, position = position_stack(vjust = 0.5)) +
    theme_minimal() +
  facet_wrap(vars(clinical)) +
  labs(title = 'Reliability score of the DECIPHER dataset and challenging scenarios',
  subtitle = 'Scenario #1 (patho (no omim) - benign (omim)\n
      Scenario #2 (patho (no omim) - benign (no omim))\n
      Scenario #3 patho & benign - (0 genes)',
    x = 'Dataset', y = 'Percentage') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



df_reli_challenges <- almost_decipher %>% 
  select(reliability_score) %>%
  mutate(dataset = 'DECIPHER') %>%
  bind_rows(
    almost_clinvar %>%
      select(reliability_score) %>%
      mutate(dataset = 'ClinVar')
  ) %>%
  bind_rows(
    almost_clinvarind %>%
      select(reliability_score) %>%
      mutate(dataset = 'ClinVar >Jan 2021')
  ) %>%
  
  bind_rows(
    decipher_setting_1 %>% 
      mutate(reliability_score = vector_quantiles_decipher_challenge1) %>%
      mutate(dataset = 'Scenario #1') %>%
      select(reliability_score, dataset)
    
  ) %>%
  bind_rows(
    decipher_setting_2 %>% 
      mutate(reliability_score = vector_quantiles_decipher_challenge2) %>%
      mutate(dataset = 'Scenario #2') %>%
      select(reliability_score, dataset)
  ) %>%
  bind_rows(
    decipher_setting_3 %>% 
      mutate(reliability_score = vector_quantiles_decipher_challenge3) %>%
      mutate(dataset = 'Scenario #3') %>%
      select(reliability_score, dataset)
  ) %>%
  bind_rows(
    decipher_setting_4 %>% 
      mutate(reliability_score = vector_quantiles_decipher_challenge4) %>%
      mutate(dataset = 'Scenario #3') %>%
      select(reliability_score, dataset)
  ) %>%
  count(reliability_score, dataset) %>%
  group_by(dataset) %>%
  mutate(perc = n / sum(n)) %>%
  mutate(reliability_score = factor(reliability_score))




plot2_reli_challenges <- df_reli_challenges %>%
  ggplot(aes(dataset, perc, group = reliability_score)) +
  scale_y_continuous(label = percent) +
  geom_col(aes(fill = reliability_score), color = 'black') +
  geom_label(aes(label = paste0(100*round(perc, 2), '%' , ' (', n, ')')), size = 3, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  labs(title = 'Reliability score of the DECIPHER dataset and challenging scenarios',
       subtitle = 'Scenario #1 (patho (no omim) - benign (omim)\n
      Scenario #2 (patho (no omim) - benign (no omim))\n
      Scenario #3 patho & benign - (0 genes)',
       x = 'Dataset', y = 'Percentage') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


plot2_reli_challenges + plot1_reli_challenges


tmp_auroc_datasets <- tibble(dataset = c("ClinVar", 'DECIPHER', 'ClinVar >Jan 2021', 'Scenario #1',
                   'Scenario #2', 'Scenario #3'),
       cnvscore = c(
         tmp_unbiased_clinvar[[3]] %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate),
         tmp_unbiased_decipher[[3]] %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate),
         tmp_unbiased_clinvarind[[3]]  %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate),
         tmp_unbiased_decipher_challenge1[[3]]  %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate),
         tmp_unbiased_decipher_challenge2[[3]]  %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate),
         tmp_unbiased_decipher_challenge3[[3]]  %>% roc_auc(clinical, .pred_pathogenic) %>% pull(.estimate)
       ))



tmp_auroc_structure <- tibble(dataset = c("ClinVar", 'DECIPHER', 'ClinVar >Jan 2021', 'Scenario #1',
                                          'Scenario #2', 'Scenario #3'),
                             strvctvre = c(
                               almost_clinvar %>% roc_auc(clinical, structure_score) %>% pull(.estimate),
                               almost_decipher %>% roc_auc(clinical, structure_score) %>% pull(.estimate),
                               almost_clinvarind  %>% roc_auc(clinical, structure_score) %>% pull(.estimate),
                               result_decipher_setting_1[[1]][result_decipher_setting_1[[1]]$tool == 'STRVCTURE',]$roc_auc,
                               result_decipher_setting_2[[1]][result_decipher_setting_2[[1]]$tool == 'STRVCTURE',]$roc_auc,
                               result_decipher_setting_3[[1]][result_decipher_setting_3[[1]]$tool == 'STRVCTURE',]$roc_auc
                             ))





plot_comparison_reliability_1 <- df_reli_challenges %>%
  filter(reliability_score == 1) %>%
  left_join(tmp_auroc_datasets, by = 'dataset') %>%
  left_join(tmp_auroc_structure, by = 'dataset') %>%
  select(-n) %>%
  pivot_longer(-c(dataset, reliability_score, perc), names_to = 'tool', values_to = 'AUROC') %>%
  ggplot(aes(perc, AUROC)) +
    geom_label(aes(label = dataset, color = tool), vjust = 0.8, nudge_y = 0.05) +
  geom_point(aes(color = tool)) +
    theme_minimal()


tmp11 <- 0.6
tmp22 <- 3
tmp33 <- 'fixed-interval'

# plot_quantile1 <- plot_rates2(almost_clinvar, tag = 'ClinVar', 
#                               tag2 = 'Unbiased features', 
#                               input_patho = 0.5, 
#                               input_benign = 0.5, 
#                               fixed_interval = tmp33,
#                               limit_y = tmp11,
#                               n_split = tmp22,
#                               tmp_sub = paste0('Score: equally-sized bins (n = ', split_score, ') ', '\n',
#                                                'Residuals: equally-sized bins (n = ', split_residuals, ')', '\n',
#                                                'Representation: ', tmp33, ' bins ', '(n = ',tmp22, ')' ))
# 
# plot_quantile2 <- plot_rates2(almost_decipher, tag = 'DECIPHER', 
#                               tag2 = 'Unbiased features', 
#                               input_patho = 0.5, 
#                               input_benign = 0.5, 
#                               fixed_interval = tmp33,
#                               limit_y = tmp11,
#                               n_split = tmp22,
#                               tmp_sub = paste0('Score: equally-sized bins (n = ', split_score, ') ', '\n',
#                                                'Residuals: equally-sized bins (n = ', split_residuals, ')', '\n',
#                                                'Representation: ', tmp33, ' bins ', '(n = ',tmp22, ')'))
# 
# plot_quantile3 <- plot_rates2(almost_clinvarind, tag = 'ClinVar independent', 
#                               tag2 = 'Unbiased features', 
#                               input_patho = 0.5, 
#                               input_benign = 0.5, 
#                               fixed_interval = tmp33,
#                               limit_y = tmp11,
#                               n_split = tmp22,
#                               tmp_sub = paste0('Score: equally-sized bins (n = ', split_score, ') ', '\n',
#                                                'Residuals: equally-sized bins (n = ', split_residuals, ')', '\n',
#                                                'Representation: ', tmp33, ' bins ', '(n = ',tmp22, ')'))
# 
# p1_rate <- plot_quantile1[[1]] + plot_quantile1[[2]] + plot_quantile2[[1]] + 
#   plot_quantile2[[2]] + plot_quantile3[[1]] + plot_quantile3[[2]] + plot_layout(nrow = 3)

almost_knowledge_based_decipher <- almost_knowledge_based_decipher %>%
  rename(cnvscore_score_knowledge_based = .pred_pathogenic) %>%
  select(id, cnvscore_score_knowledge_based, reliability_score, clinical, sd )

almost_clinvarind_knowledge_based <- almost_clinvarind_knowledge_based %>%
  rename(cnvscore_score_knowledge_based = .pred_pathogenic) %>%
  select(id, cnvscore_score_knowledge_based, reliability_score, clinical, sd )

almost_both_decipher <- almost_both_decipher %>%
  rename(cnvscore_score_both = .pred_pathogenic) %>%
  select(id, cnvscore_score_both, reliability_score, clinical, sd )

almost_clinvarind_both <- almost_clinvarind_both %>%
  rename(cnvscore_score_both = .pred_pathogenic) %>%
  select(id, cnvscore_score_both, reliability_score, clinical, sd )


plot_sampling1 <- plot_rates3(almost_clinvar, 
                              res_df_both,
                              res_df_knowledge_based,
                              tag = 'ClinVar', 
                              tag2 = '', 
                              input_patho = 0.5, 
                              input_benign = 0.5, 
                              fixed_interval = tmp33,
                              limit_y = 0.7,
                              n_split = tmp22,
                              tmp_sub = '')


plot_sampling2 <- plot_rates3(almost_decipher, 
                              almost_both_decipher,
                              almost_knowledge_based_decipher,
                              tag = 'DECIPHER', 
                              tag2 = '', 
                              input_patho = 0.5, 
                              input_benign = 0.5, 
                              fixed_interval = tmp33,
                              limit_y = 0.6,
                              n_split = tmp22,
                              tmp_sub = '')

plot_sampling3 <- plot_rates3(almost_clinvarind, 
                              almost_clinvarind_both,
                              almost_clinvarind_knowledge_based,
                              tag = 'Clinvar independent', 
                              tag2 = '', 
                              input_patho = 0.5, 
                              input_benign = 0.5, 
                              fixed_interval = tmp33,
                              limit_y = 0.6,
                              n_split = tmp22,
                              tmp_sub = '')

p1_sampling <- plot_sampling1[[1]] + plot_sampling1[[2]] + plot_sampling1[[3]] + plot_sampling2[[1]] + 
  plot_sampling2[[2]] + plot_sampling2[[3]] + plot_sampling3[[1]] + plot_sampling3[[2]] + plot_sampling3[[3]] + plot_layout(nrow = 3)


plot_prescore1 <- plot_rates4(almost_clinvar, 
                              res_df_both,
                              # res_df_combined,
                              tag = 'ClinVar', 
                              tag2 = 'Unbiased features', 
                              input_patho = 0.5, 
                              input_benign = 0.5, 
                              fixed_interval = tmp33,
                              limit_y = tmp11,
                              n_split_score = 5,
                              n_split = tmp22,
                              tmp_sub = paste0('Score: equally-sized bins (n = ', split_score, ') ', '\n',
                                               'Residuals: equally-sized bins (n = ', split_residuals, ')', '\n',
                                               'Representation: ', 'equally-sized', ' bins ', '(n = ',5, ')' ))

plot_prescore2 <- plot_rates4(almost_decipher, 
                              almost_both_decipher,
                              # almost_combined_decipher,
                              tag = 'DECIPHER', 
                              tag2 = 'Unbiased features', 
                              input_patho = 0.5, 
                              input_benign = 0.5, 
                              fixed_interval = tmp33,
                              n_split_score = 5,
                              limit_y = 0.25,
                              n_split = tmp22,
                              tmp_sub = paste0('Score: equally-sized bins (n = ', split_score, ') ', '\n',
                                               'Residuals: equally-sized bins (n = ', split_residuals, ')', '\n',
                                               'Representation: ', tmp33, ' bins ', '(n = ', 5, ')'))


plot_prescore3 <- plot_rates4(almost_clinvarind, 
                              almost_clinvarind_both,
                              # almost_combined_decipher,
                              tag = 'Clinvar independent', 
                              tag2 = '', 
                              input_patho = 0.5, 
                              input_benign = 0.5, 
                              fixed_interval = tmp33,
                              n_split_score = 5,
                              limit_y = 0.25,
                              n_split = tmp22,
                              tmp_sub = paste0('Score: equally-sized bins (n = ', split_score, ') ', '\n',
                                               'Residuals: equally-sized bins (n = ', split_residuals, ')', '\n',
                                               'Representation: ', tmp33, ' bins ', '(n = ', 5, ')'))


p1_prescore <- plot_prescore1[[1]] + plot_prescore1[[2]] + plot_prescore2[[1]] + 
  plot_prescore2[[2]] + plot_layout(nrow = 2)

plot_plot_rates4 <- plot_prescore1[[1]] + plot_prescore1[[3]] + plot_prescore1[[5]] + plot_layout(nrow = 2)

plot_plot_rates42 <- plot_prescore2[[1]] + plot_prescore2[[3]] + plot_prescore2[[5]] + plot_layout(nrow = 2)

plot_plot_rates43 <- plot_prescore3[[1]] + plot_prescore3[[3]] + plot_prescore3[[5]] + plot_layout(nrow = 2)





decipher_yhat <- predict(trendline_clinvar, 
                         plot_rates_decipher_unbiased %>%
                           rename(.pred_pathogenic = cnvscore_score))

p4_both <- plot_rates_decipher_unbiased %>%
  mutate(yhat = decipher_yhat) %>%
  mutate(gam_residuals = sd - yhat) %>% 
  mutate(source = 'decipher') %>%
  select(cnvscore_score, gam_residuals, source) %>%
  bind_rows(
    res_df %>%
      mutate(cnvscore_score = .pred_pathogenic) %>%
      mutate(source = 'clinvar') %>%
      select(cnvscore_score, gam_residuals, source)
  ) %>%
  mutate(score_interval = ntile(cnvscore_score, n = split_score)) %>%
  ggplot(aes(as.factor(score_interval), gam_residuals)) +
  geom_boxplot(aes(fill = source)) +
  labs(x = 'Equally-sized bins - CNVscore score', y = 'SD residuals',
       title = 'DECIPHER and ClinVar CNVs residuals using the linear regression built on ClinVar data ') +
  theme_minimal()


# careful if you changed the version
#there are extra four
tmp_id_decipher_clinical <- read_tsv('/data-cbl/decipher_data/decipher-cnvs-grch37-2020-12-06.txt', skip = 1) %>% 
  rename(id = `# patient_id`, clinical2 = pathogenicity) %>%
  select(id, clinical2) %>%
  filter(clinical2 %in% c('Pathogenic', 'Unknown', 'Likely pathogenic')) %>%
  mutate(id = as.character(id))


p5_both <- output_decipher_deletion %>% 
  select(id.x, id) %>%
  left_join(tmp_id_decipher_clinical, by = c('id.x' = 'id')) %>%
  distinct() %>%
  mutate(clinical = 'pathogenic') %>%
  left_join(
    plot_rates_decipher_unbiased %>%
      mutate(yhat = decipher_yhat) %>%
      mutate(gam_residuals = sd - yhat) %>% 
      mutate(source = 'decipher') %>%
      select(id, cnvscore_score, gam_residuals, source),
    by = 'id'
  ) %>%
  # there are extra four
  group_by(id) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  select(-id.x, -id) %>%
  mutate(clinical2 = if_else(is.na(clinical2), 'Benign', clinical2)) %>%
  bind_rows(
    res_df %>%
      mutate(cnvscore_score = .pred_pathogenic) %>%
      mutate(source = 'clinvar') %>%
      select(cnvscore_score, gam_residuals, source, clinical)
  ) %>% 
  mutate(clinical2 = case_when(
    is.na(clinical2) & clinical == 'pathogenic' ~ 'pathogenic',
    is.na(clinical2) & clinical == 'benign' ~ 'benign',
    TRUE ~ clinical2
  )) %>%
  mutate(tag = paste(source, '-', clinical2)) %>%
  mutate(score_interval = ntile(cnvscore_score, n = split_score)) %>%
  ggplot(aes(as.factor(score_interval), gam_residuals)) +
  geom_boxplot(aes(fill = tag)) +
  labs(x = 'Equally-sized bins - CNVscore score', y = 'SD residuals',
       title = 'DECIPHER and ClinVar CNVs residuals using the linear regression built on ClinVar data ') +
  theme_minimal()


p1_both <- almost_clinvar %>%
  select(reliability_score, cnvscore_score) %>%
  mutate(source = 'clinvar') %>%
  bind_rows(almost_decipher %>%
              select(reliability_score, cnvscore_score) %>%
              mutate(source = 'decipher')) %>%
  mutate(score_interval = ntile(cnvscore_score, n = split_score)) %>%
  count(reliability_score, source) %>%
  group_by(source) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(as.factor(source), perc)) +
  geom_col(aes(fill = as.factor(reliability_score))) +
  scale_y_continuous(label = percent) +
  geom_text(aes(label = paste0(100*round(perc, 2), '% ', '(', n, ')' )),
            size = 5, position = position_stack(vjust = 0.5)) +
  labs(x = 'Source', y = 'Percentage', fill = 'Reliability score',
       title = paste0('Score: equally-sized bins (n = ', split_score, ') ', '\n',
                      'Residuals: equally-sized bins (n = ', split_residuals, ')')) +
  theme_minimal()

p2_both <- almost_clinvar %>%
  select(reliability_score, cnvscore_score) %>%
  mutate(source = 'clinvar') %>%
  bind_rows(almost_decipher %>%
              select(reliability_score, cnvscore_score) %>%
              mutate(source = 'decipher')) %>%
  mutate(score_interval = ntile(cnvscore_score, n = split_score)) %>%
  count(reliability_score, source, score_interval) %>%
  group_by(source, score_interval) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(as.factor(score_interval), perc, group = source)) +
  geom_col(aes(fill = as.factor(reliability_score)), color = 'black') +
  scale_y_continuous(label = percent) +
  # geom_text(aes(label = paste0(100*round(perc, 2), '% ', '(', n, ')' )),
  #           size = 5, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  labs(x = 'Equally-sized bins (CNVscore score)', y = 'Percentage', 
       fill = 'Reliability score') +
  facet_wrap(vars(source))


almost_clinvar %>%
  select(reliability_score, cnvscore_score) %>%
  mutate(source = 'clinvar') %>%
  bind_rows(almost_decipher %>%
              select(reliability_score, cnvscore_score) %>%
              mutate(source = 'decipher')) %>%
  mutate(score_interval = cut_interval(cnvscore_score, n = split_score)) %>%
  count(reliability_score, source, score_interval) %>%
  group_by(source, score_interval) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(as.factor(score_interval), perc, group = source)) +
  geom_col(aes(fill = as.factor(reliability_score)), color = 'black') +
  scale_y_continuous(label = percent) +
  # geom_text(aes(label = paste0(100*round(perc, 2), '% ', '(', n, ')' )),
  #           size = 5, position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  labs(x = 'Equally-sized bins (CNVscore score)', y = 'Percentage', 
       fill = 'Reliability score') +
  facet_wrap(vars(source)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


trendline_clinvar <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
                         data= df_for_uncertainty[[3]])

trendline_decipher <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
                          data= plot_rates_decipher_unbiased %>% 
                            rename(.pred_pathogenic = cnvscore_score)
)

res_df2 <- df_for_uncertainty[[3]] %>%
  select(sd, .pred_pathogenic) %>%
  mutate(source = 'clinvar') %>%
  bind_rows(
    plot_rates_decipher_unbiased %>%
      rename(.pred_pathogenic = cnvscore_score) %>%
      select(sd, .pred_pathogenic) %>%
      mutate(source = 'decipher')
  ) %>%
  mutate(score_interval = ntile(.pred_pathogenic, n = split_score))

p3_both <- res_df2 %>%
  ggplot(aes(.pred_pathogenic, sd)) +
  geom_point(aes(color = source)) +
  geom_smooth(aes(color = source)) +
  labs(x = 'CNVscore score', y = 'SD', color = 'Source') +
  theme_minimal()

p3_both
p1_both + p2_both
p4_both + p5_both + plot_layout(nrow = 2)



# ------------------------------------------------------------------------------
# EXPORT FIGURES
# ------------------------------------------------------------------------------

current_date <- '02_04_22'





ggsave(glue("figures/{current_date}/fig_2a_length_distribution_before_matching.png"),
       figure_length_distribution_before_matching, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_2a_length_distribution_before_matching_variant_class.png"),
       figure_length_distribution_before_matching_variant_class, width = 17, height = 9.6, dpi = 300, units = "in", device='png')


ggsave(glue("figures/{current_date}/fig_1a_stacked_barplot_before_match.png"),
       figure_stacked_barplot, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_1b_stacked_barplot_after_match.png"),
       figure_stacked_barplot2, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_2b_length_distribution_after_matching.png"),
       figure_length_distribution_after_matching, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_2c_length_distribution_after_matching_simplified.png"),
       figure_length_distribution_after_matching2, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_3_qc_cnv_removed.png"),
       figure_qc_cnv_removed, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_4a_enrichment_nobias.png"),
       p_enrich1_nobias, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_4b_enrichment_bias.png"),
       p_enrich2_bias, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp1_coverage.png"),
       figure_coverage, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp2_correlation_across_models.png"),
       figure_correlation, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp3_corr_strvctvre_unbiased_cnvscore.png"),
       structure_corr1 + structure_corr2 + structure_corr3,
       width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp4_comp_gbm_bayesian_models.png"),
       p_comp_gbm_bay_clinvar / p_comp_gbm_bay_decipher, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp5_comp_gbm_bayesian_models2.png"),
       p2_comp_gbm_bay_clinvar, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp6_pca_training.png"),
       p1_pca_training + p2_pca_training, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

# ggsave(glue("figures/{current_date}/fig_supp13_analysis_ood.png"),
#        analysis_ood, width = 17, height = 9.6, dpi = 300, units = "in", device='png')
# 
# ggsave(glue("figures/{current_date}/fig_supp14_analysis2_ood.png"),
#        analysis_ood2, width = 17, height = 9.6, dpi = 300, units = "in", device='png')
# 
# 
# ggsave(glue("figures/{current_date}/fig_supp15_analysis3_ood.png"),
#        analysis_ood3, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp7_p1_uncertainty.png"),
       figure_p1_uncertainty, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp8_p2_uncertainty.png"),
       figure_p2_uncertainty, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp9_app_clinvar.png"),
       p1_app_clinvar + p2_app_clinvar + p3_app_clinvar, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp10_app_decipher.png"),
       p1_app_decipher + p2_app_decipher + p3_app_decipher, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_5a_sd_cnvscore_score.png"),
       p3_rate, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_5b_precision_npv.png"),
       p1_rate, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_5c_subsampling.png"),
       p1_sampling, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp11_reliability_score_cnvscore_clinvar.png"),
       plot_plot_rates4, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp12_reliability_score_cnvscore_decipher.png"),
       plot_plot_rates42, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_supp12_reliability_score_cnvscore_clinvarindependent.png"),
       plot_plot_rates43, width = 17, height = 9.6, dpi = 300, units = "in", device='png')

ggsave(glue("figures/{current_date}/fig_N_comparison_datasets_strvctvre_reliability1.png"),
       plot_comparison_reliability_1, width = 17, height = 9.6, dpi = 300, units = "in", device='png')



# ------------------------------------------------------------------------------
# *********************************END*****************************************
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ANTONIO'S REQUEST REMOT_GW &
# ------------------------------------------------------------------------------



# Training models-------------------------------------------------------------------



# model_bay_del_remot <- rtemis_step1(input_tbl = total_df_deletion,
#                                    tag_variant = 'deletion',
#                                    vector_features = c('max_obs_exp'),
#                                    tag_features <- 'clinical ~ Remot-GW',
#                                    input_prior = 'hs',
#                                    nc = 23,
#                                    input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
#
# model_bay_del_cadd <- rtemis_step1(input_tbl = total_df_deletion,
#                                          tag_variant = 'deletion',
#                                          vector_features = c('max_cadd'),
#                                          tag_features <- 'clinical ~ cadd',
#                                          input_prior = 'hs',
#                                          nc = 23,
#                                          input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
#
# model_bay_del_gerp <- rtemis_step1(input_tbl = total_df_deletion,
#                                    tag_variant = 'deletion',
#                                    vector_features = c('max_gerp'),
#                                    tag_features <- 'clinical ~ gerp',
#                                    input_prior = 'hs',
#                                    nc = 23,
#                                    input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
#
#
# model_bay_del_all_remot <- rtemis_step1(input_tbl = total_df_deletion,
#                                            tag_variant = 'deletion',
#                                            vector_features = vector_no_human_control[!vector_no_human_control == 'max_obs_exp'],
#                                            tag_features <- 'no_human',
#                                            input_prior = 'hs',
#                                            nc = 23,
#                                            input_hyper = tibble(trees = 500, depth = 2, min_n = 1))




#
#
#
# model_bay_del_caddsv <- rtemis_step1(input_tbl = input_caddsv,
#                                      tag_variant = 'deletion',
#                                      vector_features = c('caddsv_score'),
#                                      tag_features <- 'clinical ~ caddsv_score',
#                                      input_prior = 'hs',
#                                      nc = 23,
#                                      input_hyper = tibble(trees = 500, depth = 2, min_n = 1))
#
#
# model_bay_del_all_caddsv <- rtemis_step1(input_tbl = input_caddsv,
#                                    tag_variant = 'deletion',
#                                    vector_features = c(vector_no_human_control, 'caddsv_score'),
#                                    tag_features <- 'all_features + caddsv_score',
#                                    input_prior = 'hs',
#                                    nc = 23,
#                                    input_hyper = tibble(trees = 500, depth = 2, min_n = 1))




# Decipher total dataset-------------------------------------------------------------------


# result_rf_remot <- predict_chrom_aware(model_rf_del_remot, total_df_deletion)
# result_bay_remot <- predict_chrom_aware_rtemis(model_bay_del_remot, total_df_deletion, 'deletion', 'clinical ~ remot_gw')
# result_bay_cadd <- predict_chrom_aware_rtemis(model_bay_del_cadd, total_df_deletion, 'deletion', 'clinical ~ cadd')
# result_bay_gerp <- predict_chrom_aware_rtemis(model_bay_del_gerp, total_df_deletion, 'deletion', 'clinical ~ gerp')
# result_bay_all_remot <- predict_chrom_aware_rtemis(model_bay_del_all_remot, total_df_deletion, 'deletion', 'clinical ~ all but remot_gw')
# result_bay_all <- predict_chrom_aware_rtemis(bayesian_total_del_nohuman, total_df_deletion, 'deletion', 'clinical ~ all')
# result_bay_caddsv <- predict_chrom_aware_rtemis(model_bay_del_caddsv, input_caddsv, 'deletion', 'clinical ~ CADD-SV')
# result_bay_all_caddsv <- predict_chrom_aware_rtemis(model_bay_del_all_caddsv, input_caddsv, 'deletion', 'clinical ~ all + CADD-SV')
# result_bay_phylop100 <- predict_chrom_aware_rtemis(model_bay_del_phylop100, input_caddsv, 'deletion', 'clinical ~ Phylop100')


# Decipher noncoding dataset-------------------------------------------------------------------


# decipher_noncoding <- output_df %>% filter(n_genes == 0) %>% filter(type_variant == 'deletion')
# 
# 
# 
# 
# result_noncoding_remot <- predict_chrom_aware_rtemis(model_bay_del_remot, decipher_noncoding, 'deletion', 'clinical ~ remot_gw')
# result_noncoding_cadd <- predict_chrom_aware_rtemis(model_bay_del_cadd, decipher_noncoding, 'deletion', 'clinical ~ cadd')
# result_noncoding_gerp <- predict_chrom_aware_rtemis(model_bay_del_gerp, decipher_noncoding, 'deletion', 'clinical ~ gerp')
# result_noncoding_all_remot <- predict_chrom_aware_rtemis(model_bay_del_all_remot, decipher_noncoding, 'deletion', 'clinical ~ all but remot_gw')
# result_noncoding_all <- predict_chrom_aware_rtemis(bayesian_total_del_nohuman, decipher_noncoding, 'deletion', 'clinical ~ all')
# 
# 
# plot_test(
#   result_noncoding_remot[[1]],
#   result_noncoding_cadd[[1]],
#   result_noncoding_gerp[[1]],
#   result_noncoding_all_remot[[1]],
#   result_noncoding_all[[1]]
# ) + plot_test(
#   result_noncoding_remot[[2]],
#   result_noncoding_cadd[[2]],
#   result_noncoding_gerp[[2]],
#   result_noncoding_all_remot[[2]],
#   result_noncoding_all[[2]],
#   roc = FALSE
# )


# Decipher enhancer-omim dataset-------------------------------------------------------------------

# 
# result_enh_omim_remot <- predict_chrom_aware_rtemis(model_bay_del_remot, decipher_enh_omim, 'deletion', 'clinical ~ remot_gw')
# result_enh_omim_cadd <- predict_chrom_aware_rtemis(model_bay_del_cadd, decipher_enh_omim, 'deletion', 'clinical ~ cadd')
# result_enh_omim_gerp <- predict_chrom_aware_rtemis(model_bay_del_gerp, decipher_enh_omim, 'deletion', 'clinical ~ gerp')
# result_enh_omim_all_remot <- predict_chrom_aware_rtemis(model_bay_del_all_remot, decipher_enh_omim, 'deletion', 'clinical ~ all but remot_gw')
# result_enh_omim_all <- predict_chrom_aware_rtemis(bayesian_total_del_nohuman, decipher_enh_omim, 'deletion', 'clinical ~ all')
# 
# 
# plot_test(
#   result_enh_omim_remot[[1]],
#   result_enh_omim_cadd[[1]],
#   result_enh_omim_gerp[[1]],
#   result_enh_omim_all_remot[[1]],
#   result_enh_omim_all[[1]]
# ) + plot_test(
#   result_enh_omim_remot[[2]],
#   result_enh_omim_cadd[[2]],
#   result_enh_omim_gerp[[2]],
#   result_enh_omim_all_remot[[2]],
#   result_enh_omim_all[[2]],
#   roc = FALSE
# )


# ------------------------------------------------------------------------------
# HYPERPARAMETER RULEFIT
# ------------------------------------------------------------------------------

# https://xgboost.readthedocs.io/en/latest/parameter.html
# max_dept
# learning_rate
# minimum child weigh
# tibble(nrounds = 10, learning_rate = 0.2, max_dept = 2, min_child_weight = 1)
# nrounds learning_rate max_dept min_child_weight   auc
# 1       10           0.2        2                1 0.845
# 2       28           0.3        2                1 0.845

# 
# list_hyper <- expand.grid(list(trees = seq(100, 200, 5), learn_rate = 0.00368,
#                                tree_depth = 3))
# # list_hyper <- expand.grid(list(nrounds = c(15,25), learning_rate = 0.2, max_dept = 2, min_child_weight  = 1))
# # list_hyper <- expand.grid(list(nrounds = 15, learning_rate = c(0.25, 0.35), max_dept = 2, min_child_weight  = 1))
# # tibble(trees = 161, tree_depth = 3, learn_rate = 0.00368)
# 
# list_hyper <- list_hyper %>% mutate(auc = NA)
# tic()
# for (i in 1:nrow(list_hyper)) {
# 
# print(i)
# model_hyper <- chrom_aware(output_df_deletion_train,
#                              tag_variant = 'deletion',
#                              list_hyper = list_hyper[i,],
#                              tag_formule = 'human_no_control',
#                              model_name = 'rulefit',
#                              formule_model = human_no_control)
# 
# result_hp_tmp <- predict_chrom_aware(model_hyper, output_df_deletion_test, nc = 40)
# 
# 
# list_hyper$auc[i] <- result_hp_tmp[[1]]
# 
# }
# toc()
# 
# list_hyper %>% arrange(desc(auc))
# 
# write_tsv(x = list_hyper, path = 'cnvscore_results/hyper_rf_deletion.tsv')



# ------------------------------------------------------------------------------
# COMPARISON DECIPHER - ClinVar - DECIPHER (only patho) models
# ------------------------------------------------------------------------------

# Deletion CNVs

from_decipher_to_clinvar_human_del <- predict_chrom_aware_rtemis(bayesian_total_del_human, output_clinvar_deletion, 'Trained on DECIPHER (yes_bias) - Predict: ClinVar', '')
from_clinvar_to_decipher_human_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_human, total_df_deletion, 'Trained on ClinVar (yes_bias) - Predict: DECIPHER', '')

from_decipher_to_decipher_human_del <- predict_chrom_aware_rtemis(bayesian_total_del_human, total_df_deletion, 'Trained on DECIPHER (yes_bias) - Predict: DECIPHER', '')
from_clinvar_to_clinvar_human_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_human, output_clinvar_deletion, 'Trained on ClinVar (yes_bias) - Predict: ClinVar', '')


from_decipher_to_clinvar_del <- predict_chrom_aware_rtemis(bayesian_total_del_nohuman, output_clinvar_deletion, 'Trained on DECIPHER - Predict: ClinVar', '')
from_clinvar_to_decipher_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, total_df_deletion, 'Trained on ClinVar - Predict: DECIPHER', '')

from_clinvar_to_clinvar_del <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, output_clinvar_deletion, 'Trained on ClinVar - Predict: ClinVar', '')
from_decipher_to_decipher_del <- predict_chrom_aware_rtemis(bayesian_total_del_nohuman, total_df_deletion, 'Trained on DECIPHER - Predict: DECIPHER', '')


from_caddsv_decipher_del <- read_tsv('rival_cnvscore/cadd_sv/result.bed', skip = 1) %>%
  select(Name, `CADD-SV-score`) %>%
  rename(.pred_pathogenic = `CADD-SV-score`) %>%
  mutate(Name = str_remove(Name, 'GRCh37:')) %>%
  mutate(Name = str_replace(Name, '\\:', '-')) %>%
  separate(Name, into = c('chrom', 'start', 'end'), sep = '-') %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  mutate(start = start + 1) %>%
  rename(caddsv_score = .pred_pathogenic) %>%
  left_join(total_df_deletion, by = c('chrom', 'start', 'end')) %>%
  mutate(caddsv_score = ifelse(is.na(caddsv_score), 0, caddsv_score)) %>%
  rename(.pred_pathogenic = caddsv_score)



from_caddsv_clinvar_del <- read_tsv('rival_cnvscore/cadd_sv/result.bed', skip = 1) %>%
  select(Name, `CADD-SV-score`) %>%
  rename(.pred_pathogenic = `CADD-SV-score`) %>%
  mutate(Name = str_remove(Name, 'GRCh37:')) %>%
  mutate(Name = str_replace(Name, '\\:', '-')) %>%
  separate(Name, into = c('chrom', 'start', 'end'), sep = '-') %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  mutate(start = start + 1) %>%
  rename(caddsv_score = .pred_pathogenic) %>%
  left_join(output_clinvar_deletion, by = c('chrom', 'start', 'end')) %>%
  mutate(caddsv_score = ifelse(is.na(caddsv_score), 0, caddsv_score)) %>%
  rename(.pred_pathogenic = caddsv_score)




create_graph() %>%
  add_node(label = 'DECIPHER') %>%
  add_node(label = 'ClinVar') %>%
  add_node(label = 'CADD-SV') %>%
  add_edge(from = 1, to = 2, edge_aes = edge_aes(color = 'blue', label = from_decipher_to_clinvar_del[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 1, edge_aes = edge_aes(color = 'blue',label = from_clinvar_to_decipher_del[[3]] %>% just_auc)) %>%
  add_edge(from = 1, to = 1, edge_aes = edge_aes(color = 'blue',label = from_decipher_to_decipher_del[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 2, edge_aes = edge_aes(color = 'blue',label = from_clinvar_to_clinvar_del[[3]] %>% just_auc)) %>%
  add_edge(from = 1, to = 2, edge_aes = edge_aes(color = 'red',label = from_decipher_to_clinvar_human_del[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 1, edge_aes = edge_aes(color = 'red',label = from_clinvar_to_decipher_human_del[[3]] %>% just_auc)) %>%
  add_edge(from = 1, to = 1, edge_aes = edge_aes(color = 'red',label = from_decipher_to_decipher_human_del[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 2, edge_aes = edge_aes(color = 'red',label = from_clinvar_to_clinvar_human_del[[3]] %>% just_auc)) %>%
  add_edge(from = 3, to = 1, edge_aes = edge_aes(color = 'orange', label = from_caddsv_decipher_del %>% just_auc)) %>%
  add_edge(from = 3, to = 2, edge_aes = edge_aes(color = 'orange', label = from_caddsv_clinvar_del %>% just_auc)) %>%
  render_graph(layout = "circle", title = 'Deletion CNVs')

plot_decipher_total_deletion_roc <- plot_test(
  result_cadd_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('CADD-SV - deletion - {cadd_del_manolo_auc[1]}')),
  result_tada_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('TADA-SV - deletion - {tada_del_manolo_auc[1]}')),
  result_tadfusion_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('TADfusion - deletion - {tadfusion_del_manolo_auc[1]}')),
  result_structure_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('STRVcture - deletion - {structure_del_manolo_auc[1]}')),
  result_jarvis_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('JARVIS - deletion - {jarvis_del_manolo_auc[1]}')),
  result_gwrvis_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('GWRVIS - deletion - {gwrvis_del_manolo_auc[1]}')),
  result_classifycnv_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('ClassifyCNV - deletion - {classifycnv_del_manolo_auc[1]}')),
  result_xcnv_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('X-CNV - deletion - {xcnv_del_manolo_auc[1]}')),
  result_annotsv_manolo_deletion %>% roc_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('AnnotSV - deletion - {annotsv_del_manolo_auc[1]}')) %>%
    
)

plot_decipher_total_deletion_pr <- plot_test(
  result_cadd_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('CADD-SV - deletion - {cadd_del_manolo_auc[2]}')),
  result_tada_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('TADA-SV - deletion - {tada_del_manolo_auc[2]}')),
  result_tadfusion_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('TADfusion - deletion - {tadfusion_del_manolo_auc[2]}')),
  result_structure_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('STRVcture - deletion - {structure_del_manolo_auc[2]}')),
  result_jarvis_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('JARVIS - deletion - {jarvis_del_manolo_auc[2]}')),
  result_gwrvis_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('GWRVIS - deletion - {gwrvis_del_manolo_auc[2]}')),
  result_classifycnv_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('ClassifyCNV - deletion - {classifycnv_del_manolo_auc[2]}')),
  result_xcnv_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>%
    mutate(tag = glue('X-CNV - deletion - {xcnv_del_manolo_auc[2]}')),
  result_annotsv_manolo_deletion %>% pr_curve(clinical, .pred_pathogenic) %>% 
    mutate(tag = glue('AnnotSV - deletion - {annotsv_del_manolo_auc[2]}')),
  roc = FALSE
)

plot_decipher_total_deletion_roc + plot_decipher_total_deletion_pr



# Duplication CNVs


from_decipher_to_clinvar_human_dup <- predict_chrom_aware_rtemis(bayesian_total_dup_human, output_clinvar_duplication, 'Trained on DECIPHER (yes_bias) - Predict: ClinVar', '')
from_clinvar_to_decipher_human_dup <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_human, total_df_duplication, 'Trained on ClinVar (yes_bias) - Predict: DECIPHER', '')

from_decipher_to_decipher_human_dup <- predict_chrom_aware_rtemis(bayesian_total_dup_human, total_df_duplication, 'Trained on DECIPHER (yes_bias) - Predict: DECIPHER', '')
from_clinvar_to_clinvar_human_dup <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_human, output_clinvar_duplication, 'Trained on ClinVar (yes_bias) - Predict: ClinVar', '')


from_decipher_to_clinvar_dup <- predict_chrom_aware_rtemis(bayesian_total_dup_nohuman, output_clinvar_duplication, 'Trained on DECIPHER - Predict: ClinVar', '')
from_clinvar_to_decipher_dup <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_nohuman, total_df_duplication, 'Trained on ClinVar - Predict: DECIPHER', '')

from_clinvar_to_clinvar_dup <- predict_chrom_aware_rtemis(bayesian_clinvar_dup_nohuman, output_clinvar_duplication, 'Trained on ClinVar - Predict: ClinVar', '')
from_decipher_to_decipher_dup <- predict_chrom_aware_rtemis(bayesian_total_dup_nohuman, total_df_duplication, 'Trained on DECIPHER - Predict: DECIPHER', '')


# from_caddsv_decipher_dup <- read_tsv('rival_cnvscore/cadd_sv/result_dup2.bed', skip = 1) %>%
#   select(Name, `CADD-SV-score`) %>%
#   rename(.pred_pathogenic = `CADD-SV-score`) %>%
#   mutate(Name = str_remove(Name, 'GRCh37:')) %>%
#   mutate(Name = str_replace(Name, '\\:', '-')) %>%
#   separate(Name, into = c('chrom', 'start', 'end'), sep = '-') %>%
#   mutate(start = as.numeric(start), end = as.numeric(end)) %>%
#   mutate(start = start + 1) %>%
#   rename(caddsv_score = .pred_pathogenic) %>%
#   left_join(total_df_duplication, by = c('chrom', 'start', 'end')) %>%
#   mutate(caddsv_score = ifelse(is.na(caddsv_score), 0, caddsv_score)) %>%
#   rename(.pred_pathogenic = caddsv_score)
# 
# 
# 
# from_caddsv_clinvar_dup <- read_tsv('rival_cnvscore/cadd_sv/result_dup.bed', skip = 1) %>%
#   select(Name, `CADD-SV-score`) %>%
#   rename(.pred_pathogenic = `CADD-SV-score`) %>%
#   mutate(Name = str_remove(Name, 'GRCh37:')) %>%
#   mutate(Name = str_replace(Name, '\\:', '-')) %>%
#   separate(Name, into = c('chrom', 'start', 'end'), sep = '-') %>%
#   mutate(start = as.numeric(start), end = as.numeric(end)) %>%
#   mutate(start = start + 1) %>%
#   rename(caddsv_score = .pred_pathogenic) %>%
#   left_join(output_clinvar_duplication, by = c('chrom', 'start', 'end')) %>%
#   mutate(caddsv_score = ifelse(is.na(caddsv_score), 0, caddsv_score)) %>%
#   rename(.pred_pathogenic = caddsv_score)




create_graph() %>%
  add_node(label = 'DECIPHER') %>%
  add_node(label = 'ClinVar') %>%
  add_node(label = 'CADD-SV') %>%
  add_edge(from = 1, to = 2, edge_aes = edge_aes(color = 'blue', label = from_decipher_to_clinvar_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 1, edge_aes = edge_aes(color = 'blue',label = from_clinvar_to_decipher_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 1, to = 1, edge_aes = edge_aes(color = 'blue',label = from_decipher_to_decipher_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 2, edge_aes = edge_aes(color = 'blue',label = from_clinvar_to_clinvar_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 1, to = 2, edge_aes = edge_aes(color = 'red',label = from_decipher_to_clinvar_human_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 1, edge_aes = edge_aes(color = 'red',label = from_clinvar_to_decipher_human_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 1, to = 1, edge_aes = edge_aes(color = 'red',label = from_decipher_to_decipher_human_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 2, to = 2, edge_aes = edge_aes(color = 'red',label = from_clinvar_to_clinvar_human_dup[[3]] %>% just_auc)) %>%
  add_edge(from = 3, to = 1, edge_aes = edge_aes(color = 'orange', label = from_caddsv_decipher_dup %>% just_auc)) %>%
  add_edge(from = 3, to = 2, edge_aes = edge_aes(color = 'orange', label = from_caddsv_clinvar_dup %>% just_auc)) %>%
  render_graph(layout = "circle", title = 'Duplications CNVs')






plot_decipher_total_duplication_roc + plot_decipher_total_duplication_pr





# ------------------------------------------------------------------------------
# CONTAMINATION TRAINING
# ------------------------------------------------------------------------------

# STRVCTURE (XX IDs)------------------------------------------------------------------------------

structure_first_step <- output_clinvar_deletion %>%
  select(chrom, start, end, id, length_cnv) %>%
  bed_intersect(total_df_deletion %>% select(chrom, start, end, id, length_cnv)) %>%
  mutate(coverage = .overlap / length_cnv.x) %>%
  filter(coverage > 0.7) %>%
  rename(first_id_x = id.x, first_id_y = id.y) %>%
  select(starts_with('first'))

structure_second_step <- total_df_deletion %>%
  select(chrom, start, end, id, length_cnv) %>%
  bed_intersect(output_clinvar_deletion %>% select(chrom, start, end, id, length_cnv)) %>%
  mutate(coverage = .overlap / length_cnv.x) %>%
  filter(coverage > 0.7) %>%
  rename(second_id_x = id.x, second_id_y = id.y) %>%
  select(starts_with('second'))

remove_ids_structure <- structure_second_step %>%
  left_join(structure_first_step, by = c('second_id_x' = 'first_id_y')) %>%
  mutate(is_reciprocal = if_else(second_id_y == first_id_x, 'yes', 'no')) %>%
  filter(is_reciprocal == 'yes') %>%
  pull(second_id_x) %>%
  unique()



# X-CNV (34 IDs)------------------------------------------------------------------------------

download.file('https://ndownloader.figstatic.com/files/29234622', destfile = 'dbvar_cnvs.xlsx')


training_xcnv <- read_excel('dbvar_cnvs.xlsx', sheet = 2) %>% 
  # remove CNVs mapping X and Y
  na.omit() %>%
  rename(chrom = Chromosome, start = Start, end = End, variant_type = Type, clinical = Pathogenicity) %>%
  mutate(type_variant = if_else(variant_type == 'loss', 'deletion', 'duplication')) %>% 
  select(-variant_type) %>%
  mutate(clinical = tolower(clinical)) %>%
  mutate(source = 'dbvar') %>%
  mutate(length_cnv = end - start + 1) %>%
  mutate(id = as.integer(row_number())) %>%
  mutate(chrom = as.character(chrom)) 


xcnv_first_step <- training_xcnv %>%
  select(chrom, start, end, id, length_cnv) %>%
  bed_intersect(total_df_deletion %>% select(chrom, start, end, id, length_cnv)) %>%
  mutate(coverage = .overlap / length_cnv.x) %>%
  filter(coverage > 0.7) %>%
  rename(first_id_x = id.x, first_id_y = id.y) %>%
  select(starts_with('first'))

xcnv_second_step <- total_df_deletion %>%
  select(chrom, start, end, id, length_cnv) %>%
  bed_intersect(training_xcnv %>% select(chrom, start, end, id, length_cnv)) %>%
  mutate(coverage = .overlap / length_cnv.x) %>%
  filter(coverage > 0.7) %>%
  rename(second_id_x = id.x, second_id_y = id.y) %>%
  select(starts_with('second'))

remove_ids_xcnv <- xcnv_second_step %>%
  left_join(xcnv_first_step, by = c('second_id_x' = 'first_id_y')) %>%
  mutate(is_reciprocal = if_else(second_id_y == first_id_x, 'yes', 'no')) %>%
  filter(is_reciprocal == 'yes') %>%
  pull(second_id_x) %>%
  unique()





## JUST IN CASE



# Figure N - chromosome_distribution------------------------------------------------------------------------------

vector_clinical <- input_check_cnv %>% select(clinical) %>% distinct() %>% unique() %>% pull()


result_coverage <- tibble(chrom = c(1:22, 'X'))

for (i in 1:length(vector_clinical)) {
  
  
  result_coverage <- coord_chrom_hg19 %>% rename(end = length) %>% mutate(start = 1) %>%
    bed_coverage(input_check_cnv %>% filter(clinical == vector_clinical[i])) %>%
    select(chrom, .frac) %>%
    rename(!!vector_clinical[i] := .frac) %>%
    right_join(result_coverage, by = 'chrom')
}

result_coverage %>%
  pivot_longer(cols = -chrom, names_to = 'clinical', values_to = 'value') %>%
  mutate(chrom = factor(chrom, levels = c(1:22, 'X'))) %>%
  ggplot(aes(chrom, value)) +
  geom_col(color = 'black', aes(fill = clinical)) +
  scale_y_continuous(label = percent) +
  facet_wrap(vars(clinical)) +
  theme_minimal()


# Figure N - chromosome_distribution------------------------------------------------------------------------------

input_check_cnv %>% 
  mutate(chrom = factor(chrom, levels = c(1:22, 'X'))) %>%
  filter(variant_class == 'deletion') %>%
  count(chrom, source, variant_class) %>%
  group_by(source, variant_class) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(chrom, perc)) +
  geom_col(aes(fill = source), color = 'black') +
  scale_y_continuous(labels = percent) +
  facet_wrap(vars(source)) +
  labs(x = 'Chromosomes', y = 'Percentage', fill = 'CNVs source', title = 'Percentage of CNV deletions across chromosome and source') +
  theme_minimal() +
  input_check_cnv %>% 
  bind_rows(clinvar_cnvs_hg37 %>% select(chrom, source, variant_class)) %>%
  mutate(chrom = factor(chrom, levels = c(1:22, 'X'))) %>%
  filter(variant_class == 'duplication') %>%
  count(chrom, source, variant_class) %>%
  group_by(source, variant_class) %>%
  mutate(perc = n / sum(n)) %>%
  ggplot(aes(chrom, perc)) +
  geom_col(aes(fill = source), color = 'black') +
  scale_y_continuous(labels = percent) +
  facet_wrap(vars(source)) +
  labs(x = 'Chromosomes', y = 'Percentage', fill = 'CNVs source', title = 'Percentage of CNV duplications across chromosome and source') +
  theme_minimal()


# Figure N - cytoband distribution------------------------------------------------------------------------------

cytoband_tmp_object <- input_check_cnv %>% bind_rows(clinvar_cnvs_hg37 %>% select(chrom, start, end, source, variant_class))
cytoband_result_deletions <- coord_cytobands %>% mutate(name_cytoband = paste0(chrom, Name)) %>% select(name_cytoband)
cytoband_result_duplications <- coord_cytobands %>% mutate(name_cytoband = paste0(chrom, Name)) %>% select(name_cytoband)

name_sources <- cytoband_tmp_object %>% count(source) %>% pull(source)


# Deletions

for (i in 1:length(name_sources)) {
  
  
  cytoband_result_deletions <- coord_cytobands %>% mutate(name_cytoband = paste0(chrom, Name)) %>%
    bed_intersect(cytoband_tmp_object %>% filter(source == name_sources[i] & variant_class == 'deletion')) %>%
    count(name_cytoband.x) %>%
    rename(name_cytoband = name_cytoband.x) %>%
    rename(!!name_sources[i] := n) %>% right_join(cytoband_result_deletions, by = 'name_cytoband')
  
  cytoband_result_duplications <- coord_cytobands %>% mutate(name_cytoband = paste0(chrom, Name)) %>%
    bed_intersect(cytoband_tmp_object %>% filter(source == name_sources[i] & variant_class == 'duplication')) %>%
    count(name_cytoband.x) %>%
    rename(name_cytoband = name_cytoband.x) %>%
    rename(!!name_sources[i] := n) %>% right_join(cytoband_result_duplications, by = 'name_cytoband')
}

cytoband_result_deletions %>% 
  replace(is.na(.), 0) %>% 
  pivot_longer(cols = -name_cytoband, names_to = 'source', values_to = 'value') %>%
  group_by(source) %>%
  mutate(perc = value / sum(value)) %>%
  arrange(desc(perc), .by_group = TRUE) %>%
  mutate(cum_perc = cumsum(perc)) %>%
  mutate(n_row = row_number()) %>%
  ggplot(aes(n_row, cum_perc)) +
  geom_line(aes(color = source, group = source), size = 2) +
  scale_y_continuous(label = percent) +
  # expand_limits(x = 1) %>%
  labs(x = 'Nº of cytobands', y = 'Cumulative percentage of hits', 
       color = 'CNV sources', title = 'Cumulative percentage of hits across cytobands - Deletions') +
  theme_minimal() +
  
  cytoband_result_duplications %>% 
  replace(is.na(.), 0) %>% 
  pivot_longer(cols = -name_cytoband, names_to = 'source', values_to = 'value') %>%
  group_by(source) %>%
  mutate(perc = value / sum(value)) %>%
  arrange(desc(perc), .by_group = TRUE) %>%
  mutate(cum_perc = cumsum(perc)) %>%
  mutate(n_row = row_number()) %>%
  ggplot(aes(n_row, cum_perc)) +
  geom_line(aes(color = source, group = source), size = 2) +
  scale_y_continuous(label = percent) +
  # expand_limits(x = 1) %>%
  labs(x = 'Nº of cytobands', y = 'Cumulative percentage of hits', 
       color = 'CNV sources', title = 'Cumulative percentage of hits across cytobands - Duplications') +
  theme_minimal()
















# 
# # Figure 1 - Stacked barplot------------------------------------------------------------------------------
# 
# annotation_stacked <- bind_rows(clinvar_stacked_plot %>% select(chrom, start, end, variant_class, clinical, source),
#           decipher_stacked_plot %>% select(chrom, start, end, variant_class, clinical, source),
#           input_check_cnv %>% filter(clinical == 'benign') %>% select(chrom, start, end, variant_class, clinical, source)) %>%
#   mutate(length_cnv = end - start + 1,
#          id_tmp = row_number()) %>%
#   check_cnv_v2(factor_clinical = FALSE)
# 
# tmp_only_omim_outside <- annotation_stacked %>%
#   filter(omim == 0) %>%
#   select(id, chrom, start, end) %>%
#   bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
#   filter(.overlap > 0) %>%
#   pull(id.x) %>%
#   unique()
# 
# 
# df_genes_inside <- annotation_stacked %>%
#   filter(omim == 1) %>%
#   bed_intersect(hgcn_genes) %>%
#   select(id.x, gene.y) %>%
#   rename(id = id.x, gene_inside = gene.y)
# 
# df_genes_outside <- annotation_stacked %>%
#   filter(omim == 1) %>%
#   bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
#   select(id.x, gene.y) %>%
#   rename(id = id.x, gene_outside = gene.y)
# 
# unique_id_genes <- unique(df_genes_inside$id)
# 
# result_ids <- c()
# 
# for (i in 1:length(unique_id_genes)) {
#   
#   how_many_genes_outside_no_inside <- df_genes_outside %>% 
#     filter(id == unique_id_genes[i] ) %>%
#     filter(!gene_outside %in% df_genes_inside[df_genes_inside$id == unique_id_genes[i],]$gene_inside) %>%
#     nrow()
#   
#   if (how_many_genes_outside_no_inside > 0)  result_ids <- c(unique_id_genes[i],result_ids)
#   
#   
# }
# 
# 
# 
# tmp_inside_outside_omim <- annotation_stacked %>%
#   filter(omim == 1) %>%
#   select(id, chrom, start, end) %>%
#   bed_intersect(df_enhancers %>% filter(gene %in% hgcn_genes$gene[hgcn_genes$omim == 'Yes'])) %>%
#   filter(.overlap > 0) %>%
#   pull(id.x) %>%
#   unique()
# 
# 
# remove_small_clinicals <- annotation_stacked %>% 
#   mutate(tag2 = paste(source, '-', clinical)) %>% 
#   count(tag2) %>% 
#   arrange(desc(n)) %>%
#   filter(n < 22) %>%
#   pull(tag2)
# 
# show_n <- annotation_stacked %>% 
#   mutate(tag2 = paste(source, '-', clinical)) %>% 
#   count(tag2)
# 
# 
# 
# 
# figure_stacked_barplot <- annotation_stacked %>% 
#   mutate(target_omim_only_outside = if_else(id %in% tmp_only_omim_outside, 'yes', 'no')) %>%
#   mutate(tag = case_when(
#     id %in% result_ids ~ 'OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV',
#     omim == 0 & target_omim_only_outside == 'yes' ~ 'Enhancers targeted of OMIM genes outside the CNV',
#     n_genes > 0 & omim == 0 ~ 'Other protein-coding genes targeted',
#     n_genes == 0 ~ 'No protein-coding genes targeted',
#     omim > 0 ~ 'OMIM genes targeted'
#   )) %>% 
#   # mutate(source = if_else(source == 'clinvar', 'ClinVar', source)) %>%
#   # mutate(source = if_else(source == 'decipher', 'DECIPHER', source)) %>%
#   mutate(tag2 = paste(source, '-', clinical)) %>%
#   filter(!tag2 %in% remove_small_clinicals) %>%
#   left_join(show_n, by = 'tag2') %>%
#   mutate(tag2 = glue('{tag2} ({n})')) %>%
#   count(tag, tag2, source) %>%
#   group_by(tag2, source) %>%
#   mutate(perc = n / sum(n)) %>%
#   mutate(tag = factor(tag, c('No protein-coding genes targeted', 'Other protein-coding genes targeted', 'Enhancers targeted of OMIM genes outside the CNV','OMIM genes targeted + enhancers of \nadditional OMIM genes outside the CNV', 'OMIM genes targeted'))) %>%
#   mutate(tag3 = case_when(
#     str_detect(source, 'decipher_control') ~ 'General population',
#     str_detect(source, 'decipher') ~ 'DECIPHER',
#     str_detect(source,'clinvar') ~ 'ClinVar',
#     TRUE ~ 'General population'
#   )) %>%
#   # mutate(tag2 = str_remove(tag2, 'clinvar - ')) %>%
#   # mutate(tag2 = str_remove(tag2, 'decipher - ')) %>%
#   mutate(tag2 = str_remove(tag2, ' - benign')) %>%
#   mutate(tag2 = str_to_title(tag2)) %>%
#   filter((tag3 == 'ClinVar' & str_detect(tag2, 'Pathogenic')) | 
#            (tag3 == 'DECIPHER' & str_detect(tag2, 'Pathogenic |Unknown')) |
#               (tag3 == 'General population')) %>%
#   mutate(tag4 = if_else(tag3 == 'General population', 'Benign CNVs', 'Pathogenic CNVs')) %>%
#   mutate(tag4 = factor(tag4, levels = c('Pathogenic CNVs', 'Benign CNVs'))) %>%
#   mutate(tag5 = case_when(
#     str_detect(tag2, 'Clinvar - Pathogenic') ~ 1,
#     str_detect(tag2, 'Clinvar - Likely Pathogenic') ~ 3,
#     str_detect(tag2, 'Decipher - Pathogenic') ~ 2,
#     str_detect(tag2, 'Decipher - Likely Pathogenic') ~ 3,
#     str_detect(tag2, 'Decipher - Unknown') ~ 4,
#     str_detect(tag2, 'Decipher_control') ~ 5,
#     str_detect(tag2, 'Dgv') ~ 6,
#     str_detect(tag2, 'Gnomad_v2.1') ~ 9,
#     str_detect(tag2, 'Chaisson') ~ 10,
#     str_detect(tag2, 'Audano') ~ 8,
#     str_detect(tag2, 'Beyter') ~ 7
#   )) %>%
#   ungroup() %>%
#   mutate(tag2 = factor(tag2)) %>%
#   mutate(tag2 = fct_reorder(tag2, tag5)) %>%
#   ggplot(aes(tag2, perc, group = tag)) +
#     geom_col(aes(fill = tag), color = 'black') +
#   scale_y_continuous(label = percent) +
#   geom_label(aes(label = paste0(100*round(perc, 2), '%' , ' (', n, ')')), size = 3, position = position_stack(vjust = 0.5)) +
#   scale_fill_manual(values = c(hue_pal()(4)[2], '#Ecd26d', '#Eca9a9', hue_pal()(1), '#DD2E44')) +
#   facet_wrap(vars(tag4), scale = 'free_x') +
#     theme_minimal() +
#   labs(fill = 'Category', y = 'Percentage') +
#   theme(axis.text.x=element_text(angle = -45, hjust = 0))
# 
# 
# 
# (glue("figures/{current_date}/stacked_barplot.png"),
#        figure_stacked_barplot, width = 17, height = 9.6, dpi = 300, units = "in", device='png')
# 
# 
#   


## DECIPHER

# Unbiased model

# distances_tbl <- get_distances(
#   training_tbl = output_clinvar_deletion,
#   input_tbl = output_decipher_deletion,
#   formule_input = human_no_control, 
#   threshold = 0.95
# )
# 
# 
# comp_structure_decipher2 <- comp_structure_decipher1 %>%
#   left_join(distances_tbl %>% select(id, dist_pathogenic, dist_benign), by = 'id')
# 
# 
# tmp1 <- read_tsv('structure_tmp/output_decipher.bed', 
#                  col_names = c('chrom', 'start', 'end')) %>%
#   mutate(start = start + 1) %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   rename(.pred_pathogenic = X5) %>%
#   mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
#   mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
#   rename(structure_score = .pred_pathogenic) %>%
#   select(chrom, start, end, structure_score) %>%
#   GRanges()
# 
# seqlevelsStyle(tmp1) = "UCSC" 
# 
# tmp2 <- liftOver(tmp1, from_hg38_to_hg19) %>%
#   as_tibble() %>%
#   rename(chrom = seqnames) %>%
#   select(chrom, start, end, structure_score) %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   left_join(comp_structure_decipher2) %>%
#   na.omit() %>%
#   select(-c('chrom', 'start', 'end'))
# 
# 
# 
# tpr_decipher_dist <- tmp2 %>% 
#   filter(clinical == 'pathogenic') %>%
#   mutate(tile_dist = ntile(dist_pathogenic, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TPR) - Unbiased model',
#        x = 'Pathogenic distance (Percentile 1-10)', y = 'TPR (%)')
# 
# 
# tnr_decipher_dist <- tmp2 %>% 
#   filter(clinical == 'benign') %>%
#   mutate(tile_dist = ntile(dist_benign, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TNR) - Unbiased model',
#        x = 'Benign distance (Percentile 1-10)', y = 'TNR (%)')
# 
# 
# tpr_decipher_sd <- tmp2 %>% 
#   filter(clinical == 'pathogenic') %>%
#   mutate(tile_dist = ntile(sd, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TPR) - Unbiased model',
#        x = 'SD (Percentile 1-10)', y = 'TPR (%)')
# 
# 
# tnr_decipher_sd <- tmp2 %>% 
#   filter(clinical == 'benign') %>%
#   mutate(tile_dist = ntile(sd, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TNR) - Unbiased model',
#        x = 'SD (Percentile 1-10)', y = 'TNR (%)')
# 
# 
# # Total model 
# 
# distances_tbl_both <- get_distances(
#   training_tbl = output_clinvar_deletion,
#   input_tbl = output_decipher_deletion,
#   formule_input = formule_total, 
#   threshold = 0.95
# )
# 
# 
# comp_structure_both_decipher <- comp_structure_both_decipher1 %>%
#   left_join(distances_tbl_both %>% select(id, dist_pathogenic, dist_benign), by = 'id')
# 
# 
# tmp1 <- read_tsv('structure_tmp/output_decipher.bed', 
#                  col_names = c('chrom', 'start', 'end')) %>%
#   mutate(start = start + 1) %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   rename(.pred_pathogenic = X5) %>%
#   mutate(.pred_pathogenic = ifelse(.pred_pathogenic == 'not_exonic', 0, .pred_pathogenic)) %>%
#   mutate(.pred_pathogenic = as.numeric(.pred_pathogenic)) %>%
#   rename(structure_score = .pred_pathogenic) %>%
#   select(chrom, start, end, structure_score) %>%
#   GRanges()
# 
# seqlevelsStyle(tmp1) = "UCSC" 
# 
# tmp2 <- liftOver(tmp1, from_hg38_to_hg19) %>%
#   as_tibble() %>%
#   rename(chrom = seqnames) %>%
#   select(chrom, start, end, structure_score) %>%
#   mutate(chrom = str_remove(chrom, 'chr')) %>%
#   left_join(comp_structure_both_decipher) %>%
#   na.omit() %>%
#   select(-c('chrom', 'start', 'end'))
# 
# 
# tpr_decipher_both_dist <- tmp2 %>% 
#   filter(clinical == 'pathogenic') %>%
#   mutate(tile_dist = ntile(dist_pathogenic, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TPR) - Total model',
#        x = 'Pathogenic distance (Percentile 1-10)', y = 'TPR (%)')
# 
# tnr_decipher_both_dist <- tmp2 %>% 
#   filter(clinical == 'benign') %>%
#   mutate(tile_dist = ntile(dist_benign, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TNR) - Total model',
#        x = 'Benign distance (Percentile 1-10)', y = 'TNR (%)')
# 
# 
# tpr_decipher_both_sd <- tmp2 %>% 
#   filter(clinical == 'pathogenic') %>%
#   mutate(tile_dist = ntile(sd, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TPR) - Total model',
#        x = 'SD (Percentile 1-10)', y = 'TPR (%)')
# 
# tmp2 %>% 
#   filter(clinical == 'benign') %>%
#   mutate(tile_dist = ntile(sd, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(tile_dist, tool) %>%
#   count(value) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   ggplot(aes(tile_dist, perc)) +
#   geom_line(aes(group = tool)) +
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   scale_y_continuous(label = percent) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset (TNR) - Total model',
#        x = 'SD (Percentile 1-10)', y = 'TNR (%)')

# 
# tpr_decipher_dist + 
#   tnr_decipher_dist + 
#   tpr_decipher_both_dist + 
#   tnr_decipher_both_dist + plot_layout(nrow = 2)
# 
# tpr_decipher_sd + 
#   tnr_decipher_sd + 
#   tpr_decipher_both_sd + 
#   tnr_decipher_both_sd + plot_layout(nrow = 2)

# 
# 
# tmp2 %>%
#   filter(clinical == 'pathogenic') %>%
#   mutate(tile_dist = ntile(dist_pathogenic, 10)) %>%
#   mutate(structure_score = if_else(structure_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score >= 0.5, 'pathogenic', 'benign')) %>%
#   mutate(structure_score = if_else(structure_score == clinical, 'ok', 'wrong')) %>%
#   mutate(cnvscore_score = if_else(cnvscore_score == clinical, 'ok', 'wrong')) %>%
#   select(structure_score, cnvscore_score, tile_dist) %>%
#   pivot_longer(-c(tile_dist), names_to = 'tool', values_to = 'value') %>%
#   group_by(value, tool) %>%
#   count(tile_dist, value, tool) %>%
#   mutate(perc = n / sum(n)) %>%
#   ungroup() %>%
#   filter(value == 'ok') %>%
#   group_by(tool) %>%
#   mutate(total_n = cumsum(n)) %>%
#   mutate(total_n / sum(total_n)) %>%
#   select(tile_dist, tool, total_n) %>%
#   ggplot(aes(tile_dist, total_n)) +
#   geom_line(aes(group = tool)) +
# 
#   geom_segment(aes(x=0,
#                    xend=10,
#                    y= nrow(tmp2[tmp2$clinical == 'pathogenic',]),
#                    yend= nrow(tmp2[tmp2$clinical == 'pathogenic',])
#                    ),linetype='dotted', col = 'red', show.legend = FALSE) +
# 
# 
#   annotate("text", x = 4, y = nrow(tmp2[tmp2$clinical == 'pathogenic',]),
#            label = paste(nrow(tmp2[tmp2$clinical == 'pathogenic',]), ' Pathogenic CNVs'), vjust = -0.2) +
# 
#   geom_point(aes(fill = tool), shape = 21, color = 'black') +
#   coord_cartesian(ylim = c(0,nrow(tmp2[tmp2$clinical == 'pathogenic',]) + 200)) +
#   # scale_y_continuous(label = percent) +
#   scale_x_continuous(breaks = pretty_breaks()) +
#   theme_minimal() +
#   labs(title = 'ClinVar dataset (TPR) - Unbiased model',
#        x = 'Pathogenic distance (Percentile 1-10)', y = 'Cumulative Pathogenic CNVs (%)')





# res_df_human %>% select(sd) %>% mutate(model = 'human') %>% bind_rows(
# res_df %>% select(sd) %>% mutate(model = 'unbiased')
# ) %>% bind_rows(
# res_df_both %>% select(sd) %>% mutate(model = 'both_unbiased_human')
# ) %>%
#   ggplot(aes(sd, model)) +
#   geom_density_ridges(aes(fill = model), alpha = 0.4) + theme_minimal()

# check_noprotein_match_del <- input_check_cnv_del_pathogenic_clinvar %>%
#   bind_rows(input_check_cnv_del_benign) %>%
#   bed_intersect(hgcn_genes %>% 
#                   select(chrom, start, end),
#                 suffix = c('', '.y')) %>%
#   pull(id_tmp) %>% unique() 
# 
# check_noprotein_match_del <- input_check_cnv_del_pathogenic_clinvar %>%
#   bind_rows(input_check_cnv_del_benign) %>%
#   filter(!id_tmp %in% check_noprotein_match_del)
# 
# check_noprotein_matched_del <- matching_length(bin_length = 100, check_noprotein_match_del)
# Reliability score ------------------------------------------------------------------------------


# result_clinvar_del_realscore <- result_clinvar_del %>%
#   left_join(res_df %>% select(id, reliability_score, clinical), by = 'id') %>%
#   filter(str_detect(tag, 'xcnv|tada|strvctvre|cadd|bayesian_unbiased')) %>%
#   group_by(tag, reliability_score) %>%
#   roc_auc(clinical, .pred_pathogenic) %>%
#   rename(model = tag, auroc = .estimate) %>%
#   select(model, reliability_score, auroc)
# 
# realscore_decipher_del <- get_reliability_score(ref_quantiles, ref_sd_decipher_del[[3]])
# 
# result_decipher_del_realscore <- result_decipher_del %>%
#   left_join(realscore_decipher_del %>% 
#               select(id, reliability_score, clinical), by = 'id') %>%
#   filter(str_detect(tag, 'xcnv|tada|jarvis|gwrvis|strvctvre|cadd|bayesian_unbiased')) %>%
#   group_by(tag, reliability_score) %>%
#   roc_auc(clinical, .pred_pathogenic) %>%
#   rename(model = tag, auroc = .estimate) %>%
#   select(model, reliability_score, auroc)
# 
# p1_realscore <- result_clinvar_del_realscore %>%
#   ggplot(aes(reliability_score, auroc)) +
#   geom_point() +
#   geom_path(aes(group = model, color = model)) +
#   theme_minimal() +
#   labs(title = 'ClinVar dataset')
# 
# p2_realscore <- result_decipher_del_realscore %>%
#   ggplot(aes(reliability_score, auroc)) +
#   geom_point() +
#   geom_path(aes(group = model, color = model)) +
#   theme_minimal() +
#   labs(title = 'DECIPHER dataset')
# 
# p1_realscore + p2_realscore

# 
# 
# trendline_clinvar_human <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
#                               data= tmp_human_clinvar_prescore[[3]])
# 
# 
# 
# res_df_knowledge_based <- tmp_human_clinvar_prescore[[3]] %>%
#   mutate(gam_residuals = residuals(trendline_clinvar_human)) %>%
#   mutate(score_interval = ntile(.pred_pathogenic, n = split_score)) %>%
#   group_by(score_interval) %>%
#   mutate(reliability_score = ntile(gam_residuals, split_residuals)) %>%
#   ungroup() %>%
#   rename(cnvscore_score_knowledge_based = .pred_pathogenic) %>%
#   select(id, cnvscore_score_knowledge_based, reliability_score, clinical, sd, score_interval)
# 
# trendline_clinvar_both <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
#                               data= tmp_both_clinvar_prescore[[3]])
# 
# res_df_both <- tmp_both_clinvar_prescore[[3]] %>%
#   mutate(gam_residuals = residuals(trendline_clinvar_both)) %>%
#   mutate(score_interval = ntile(.pred_pathogenic, n = split_score)) %>%
#   group_by(score_interval) %>%
#   mutate(reliability_score = ntile(gam_residuals, split_residuals)) %>%
#   ungroup() %>%
#   rename(cnvscore_score_both = .pred_pathogenic) %>%
#   select(id, cnvscore_score_both, reliability_score, clinical, sd, score_interval)

# res_df_combined <- df_for_uncertainty[[3]] %>%
#   mutate(cnvscore_score_combined = .pred_pathogenic) %>%
#   mutate(gam_residuals = residuals(trendline_clinvar)) %>%
#   mutate(score_interval = ntile(.pred_pathogenic, n = split_score)) %>%
#   group_by(score_interval) %>%
#   mutate(sd_residuals = ntile(-gam_residuals, 10) / 10) %>%
#   ungroup() %>%
#   mutate(rank_score = ntile(.pred_pathogenic, n = 3) / 3) %>%
#   mutate(reliability_score = if_else(.pred_pathogenic >= 0.5,
#                                      rank_score + sd_residuals,
#                                      rank_score - sd_residuals)) 
# ggplot(aes(as.factor(sd_residuals), combined_score)) + geom_boxplot()
# ggplot(aes(as.factor(rank_score), combined_score)) + geom_boxplot()
# res_df_combined %>% roc_auc(clinical, reliability_score)
# ClinVar special V ------------------------------------------------------------------------------

# result_clinvar_setting_5 <- get_results(clinvar_setting_5, 'DEL', 'ClinVar - patho (>0 omim) - benign (>0 omim)')
# ClinVar OMIM------------------------------------------------------------------------------

# result_clinvar_del_omim <- get_results(output_clinvar_omim_deletion, 'DEL', 'ClinVar - OMIM genes')

# tibble(name = c('Supp Table 1. Reference variants databases',
#                 'Supp Table 2. Summary of the tools used in benchmarks',
#                 'Supp Table 3. List of unbiased featues',
#                 'Supp Table 4. List of knowledge-based featues',
#                 # 'Supp Table 5. Nº CNVs across chromosomes - Training dataset (duplications)',
#                 'Supp Table 5. Nº pathogenic and benign CNVs across benchmark datasets',
#                 'Supp Table 6. Results - ClinVar DEL',
#                 'Supp Table 7. Results - ClinVar DEL - Reliability scores',
#                 'Supp Table 8. Results - ClinVar DUP',
#                 'Supp Table 9. Results - ClinVar Ind. (First evaluation since Jan 2021)',
#                 'Supp Table 10. Results - DECIPHER DEL',
#                 'Supp Table 11. Results - DECIPHER DEL - Reliability scores',
#                 'Supp Table 12. Results - DECIPHER DUP',
#                 'Supp Table 13. Results - DECIPHER Scenario #1 (Pathogenic CNVs mapping 0 OMIM genes - Benign CNVs mapping >0 OMIM genes)',
#                 'Supp Table 14. Results - DECIPHER Scenario #2 (Pathogenic CNVs mapping 0 OMIM genes - Benign CNVs mapping 0 OMIM genes))',
#                 'Supp Table 15. Results - DECIPHER Scenario #3 (Pathogenic CNVs mapping 0 protein-coding genes - Benign CNVs mapping 0 protein-coding genes)',
#                 'Supp Table 16. Results - DECIPHER Scenario #4 (Pathogenic CNVs mapping >0 OMIM genes - Benign CNVs mapping >0 OMIM genes)',
#                 'Supp Table 18')) %>%
#   sheet_write(google_calc_results, sheet = "index")
# Figure VUS-DECIPHER and VUS-ClinVar ------------------------------------------------------------------------------

# vus_decipher <- read_tsv('/data-cbl/decipher_data/decipher-cnvs-grch37-2020-12-06.txt', skip = 1)
# 
# vus_decipher <- vus_decipher %>% 
#   rename(id = `# patient_id`, clinical = pathogenicity) %>%
#   # select(id, clinical2, contribution) %>%
#   filter(clinical == 'Uncertain') %>%
#   # filter(clinical %in% c('Pathogenic', 'Unknown', 'Likely pathogenic')) %>%
#   mutate(id = as.character(id)) %>%
#   rename(chrom = chr, id_tmp = id) %>%
#   filter(variant_class == 'Deletion') %>%
#   filter(inheritance %in% 'De novo') %>%
#   filter(genotype == 'Heterozygous') %>%
#   mutate(variant_class = tolower(variant_class)) %>%
#   mutate(source = 'clinvar') %>%
#   mutate(length_cnv = end - start + 1) %>%
#   select(chrom, start, end, variant_class, clinical, source, length_cnv) %>%
#   mutate(id_tmp = row_number())
# 
# 
# # 0
# vus_decipher_remove_overlap_ids <- vus_decipher %>% 
#   bed_coverage(problematic_regions) %>% 
#   select(id_tmp, .frac) %>% filter(.frac >= 0.3) %>% pull(id_tmp)
# 
# vus_decipher <- vus_decipher %>% filter(!id_tmp %in% vus_decipher_remove_overlap_ids)
# 
# # 1
# vus_decipher_remove_split_ids <- report_split_cnvs(vus_decipher)
# 
# vus_decipher <- vus_decipher %>% filter(id_tmp %in% vus_decipher_remove_split_ids)
# # 2. Remove identical CNVs------------------------------------------------------------------------------
# 
# vus_decipher <- vus_decipher %>%
#   distinct(chrom, start, end, clinical, .keep_all = TRUE)
# # 3. Reciprocal overlap------------------------------------------------------------------------------
# 
# 
# vus_decipher_remove_reciprocal_ids <- reciprocal_overlap(vus_decipher)
# 
# vus_decipher <- vus_decipher %>% filter(!id_tmp %in% vus_decipher_remove_reciprocal_ids)
# 
# # Check pathogenic - benign overlap ------------------------------------------------------------------------------
# 
# vus_decipher_remove_second_overlap_ids <- overlap_benign_pathogenic(vus_decipher, input_check_cnv_del_benign)
# 
# vus_decipher <-  vus_decipher %>% filter(!id_tmp %in% vus_decipher_remove_second_overlap_ids$id_tmp_patho)
# 
# # Matching by length------------------------------------------------------------------------------
# 
# vus_decipher_before_match <- vus_decipher %>%
#   mutate(clinical = 'pathogenic') %>%
#   bind_rows(input_check_cnv_del_benign %>% 
#               filter(!id_tmp %in% vus_decipher_remove_second_overlap_ids$id_tmp_benign) %>%
#               select(-LastEvaluated, -NumberSubmitters))
# 
# vus_decipher_match_deletion <- matching_length(bin_length = 100, vus_decipher_before_match)
# 
# vus_decipher_annot <- check_cnv_v2(vus_decipher_match_deletion %>% select(-id))
# 
# 
# vus_decipher_scores <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, vus_decipher_annot, 'deletion', 'unbiased approach', 
#                                                   only_table = TRUE)
# 
# 
# ref_sd_decipher_uncertain <- get_reliability_score_mid(ref_quantiles, vus_decipher_scores %>% filter(clinical == 'pathogenic'))
# 
# vus_decipher_results <- ref_sd_decipher_uncertain %>% count(reliability_score) %>%
#   # group_by(clinical) %>%
#   # mutate(perc = n / sum(n)) %>%
#   mutate(clinical2 = 'VOUS') %>%
#   select(reliability_score, clinical2, n)
# 
# 
# # ClinVar
# 
# vus_clinvar <- clinvar_stacked_plot %>%
#   filter(clinical %in% c('uncertain significance')) %>%
#   rename(clinical2 = clinical) %>%
#   mutate(clinical = 'pathogenic') %>%
#   mutate(id_tmp = row_number())
# 
# 
# # 0
# vus_clinvar_remove_overlap_ids <- vus_clinvar %>% 
#   bed_coverage(problematic_regions) %>% 
#   select(id_tmp, .frac) %>% filter(.frac >= 0.3) %>% pull(id_tmp)
# 
# vus_clinvar <- vus_clinvar %>% filter(!id_tmp %in% vus_decipher_remove_overlap_ids)
# 
# # 1
# vus_clinvar_remove_split_ids <- report_split_cnvs(vus_clinvar)
# 
# vus_clinvar <- vus_clinvar %>% filter(id_tmp %in% vus_clinvar_remove_split_ids)
# 
# # 2
# vus_clinvar <- vus_clinvar %>%
#   distinct(chrom, start, end, clinical, .keep_all = TRUE)
# 
# # 3
# 
# vus_clinvar_remove_reciprocal_ids <- reciprocal_overlap(vus_clinvar)
# 
# vus_clinvar <- vus_clinvar %>% filter(!id_tmp %in% vus_clinvar_remove_reciprocal_ids)
# 
# 
# 
# # 4
# 
# vus_clinvar_remove_second_overlap_ids <- overlap_benign_pathogenic(vus_clinvar, input_check_cnv_del_benign)
# 
# vus_clinvar <-  vus_clinvar %>% filter(!id_tmp %in% vus_clinvar_remove_second_overlap_ids$id_tmp_patho)
# 
# # 5
# 
# vus_clinvar_before_match <- vus_clinvar %>%
#   mutate(clinical = 'pathogenic') %>%
#   bind_rows(input_check_cnv_del_benign %>% 
#               filter(!id_tmp %in% vus_clinvar_remove_second_overlap_ids$id_tmp_benign) %>%
#               select(-LastEvaluated, -NumberSubmitters))
# 
# # 6
# 
# vus_clinvar_match_deletion <- matching_length(bin_length = 100, vus_clinvar_before_match %>% mutate(id_tmp = row_number()))
# 
# vus_clinvar_annot <- check_cnv_v2(vus_clinvar_match_deletion %>% select(-id) %>%
#                                     filter(clinical == 'pathogenic'))
# 
# # 7
# 
# 
# vus_clinvar_scores <- predict_chrom_aware_rtemis(bayesian_clinvar_del_nohuman, vus_clinvar_annot, 'deletion', 'unbiased approach', 
#                                                  only_table = TRUE)
# 
# 
# ref_sd_clinvar_uncertain <- get_reliability_score_mid(ref_quantiles, vus_clinvar_scores)
# 
# vus_clinvar_results <- ref_sd_clinvar_uncertain %>% count(reliability_score) %>%
#   # group_by(clinical) %>%
#   # mutate(perc = n / sum(n)) %>%
#   mutate(clinical2 = 'VOUS') %>%
#   select(reliability_score, clinical2, n)
# 
# # Final merge
# 
# bind_rows(
#   vus_decipher_results %>% rename(clinical = clinical2) %>% mutate(clinical = 'DECIPHER - VOUS') %>% mutate(tag = 'DECIPHER'),
#   vus_clinvar_results %>% rename(clinical = clinical2) %>% mutate(clinical = 'ClinVar - VOUS') %>% mutate(tag = 'ClinVar'),
#   ref_sd_clinvar_del_real %>% count(clinical, reliability_score) %>% mutate(clinical = paste('ClinVar -', clinical)) %>% mutate(tag = 'ClinVar'), 
#   ref_sd_decipher_del_real %>% count(clinical, reliability_score) %>% mutate(clinical = paste('DECIPHER -', clinical)) %>% mutate(tag = 'DECIPHER')
# ) %>%
#   group_by(clinical) %>%
#   mutate(perc = n / sum(n)) %>%
#   mutate(reliability_score = factor(reliability_score, levels = c('3', '2', '1'))) %>%
#   ggplot(aes(clinical, perc, group = reliability_score)) +
#   geom_col(aes(fill = as.factor(reliability_score)), color = 'black') +
#   scale_y_continuous(label = percent) +
#   labs(x = 'Clinical assessment', y = 'Percentage', fill = 'Uncertainty level') +
#   facet_wrap(vars(tag), scales = 'free') +
#   geom_text(aes(label = paste0(100*round(perc, 2), '% ', '(', n, ')' )),
#             size = 3, position = position_stack(vjust = 0.5)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.8))
# 
# tmp_id_decipher_clinical <- read_tsv('/data-cbl/decipher_data/decipher-cnvs-grch37-2020-12-06.txt', skip = 1)
# 
# tmp_id_decipher_clinical <- tmp_id_decipher_clinical %>% 
#   rename(id = `# patient_id`, clinical2 = pathogenicity) %>%
#   select(id, clinical2, contribution) %>%
#   filter(clinical2 %in% c('Pathogenic', 'Unknown', 'Likely pathogenic')) %>%
#   mutate(id = as.character(id))
# 
# 
# output_decipher_deletion %>% select(id.x, id) %>% 
#   left_join(tmp_id_decipher_clinical, by = c('id.x' = 'id')) %>%
#   distinct() %>%
#   right_join(ref_sd_decipher_del_real %>% select(id, reliability_score), by = 'id') %>%
#   mutate(clinical2 = if_else(is.na(clinical2), 'Benign', clinical2)) %>%
#   mutate(clinical2 = factor(clinical2, levels = c('Benign', 'Pathogenic', 'Likely pathogenic', 'Unknown'))) %>%
#   count(reliability_score, clinical2) %>%
#   bind_rows(vus_decipher_results) %>%
#   group_by(clinical2) %>%
#   mutate(perc = n / sum(n)) %>%
#   mutate(clinical2 = factor(clinical2, levels = c('Pathogenic', 'Likely pathogenic', 'Unknown', 'VOUS', 'Benign'))) %>%
#   mutate(reliability_score = factor(reliability_score, levels = c('3', '2', '1'))) %>%
#   ggplot(aes(clinical2, perc)) +
#   geom_col(aes(fill = as.factor(reliability_score)), color = 'black') +
#   scale_y_continuous(label = percent) +
#   labs(x = 'Clinical assessment', y = 'Percentage', fill = 'Uncertainty level') +
#   theme_minimal()
# 
# output_decipher_deletion %>% select(id.x, id) %>% 
#   left_join(tmp_id_decipher_clinical, by = c('id.x' = 'id')) %>%
#   distinct() %>%
#   right_join(ref_sd_decipher_del_real %>% select(id, reliability_score), by = 'id') %>%
#   mutate(clinical2 = if_else(is.na(clinical2), 'Benign', clinical2)) %>%
#   mutate(clinical2 = factor(clinical2, levels = c('Benign', 'Pathogenic', 'Likely pathogenic', 'Unknown'))) %>%
#   filter(!is.na(contribution)) %>%
#   count(reliability_score, contribution) %>%
#   group_by(contribution) %>%
#   mutate(perc = n / sum(n)) %>%
#   ggplot(aes(contribution, perc)) +
#   geom_col(aes(fill = as.factor(reliability_score))) +
#   theme_minimal()

# Test - CNVs mapping ------------------------------------------------------------------------------


# check_omim_match_del <- input_check_cnv_del_pathogenic_clinvar %>%
#   bind_rows(input_check_cnv_del_benign) %>%
#   bed_intersect(hgcn_genes %>% 
#                   filter(omim == 'Yes') %>% 
#                   select(chrom, start, end),
#                 suffix = c('', '.y')) %>%
#   select(-c(start.y, end.y, .source, .overlap)) %>%
#   distinct(chrom, start, end, .keep_all = TRUE) 
# 
# 
# check_omim_matched_del <- matching_length(bin_length = 100, check_omim_match_del)
# 
# 
# check_protein_match_del <- input_check_cnv_del_pathogenic_clinvar %>%
#   bind_rows(input_check_cnv_del_benign) %>%
#   bed_intersect(hgcn_genes  %>% 
#                   select(chrom, start, end),
#                 suffix = c('', '.y')) %>%
#   select(-c(start.y, end.y, .source, .overlap)) %>%
#   distinct(chrom, start, end, .keep_all = TRUE) 
# 
# 
# check_protein_matched_del <- matching_length(bin_length = 100, check_protein_match_del)
# 
# 
# check_omim_matched_del %>% count(clinical)
# check_protein_matched_del  %>% count(clinical)
# clinvar_match_deletion  %>% count(clinical)


# result_bancco <-  bancco_del_decipher_rel %>%
#   mutate(tag2 = 'DECIPHER DEL') %>%
#   bind_rows(
# bancco_dup_decipher_rel %>%
#   mutate(tag2 = 'DECIPHER DUP') ,
# bancco_del_rel %>%
#   mutate(tag2 = 'ClinVar DEL'),
# bancco_dup_rel %>%
#   mutate(tag2 = 'ClinVar DUP')
# ) %>%
#   # group_by(tag, tag2) %>%
#   mutate(clinical = factor(clinical, levels = c('pathogenic', 'benign')))
# 
# result_bancco_auroc <- result_bancco %>%
#   group_by(tag, tag2) %>%
#   roc_auc(clinical, .pred_pathogenic) %>%
#   rename(auroc_total = .estimate) %>%
#   select(tag2, auroc_total)

# result_bancco_prauc <- result_bancco %>%
#   group_by(tag, tag2) %>%
#   pr_auc(clinical, .pred_pathogenic) %>%
#   rename(aupr_total = .estimate) %>%
#   select(tag2, aupr_total)
# 
# result_bancco_auroc_split <- result_bancco %>%
#   group_by(tag, reliability_score, tag2) %>%
#   roc_auc(clinical, .pred_pathogenic) %>%
#   rename(auroc_split = .estimate) %>%
#   select(tag2, reliability_score, auroc_split)
# 
# 

# 
# result_bancco_prauc_split <- result_bancco %>%
#   group_by(tag, reliability_score, tag2) %>%
#   pr_auc(clinical, .pred_pathogenic) %>%
#   rename(prauc_split = .estimate) %>%
#   select(tag2, reliability_score, prauc_split)
# 
# result_bancco_count <- result_bancco %>%
#   count(reliability_score, tag2) %>%
#   group_by(tag2) %>%
#   mutate(perc = 100*(n / sum(n))) %>%
#   select(-n)
# 
# result_bancco_auroc %>%
#   left_join(result_bancco_prauc, by = 'tag2') %>%
#   left_join(result_bancco_count, by = 'tag2') %>%
#   left_join(result_bancco_auroc_split, by = c('tag2', 'reliability_score')) %>%
#   left_join(result_bancco_prauc_split, by = c('tag2', 'reliability_score')) %>%
#   rename(uncertainty_level = reliability_score) %>%
#   mutate(across(where(is.numeric), round, digits= 2)) %>%
#   gt()

# let's do a 6 panel figure, 2 rows x 3 columns, first row for deletions, second row for duplications:

# ROC (left panel) and PR (middle panel) curves on the TOTAL Bancco set, for the model trained on Clinvar together with the model trained on Decipher.  And label them as "pathogenicity CNVScore - trained on Clinvar" and "pathogenicity CNVScore - trained on Decipher", respectively
# Right panel : Barplots representing the distribution across uncertainty bins of the Bancco dataset. Label them as "uncertainty CNVScore - trained on Clinvar" and "unreliability CNVScore - trained on Decipher". Rather than using "1" , "2", "3", use "lowly uncertain" , "moderately uncertain" and "highly uncertain" (do so across all figures throughout the mansucript, as "1", "2", "3" are difficult to interpret and in addition it can be confused with the Scenarios, numbering)
# When putting headers, indicate clearly that this is the Bancco dataset, so that it is not confused with the evaluation on Clinvar or Decipher

#######



# 
# # Deletion CNVs----------
# 
# 
# decipher_by_decipher_del <- predict_chrom_aware_rtemis(bayesian_decipher_del_nohuman, 
#                                                    output_decipher_deletion, '', '')
# 
# trendline_decipher <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
#                          data= decipher_by_decipher_del[[3]])
# 
# 
# res_df_decipher_del <- decipher_by_decipher_del[[3]] %>%
#   mutate(gam_residuals = residuals(trendline_decipher)) %>%
#   mutate(score_interval = ntile(.pred_pathogenic, n = 3)) %>%
#   group_by(score_interval) %>%
#   mutate(reliability_score = ntile(gam_residuals, split_residuals)) %>%
#   ungroup()
# 
# 
# score_intervals_df_decipher_del <- res_df_decipher_del %>% 
#   group_by(score_interval) %>% 
#   summarise(max_score = max(.pred_pathogenic), 
#             min_score = min(.pred_pathogenic)) %>%
#   arrange(score_interval)
# 
# ref_quantiles_decipher_del <- res_df_decipher_del %>%
#   left_join(score_intervals_df_decipher_del, by = 'score_interval') %>%
#   mutate(score_interval = as.character(score_interval)) %>%
#   select(score_interval, min_score, max_score, reliability_score, sd, gam_residuals)
# 
# 
# # Duplication CNVs----------
# 
# 
# decipher_by_decipher_dup <- predict_chrom_aware_rtemis(bayesian_decipher_dup_nohuman, 
#                                                        output_decipher_duplication, '', '')
# 
# trendline_decipher_dup <- gam(sd ~ s(.pred_pathogenic, bs="cr"), 
#                           data= decipher_by_decipher_dup[[3]])
# 
# 
# res_df_decipher_dup <- decipher_by_decipher_dup[[3]] %>%
#   mutate(gam_residuals = residuals(trendline_decipher_dup)) %>%
#   mutate(score_interval = ntile(.pred_pathogenic, n = 3)) %>%
#   group_by(score_interval) %>%
#   mutate(reliability_score = ntile(gam_residuals, split_residuals)) %>%
#   ungroup()
# 
# 
# score_intervals_df_decipher_dup <- res_df_decipher_dup %>% 
#   group_by(score_interval) %>% 
#   summarise(max_score = max(.pred_pathogenic), 
#             min_score = min(.pred_pathogenic)) %>%
#   arrange(score_interval)
# 
# ref_quantiles_decipher_dup <- res_df_decipher_dup %>%
#   left_join(score_intervals_df_decipher_dup, by = 'score_interval') %>%
#   mutate(score_interval = as.character(score_interval)) %>%
#   select(score_interval, min_score, max_score, reliability_score, sd, gam_residuals)

