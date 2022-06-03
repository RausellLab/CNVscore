# Load packages
library(readxl)
library(liftOver)
library(rtracklayer)
library(glue)
library(rentrez)
library(biomaRt) 
library(data.table)
library(tidyverse)
library(clusterProfiler)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicScores)
library(Biostrings)
select <- dplyr::select
rename <- dplyr::rename




# ------------------------------------------------------------------------------
# Centromeric and telomeric regions
# Source: UCSC
# ------------------------------------------------------------------------------  

# download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz',
#               destfile = 'gap.txt.gz')
# 
# read_tsv('gap.txt.gz', col_names = FALSE) %>%
#   filter(X2 == 'chr17') %>%
#   filter(X8 == 'telomere')

session <- browserSession("UCSC")
genome(session) <- "hg19"


region_gaps <- getTable(ucscTableQuery(session, "Gap"))


region_gaps <- region_gaps %>%
  as_tibble() %>%
  filter(type %in% c('telomere', 'centromere')) %>%
  mutate(chrom = str_remove(as.character(chrom), 'chr')) %>%
  rename(start = chromStart, end = chromEnd ) %>%
  select(chrom, start, end, type) %>%
  # In this version (hg19), there are not chrom 17 telomeres
  # I checked length in hg38 and included
  # 81195210 is the length of the chrom 
  bind_rows(tibble(chrom = '17', 
                   start = c(1, 81195210 - 1e4), 
                   end = c(1e4,81195210 ),
                   type = 'telomere'))


# ------------------------------------------------------------------------------
# LiftOver
# Instructions: https://www.bioconductor.org/packages/release/workflows/vignettes/liftOver/inst/doc/liftov.html
# ------------------------------------------------------------------------------

path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
from_hg38_to_hg19 = import.chain(path)

download.file('http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz',
              destfile = 'hg18ToHg19.over.chain.gz')

system('gunzip hg18ToHg19.over.chain.gz')

from_hg18_to_hg19 = import.chain('hg18ToHg19.over.chain')

file.remove('hg18ToHg19.over.chain')


download.file('http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz',
              destfile = 'hg19ToHg38.over.chain.gz')

system('gunzip hg19ToHg38.over.chain.gz')

from_hg19_to_hg38 = import.chain('hg19ToHg38.over.chain')



# ------------------------------------------------------------------------------
# Dataset: STRING
# Source: https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz
# Version: 11.5
# Threshold score (700) based on STRING FAQs (https://string-db.org/cgi/help.pl?sessionId=H0T1nAmrXuLJ)
# combine_subscores.py
# ------------------------------------------------------------------------------


library(tidygraph)
library(igraph)
library(biomaRt)



# system('gunzip 9606.protein.links.v11.0.txt.gz')

string_db_pre <- read_delim('9606.protein.links.full.v11.5.txt.gz', delim = ' ')

string_db_pre <- string_db_pre %>%
  mutate_if(is.double, function(x) x / 1e3) %>%
  mutate_if(is.double, function(x) if_else(x < 0.041, 0.041, (x - 0.041) / (1 - 0.041))) %>%
  mutate(neighborhood_both = 1.0 - (1.0 - neighborhood) * (1.0 - neighborhood_transferred)) %>%
  mutate(coexpression_both = 1.0 - (1.0 - coexpression) * (1.0 - coexpression_transferred)) %>%
  # mutate(experiments_both = 1.0 - (1.0 - experiments) * (1.0 - experiments_transferred)) %>%
  # mutate(database_both = 1.0 - (1.0 - database) * (1.0 - database_transferred)) %>%
  # mutate(textmining_both = 1.0 - (1.0 - textmining) * (1.0 - textmining_transferred)) %>%
  mutate(cooccurence_corrected = cooccurence * (1.0 - homology)) %>%
  # mutate(textmining_corrected = textmining_both * (1.0 - homology)) %>%
  select(!matches('experiments|database|textmining|combined_score')) %>%
  select(protein1, protein2, neighborhood_both, fusion, cooccurence_corrected, coexpression_both) %>%
  rename(neighborhood = neighborhood_both, cooccurence = cooccurence_corrected, coexpression = coexpression_both) %>%
  mutate(combined_score_one_minus = (1 - neighborhood) * (1 - cooccurence) * (1 - coexpression) * (1 - fusion)) %>%
  mutate(combined_score = ((1 - combined_score_one_minus) * (1 - 0.041)) + 0.041) %>%
  mutate(combined_score = combined_score * 1e3) %>%
  select(protein1, protein2, combined_score) %>%
  filter(combined_score >= 700) %>%
  mutate(protein1 = str_remove(protein1, '9606.'),
         protein2 = str_remove(protein2, '9606.')) %>%
  select(-combined_score) %>%
  rename(p1 = protein1, p2 = protein2)
  



human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")

interval_genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_peptide_id'),
                        mart = human ) %>% 
  as_tibble()

string_db_pre <- string_db_pre %>% 
  left_join(interval_genes, by = c('p1' = 'ensembl_peptide_id')) %>%
  left_join(interval_genes, by = c('p2' = 'ensembl_peptide_id')) %>%
  rename(from = hgnc_symbol.x, to = hgnc_symbol.y) %>%
  select(from, to) %>%
  na.omit() %>%
  filter(!(from == '' | to == ''))


degree_haplo_genes <- string_db_pre %>% filter(to %in% (haplo_triplo_genes %>% filter(haplo == 'yes') %>% pull(gene))) %>% 
  count(from)
degree_triplo_genes <- string_db_pre %>% filter(to %in% (haplo_triplo_genes %>% filter(triplo == 'yes') %>% pull(gene))) %>% 
  count(from)

string_db_pre <- string_db_pre %>% as_tbl_graph(string_db_pre, directed = FALSE)



string_db <- string_db_pre %>%
  mutate(page_rank = centrality_pagerank(),
         degree = centrality_degree()) %>%
  activate(nodes) %>%
  as_tibble()


vector_haplo_genes <- haplo_triplo_genes %>% filter(haplo == 'yes') %>% pull(gene)
vector_haplo_genes <- vector_haplo_genes[vector_haplo_genes %in% string_db$name]

vector_triplo_genes <- haplo_triplo_genes %>% filter(triplo == 'yes') %>% pull(gene)
vector_triplo_genes <- vector_triplo_genes[vector_triplo_genes %in% string_db$name]

get_shortest_path <- function(input_name) {
  
  # print(glue('{input_name}'))
  
  # tmp_id <- string_db$name[3]
  
  tmp_id <- input_name
  
  # string_db_pre %>% activate(nodes) %>% filter(name == 'ARF5')
  
  # Nº interactions haplo
  # string_db_pre %>% filter(id %in% tmp_id) %>% filter(id2 %in% vector_haplo_genes)
  
  # haplo
  if (tmp_id %in% vector_haplo_genes) {
    
    tmp_result_haplo <- 0
  } else {
    tmp_collector <- c()
    for (i in 1:length(vector_haplo_genes)) {
      
      tmp_value <- length(shortest_paths(string_db_pre, tmp_id, vector_haplo_genes[i])$vpath[[1]])
      tmp_collector <- c(tmp_collector, tmp_value)
    }
    tmp_result_haplo <- min(tmp_collector[tmp_collector != 0])
  }
  
  # triplo
  if (tmp_id %in% vector_triplo_genes) {
    
    tmp_result_triplo <- 0
  } else {
    tmp_collector <- c()
    for (i in 1:length(vector_triplo_genes)) {
      
      tmp_value <- length(shortest_paths(string_db_pre, tmp_id, vector_triplo_genes[i])$vpath[[1]])
      tmp_collector <- c(tmp_collector, tmp_value)
    }
    tmp_result_triplo <- min(tmp_collector[tmp_collector != 0])
  }
  
  
  
  return(tibble('gene' = input_name, 'haplo' = tmp_result_haplo, 'triplo' = tmp_result_triplo))
}



plan("multiprocess", workers = 60)

haplo_triplo_shortest_path_tibble <- string_db$name %>% future_map_dfr(~ get_shortest_path(.), .progress = TRUE)

string_db <- haplo_triplo_shortest_path_tibble %>%
  mutate(haplo = if_else(haplo == 0, haplo, haplo - 1)) %>%
  mutate(triplo = if_else(triplo == 0, triplo, triplo - 1)) %>%
  mutate(haplo = ifelse(haplo == Inf, 3, haplo)) %>% 
  mutate(triplo = ifelse(triplo == Inf, 3, triplo)) %>%
  rename(haplo_short = haplo, triplo_short = triplo) %>%
  right_join(string_db, by = c('gene' = 'name'))


string_db <- string_db %>%
  left_join(degree_haplo_genes %>% rename(degree_to_haplo = n), by = c('gene'= 'from')) %>%
  left_join(degree_triplo_genes %>% rename(degree_to_triplo = n), by = c('gene'= 'from')) %>%
  replace_na(list(degree_to_haplo = 0, degree_to_triplo = 0))


# ------------------------------------------------------------------------------
# pNull
# ------------------------------------------------------------------------------


url <- 'https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz'
download.file(url, destfile = basename(url))

system('mv gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz gnomad.v2.1.1.lof_metrics.by_transcript.txt.gz &&
gunzip gnomad.v2.1.1.lof_metrics.by_transcript.txt.gz')

pnull <- read_tsv('gnomad.v2.1.1.lof_metrics.by_transcript.txt') %>%
  filter(canonical == 'TRUE') %>% 
  # filter(exp_lof < 10) %>% # take into account - important for remot gw project
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  select(gene, pNull) %>%
  mutate(pNull = round(pNull, 2)) %>%
  rename(pnull = pNull) %>%
  mutate(pnull = ntile(-pnull, 100))


file.remove('gnomad.v2.1.1.lof_metrics.by_transcript.txt')

# ------------------------------------------------------------------------------
# GTEx - expression features
# ------------------------------------------------------------------------------

url <- 'https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'

download.file(url, 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz')

mean_median_expression <- read_tsv('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz', skip = 2)


expression_features <- mean_median_expression %>%
  select(-Name) %>%
  pivot_longer(-Description, names_to = 'tissue', values_to = 'median') %>%
  group_by(Description) %>%
  summarise(mean_expression = mean(median), .groups = 'keep', min_expression = min(median)) %>%
  ungroup() %>%
  rename(gene = Description) %>%
  mutate(mean_expression = ntile(mean_expression, 100)) %>%
  mutate(min_expression = ntile(-min_expression, 100)) 

file.remove('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz')



# ------------------------------------------------------------------------------
# Dataset: Recombination rates
# Source: https://science.sciencemag.org/content/suppl/2019/01/23/363.6425.eaau1043.DC1
# Genome assembly: hg38
# https://pubmed.ncbi.nlm.nih.gov/30679340/
# ------------------------------------------------------------------------------
url <- "https://science.sciencemag.org/highwire/filestream/721792/field_highwire_adjunct_files/4/aau1043_DataS3.gz"
download.file(url, 'aau1043_DataS3.gz')
system("gunzip  aau1043_DataS3.gz")  
recomb <- read_tsv('aau1043_DataS3', skip = 8, col_names = c('chrom', 'start', 'end', 'cm_mb', 'cm'))

file.remove('aau1043_DataS3')

tmp_granges_recomb <- recomb %>% GRanges()

seqlevelsStyle(tmp_granges_recomb) = "UCSC"  # necessary
tmp_granges_recomb = liftOver(tmp_granges_recomb, from_hg38_to_hg19)

recomb <- tmp_granges_recomb %>% as_tibble() %>% select(seqnames, start, end, cm_mb, cm) %>%
  rename(chrom = seqnames) %>% mutate(chrom = str_remove(as.character(chrom), 'chr'))




# ------------------------------------------------------------------------------
# LADs
# ------------------------------------------------------------------------------

lads <- read_tsv('https://static-content.springer.com/esm/art%3A10.1038%2Fnature06947/MediaObjects/41586_2008_BFnature06947_MOESM252_ESM.txt', col_names = FALSE)

lads_hg18 <- lads %>% 
  rename(chrom = X1, start = X4, end = X5) %>%
  select(chrom, start, end) %>%
  GRanges()


seqlevelsStyle(lads_hg18) = "UCSC"

lads = liftOver(lads_hg18, from_hg18_to_hg19)

lads <- lads %>%
  as_tibble() %>%
  rename(chrom = seqnames) %>%
  select(chrom, start, end) %>%
  mutate(chrom = str_remove(chrom, 'chr'))

# ------------------------------------------------------------------------------
# LOEUF
# ------------------------------------------------------------------------------

download.file('https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', 
              destfile = 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')

loeuf_score <- read_tsv('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz') %>% select(gene, oe_lof_upper) %>%
  rename(loeuf = oe_lof_upper) %>%
  mutate(loeuf = 1 - loeuf)

file.remove('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')

# ------------------------------------------------------------------------------
# Dataset: Pubmed articles associated with G-bands and "deletion" "duplication" keywords
# ------------------------------------------------------------------------------


pubmed_bands <- chromPlot::hg_cytoBandIdeo %>% 
  as_tibble() %>%
  select(Name, Chrom)

get_band <- function(band_input, chrom_input) {
  
  
  print(glue('{chrom_input} - {band_input}'))
  
  # band_input <- 'q24.3'
  # chrom_input <- '5'
  Sys.sleep(1)
  
  band_tmp <- band_input
  band2_tmp <- paste0(chrom_input, band_input)
  chrom_tmp <- paste('chromosome', chrom_input)
  
  result_del <- length(entrez_search(db="pubmed", 
                                     term= paste('(','(', chrom_tmp,'AND', band_tmp, ')','OR', band2_tmp,')', 'AND (deletion OR microdeletion) AND homo sapiens'), retmax = 1000 )$ids)
  result_dup <- length(entrez_search(db="pubmed", 
                                     term= paste('(','(', chrom_tmp,'AND', band_tmp, ')','OR', band2_tmp,')', 'AND duplication AND homo sapiens'), retmax = 1000 )$ids)
  
  result <- tibble(band = band_tmp, chrom = chrom_input, hits_del = result_del, hits_dup = result_dup)
  
  return(result)
}


# plan("multiprocess", workers = 2)


pubmed_lists <- pmap_dfr(list(pubmed_bands$Name, 
                              pubmed_bands$Chrom),
                         get_band)

pubmed_df <- pubmed_lists

pubmed_df <- pubmed_df %>% left_join( chromPlot::hg_cytoBandIdeo %>% select(Chrom, Start, End, Name), by =
                                        c('band' = 'Name', 'chrom' = 'Chrom')) %>%
  rename(start = Start, end = End) %>% 
  mutate(start = start + 1)


# ------------------------------------------------------------------------------
# Annotation promoter region (-2kbs - TSS - + 2kbs)
# ------------------------------------------------------------------------------


phast46pla <- getGScores("phastCons46wayPlacental.UCSC.hg19")
plan("multiprocess", workers = 40)



get_annot_promoter <- function(gene, chrom, tss) {
  
  # gene <- 'PRY'
  # chrom <- hgcn_genes[hgcn_genes$gene == 'PRY',]$chrom
  # tss <- hgcn_genes[hgcn_genes$gene == 'PRY',]$start
  
  tmp_granges <- GRanges(seqnames= paste0('chr', chrom), 
                         IRanges(start=(tss-2000):(tss+2000), width=1))
  
  # Calculation mean PhastCons46
  phast_temp <- gscores(phast46pla, tmp_granges)
  phast_temp <- phast_temp$default %>% mean(na.rm = TRUE)
  
  
  seqs <- getSeq(Hsapiens, paste0('chr', chrom), (tss-2000), (tss+2000))
  
  #maybe remove length(seqs)
  cpg_temp <- length(seqs)*dinucleotideFrequency(seqs)["CG"]/
    (letterFrequency(seqs, "C")*letterFrequency(seqs, "G"))
  
  
  result_temp <- tibble(gene = gene, mean_phast = phast_temp, cpg_density = cpg_temp )
  return(result_temp)
}


tic()

output_temp <- future_pmap_dfr(list(hgcn_genes$gene, 
                                    hgcn_genes$chrom, 
                                    hgcn_genes$start), 
                               get_annot_promoter, .progress = TRUE)

genes_promoter <- output_temp

toc()

genes_promoter <- genes_promoter %>%
  mutate(mean_phast = ntile(mean_phast, 10)) %>%
  mutate(cpg_density = ntile(cpg_density, 10))

# ------------------------------------------------------------------------------
# Ensembl Regulatory Build
# Aggregation from ENCODE, Roadmap Epigenomics and Blueprint
# Version: Ensembl 103
# ------------------------------------------------------------------------------

download.file('http://ftp.ensembl.org/pub/grch37/release-103/regulation/homo_sapiens/RegulatoryFeatureActivity/H1_hESC_3/homo_sapiens.GRCh37.H1_hESC_3.Regulatory_Build.regulatory_activity.20191101.gff.gz',
              'homo_sapiens.GRCh37.H1_hESC_3.Regulatory_Build.regulatory_activity.20191101.gff.gz')

ensembl_reg <- read_tsv('homo_sapiens.GRCh37.H1_hESC_3.Regulatory_Build.regulatory_activity.20191101.gff.gz', col_names = FALSE)

file.remove('homo_sapiens.GRCh37.H1_hESC_3.Regulatory_Build.regulatory_activity.20191101.gff.gz')  

ensembl_reg <- ensembl_reg %>% 
  filter(!str_detect(X9, 'INACTIVE')) %>%
  select(X1, X3, X4, X5) %>%
  filter(nchar(X1) <= 2) %>%
  rename(chrom = X1,
         start = X4,
         end = X5,
         type = X3) %>%
  relocate(chrom, start, end, type)

# ------------------------------------------------------------------------------
# CRISPR scores (CS)
# More info: https://science.sciencemag.org/content/suppl/2015/10/14/science.aac7041.DC1 Table S3
# ------------------------------------------------------------------------------
library(readxl)
download.file('https://science.sciencemag.org/highwire/filestream/637113/field_highwire_adjunct_files/3/aac7041_SM_Table_S3.xlsx',
              'aac7041_SM_Table_S3.xlsx')

crispr_score_df <- read_xlsx('aac7041_SM_Table_S3.xlsx') %>% select(Gene, `KBM7 CS`) %>% rename(gene = Gene, score = `KBM7 CS`) %>%
  rename(crispr_score = score) %>%
  mutate(crispr_score = ntile(-crispr_score, 100))


# hgcn_genes %>% 
#   select(gene, fusil) %>% 
#   left_join(crispr_score, by = 'gene') %>% 
#   na.omit() %>% 
#   mutate(score = ntile(score, 10)) %>% 
#   count(fusil, score) %>% 
#   group_by(score) %>% 
#   mutate(perc = n / sum(n)) %>% 
#   ggplot(aes(score, perc)) + 
#   geom_col(aes(fill = fusil))

# ------------------------------------------------------------------------------
# UCNEbase
# More info: https://ccg.epfl.ch//UCNEbase/
# ------------------------------------------------------------------------------

ucne <- read_tsv('https://ccg.epfl.ch//UCNEbase/data/download/ucnes/hg19_UCNE_coord.bed', col_names = FALSE)
colnames(ucne) <- c('chrom', 'start', 'end', 'delete1', 'delete2')
ucne <- ucne %>% 
  select(-contains('delete')) %>%
  mutate(start = start + 1) %>% # 0-based -> 1-based
  mutate(chrom = str_remove(chrom, 'chr'))

# ------------------------------------------------------------------------------
# Human Accelerated Regions (HARs)
# Three sources extracted from Supp. file gnomad SVs (Collins et al. 2020)
# More info: https://ccg.epfl.ch//UCNEbase/
# ------------------------------------------------------------------------------

# source one:  https://www.nature.com/articles/nature10530
# source two:
# source three:

hars <- read_tsv('https://docpollard.org/wordpress/wp-content/research/hars_merged_hg19.bed', 
                 col_names = FALSE) %>%
  rename(chrom = X1, start = X2, end = X3, name = X4) %>%
  mutate(chrom = str_remove(chrom, 'chr')) %>%
  mutate(start = start + 1) # 0-based -> 1-based


# ------------------------------------------------------------------------------
# Dataset: SV Hotspot regions
# Link: https://science.sciencemag.org/highwire/filestream/757669/field_highwire_adjunct_files/3/abf7117_Ebert_Tables-S1-S56.xlsx
# Reference genome: hg38
# ------------------------------------------------------------------------------

download.file('https://science.sciencemag.org/highwire/filestream/757669/field_highwire_adjunct_files/3/abf7117_Ebert_Tables-S1-S56.xlsx',
              destfile = 'abf7117_Ebert_Tables-S1-S56.xlsx')


test_hotspot <- read_xlsx('abf7117_Ebert_Tables-S1-S56.xlsx', sheet = 15, skip = 1) %>%
  rename(chrom = Chr, start = Start, end = End) %>%
  select(chrom, start, end) %>% 
  na.omit() %>%
  GRanges()

file.remove('abf7117_Ebert_Tables-S1-S56.xlsx')

seqlevelsStyle(test_hotspot) = "UCSC"

test_hotspot_hg19 = liftOver(test_hotspot, from_hg38_to_hg19)

count_n_group <- test_hotspot_hg19 %>% as_tibble() %>% count(group) %>%
  filter(n == 1) %>%
  pull(group)

hotspot <- test_hotspot_hg19 %>% as_tibble() %>% filter(group %in% count_n_group) %>%
  select(seqnames, start, end) %>%
  rename(chrom = seqnames) %>%
  mutate(id = row_number()) %>%
  mutate(chrom = str_remove(chrom, 'chr'))


# ------------------------------------------------------------------------------
# EDS
# ------------------------------------------------------------------------------


download.file('https://www.cell.com/cms/10.1016/j.ajhg.2020.01.012/attachment/16ff2f31-e4aa-46b0-a384-a242843ac763/mmc2.xlsx',
              destfile = 'mmc2.xlsx')

eds <- read_xlsx('mmc2.xlsx')

eds <- eds %>%
  select(GeneSymbol, EDS) %>%
  left_join(hgcn_genes %>% select(ensembl_gene_id, gene), by = c('GeneSymbol' = 'ensembl_gene_id')) %>%
  select(gene, EDS) %>%
  rename(eds = EDS) %>%
  mutate(eds = ntile(eds, 100))


file.remove('mmc2.xlsx')

# ------------------------------------------------------------------------------
# Protein complex genes
# CORUM - Mammals Protein complexes 
# ------------------------------------------------------------------------------


download.file('http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip', 'allComplexes.txt.zip')
prot_complex <- read_tsv('allComplexes.txt.zip')

prot_complex <- prot_complex %>% 
  filter(Organism == 'Human') %>%
  select(ComplexID, ComplexName, `subunits(Gene name)`) %>%
  rename(id = ComplexID, name = ComplexName, gene =  `subunits(Gene name)`) %>%
  separate_rows(gene, sep = ';')


file.remove('allComplexes.txt.zip')


# ------------------------------------------------------------------------------
# Protein complex genes
# Ref: http://humap2.proteincomplexes.org
# ------------------------------------------------------------------------------
#  1=Extremely High, 2=Very High, 3=High, 4=Medium High, 5=Medium

hu_map <- read_csv('http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt')
hu_map <- hu_map %>% filter(Confidence == 1) %>% separate_rows(genenames, sep = ' ') %>% rename(gene = genenames) %>% select(gene)


# just checking enrichment
# hgcn_genes %>% 
#   mutate(prot_complex_nohuman = if_else(gene %in% hu_map$gene, 'Yes', 'No')) %>%
#   mutate(pLI = if_else(pLI >= 90, 'Yes', 'No')) %>%
#   count(pLI, prot_complex_nohuman) %>%
#   group_by(pLI) %>%
#   mutate(perc = n / sum(n)) %>%
#   na.omit()




# ------------------------------------------------------------------------------
# Generate sliding windows 100 b.p
# ------------------------------------------------------------------------------



windows_gw <- c(1:22, 'X') %>% map_dfr(function(x) {
  
  
  tmp_length <- coord_chrom_hg19 %>% filter(chrom == x) %>% pull(length)
  
  tmp_chunks <- tibble('chrom' = x, 'start' = seq(1, tmp_length, by = 100), 'end' = seq(1, tmp_length, by = 100) + 99)
  # We do this because it's out of the genomic boundaries
  tmp_chunks <- tmp_chunks %>% head(-1)
  
  tmp_chunks <- tmp_chunks %>% 
    mutate(start = start - 1) 
  
  return(tmp_chunks)
})



c(1:22, 'X') %>% map(function(x) {
  
  windows_gw %>%
    filter(chrom == x) %>%
    write_tsv(glue('/data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{x}.bed'), col_names = FALSE)
  
})

# ------------------------------------------------------------------------------
# conda activate giggle
# parallel --jobs 23 'bedtools intersect -sorted -wa -wb -a /data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{}.bed  -b /data/non-coding/CADD_v.1.6/whole_genome_SNVs.bed > /data-cbl/frequena_data/cnvscore/cadd_feature/output_chrom_{}.bed' ::: {1..22} X ; parallel --jobs 23 'bedtools intersect -sorted -wa -wb -a /data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{}.bed  -b /data/non-coding/gerp/gerp_scores.bed > /data-cbl/frequena_data/cnvscore/gerp_feature/output_chrom_{}.bed' ::: {1..22} X
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# CADD
# From TSV to BED format:
# gunzip whole_genome_SNVs.tsv.gz && cat whole_genome_SNVs.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' | awk -v n=1 '{print $1, $2-n, $3, $7}' OFS='\t' | tail -n +3 > whole_genome_SNVs.bed
# ------------------------------------------------------------------------------
plan('multiprocess', workers = 5)

tic()
cadd_max <- c(1:22, 'X') %>% map_dfr(function(x) {
  
  print(x)
  tmp_df <- read_tsv(glue('cadd_feature/output_chrom_{x}.bed'), col_names = FALSE)
  tmp_df2 <- tmp_df %>%
    mutate(X1 = as.character(X1)) %>%
    select(-X4, -X5, -X6) %>%
    group_by(X1, X2, X3) %>%
    summarise(max_cadd = max(X7)) %>%
    ungroup()
  
  return(tmp_df2)
  
})
toc()

cadd_max <- cadd_max %>% rename(chrom = X1, start = X2, end = X3) %>% mutate(start = start + 1)

write_tsv(cadd_max, 'cadd_feature/result_cadd.tsv')



# ------------------------------------------------------------------------------
# GERP
# gunzip /data/non-coding/CADD_v1.3/annotations/gerp/gerp_scores.tsv.gz
# cat gerp_scores.tsv | awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' | awk -v n=1 '{print $1, $2-n, $3, $5}' OFS='\t' > gerp_scores.bed
# Find intersects: 
# parallel --jobs 22 'bedtools intersect -sorted -wa -wb -a /data-cbl/frequena_data/cnvscore/cadd_feature/input_chrom_{}.bed  -b /data/non-coding/gerp/gerp_scores.bed > /data-cbl/frequena_data/cnvscore/gerp_feature/output_chrom_{}.bed' ::: {1..22}
# ------------------------------------------------------------------------------

gerp_raw <- c(1:22, 'X') %>% map_dfr(function(x) {
  
  print(x)
  tmp_df <- read_tsv(glue('gerp_feature/output_chrom_{x}.bed'), col_names = FALSE)
  tmp_df2 <- tmp_df %>%
    mutate(X1 = as.character(X1)) %>%
    select(-X4, -X5, -X6) %>%
    group_by(X1, X2, X3) %>%
    summarise(max_gerp = max(X7)) %>%
    ungroup()
  
  return(tmp_df2)
  
})

gerp_max <- gerp_raw %>% 
  rename(chrom = X1, start = X2, end = X3) %>% 
  mutate(chrom = as.character(chrom)) %>%
  mutate(start = start + 1)

write_tsv(gerp_max, 'gerp_feature/result_gerp.tsv')


# ------------------------------------------------------------------------------
# PhastCons scores: 
# ------------------------------------------------------------------------------



windows_gw2 <- windows_gw %>%
  mutate(tag = sample(1:1000, size = nrow(.), replace = TRUE)) %>%
  mutate(start = start + 1)

phast46pla <- getGScores("phastCons46wayPlacental.UCSC.hg19")
phast46pri <- getGScores("phastCons46wayPrimates.UCSC.hg19")
phylop100 <- getGScores("phyloP100way.UCSC.hg19")

library(phastCons100way.UCSC.hg19)

phast100 <- phastCons100way.UCSC.hg19



get_scores_max <- function(df, tag_input) {
  
  df_tmp <- df %>%
    # slice_sample(n = 100)
    filter(tag == tag_input)
  
  # df_tmp_phast100 <- gscores(phast100, df_tmp %>% GRanges(), summaryFun = max)
  # df_tmp_phast46pri <- gscores(phast46pri, df_tmp %>% GRanges(), summaryFun = max)
  df_tmp_phylop100 <- gscores(phylop100, df_tmp %>% GRanges(), summaryFun = max)
  
  # df_tmp_phast100 <- df_tmp_phast100 %>% as_tibble() %>% 
  #   rename(chrom = seqnames, max_phast100 = default) %>% 
  #   mutate(chrom = as.integer(str_remove(chrom, 'chr'))) %>%
  #   mutate(max_phast100 = if_else(is.na(max_phast100), 0, max_phast100))
  # 
  # df_tmp_phast46pri <- df_tmp_phast46pri %>% as_tibble() %>% rename(chrom = seqnames, max_phast46pri = default) %>%
  #   mutate(chrom = as.integer(str_remove(chrom, 'chr'))) %>%
  #   mutate(max_phast46pri = if_else(is.na(max_phast46pri), 0, max_phast46pri))
  
  df_tmp_phylop100 <- df_tmp_phylop100 %>% as_tibble() %>% 
    rename(chrom = seqnames, max_phylop100 = default) %>%
    mutate(chrom = as.integer(str_remove(chrom, 'chr'))) %>%
    mutate(max_phylop100 = if_else(is.na(max_phylop100), 0, max_phylop100))
  
  df_tmp %>%
    # left_join(df_tmp_phast100, by = c('chrom', 'start', 'end')) %>%
    # left_join(df_tmp_phast46pri, by = c('chrom', 'start', 'end')) %>%
    left_join(df_tmp_phylop100, by = c('chrom', 'start', 'end')) %>%
    select(chrom, start, end, max_phylop100)
}
  
  

plan("multiprocess", workers = 60)

tic()

output_temp <- windows_gw2 %>% 
  pull(tag) %>% 
  # .[1:40] %>%
  unique() %>% 
  sort() %>%
  future_map_dfr(function(x) get_scores_max(windows_gw2, tag_input = x ), .progress = TRUE)

phylop100_max <- output_temp %>% mutate(chrom = as.character(chrom))

toc()










gscores(phast46pla, df_tmp %>% GRanges(), summaryFun = max)



# phyloP100way.UCSC.hg19


# ------------------------------------------------------------------------------
# Dataset: List TFs human
# Source: https://www.nature.com/articles/nature13182
# ------------------------------------------------------------------------------
library(httr)

url1 <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fnature13182/MediaObjects/41586_2014_BFnature13182_MOESM86_ESM.xlsx'

GET(url1, write_disk(tf <- tempfile(fileext = ".xlsx")))

tf_genes <- read_excel(tf, sheet = 8) %>% select(SYMBOL) %>% distinct() %>% pull()



# ------------------------------------------------------------------------------
# Dataset: Paralogous genes
# Source: biomart
# ------------------------------------------------------------------------------


human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                  host    = "grch37.ensembl.org",
                  path    = "/biomart/martservice")


para_genes <- getBM(attributes = c('external_gene_name', 'hsapiens_paralog_ensembl_gene'),
                    mart = human )

para_genes <- para_genes %>% 
  as_tibble() %>%
  rename(gene = external_gene_name, para = hsapiens_paralog_ensembl_gene) %>%
  filter(para != '') %>%
  count(gene)


# ------------------------------------------------------------------------------
# Gene density
# ------------------------------------------------------------------------------  

coord_chrom_hg19 <- read_tsv('https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',
                             col_names = c('chrom', 'length')) %>%
  filter(nchar(chrom) < 6) %>% filter(!str_detect(chrom, 'chrM')) %>%
  mutate(chrom = str_remove(chrom, 'chr'))

result_tbl <- tibble()

for (i in 1:nrow(coord_chrom_hg19)) {
  
  last_nt <- coord_chrom_hg19 %>% slice(i) %>% pull(length)
  
  tmp_tbl <- tibble('chrom' = coord_chrom_hg19 %>% slice(i) %>% pull(chrom), 
                    'start' = seq(1, (last_nt-10**6), 10**6 ), 
                    'end' = seq(1000000, last_nt, 10**6,  ))
  
  tmp2_tbl <- tibble('chrom' = coord_chrom_hg19 %>% slice(i) %>% pull(chrom),
                     'start' = tmp_tbl %>% tail(1) %>% pull(end) + 1,
                     'end' = last_nt)
  
  
  result_tbl <- result_tbl %>% bind_rows(tmp_tbl, tmp2_tbl)
}

# Extra step: no  final regions have length = 1e6 so we create them and remove
# previous one
result_tbl <- result_tbl %>% 
  filter(!end %in% coord_chrom_hg19$length) %>%
  bind_rows(coord_chrom_hg19 %>% rename(end = length) %>%
              mutate(start = end - 1e6))

count_genes <- function(chrom, start, end) {
  
  result_tmp <- hgcn_genes %>% 
    bed_intersect(tibble('chrom' = chrom, 'start' = start, 'end' = end)) %>%
    nrow()
  
  return(result_tmp)
}


gene_density_tbl <- result_tbl %>% 
  mutate(gene_density = pmap_dbl(list(chrom, start, end), count_genes))
