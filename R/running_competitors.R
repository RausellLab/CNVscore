




# ------------------------------------------------------------------------------
# CHOOSE DATASET
# ------------------------------------------------------------------------------

to_write <- df_manolo
which_class <- 'DEL'

# ------------------------------------------------------------------------------
# CADD-SV
# ------------------------------------------------------------------------------


write_tsv(to_write %>%
            mutate(start = start - 1) %>%
            mutate(type = which_class) %>%
            select(chrom, start, end, type), 'rival_cnvscore/cadd_sv/input.bed', col_names = FALSE)



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


tmp_grange_structure_del <- to_write %>% 
  select(chrom, start, end, id) %>% 
  GRanges()

seqlevelsStyle(tmp_grange_structure_del) = "UCSC"  # necessary
tmp_grange_structure_del = liftOver(tmp_grange_structure_del, from_hg19_to_hg38)  %>% as_tibble()

tmp_grange_structure_del %>%
  rename(chrom = seqnames) %>%
  mutate(type = which_class) %>%
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



# ------------------------------------------------------------------------------
# SVSCORE
# ./generateannotations.pl (generate refGenes* files)
#  mv refGene.* tests/
# conda activate svscore
# parallel --jobs 22 'perl svscore.pl -i input_{}_del.vcf -o max,mean -e tests/refGene.exons.bed -c /data/non-coding/CADD_v.1.6/whole_genome_SNVs.tsv.gz -f tests/refGene.introns.bed > output_{}.del.vcf' ::: {1..22} && parallel --jobs 22 'perl svscore.pl -i input_{}_dup.vcf -o max,mean -e tests/refGene.exons.bed -c /data/non-coding/CADD_v.1.6/whole_genome_SNVs.tsv.gz -f tests/refGene.introns.bed > output_{}_dup.vcf' ::: {1..22}
#  
# ------------------------------------------------------------------------------

# 
# 1:22 %>% map(function(x) {
#   
#   write(header_file, glue('rival_cnvscore/svscore/input_{x}_del.vcf'))
#   
#   to_write %>% 
#     filter(chrom == x) %>%
#     arrange(start) %>%
#     mutate(ref = 'N', 
#            alt = '<DEL>', 
#            qual = 841.80, 
#            filter = '.', 
#            info = paste0('SVTYPE=DEL;END=', end), 
#            gt = 'GT', 
#            ratio = '1/1') %>%
#     select(chrom, start, id, ref, alt, qual, filter, info, gt, ratio) %>%
#     write_tsv(glue('rival_cnvscore/svscore/input_{x}_del.vcf'), col_names = FALSE, append = TRUE)
#   
# })
# 
# 
# 1:22 %>% map(function(x) {
#   
#   write(header_file, glue('rival_cnvscore/svscore/input_{x}_dup.vcf'))
#   
#   output_clinvar_duplication %>% 
#     filter(chrom == x) %>%
#     arrange(start) %>%
#     mutate(ref = 'N', 
#            alt = '<DUP>', 
#            qual = 841.80, 
#            filter = '.', 
#            info = paste0('SVTYPE=DUP;END=', end), 
#            gt = 'GT', 
#            ratio = '1/1') %>%
#     select(chrom, start, id, ref, alt, qual, filter, info, gt, ratio) %>%
#     write_tsv(glue('rival_cnvscore/svscore/input_{x}_dup.vcf'), col_names = FALSE, append = TRUE)
#   
# })


# ------------------------------------------------------------------------------
# HEADER FOR SVTOOL
# ------------------------------------------------------------------------------


header_file <- '##fileformat=VCFv4.2
##fileDate=20151026
##reference=
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description="Secondary breakend in a multi-line variants">
##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=EV,Number=.,Type=String,Description="Type of LUMPY evidence contributing to the variant call">
##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">
##INFO=<ID=PREND,Number=.,Type=String,Description="LUMPY probability curve of the END breakend">
##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=ALG,Number=1,Type=String,Description="Evidence PDF aggregation algorithm">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=PRPOS,Number=.,Type=String,Description="LUMPY probability curve of the POS breakend">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=AS,Number=A,Type=Integer,Description="Alternate allele split-read observation count, with partial observations recorded fractionally">
##FORMAT=<ID=RP,Number=1,Type=Integer,Description="Reference allele paired-end observation count, with partial observations recorded fractionally">
##FORMAT=<ID=AP,Number=A,Type=Integer,Description="Alternate allele paired-end observation count, with partial observations recorded fractionally">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count, with partial observations recorded fractionally">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of reference observations">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of alternate observations">
##FORMAT=<ID=AB,Number=A,Type=Float,Description="Allele balance, fraction of observations from alternate allele, QA/(QR+QA)">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
##FORMAT=<ID=SQ,Number=1,Type=Float,Description="Phred-scaled probability that this site is variant (non-reference in this sample">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">
##FORMAT=<ID=BD,Number=1,Type=Integer,Description="Amount of BED evidence supporting the variant">
##FORMAT=<ID=RS,Number=1,Type=Integer,Description="Reference allele split-read observation count, with partial observations recorded fractionally">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number of structural variant segment.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO'



