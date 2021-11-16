
library(magrittr)

args = commandArgs(trailingOnly = TRUE)
### Arguments script neeeds
# Rscript chip_vcf_to_tsv.R vcf_file quality read_depth filename_to_write

vroom::vroom(args[1],
             delim = '\t', comment = '#',
             col_names = c('chr', 'pos', 'id', 'ref', 'alt', 'qual',
                           'filter', 'info', 'format', 'sample')) %>%
  tidyr::separate(sample, into = c('genotype', 'allele_depth', 'read_depth', 
                                   'genotype_quality', 
                                   'phred-scaled_likelihood', 'statistics'), 
                  sep = ':', convert = T) %>%
  dplyr::mutate(qual = as.numeric(qual), 
                read_depth = as.numeric(read_depth)) %>%
  dplyr::filter(qual >= as.integer(args[2]), 
                read_depth >= as.integer(args[3])) %>%
  dplyr::mutate(alt = stringr::str_remove(alt, ',<NON_REF>'),
                ann = stringr::str_remove_all(stringr::str_extract(info, 'ANN=(.*?)\\|\\|\\|'),
                                              '^ANN=|\\|\\|\\|$')) %>%
  tidyr::separate(ann, sep = '\\|',
                  into = c('allele', 'annotation', 'annotation_impact', 
                           'gene_name', 'gene_id', 'feature_type', 'feature_id', 
                           'transcript_biotype', 'rank', 'hgvs.c', 'hgvs.p', 
                           'cDNA', 'cds', 'aa', 'distance', 'other')) %>%
  dplyr::select(chr, pos, ref, alt, rsid = id, annotation, annotation_impact, 
                gene = gene_name, aa_change = `hgvs.c`, block = rank) -> filtered_data

readr::write_tsv(filtered_data, args[4])