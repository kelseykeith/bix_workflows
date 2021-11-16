
library(magrittr)
library(tidyverse)
library(vroom)


args = commandArgs(trailingOnly = TRUE)
### Arguments script neeeds
# Rscript get_meth_cpg_counts.R CpG_OT_filename CpG_OB_filename 
# cpg_whitelist.tsv filename_to_write_summarized_counts.tsv


### Read in files and whitelisted CpGs
vroom(args[1], delim = '\t', skip = 1, 
      col_names = c('read_id', 'strand', 'chr', 'position', 'meth_status')) %>%
  select(-strand) %>%
  distinct() -> cpg_ot

vroom(args[2], 
      delim = '\t', skip = 1,
      col_names = c('read_id', 'strand', 'chr', 'position', 'meth_status')) %>%
  select(-strand) %>%
  distinct() -> cpg_ob

vroom(args[3]) -> cpg_whitelist

### Count the number of methylated CpGs per read
rbind(cpg_ot, cpg_ob) %>%
  unite(chr_base, c(chr, position), sep = '_') %>%
  left_join(cpg_whitelist, by = 'chr_base') %>%
  na.omit() %>%
  add_count(read_id, target) %>%
  filter(n == target_cpg_count) %>%
  select(-n) %>%
  count(read_id, target, target_cpg_count, meth_status) %>%
  pivot_wider(names_from = meth_status, values_from = n) %>%
  replace_na(list(z = 0, Z = 0)) %>%
  select(read_id, target, target_cpg_count, 
         count_meth_cpgs = Z) -> read_counts
#write_tsv(read_counts, args[4])

### Summarize counts by target and number of methylated CpGs
tbl <- NULL
for (i in unique(cpg_whitelist$target)) {
  rbind(tibble(target = i,
               count_meth_cpgs = 0:unique(filter(cpg_whitelist, 
                                                 target == i)$target_cpg_count)),
        tbl) -> tbl
}

read_counts %>%
  add_count(target, name = 'target_reads') %>% 
  count(target, target_reads, count_meth_cpgs, 
        name = 'read_count') %>%
  mutate(probability = read_count / target_reads) %>%
  full_join(tbl, by = c('target', 'count_meth_cpgs')) %>%
  arrange(target, count_meth_cpgs) %>%
  fill(target_reads) %>%
  replace_na(list(read_count = 0, probability = 0)) -> counts
write_tsv(counts, args[4])

