

library(magrittr)
library(tidyverse)
library(vroom)

args = commandArgs(trailingOnly = TRUE)
### Arguments script neeeds
# Rscript calc_avg_meth_jsd.R reference_dist_filename allele_counts_filename
# file_to_save


### Read in files and whitelisted CpGs
# reference distribution for calculating JDS
vroom::vroom(args[1]) -> cb_ref
# allele count file
vroom::vroom(args[2]) %>%
  group_by(target) %>%
  # Add in the number of CpGs in the target for average percent meth calculation
  # Also need to recalculated total reads per target since that column is
  # inaccurate in multiple places, but the read counts for the allele frequencies
  # are correct and I can just sum them; can be removed in future
  mutate(target_cpg_count = max(count_meth_cpgs),
         target_reads = sum(read_count)) %>%
  ungroup() %>% 
  # drop targets missing reads
  na.omit() %>%
  # drop targets without reads
  filter(target_reads > 0) -> data

### Calculate average percent methylation over the amplicon
data %>% 
  # calculate average percent methylation
  group_by(target, target_reads) %>%
  summarize(avg_perc_meth = (sum((count_meth_cpgs / target_cpg_count) * read_count) / target_reads) * 100) %>%
  ungroup() %>%
  distinct() -> perc_meth

### Calculate JSD for the amplicon
data %>% 
  ### add in the cord blood reference information
  left_join(cb_ref, by = c("target", "count_meth_cpgs")) %>%
  ### reformat for JSD function
  select(target, count_meth_cpgs, test = probability, ref_cb) %>%
  pivot_longer(c(test, ref_cb), names_to = 'set', values_to = 'prob') %>%
  pivot_wider(names_from = count_meth_cpgs, values_from = prob) %>%
  select(-set) %>%
  group_by(target) %>%
  nest() %>%
  ungroup() %>%
  ### reshape and use philentropy::JSD() to calculate the JSD
  mutate(data = map(data, ~ t(as.matrix(.))),
         data = map(data, ~ t(na.omit(.))),
         jsd = map(data, ~ sqrt(philentropy::JSD(., unit = "log2")))) %>%
  unnest(jsd) %>%
  select(-data) -> jsd

### Join together and save
left_join(perc_meth, jsd, by = "target") %>%
  write_tsv(args[3])