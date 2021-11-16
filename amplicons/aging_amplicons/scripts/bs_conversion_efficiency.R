
library(magrittr)

args = commandArgs(trailingOnly = TRUE)
### Arguments script neeeds
# Rscript bs_conversion_efficiency.R bismark_splitting_report.txt file_to_save

vroom::vroom(args[1], id = 'file_path',
      delim = '\t', skip = 13, col_names = c('parameter', 'percentage')) %>%
  dplyr::filter(stringr::str_detect(percentage, '%'), 
                !stringr::str_detect(parameter, 'CpG')) %>%
  dplyr::mutate(context = stringr::str_extract(parameter, 'CH[GH]'),
                percentage = as.numeric(stringr::str_remove(percentage, '%'))) %>%
  dplyr::select(context, percentage) -> non_cg_meth

# save
readr::write_tsv(non_cg_meth, args[2])