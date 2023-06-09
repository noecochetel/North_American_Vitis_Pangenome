library(tidyverse)
library(here)

# summarize ---------------------------------------------------------------
var_vec <- c("snp", "ins", "del", "mnp")

for (var in var_vec){
  
  snpeff_res_file <- list.files(path = "snpEff/results", pattern = paste0(var, ".ann.ez.txt"), full.names = T)
  snpeff_res_name <- gsub("^.*seqwish_k49.smooth.on.", "", snpeff_res_file)
  snpeff_res_name <- gsub("_ae.lv0.*$", "", snpeff_res_name)
  names(snpeff_res_file) <- snpeff_res_name
  
  map(snpeff_res_name, function(x) {
    
    tbx <- read_delim(snpeff_res_file[[x]], comment = '#', delim = '|', col_names = c("type", "impact", "locus"))
    
    # site level (comparable with the snpEff output)
    stats_tb <- tbx %>% 
      group_by(type, impact) %>% 
      summarize(nb = length(locus))
    
    # gene level
    gene_lvl_tb <- tbx %>% 
      separate_rows(locus, sep = '-') %>%
      separate_rows(locus, sep = '&') %>%
      separate_rows(type, sep = '&')
    
    gene_lvl_stat_tb <- gene_lvl_tb %>% 
      group_by(locus, type, impact) %>%
      summarize(nb = length(locus))
    
    # save
    save(stats_tb, file = here("snpEff", "results", paste0(x, "_ae.lv0.norm.",var,".ann.ez.stats_per_site.rda")))
    save(gene_lvl_stat_tb, file = here("snpEff", "results", paste0(x, "_ae.lv0.norm.",var,".ann.ez.stats_per_locus.rda")))
    
  })
}

# other
snpeff_res_file <- list.files(path = "snpEff/results", pattern = "other.ann.ez.txt", full.names = T)
snpeff_res_name <- gsub("^.*seqwish_k49.smooth.on.", "", snpeff_res_file)
snpeff_res_name <- gsub("_ae.lv0.*$", "", snpeff_res_name)
names(snpeff_res_file) <- snpeff_res_name

map(snpeff_res_name, function(x) {
  
  tbx <- read_delim(snpeff_res_file[[x]], delim = '|', col_names = c("type", "impact", "locus"))
  
  # site level (comparable with the snpEff output)
  stats_tb <- tbx %>% 
    group_by(type, impact) %>% 
    summarize(nb = length(locus))
  
  # gene level
  gene_lvl_tb <- tbx %>% 
    separate_rows(locus, sep = '-') %>%
    separate_rows(locus, sep = '&') %>%
    separate_rows(type, sep = '&')
  
  gene_lvl_stat_tb <- gene_lvl_tb %>% 
    group_by(locus, type, impact) %>%
    summarize(nb = length(locus))
  
  # save
  save(stats_tb, file = here("snpEff", "results", paste0(x, "_ae.lv0.norm.other.ann.ez.stats_per_site.rda")))
  save(gene_lvl_stat_tb, file = here("snpEff", "results", paste0(x, "_ae.lv0.norm.other.ann.ez.stats_per_locus.rda")))
  
})


# parse -------------------------------------------------------------------
# per site for basic stats
res_per_site_file <- list.files(path = "snpEff/results", pattern = ".ann.ez.stats_per_site.rda", full.names = T)
res_per_site_name <- gsub("^.*/", "", res_per_site_file)
res_per_site_name <- gsub(".ann.ez.stats_per_site.rda", "", res_per_site_name)
res_per_site_name <- gsub("_ae.lv0.norm", "", res_per_site_name)
names(res_per_site_file) <- res_per_site_name

res_per_site_list <- map(res_per_site_file, function(x){load(x); return(stats_tb)})

res_per_site_tb <- bind_rows(res_per_site_list, .id = 'name')

res_per_site_tb2 <- res_per_site_tb %>%
  mutate(chr = gsub("^.*chr", "chr", name),
         chr = gsub("[.].*$", "", chr),
         genome = gsub(".chr.*$", "", name),
         var_type = gsub("^.*[.]", "", name))

# from here we can know the number of sites per var type per genome
# let's save it
write_tsv(res_per_site_tb2, file = here("snpEff", "results", "summary_stats_per_site.txt"))

# per gene and per type decomposed for number/class of genes highly impacted
res_per_locus_file <- list.files(path = "snpEff/results", pattern = ".ann.ez.stats_per_locus.rda", full.names = T)
res_per_locus_name <- gsub("^.*/", "", res_per_locus_file)
res_per_locus_name <- gsub(".ann.ez.stats_per_locus.rda", "", res_per_locus_name)
res_per_locus_name <- gsub("_ae.lv0.norm", "", res_per_locus_name)
names(res_per_locus_file) <- res_per_locus_name

res_per_locus_list <- map(res_per_locus_file, function(x){load(x); return(gene_lvl_stat_tb)})

res_per_locus_tb <- bind_rows(res_per_locus_list, .id = 'name')

res_per_locus_tb2 <- res_per_locus_tb %>%
  mutate(chr = gsub("^.*chr", "chr", name),
         chr = gsub("[.].*$", "", chr),
         genome = gsub(".chr.*$", "", name),
         var_type = gsub("^.*[.]", "", name))

# let's save it
write_tsv(res_per_locus_tb2, file = here("snpEff", "results", "summary_stats_per_locus.txt"))


