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


# plot --------------------------------------------------------------------
# prepare variables -------------------------------------------------------
core_gene_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "_reclass_core.id", full.names = T)
core_gene_name <- gsub("^.*/", "", core_gene_file)
core_gene_name <- gsub("_reclass_core.id", "", core_gene_name)
names(core_gene_file) <- core_gene_name

disp_gene_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "_reclass_disp.id", full.names = T)
disp_gene_name <- gsub("^.*/", "", disp_gene_file)
disp_gene_name <- gsub("_reclass_disp.id", "", disp_gene_name)
names(disp_gene_file) <- disp_gene_name

priv_gene_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "_reclass_priv.id", full.names = T)
priv_gene_name <- gsub("^.*/", "", priv_gene_file)
priv_gene_name <- gsub("_reclass_priv.id", "", priv_gene_name)
names(priv_gene_file) <- priv_gene_name

ambi_gene_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "_reclass_ambi.id", full.names = T)
ambi_gene_name <- gsub("^.*/", "", ambi_gene_file)
ambi_gene_name <- gsub("_reclass_ambi.id", "", ambi_gene_name)
names(ambi_gene_file) <- ambi_gene_name


# load data ---------------------------------------------------------------
res_per_site_tb <- read_tsv(here("assembly", "pangenome", "snpEff", "summary_stats_per_site.txt"))
res_per_locus_tb <- read_tsv(here("assembly", "pangenome", "snpEff", "summary_stats_per_locus.txt")) %>%
  filter(! locus %in% c("CHR_START", "CHR_END"), type!="chromosome_number_variation")

core_gene_list <- map(core_gene_file, function(x) read_tsv(x, col_names = "gene_name"))
disp_gene_list <- map(disp_gene_file, function(x) read_tsv(x, col_names = "gene_name"))
priv_gene_list <- map(priv_gene_file, function(x) read_tsv(x, col_names = "gene_name"))
ambi_gene_list <- map(ambi_gene_file, function(x) read_tsv(x, col_names = "gene_name"))

core_gene_tb <- bind_rows(core_gene_list) %>% mutate(class = "core")
disp_gene_tb <- bind_rows(disp_gene_list) %>% mutate(class = "disp")
priv_gene_tb <- bind_rows(priv_gene_list) %>% mutate(class = "priv")
ambi_gene_tb <- bind_rows(ambi_gene_list) %>% mutate(class = "ambi")

gene_tb <- bind_rows(core_gene_tb, disp_gene_tb, priv_gene_tb, ambi_gene_tb)


# per locus ---------------------------------------------------------------
# add the gene classes
res_per_locus_class_tb <- left_join(res_per_locus_tb, gene_tb, by = c("locus" = "gene_name")) %>%
  mutate(group = paste(genome, class, sep = "_"))

# build a summary of the gene content per class to have a total for the perc
gene_tb2 <- gene_tb %>%
  mutate(genome = gsub(".chr.*$", "", gene_name)) %>%
  group_by(genome, class) %>%
  summarize(n = length(gene_name)) %>%
  mutate(group = paste(genome, class, sep = "_")) %>%
  ungroup() %>%
  dplyr::select(group, n)

# summarize
res_per_locus_class_tb2 <- res_per_locus_class_tb %>%
  filter(class != "ambi") %>%
  group_by(genome, class, var_type, type, impact) %>%
  summarize(locus_nb = length(unique(locus)), site_nb = sum(nb)) %>%
  mutate(group = paste(genome, class, sep = "_"))

# get site_per_gene
res_per_locus_class_tb3 <- left_join(res_per_locus_class_tb2, gene_tb2, by = "group") %>%
  mutate(perc = locus_nb * 100 / n, site_per_gene = site_nb/n)


# top5 --------------------------------------------------------------------
## let's select the top 5 categories overall (all var_types) per impact
site_per_gene_overall <- res_per_locus_class_tb3 %>%
  group_by(genome, impact, type) %>%
  summarize(site_per_gene = sum(site_per_gene)) %>%
  ungroup() %>%
  group_by(impact, type) %>%
  summarize(site_per_gene = mean(site_per_gene)) %>%
  ungroup() %>%
  group_by(impact) %>%
  arrange(-site_per_gene) %>%
  dplyr::slice(1:5) %>% # select the top 5 per impact
  mutate(top5var = paste(impact, type, sep = "_"))

## top5 variant types were defined per impact
## let's filter the data now
res_per_locus_class_4top5 <- res_per_locus_class_tb3 %>%
  mutate(top5var = paste(impact, type, sep = "_"))

res_per_locus_class_top5 <- res_per_locus_class_4top5 %>%
  filter(top5var %in% site_per_gene_overall$top5var) %>%
  mutate(impact = factor(impact, levels = c("LOW", "MODERATE", "HIGH", "MODIFIER"), labels = c("Low", "Moderate", "High", "Modifier")),
         var_type = factor(var_type, levels = c("snp", "ins", "del", "mnp", "other"), labels = c("SNP", "INS", "DEL", "MNP", "Other")),
         class = factor(class, levels = c("priv", "disp", "core"))) %>%
  dplyr::select(genome, class, var_type, type, impact, site_per_gene) %>%
  ungroup %>%
  group_by(impact) %>%
  complete(genome, class, var_type, type, fill = list(site_per_gene = 0)) %>%
  mutate(type = gsub("5_prime_UTR", "5'UTR", type),
         type = gsub("3_prime_UTR", "3'UTR", type),
         type = gsub("premature", "prem", type),
         type = gsub("_variant", "", type),
         type = gsub("_", " ", type))


