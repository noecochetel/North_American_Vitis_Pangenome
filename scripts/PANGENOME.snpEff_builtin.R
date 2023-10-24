library(tidyverse)
library(here)


# prepare variables -------------------------------------------------------
snpeff_file <- list.files(path = here("snpeff_res_dir"), pattern = "snpEff_summary.genes.txt", full.names = T)
snpeff_name <- gsub(".snpEff_summary.genes.txt", "", basename(snpeff_file))
names(snpeff_file) <- snpeff_name


# load data ---------------------------------------------------------------
snpeff_list <- map(snpeff_file, function(x) read_tsv(x, skip = 1))


# prepare data ------------------------------------------------------------
snpeff_tb <- bind_rows(snpeff_list, .id = "name")

# create a column with the variant type
snpeff_tb2 <- snpeff_tb %>% mutate(var = gsub("^.*[.]", "", name)) 

# separate impacts and effects
var_impact_tb <- snpeff_tb2[,c("GeneId", colnames(snpeff_tb2)[grep("variants_impact", colnames(snpeff_tb2))], "var")]
var_effect_tb <- snpeff_tb2[,c("GeneId", colnames(snpeff_tb2)[grep("variants_effect", colnames(snpeff_tb2))], "var")]

# prepare impacts
var_impact_tb2 <- var_impact_tb %>%
  pivot_longer(cols = -c(GeneId,var)) %>%
  mutate(name = gsub("variants_impact_", "", name))

var_impact_tb3 <- var_impact_tb2 %>%
  group_by(GeneId, var, name) %>%
  summarize(value = mean(value))

# prepare effects
var_effect_tb2 <- var_effect_tb %>%
  pivot_longer(cols = -c(GeneId,var)) %>%
  mutate(name = gsub("variants_effect_", "", name))

var_effect_tb3 <- var_effect_tb2 %>%
  group_by(GeneId, var, name) %>%
  summarize(value = mean(value))


# export parsed results ---------------------------------------------------
save(var_impact_tb3, file = here("snpeff_res_dir", "var_impact_builtin.rda"))
save(var_effect_tb3, file = here("snpeff_res_dir", "var_effect_builtin.rda"))


# calculate site frequency per gene locus ---------------------------------
load(here("snpeff_res_dir", "var_impact_builtin.rda"))
load(here("snpeff_res_dir", "var_effect_builtin.rda"))

# upload all the gene loci from the different pangenome classes
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

core_gene_list <- map(core_gene_file, function(x) read_tsv(x, col_names = "gene_name"))
disp_gene_list <- map(disp_gene_file, function(x) read_tsv(x, col_names = "gene_name"))
priv_gene_list <- map(priv_gene_file, function(x) read_tsv(x, col_names = "gene_name"))
ambi_gene_list <- map(ambi_gene_file, function(x) read_tsv(x, col_names = "gene_name"))

core_gene_tb <- bind_rows(core_gene_list) %>% mutate(class = "core")
disp_gene_tb <- bind_rows(disp_gene_list) %>% mutate(class = "disp")
priv_gene_tb <- bind_rows(priv_gene_list) %>% mutate(class = "priv")
ambi_gene_tb <- bind_rows(ambi_gene_list) %>% mutate(class = "ambi")

gene_tb <- bind_rows(core_gene_tb, disp_gene_tb, priv_gene_tb, ambi_gene_tb)

# build a summary of the gene content per class to have a total for the perc
gene_tb2 <- gene_tb %>%
  mutate(genome = gsub(".chr.*$", "", gene_name)) %>%
  group_by(genome, class) %>%
  summarize(n = length(gene_name)) %>%
  mutate(group = paste(genome, class, sep = "_")) %>%
  ungroup() %>%
  dplyr::select(group, n)

# class the impact
## add the gene classes
var_impact_gene_tb <- left_join(var_impact_tb3, gene_tb, by = c("GeneId" = "gene_name")) %>%
  mutate(genome = gsub(".chr.*$", "", GeneId))

## summarize
var_impact_gene_tb2 <- var_impact_gene_tb %>%
  filter(class != "ambi") %>%
  group_by(genome, class, var, name) %>%
  summarize(locus_nb = length(unique(GeneId)), site_nb = sum(value)) %>%
  mutate(group = paste(genome, class, sep = "_"))

## get site_per_gene
var_impact_gene_tb3 <- left_join(var_impact_gene_tb2, gene_tb2, by = "group") %>%
  mutate(site_per_gene = site_nb/n)

save(var_impact_gene_tb3, file = here("snpeff_res_dir", "var_impact_builtin_classed.rda"))


# class the effect
## add the gene classes
var_effect_gene_tb <- left_join(var_effect_tb3, gene_tb, by = c("GeneId" = "gene_name")) %>%
  mutate(genome = gsub(".chr.*$", "", GeneId),
         value = replace_na(value, 0))

## summarize
var_effect_gene_tb2 <- var_effect_gene_tb %>%
  filter(class != "ambi") %>%
  group_by(genome, class, var, name) %>%
  summarize(locus_nb = length(unique(GeneId)), site_nb = sum(value)) %>%
  mutate(group = paste(genome, class, sep = "_"))

## get site_per_gene
var_effect_gene_tb3 <- left_join(var_effect_gene_tb2, gene_tb2, by = "group") %>%
  mutate(site_per_gene = site_nb/n)
  

save(var_effect_gene_tb3, file = here("snpeff_res_dir", "var_effect_builtin_classed.rda"))


