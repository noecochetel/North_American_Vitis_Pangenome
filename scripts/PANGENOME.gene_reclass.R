library(tidyverse)
library(here)


# prepare variable --------------------------------------------------------
# intersect
intersect_list_file <- list.files(path = here("intersect"), pattern = "thru.all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_gene_ez.txt", full.names = T)
intersect_list_name <- gsub("^.*/", "", intersect_list_file)
intersect_list_name <- gsub(".thru.all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_gene_ez.txt", "", intersect_list_name)
names(intersect_list_file) <- intersect_list_name

# seq type
core_seq_list_file <- list.files(path = here("fasta"), pattern = "all.on.all.wfmash_s10000p85n1.*.seqwish_k49.smooth.join.core.fasta.len", full.names = T)
disp_seq_list_file <- list.files(path = here("fasta"), pattern = "all.on.all.wfmash_s10000p85n1.*.seqwish_k49.smooth.join.disp.fasta.len", full.names = T)
priv_seq_list_file <- list.files(path = here("fasta"), pattern = "all.on.all.wfmash_s10000p85n1.*.seqwish_k49.smooth.join.priv.fasta.len", full.names = T)


# load data ---------------------------------------------------------------
# intersect
intersect_list <- map(intersect_list_file, function(x){read_tsv(x, col_names = c("gene_name", "seg_name", "int_len"))})

# seq type
core_seq_list <- map(core_seq_list_file, function(x){read_tsv(x, col_names = c("seg_name", "seg_len"))})
core_seq_tb <- bind_rows(core_seq_list)
disp_seq_list <- map(disp_seq_list_file, function(x){read_tsv(x, col_names = c("seg_name", "seg_len"))})
disp_seq_tb <- bind_rows(disp_seq_list)
priv_seq_list <- map(priv_seq_list_file, function(x){read_tsv(x, col_names = c("seg_name", "seg_len"))})
priv_seq_tb <- bind_rows(priv_seq_list)


# prepare data ------------------------------------------------------------
# attribute a genome type and a genic type
intersect_list2 <- map(names(intersect_list), function(x){
  
  tbx <- intersect_list[[x]]
  
  tbx[tbx$seg_name %in% core_seq_tb$seg_name, "seq_type"] <- "core"
  tbx[tbx$seg_name %in% disp_seq_tb$seg_name, "seq_type"] <- "disp"
  tbx[tbx$seg_name %in% priv_seq_tb$seg_name, "seq_type"] <- "priv"
  
  return(tbx)
  
})

names(intersect_list2) <- names(intersect_list)

# all SVs
## calculate perc
intersect_list3 <- map(intersect_list2, function(x){
  
  x %>%
    complete(gene_name, seq_type, fill = list(seg_name = "segadd", int_len = 0)) %>%
    group_by(gene_name, seq_type) %>%
    mutate(sum_len = sum(int_len)) %>%
    ungroup() %>%
    group_by(gene_name) %>%
    mutate(tot_len = sum(int_len)) %>%
    ungroup() %>%
    mutate(perc = sum_len * 100 / tot_len) %>%
    select(gene_name, seq_type, perc) %>%
    distinct()
  
})

## summarize
intersect_stat_list <- map(intersect_list3, function(x){
  
  x %>%
    group_by(seq_type) %>%
    summarize(mean = mean(perc), median = median(perc), sd = sd(perc), n = length(gene_name), se = sd/sqrt(n)) 
  
})

## merge
intersect_stat_tb <- bind_rows(intersect_stat_list, .id = "name")

## export
save(intersect_list3, file = here("intersect", "all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_list.rda"))
save(intersect_stat_tb, file = here("intersect", "all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_stat_tb.rda"))



# class the genes based on the seq pangenome ------------------------------
# load data
load(file = here("intersect", "all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_list.rda"))

# prepare data
ltb <- map(intersect_list3, function(x){
  
  x %>%
    pivot_wider(names_from = seq_type, values_from = perc)
  
})

# class
ltb2 <- map(ltb, function(x){
  
  x %>%
    rowwise() %>%
    mutate(class = if(core > 80){"core"}else if(core + disp > 80){"disp"}else if(disp > 80){"disp"}else if(priv > 80){"priv"}else{"ambi"})
  
})

map(ltb2, function(x){table(x$class)})

# export
map(names(ltb2), function(x){
  
  tbx <- ltb2[[x]]
  core_id <- tbx %>% filter(class == "core") %>% dplyr::select(gene_name)
  disp_id <- tbx %>% filter(class == "disp") %>% dplyr::select(gene_name)
  priv_id <- tbx %>% filter(class == "priv") %>% dplyr::select(gene_name)
  ambi_id <- tbx %>% filter(class == "ambi") %>% dplyr::select(gene_name)
  
  write_tsv(core_id, here("intersect", paste0(x, "_reclass_core.id")), col_names = F)
  write_tsv(disp_id, here("intersect", paste0(x, "_reclass_disp.id")), col_names = F)
  write_tsv(priv_id, here("intersect", paste0(x, "_reclass_priv.id")), col_names = F)
  write_tsv(priv_id, here("intersect", paste0(x, "_reclass_ambi.id")), col_names = F)
  
})

