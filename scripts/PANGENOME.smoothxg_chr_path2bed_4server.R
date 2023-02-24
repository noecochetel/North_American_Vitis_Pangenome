library(tidyverse)
library(here)

# prepare variable --------------------------------------------------------
path_list_file <- list.files(path = here("path"), pattern = "chr.*.thru.all.on.all.wfmash_s10000p85n1.chr.*.seqwish_k49.smooth.join.path", full.names = T)
path_list_name <- gsub("^.*/", "", path_list_file)
path_list_name <- gsub(".thru.all.on.all.wfmash_s10000p85n1.chr.*.seqwish_k49.smooth.join.path", "", path_list_name)
names(path_list_file) <- path_list_name

len_list_file <- list.files(path = here("fasta"), pattern = "all.on.all.wfmash_s10000p85n1.chr.*.seqwish_k49.smooth.join.fasta.len", full.names = T)
len_list_name <- gsub("^.*/", "", len_list_file)
len_list_name <- gsub(".seqwish_k49.smooth.join.fasta.len", "", len_list_name)
names(len_list_file) <- len_list_name

# load data ---------------------------------------------------------------
path_list <- map(path_list_file, function(x){read_tsv(x, col_names = c("seg_name", "strand"))})
len_list <- map(len_list_file, function(x){read_tsv(x, col_names = c("seg_name", "len"))})


# create bed files of the paths -------------------------------------------
map(names(path_list), function(x){

  chrx <- gsub("^.*chr", "chr", x)
  path_listx <- path_list[[x]]
  len_tbx <- len_list[[paste0("all.on.all.wfmash_s10000p85n1.", chrx)]]
  tbx <- left_join(path_listx, len_tbx, by = "seg_name")
  tbx2 <- tbx %>%
    mutate(cumsum = cumsum(len),
           start = cumsum - len,
           seq_name = x,
           score = ".")
  
  tbx3 <- tbx2[,c("seq_name", "start", "cumsum", "seg_name", "score", "strand")]
  
  write_tsv(tbx3, here("bed", paste0(x, ".thru.all.on.all.wfmash_s10000p85n1.", chrx, ".seqwish_k49.smooth.join.bed")), col_names = F)

  })