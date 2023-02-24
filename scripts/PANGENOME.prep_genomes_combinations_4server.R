library(tidyverse)
library(here)

tb <- read_tsv("genome_loc.txt", col_names = F)
tb2 <- tb %>%
  mutate(X1 = gsub("^.*/", "", X1))

genome_names_vec <- tb2 %>% select(X1) %>% pull()

map(seq(length(genome_names_vec)), function(genome_nb){
  
  comb_tb <- as_tibble(combn(genome_names_vec, genome_nb))
  
  map(seq(ncol(comb_tb)), function(column_x){
    
    combx <- comb_tb %>% select(all_of(column_x)) %>% pull()
    write_tsv(as_tibble(combx), here("path", paste0("genome_nb", genome_nb, ".", "comb", column_x, ".txt")), col_names = F)
    
  })
})