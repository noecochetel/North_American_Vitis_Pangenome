library(tidyverse)
library(here)

# cd /DATA13/Projects/North_American_Species/assembly/pangenome


# calculate the seq class composition for each gene per combination --------------------------
# prep intersect var
intersect_list_file <- list.files(path = here("/DATA8/Projects/North_American_Species/assembly/pangenome/intersect"), pattern = "thru.all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_gene_ez.txt.gz", full.names = T)
intersect_list_name <- gsub("^.*/", "", intersect_list_file)
intersect_list_name <- gsub(".thru.all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_gene_ez.txt.gz", "", intersect_list_name)
names(intersect_list_file) <- intersect_list_name

# load intersect
# genes are composed of nodes from the final pangenome
# each node has a specific length
# each node is classed based on the genomes in the comb
intersect_list <- map(intersect_list_file, function(x){read.table(x, header = F, col.names = c("gene_name", "seg_name", "int_len"))})

# prep comb var
comb_list <- list.files(path = here("/DATA8/Projects/North_American_Species/assembly/pangenome/path"), pattern = "genome_nb", full.names = T)

comb_list <- comb_list[-grep("[.]path", comb_list)]
comb_list <- comb_list[-grep("[.]freq", comb_list)]
comb_list <- comb_list[-grep("[.]len", comb_list)]
comb_list <- comb_list[-grep("[.]stats", comb_list)]
comb_list <- comb_list[-grep("pdf", comb_list)]
comb_list <- comb_list[-grep("_freq", comb_list)]
comb_name <- gsub("/DATA8/Projects/North_American_Species/assembly/pangenome/path/", "", comb_list)
comb_name <- gsub(".txt", "", comb_name)

# for modeling
comb_name <- comb_name[grep("genome_nb2.comb", comb_name)]

for (combx in comb_name){
  
  # prep seq type var
  core_seq_file <- paste0("/DATA13/Projects/ncochtl_sva/path/", combx, ".core.len.gz")
  disp_seq_file <- paste0("/DATA13/Projects/ncochtl_sva/path/", combx, ".disp.len.gz")
  priv_seq_file <- paste0("/DATA13/Projects/ncochtl_sva/path/", combx, ".priv.len.gz")
  
  # prep genomes in comb var
  genome_vec_file <- paste0("/DATA8/Projects/North_American_Species/assembly/pangenome/path/", combx, ".txt")
  
  # load seq type
  core_seq_tb <- read.table(core_seq_file, header = F, col.names = c("seg_name", "seg_len"))
  disp_seq_tb <- read.table(disp_seq_file, header = F, col.names = c("seg_name", "seg_len"))
  priv_seq_tb <- read.table(priv_seq_file, header = F, col.names = c("seg_name", "seg_len"))
  
  # load genomes in comb
  genome_vec <- read_tsv(genome_vec_file, col_names = "name")
  
  # attribute a seq type to each path
  intersect_list2 <- map(genome_vec$name, function(x){
    
    tbx <- as_tibble(intersect_list[[x]])
    
    tbx[tbx$seg_name %in% core_seq_tb$seg_name, "seq_type"] <- "core"
    tbx[tbx$seg_name %in% disp_seq_tb$seg_name, "seq_type"] <- "disp"
    tbx[tbx$seg_name %in% priv_seq_tb$seg_name, "seq_type"] <- "priv"
    
    return(tbx)
    
  })
  
  names(intersect_list2) <- genome_vec$name
  
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
  
  
  ## export
  save(intersect_list3, file = here("intersect", paste0(combx, "seq_type_perGene.rda")))
}


# # class the genes based on the seq pangenome ------------------------------
# comb_list_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "seq_type_perGene.rda", full.names = T)
# comb_list_name <- gsub("^.*/", "", comb_list_file)
# comb_list_name <- gsub("seq_type_perGene.rda", "", comb_list_name)
# names(comb_list_file) <- comb_list_name
# 
# comb_vec <- unique(gsub("[.]comb.*$", "", comb_list_name))
# 
# for (comb_vecx in comb_vec[6:9]){
#   
#   comb_list_filex <- comb_list_file[grep(paste0(comb_vecx, ".comb"), comb_list_name)]
#   
#   # load data
#   comb_listx <- map(comb_list_filex, function(x) {
#     
#     load(x)
#     tbx <- bind_rows(intersect_list3, .id = "name")
#     
#   })
#   
#   # prepare data
#   ltbx <- map(comb_listx, function(x){
#     
#     x %>%
#       pivot_wider(names_from = seq_type, values_from = perc)
#     
#   })
#   
#   # class
#   if(comb_vecx == "genome_nb1"){
#     
#     ltb2x <- map(ltbx, function(x){
#       
#       x %>%
#         rowwise() %>%
#         mutate(class = "priv")
#       
#     })
#     
#   } else if (comb_vecx == "genome_nb2"){
#     
#     ltb2x <- map(ltbx, function(x){
#     
#     x %>%
#       rowwise() %>%
#       mutate(class = if(core > 80){"core"}else if(priv > 80){"priv"}else{"ambi"})
#       
#     })
#     
#   } else {
#     
#     ltb2x <- map(ltbx, function(x){
#     
#     x %>%
#       rowwise() %>%
#       mutate(class = if(core > 80){"core"}else if(core + disp > 80){"disp"}else if(disp > 80){"disp"}else if(priv > 80){"priv"}else{"ambi"})
#       
#     })
#     
#   }
#   
#   # export
#   map(names(ltb2x), function(x){
#     
#     genome_nbx <- gsub("[.]comb.*$", "", x)
#     genome_nbx <- gsub("genome_nb", "", genome_nbx)
#     combx <- gsub("^.*[.]comb", "", x)
#     
#     tbx <- ltb2x[[x]] %>%
#       group_by(name, class) %>%
#       summarize(n = length(unique(gene_name))) %>%
#       mutate(genome_nb = genome_nbx,
#              comb = combx) %>%
#       ungroup()
#     
#     write_tsv(tbx, here("assembly", "pangenome", "intersect", "modeling", paste0(x, "_reclass_stats.txt")))
#     
#   })
#   
# }


# class the genes based on the perc ------------------------------------------------------------
comb_list_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "seq_type_perGene.rda", full.names = T)
comb_list_name <- gsub("^.*/", "", comb_list_file)
comb_list_name <- gsub("seq_type_perGene.rda", "", comb_list_name)
names(comb_list_file) <- comb_list_name

comb_vec <- unique(gsub("[.]comb.*$", "", comb_list_name))

# selection of genome_nb2
comb_vecx <- comb_vec[2]

comb_list_filex <- comb_list_file[grep(paste0(comb_vecx, ".comb"), comb_list_name)]

# load data
comb_listx <- map(comb_list_filex, function(x) {
  
  load(x)
  tbx <- bind_rows(intersect_list3, .id = "name")
  
})

# prepare data
ltbx <- map(comb_listx, function(x){
  
  x %>%
    pivot_wider(names_from = seq_type, values_from = perc)
  
})

# class
ltb2x <- map(ltbx, function(x){
  
  x %>%
    rowwise() %>%
    mutate(class = if(core > 80){1}else{0})
  
})
    
# parse
ltb3x <- map(ltb2x, function(x){
  
  genome_vec <- unique(x$name)
  
  x %>%
    rowwise() %>%
    mutate(other_genome = genome_vec[genome_vec != name]) %>%
    dplyr::select(name, gene_name, class, other_genome)
  
}) 


bin_tb <- bind_rows(ltb3x) %>%
  pivot_wider(names_from = other_genome, values_from = class) %>%
  replace(is.na(.), 1)


# export the binary matrix
write_tsv(bin_tb, here("assembly", "pangenome", "intersect", "modeling", "binary_matrix.txt"))


# model -------------------------------------------------------------------
# list comb
comb_list_file <- list.files(path = here("assembly", "pangenome", "comb"), pattern = "comb", full.names = T)
comb_list_name <- gsub("^.*/", "", comb_list_file)
comb_list_name <- gsub(".txt", "", comb_list_name)
names(comb_list_file) <- comb_list_name

comb_list <- map(comb_list_file, function(x){read_tsv(x, col_names = F)})

# re-open the matrix
bin_tb <- read_tsv(here("assembly", "pangenome", "intersect", "modeling", "binary_matrix.txt"))

lpan_stats <- map(comb_list, function(combx){
  
  genome_vec <- combx$X1
  bin_tbx <- bin_tb %>%
    filter(name %in% genome_vec) %>%
    dplyr::select(gene_name, all_of(genome_vec)) %>%
    mutate(sum = rowSums(across(where(is.numeric)))) %>%
    rowwise() %>%
    mutate(class = if(sum == 1){"priv"}else if (sum>1 & sum<length(genome_vec)){"disp"}else{"core"})
  stats_tbx <- bin_tbx %>%
    rowwise() %>%
    mutate(sum_norm = if(class == "core"){1 / length(genome_vec)}else if(class == "disp"){1 / sum} else {sum}) %>%
    ungroup() %>%
    group_by(class) %>%
    summarize(n = length(class),
              n_norm = round(sum(sum_norm)))
  stats_tbx2 <- stats_tbx %>%
    add_row(class = "pan", n = sum(stats_tbx$n), n_norm = sum(stats_tbx$n_norm))
  
  return(stats_tbx2)
  
})

pan_stats_tb <- bind_rows(lpan_stats, .id = "comb") %>%
  mutate(genome_nb = gsub("[.]comb.*$", "", comb),
         genome_nb = gsub("genome_nb", "", genome_nb),
         comb_nb = gsub("^.*comb", "", comb))

write_tsv(pan_stats_tb, here("assembly", "pangenome", "intersect", "modeling", "pangenome_modeling_stats.txt"))


# # plot --------------------------------------------------------------------
# # load data
# tb_list_file <- list.files(path = here("assembly", "pangenome", "intersect", "modeling"), pattern = "_reclass_stats.txt", full.names = T)
# 
# ltb <- map(tb_list_file, function(x) read_tsv(x))
# 
# # prepare data
# tb <- bind_rows(ltb)
# 
# tb2 <- tb %>%
#   filter(class != "ambi")
# 
# tb_priv1 <- tb2 %>%
#   filter(genome_nb == 1) %>%
#   mutate(n = 0)
# 
# tb_core1 <- tb2 %>%
#   filter(genome_nb == 1) %>%
#   mutate(class = "core")
# 
# tb_disp1 <- tb2 %>%
#   filter(genome_nb == 1) %>%
#   mutate(n = 0, class = "disp")
# 
# tb_pan1 <- bind_rows(tb_priv1, tb_core1, tb_disp1)
# 
# tb_disp2 <- tb2 %>%
#   filter(genome_nb == 2, class == "core") %>%
#   mutate(n = 0, class = "disp")
# 
# tb_notdisp2 <- tb2 %>%
#   filter(genome_nb == 2)
# 
# tb_pan2 <- bind_rows(tb_disp2, tb_notdisp2)
# 
# tb_gt2 <- tb2 %>%
#   filter(! genome_nb %in% c(1:2))
# 
# full_tb <- bind_rows(tb_pan1, tb_pan2, tb_gt2)
# 
# # understand the calculations
# demo_tb <- tibble(genome = rep(letters[1:3], each = 4),
#                   gene = rep(paste0("gene_", seq(4)), 3),
#                   class = c("core", "dispab", "dispac", "priv", "core", "dispab", "dispbc", "priv", "core", "dispac", "dispbc", "priv"))
# 
# 
# 
# ## prepare
# ### for genome number = 1, the pangenes and coregenes are identical
# res_comb1_tb <- full_tb %>%
#   filter(genome_nb == 1) %>%
#   pivot_wider(names_from = class, values_from = n) %>%
#   group_by(name, comb, genome_nb) %>%
#   summarise(pangene = core, coregene = core, dispgene = disp, privgene = priv) %>%
#   ungroup()
# 
# ### for each comp and genome, the sum is the genome coregenes + the uniq of the genomes compared
# res_comb2_tb <- full_tb %>%
#   filter(genome_nb == 2) %>%
#   pivot_wider(names_from = class, values_from = n) %>%
#   group_by(comb, genome_nb) %>%
#   mutate(priv_sum = sum(priv)) %>%
#   ungroup() %>%
#   group_by(name, comb, genome_nb) %>%
#   summarise(pangene = sum(core, priv_sum), coregene = core, dispgene = disp, privgene = priv) %>%
#   ungroup()
# 
# ### for each comp and genome, the sum is the genome coregenes + the uniq of the genomes compared
# res_combgt2_tb <- full_tb %>%
#   filter(genome_nb > 2) %>%
#   pivot_wider(names_from = class, values_from = n) %>%
#   group_by(comb, genome_nb) %>%
#   mutate(priv_sum = sum(priv)) %>%
#   ungroup() %>%
#   group_by(name, comb, genome_nb) %>%
#   summarise(pangene = sum(core, disp, priv_sum), coregene = core, dispgene = disp, privgene = priv) %>%
#   ungroup()
# 
# ## merge
# res_tb <- bind_rows(res_comb1_tb, res_comb2_tb, res_combgt2_tb) %>%
#   pivot_longer(cols = c(pangene, coregene, dispgene, privgene), names_to = "class") %>%
#   mutate(class = factor(class, levels = c("pangene", "coregene", "dispgene", "privgene")))
# 
# avg_tb <- res_tb %>%
#   group_by(name, genome_nb, class) %>%
#   summarise(mean = mean(value), sd = sd(value)) %>%
#   ungroup()
# 
# 
# 
# # plot
# p <- ggplot() +
#   geom_point(data = avg_tb, mapping = aes(x = genome_nb, y = mean, color = class), alpha = 0.5, size = 1) +
#   geom_smooth(data = res_tb, mapping = aes(x = genome_nb, y = value, color = class), method = "loess", show.legend = F, size = 0.5) +
#   # scale_x_continuous(breaks = seq(9)) +
#   # scale_y_continuous(breaks = seq(0, 80000, 10000)) +
#   labs(y = "Number of genes", x = "Number of genomes") +
#   scale_color_manual(name = "", labels = c("Pan-genome", "Core genome", "Dispensable genome", "Private genome"), values=c("#F01579", "#90E711", "#F4D213", "#A657F5")) +
#   guides(color = guide_legend(override.aes = list(alpha = 1)))
# 
# p <- p + theme(plot.title = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
#                axis.title.x = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
#                axis.title.y = element_text(color="black", size=7, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
#                axis.text.x = element_text(color="black", size=5, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
#                axis.text.y = element_text(color="black", size=5, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
#                axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
#                axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
#                axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
#                axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
#                plot.background = element_blank(),
#                panel.background = element_blank(),
#                panel.grid.major.y = element_blank(),
#                panel.grid.major.x = element_blank(),
#                panel.grid.minor.y = element_blank(),
#                panel.grid.minor.x = element_blank(),
#                plot.margin=unit(c(1,1,1,1),"mm"),
#                legend.title=element_text(size=8), 
#                legend.text=element_text(size=6),
#                legend.background = element_blank(),
#                legend.key.size = unit(2, "mm"),
#                legend.key=element_blank(),
#                legend.position = c(1, 0.3),
#                legend.justification = c("right", "bottom"))
# 
# p
# 
# ggsave(p, filename = here("plots", "assembly", "pangenome", "intersect", "gene_reclass_pangenome_modeling.pdf"), width = 3.5, height = 2.5)


# plot model --------------------------------------------------------------------
# load data
tb <- read_tsv(here("assembly", "pangenome", "intersect", "modeling", "pangenome_modeling_stats.txt"))

# prepare data
tb_nb1_core <- tb %>%
  filter(genome_nb == 1) %>%
  mutate(class = "core") %>%
  distinct()

tb_nb1_pan <- tb_nb1_core %>%
  mutate(class = "pan")

tb_nb1_disp <- tb %>%
  filter(genome_nb == 1) %>%
  mutate(class = "disp", n=0, n_norm=0) %>%
  distinct()

tb_nb1_priv <- tb %>%
  filter(genome_nb == 1) %>%
  mutate(class = "priv", n=0, n_norm=0) %>%
  distinct()

tb_nb1 <- bind_rows(tb_nb1_core, tb_nb1_pan, tb_nb1_disp, tb_nb1_priv)

tb_nb2_notdisp <- tb %>%
  filter(genome_nb == 2) %>%
  distinct()

tb_nb2_disp <- tb %>%
  filter(genome_nb == 2, class=="core") %>%
  mutate(class = "disp", n =0, n_norm=0) %>%
  distinct()

tb_nb2 <- bind_rows(tb_nb2_notdisp, tb_nb2_disp)

tb_not_nb1nb2 <- tb %>%
  filter(!genome_nb %in% c(1,2))

tb2 <- bind_rows(tb_nb1, tb_nb2, tb_not_nb1nb2) %>%
  mutate(class = factor(class, levels = c("pan", "core", "disp", "priv"), labels = c("Pan-genome", "Core genome", "Dispensable genome", "Private genome")))

avg_tb <- tb2 %>%
  group_by(genome_nb, class) %>%
  summarise(mean = mean(n_norm), sd = sd(n_norm)) %>%
  ungroup()

# for in-text stats
avg_tb %>%
  filter(genome_nb %in% c(8,9))

# plot
p <- ggplot() +
  geom_point(data = tb2, mapping = aes(x = genome_nb, y = n_norm, color = class), alpha = 0.5, size = 1) +
  geom_smooth(data = tb2, mapping = aes(x = genome_nb, y = n_norm, color = class), method = "loess", show.legend = F, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(9)) +
  scale_y_continuous(breaks = seq(0, 125000, 20000)) +
  labs(y = "Number of genes", x = "Number of genomes") +
  scale_color_manual(values=c("#E7298A", "#66A61E", "#E6AB02", "#7570B3"), name = "") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

p <- p + theme(plot.title = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
               axis.title.x = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
               axis.title.y = element_text(color="black", size=7, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
               axis.text.x = element_text(color="black", size=5, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
               axis.text.y = element_text(color="black", size=5, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
               axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
               axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
               axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
               axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
               plot.background = element_blank(),
               panel.background = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.grid.minor.x = element_blank(),
               plot.margin=unit(c(1,1,1,1),"mm"),
               legend.title=element_text(size=8), 
               legend.text=element_text(size=6),
               legend.background = element_blank(),
               legend.key.size = unit(2, "mm"),
               legend.key=element_blank(),
               legend.position = c(1, 0.5),
               legend.justification = c("right", "bottom"))

p

ggsave(p, filename = here("plots", "assembly", "pangenome", "intersect", "gene_reclass_pangenome_new_modeling.pdf"), width = 3.5, height = 2.5)

save(p, file = here("plots", "assembly", "pangenome", "intersect", "gene_reclass_pangenome_new_modeling.rda"))


