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
# attribute a genome type and a genic type to each path
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
load(file = here("assembly", "pangenome", "intersect", "all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_list.rda"))

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
  
  write_tsv(core_id, here("assembly", "pangenome", "intersect", paste0(x, "_reclass_core.id")), col_names = F)
  write_tsv(disp_id, here("assembly", "pangenome", "intersect", paste0(x, "_reclass_disp.id")), col_names = F)
  write_tsv(priv_id, here("assembly", "pangenome", "intersect", paste0(x, "_reclass_priv.id")), col_names = F)
  
})

# try to class the ambi
lambi_tb <- map(ltb2, function(x){
  
  x %>%
    filter(class == "ambi")
  
})

ambi_tb <- bind_rows(lambi_tb, .id = "name")

ambi_tb2 <- ambi_tb %>%
  rowwise() %>%
  mutate(ambi_class = if(core > 60) {"ambi_core"} else if (priv > 60){"ambi_priv"} else if(core + priv > 80){"ambi_corepriv"} else if (disp > 60){"ambi_disp"} else if (disp + priv > 80) {"ambi_dispriv"} else {"ambi"}) %>%
  ungroup()

# export
map(unique(ambi_tb2$name), function(x){
  
  tbx <- ambi_tb2 %>% filter(name == x)
  
  map(unique(tbx$ambi_class), function(classx){
    
    classx_id <- tbx %>% filter(ambi_class == classx) %>% dplyr::select(gene_name)
    write_tsv(classx_id, here("assembly", "pangenome", "intersect", paste0(x, "_reclass_", classx, ".id")), col_names = F)
    
  }) 
})


# annotation per reclass --------------------------------------------------
# prepare variables
id_core_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "_reclass_core.id", full.names = T)
id_core_name <- gsub("^.*/", "", id_core_file)
id_core_name <- gsub("_reclass_core.id", "", id_core_name)
names(id_core_file) <- id_core_name

id_disp_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "_reclass_disp.id", full.names = T)
id_disp_name <- gsub("^.*/", "", id_disp_file)
id_disp_name <- gsub("_reclass_disp.id", "", id_disp_name)
names(id_disp_file) <- id_disp_name

id_priv_file <- list.files(path = here("assembly", "pangenome", "intersect"), pattern = "_reclass_priv.id", full.names = T)
id_priv_name <- gsub("^.*/", "", id_priv_file)
id_priv_name <- gsub("_reclass_priv.id", "", id_priv_name)
names(id_priv_file) <- id_priv_name

annot_list_file <- list.files(path = here("annotation"), pattern = "[.]annotation.txt", full.names = T)
annot_list_name <- gsub("^.*/", "", annot_list_file)
annot_list_name <- gsub(".annotation.txt", "", annot_list_name)
names(annot_list_file) <- annot_list_name

# load data
id_core_list <- map(id_core_file, function(x){read_tsv(x, col_names = "gene_name")})
id_disp_list <- map(id_disp_file, function(x){read_tsv(x, col_names = "gene_name")})
id_priv_list <- map(id_priv_file, function(x){read_tsv(x, col_names = "gene_name")})

annot_list <- map(annot_list_file, function(x){read_tsv(x, col_names = F)})

# prepare data
annot_list2 <- map(annot_list, function(x){x %>%
  mutate(gene_name = gsub("[.]t.*$", "", X2)) %>%
    dplyr::select(gene_name, X1, X6, X12)})

# merge
annot_core <- map(names(id_core_list), function(x){
  
  left_join(id_core_list[[x]], annot_list2[[x]], by = "gene_name")
  
})
names(annot_core) <- names(id_core_list)

annot_disp <- map(names(id_disp_list), function(x){
  
  left_join(id_disp_list[[x]], annot_list2[[x]], by = "gene_name")
  
})
names(annot_disp) <- names(id_disp_list)

annot_priv <- map(names(id_priv_list), function(x){
  
  left_join(id_priv_list[[x]], annot_list2[[x]], by = "gene_name")
  
})
names(annot_priv) <- names(id_priv_list)

# export
map(names(annot_core), function(x){write_csv(annot_core[[x]], here("assembly", "pangenome", "intersect", paste0(x, "_reclass_core.annotation.csv")))})
map(names(annot_disp), function(x){write_csv(annot_disp[[x]], here("assembly", "pangenome", "intersect", paste0(x, "_reclass_disp.annotation.csv")))})
map(names(annot_priv), function(x){write_csv(annot_priv[[x]], here("assembly", "pangenome", "intersect", paste0(x, "_reclass_priv.annotation.csv")))})


# # plot ------------------------------------------------------------
# load(file = here("assembly", "pangenome", "intersect", "all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_list.rda"))
# load(file = here("assembly", "pangenome", "intersect", "all.on.all.wfmash_s10000p85n1.seqwish_k49.smooth.join.intersect_stat_tb.rda"))
# 
# intersect_tb <- bind_rows(intersect_list3)
# 
# intersect_tb2 <- intersect_tb %>%
#   mutate(name = gsub(".hap.*$", "", gene_name))
# 
# intersect_stat_tb2 <- intersect_stat_tb %>%
#   mutate(gene_type = factor(gene_type, levels = c("core", "disp", "priv"), labels = c("core genes", "dispensable genes", "private genes")),
#          seq_type = factor(seq_type, levels = c("core", "disp", "priv"), labels = c("core sequence", "dispensable sequence", "private sequence")))
# 
# p <- ggplot(data = intersect_stat_tb2, mapping = aes(x = seq_type, y = mean, fill = seq_type)) +
#   geom_boxplot(position = position_dodge(0.75), width = 0.7) +
#   labs(x = "", y = "Contribution of sequence type (% of bp)", title = "Sequence pangenome classes per genic pangenome classes", subtitle = "Genic pangenome classes are represented in column facets.\nPer facet, the contribution of each class of the sequence pangenome is represented as a percentage of the gene sequence.\nWe can see that each class of gene is, in average, well represented by its class of sequence,\ni.e. a private gene (on the right side column facet) is mostly composed of bases classified as private in the sequence pangenome.") +
#   facet_grid(cols = vars(gene_type)) +
#   scale_fill_manual(name = "", values=c("#90E711", "#F4D213", "#A657F5"))
# 
# # ggplot(data = intersect_tb2, mapping = aes(x = seq_type, y = perc, fill = seq_type)) +
# #   geom_boxplot(position = position_dodge(0.75), width = 0.7) +
# #   labs(x = "Genic pangenome class", y = "Sequence pangenome class (bp %)") +
# #   facet_grid(cols = vars(gene_type)) +
# #   scale_fill_manual(name = "", values=c("#90E711", "#F4D213", "#A657F5")) +
# #   stat_summary(fun=mean, geom="point", size=4, color = "red", show.legend = F) 
# 
# p <- p + theme(plot.title = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
#                plot.subtitle = element_text(color="grey30", size=6, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
#                axis.title.x = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
#                axis.title.y = element_text(color="black", size=10, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
#                axis.text.x = element_blank(),
#                axis.text.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
#                axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
#                axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
#                axis.ticks.x = element_blank(),
#                axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
#                plot.background = element_blank(),
#                panel.background = element_blank(),
#                panel.grid.major.y = element_blank(),
#                panel.grid.major.x = element_blank(),
#                panel.grid.minor.y = element_blank(),
#                panel.grid.minor.x = element_blank(),
#                plot.margin=unit(c(1,1,1,1),"mm"),
#                legend.position="bottom")
# 
# 
# p
# 
# ggsave(p, filename = here("plots", "assembly", "pangenome", "smoothxg", "genicVSsequence_pangenome.pdf"), width = 6, height = 4)
# 
# 
# # summarized plot ---------------------------------------------------------
# intersect_stat_tb3 <- intersect_stat_tb2 %>%
#   group_by(gene_type, seq_type) %>%
#   summarize(avg_perc = mean(mean)) %>%
#   ungroup()
# 
# intersect_stat_tb4 <- intersect_stat_tb3 %>%
#   group_by(gene_type) %>%
#   mutate(max_perc = max(avg_perc), min_perc = min(avg_perc)) %>%
#   ungroup() %>%
#   rowwise() %>%
#   mutate(col_val = if(avg_perc == max_perc){"A"} else if(avg_perc == min_perc){"C"} else {"B"}) %>%
#   ungroup() %>%
#   mutate(col_val = factor(col_val, levels = LETTERS[c(3,2,1)]),
#          gene_type = factor(gene_type, levels = c("private genes", "dispensable genes", "core genes")))
# 
# p <- ggplot() +
#   geom_bar(data = intersect_stat_tb4, mapping = aes(x = gene_type, y = avg_perc, group = col_val, fill = seq_type), stat = "identity", position = "stack", show.legend = F, width = 0.75) +
#   coord_flip() +
#   labs(x = "", y = "Pangenome sequence class (%)") +
#   scale_fill_manual(values = c("#90E711", "#F4D213", "#A657F5"))
# 
# p <- p + theme(plot.title = element_blank(),
#                axis.title.x = element_text(color="black", size=7, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
#                axis.title.y = element_blank(),
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
#                legend.position="none")
# 
# p
# 
# 
# ggsave(p, filename = here("plots", "assembly", "pangenome", "smoothxg", "genicVSsequence_pangenome_summarized.pdf"), width = 7/3, height = 1.5)
# 
# 
# 
# 
# 
