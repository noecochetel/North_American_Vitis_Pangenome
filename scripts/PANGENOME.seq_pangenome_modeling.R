library(tidyverse)
library(here)


# prepare variable --------------------------------------------------------
seq_pan_stats_file <- here("assembly", "pangenome", "path", "all_stats.txt")


# load data ---------------------------------------------------------------
seq_pan_stats_tb <- read_tsv(seq_pan_stats_file, col_names = c("genome_nb", "comb", "type", "n", "len"))


# prepare data ------------------------------------------------------------
seq_pan_stats_tb_priv1 <- seq_pan_stats_tb %>%
  filter(genome_nb == 1, type == "priv") %>%
  mutate(n = 0, len = 0)

seq_pan_stats_tb_notpriv1 <- seq_pan_stats_tb %>%
  filter((genome_nb == 1 & type != "priv"))

seq_pan_stats_tb_not1 <- seq_pan_stats_tb %>%
  filter((genome_nb != 1))

seq_pan_stats_tb2 <- bind_rows(seq_pan_stats_tb_priv1, seq_pan_stats_tb_notpriv1, seq_pan_stats_tb_not1) %>%
  replace_na(replace = list(len = 0)) %>%
  mutate(n = n / 1000000, len = len / 1000000000)

n_stats_tb <- seq_pan_stats_tb2 %>%
  dplyr::select(-len) %>%
  pivot_wider(names_from = type, values_from = n) %>%
  mutate(pan = priv + disp + core) %>%
  pivot_longer(cols = priv:pan) %>%
  mutate(name = factor(name, levels = c("pan", "core", "disp", "priv")))

n_stats_tb2 <- n_stats_tb %>%
  group_by(genome_nb, name) %>%
  summarize(n_mean = mean(value))

len_stats_tb <- seq_pan_stats_tb2 %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = type, values_from = len) %>%
  mutate(pan = priv + disp + core) %>%
  pivot_longer(cols = priv:pan) %>%
  mutate(name = factor(name, levels = c("pan", "core", "disp", "priv")))

len_stats_tb2 <- len_stats_tb %>%
  group_by(genome_nb, name) %>%
  summarize(len_mean = mean(value))




# plot --------------------------------------------------------------------
n_p <- ggplot() +
  geom_point(data = n_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), alpha = 0.5, size = 1) +
  geom_smooth(data = n_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), method = "loess", show.legend = F, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(9)) +
  # scale_y_continuous(breaks = seq(0, 60, 10)) +
  labs(y = "Genome segments composition (Millions)", x = "Number of genomes") +
  scale_color_manual(name = "", labels = c("Pan-genome", "Core genome", "Dispensable genome", "Private genome"), values=c("#F01579", "#90E711", "#F4D213", "#A657F5")) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

n_p <- n_p + theme(plot.title = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
               axis.title.x = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
               axis.title.y = element_text(color="black", size=10, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
               axis.text.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
               axis.text.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
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
               legend.position = c(0, 1),
               legend.justification = c("left", "top"))

n_p

ggsave(n_p, filename = here("plots", "assembly", "pangenome", "smoothxg", "sequence_pangenome_modeling_n.pdf"), width = 3.5, height = 3)

save(n_p, file = here("plots", "assembly", "pangenome", "smoothxg", "sequence_pangenome_modeling_n.rda"))

len_p <- ggplot() +
  geom_point(data = len_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), alpha = 0.5, size = 1) +
  geom_smooth(data = len_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), method = "loess", show.legend = F, linewidth = 0.5) +
  scale_x_continuous(breaks = seq(9)) +
  scale_y_continuous(breaks = seq(0,2,0.5), limits = c(0,2)) +
  # scale_y_continuous(breaks = seq(0, 60, 10)) +
  labs(y = "Genome length (Gb)", x = "Number of genomes") +
  scale_color_manual(name = "", labels = c("Pan-genome", "Core genome", "Dispensable genome", "Private genome"), values=c("#F01579", "#90E711", "#F4D213", "#A657F5")) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

len_p <- len_p + theme(plot.title = element_blank(),
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
                   legend.position = c(0, 1),
                   legend.justification = c("left", "top"))

len_p

ggsave(len_p, filename = here("plots", "assembly", "pangenome", "smoothxg", "sequence_pangenome_modeling_len.pdf"), width = 3.5, height = 2.5)

save(len_p, file = here("plots", "assembly", "pangenome", "smoothxg", "sequence_pangenome_modeling_len.rda"))

library(patchwork)

p <- n_p | len_p
p

ggsave(p, filename = here("plots", "assembly", "pangenome", "smoothxg", "sequence_pangenome_modeling.pdf"), width = 7, height = 3)



# # post-filtering ----------------------------------------------------------
# # prepare variable --------------------------------------------------------
# seq_pan_stats_file <- here("assembly", "pangenome", "path", "all_stats.gt10k.txt")
# 
# 
# # load data ---------------------------------------------------------------
# seq_pan_stats_tb <- read_tsv(seq_pan_stats_file, col_names = c("genome_nb", "comb", "type", "n", "len"))
# 
# 
# # prepare data ------------------------------------------------------------
# seq_pan_stats_tb_priv1 <- seq_pan_stats_tb %>%
#   filter(genome_nb == 1, type == "priv") %>%
#   mutate(n = 0, len = 0)
# 
# seq_pan_stats_tb_notpriv1 <- seq_pan_stats_tb %>%
#   filter((genome_nb == 1 & type != "priv"))
# 
# seq_pan_stats_tb_not1 <- seq_pan_stats_tb %>%
#   filter((genome_nb != 1))
# 
# seq_pan_stats_tb2 <- bind_rows(seq_pan_stats_tb_priv1, seq_pan_stats_tb_notpriv1, seq_pan_stats_tb_not1) %>%
#   replace_na(replace = list(len = 0)) %>%
#   mutate(n = n / 1000000, len = len / 1000000000)
# 
# n_stats_tb <- seq_pan_stats_tb2 %>%
#   select(-len) %>%
#   pivot_wider(names_from = type, values_from = n) %>%
#   mutate(pan = priv + disp + core) %>%
#   pivot_longer(cols = priv:pan) %>%
#   mutate(name = factor(name, levels = c("pan", "core", "disp", "priv")))
# 
# n_stats_tb2 <- n_stats_tb %>%
#   group_by(genome_nb, name) %>%
#   summarize(n_mean = mean(value))
# 
# len_stats_tb <- seq_pan_stats_tb2 %>%
#   select(-n) %>%
#   pivot_wider(names_from = type, values_from = len) %>%
#   mutate(pan = priv + disp + core) %>%
#   pivot_longer(cols = priv:pan) %>%
#   mutate(name = factor(name, levels = c("pan", "core", "disp", "priv")))
# 
# len_stats_tb2 <- len_stats_tb %>%
#   group_by(genome_nb, name) %>%
#   summarize(len_mean = mean(value))
# 
# 
# 
# 
# # plot --------------------------------------------------------------------
# n_p <- ggplot() +
#   geom_point(data = n_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), alpha = 0.5, size = 1) +
#   geom_smooth(data = n_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), method = "loess", show.legend = F, size = 0.5) +
#   scale_x_continuous(breaks = seq(9)) +
#   # scale_y_continuous(breaks = seq(0, 60, 10)) +
#   labs(y = "Genome segments composition (Millions)", x = "Number of genomes") +
#   scale_color_manual(name = "", labels = c("Pan-genome", "Core genome", "Dispensable genome", "Private genome"), values=c("#F01579", "#90E711", "#F4D213", "#A657F5")) +
#   guides(color = guide_legend(override.aes = list(alpha = 1)))
# 
# n_p <- n_p + theme(plot.title = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
#                    axis.title.x = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
#                    axis.title.y = element_text(color="black", size=10, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
#                    axis.text.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
#                    axis.text.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
#                    axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
#                    axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
#                    axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
#                    axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
#                    plot.background = element_blank(),
#                    panel.background = element_blank(),
#                    panel.grid.major.y = element_blank(),
#                    panel.grid.major.x = element_blank(),
#                    panel.grid.minor.y = element_blank(),
#                    panel.grid.minor.x = element_blank(),
#                    plot.margin=unit(c(1,1,1,1),"mm"),
#                    legend.title=element_text(size=8), 
#                    legend.text=element_text(size=6),
#                    legend.background = element_blank(),
#                    legend.key.size = unit(2, "mm"),
#                    legend.key=element_blank(),
#                    legend.position = c(0, 1),
#                    legend.justification = c("left", "top"))
# 
# n_p
# 
# len_p <- ggplot() +
#   geom_point(data = len_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), alpha = 0.5, size = 1) +
#   geom_smooth(data = len_stats_tb, mapping = aes(x = genome_nb, y = value, color = name), method = "loess", show.legend = F, size = 0.5) +
#   scale_x_continuous(breaks = seq(9)) +
#   # scale_y_continuous(breaks = seq(0, 60, 10)) +
#   labs(y = "Genome length (Gb)", x = "Number of genomes") +
#   scale_color_manual(name = "", labels = c("Pan-genome", "Core genome", "Dispensable genome", "Private genome"), values=c("#F01579", "#90E711", "#F4D213", "#A657F5")) +
#   guides(color = guide_legend(override.aes = list(alpha = 1)))
# 
# len_p <- len_p + theme(plot.title = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=1,l=0,"mm")),
#                        axis.title.x = element_text(color="black", size=10, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm")),
#                        axis.title.y = element_text(color="black", size=10, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm")),
#                        axis.text.x = element_text(color="black", size=8, face="plain", margin=margin(t=1,r=0,b=0,l=0,"mm"), hjust = 0.5),
#                        axis.text.y = element_text(color="black", size=8, face="plain", margin=margin(t=0,r=1,b=0,l=0,"mm"), vjust = 0.5),
#                        axis.line.x = element_line(color="black", size=0.5, linetype="solid"),
#                        axis.line.y = element_line(color="black", size=0.5, linetype="solid"),
#                        axis.ticks.x = element_line(color="black", size=0.5, linetype="solid"),
#                        axis.ticks.y = element_line(color="black", size=0.5, linetype="solid"),
#                        plot.background = element_blank(),
#                        panel.background = element_blank(),
#                        panel.grid.major.y = element_blank(),
#                        panel.grid.major.x = element_blank(),
#                        panel.grid.minor.y = element_blank(),
#                        panel.grid.minor.x = element_blank(),
#                        plot.margin=unit(c(1,1,1,1),"mm"),
#                        legend.title=element_text(size=8), 
#                        legend.text=element_text(size=6),
#                        legend.background = element_blank(),
#                        legend.key.size = unit(2, "mm"),
#                        legend.key=element_blank(),
#                        legend.position = c(0, 1),
#                        legend.justification = c("left", "top"))
# 
# len_p
# 
# library(patchwork)
# 
# p <- n_p | len_p
# p
# 
# 
# ggsave(p, filename = here("plots", "assembly", "pangenome", "smoothxg", "sequence_pangenome_modeling.gt10k.pdf"), width = 7, height = 3)
# 
# 
