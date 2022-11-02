# ABSTRACT
# Script for calculating alpha- and beta
# diversity from the kraken data

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(impoRt)
library(readxl)
library(tibble)
library(vegan)
library(phyloseq)
library(here)

# 01. Alpha diversity ----
group_info <- read_delim(
  here(
    "results",
    "tables",
    "02_group_info_metagenomics_detailed.txt"
  ),
  delim = "\t"
)

biom_file <- import_biom(
  here(
    "data",
    "bracken",
    "biom",
    "bracken_species.biom"
    )
  )

## Calculate overall alpha diversity
alpha_div <- estimate_richness(biom_file) %>%
  rownames_to_column("ref") %>%
  mutate(ref = sub("X.+(2020.01.....)_S.+", "\\1", ref),
         ref = gsub("\\.", "-", ref)) %>%
  left_join(group_info, by = c("ref" = "saksnr")) %>%
  select(-c(se.chao1, ACE, se.ACE, InvSimpson)) %>%
  pivot_longer(cols = c("Observed",
                        "Chao1",
                        "Shannon",
                        "Simpson",
                        "Fisher"),
               names_to = "measure",
               values_to = "value")

## Calculate bacterial alpha diversity
biom_bacteria <- subset_taxa(biom_file, Rank1 == "k__Bacteria")

alpha_bacteria <- estimate_richness(biom_bacteria) %>%
  rownames_to_column("ref") %>%
  mutate(ref = sub("X.+(2020.01.....)_S.+", "\\1", ref),
         ref = gsub("\\.", "-", ref)) %>%
  left_join(group_info, by = c("ref" = "saksnr")) %>%
  select(-c(se.chao1, ACE, se.ACE, InvSimpson)) %>%
  pivot_longer(cols = c("Observed",
                      "Chao1",
                      "Shannon",
                      "Simpson",
                      "Fisher"),
             names_to = "measure",
             values_to = "value")

## Plot figures
palette <- c("Broiler" = colorspace::lighten("#006c89", 
                                             amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", 
                                            amount = 0.5))

p_alpha <- ggplot(alpha_div, aes(origin, value)) +
  geom_violin(aes(fill = origin),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  scale_fill_manual(values = palette) +
  labs(y = "Distance measure value",
       fill = NULL) +
  facet_wrap(~measure, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.87, 0.3))

p_alpha_bacteria <- ggplot(alpha_bacteria, aes(origin, value)) +
  geom_violin(aes(fill = origin),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  scale_fill_manual(values = palette) +
  labs(y = "Distance measure value",
       fill = NULL) +
  facet_wrap(~measure, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.87, 0.3))

ggsave(
  here(
    "results",
    "figures",
    "04_alpha_diversity.png"
  ),
  p_alpha,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 15
)

ggsave(
  here(
    "results",
    "figures",
    "04_alpha_diversity_bacteria.png"
  ),
  p_alpha_bacteria,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 15
)

# 02. Beta diversity ----
## Transform to percentages
biom_percent_bacteria <- transform_sample_counts(
  biom_bacteria, function(x) x*100 / sum(x)
  )

beta_ordination_bray <- ordinate(
  physeq = biom_percent_bacteria,
  method = "NMDS",
  distance = "bray"
)

beta_plot_data_bray <- as.data.frame(beta_ordination_bray$points) %>%
  rownames_to_column("ref") %>%
  mutate(ref = sub(".+-(2020.+)_S.+", "\\1", ref)) %>%
  left_join(group_info, by = c("ref" = "saksnr"))

## Run stats
otu_table <- psmelt(biom_percent) %>%
  select(-contains("Rank")) %>%
  pivot_wider(names_from = "OTU",
              values_from = "Abundance") %>%
  mutate(Sample = sub(".+-(2020.+)_S.+", "\\1", Sample)) %>%
  column_to_rownames("Sample")

stats_metadata <- group_info %>%
  select(-cfu_g_total) %>%
  filter(saksnr %in% rownames(otu_table)) %>%
  column_to_rownames("saksnr")

bray_dist <- vegdist(otu_table, method = "bray")

adonis2(bray_dist ~ origin, data = stats_metadata, permutations = 9999, method = "bray")

## Plot figure

label <- "<b><i>F</i></b>(1,42) = 2.6907, <b><i>p</i></b> = 0.0156"

p_beta <- ggplot(beta_plot_data_bray, aes(MDS1, MDS2, color = origin)) +
  stat_ellipse(linetype = "dotted",
               size = 1) +
  geom_point(aes(shape = result),
             size = 2.5) +
  annotation_custom(richtext_grob(label),
                    xmin = -0.13,
                    xmax = -0.01,
                    ymin = -0.21,
                    ymax = -0.15) +
  labs(shape = NULL,
       color = NULL,
       title = "Bray-Curtis Beta Diversity") +
  scale_color_manual(values = palette) +
  scale_y_continuous(limits = c(-0.17, 0.2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9,0.83))
  


ggsave(
  here(
    "results",
    "figures",
    "04_beta_diversity_bacteria.png"
  ),
  p_beta,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 15
)



