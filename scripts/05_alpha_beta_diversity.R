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
library(stringr)
library(vegan)
library(phyloseq)
library(here)
library(microbial)
library(rstudioapi)
library(microbiome)

# 01. Alpha diversity ----
group_info <- read_delim(
  here(
    "results",
    "tables",
    "02_group_info_metagenomics_detailed.txt"
  ),
  delim = "\t"
)

## Normalize data
jobRunScript(
  path = here(
    "scripts",
    "04_normalize_and_diff_abund.R"
  ),
  exportEnv = "R_GlobalEnv"
)

## Merge data
bracken_clean_tax <- psmelt(bracken_data) %>%
  select(OTU, contains("Rank")) %>%
  group_by(OTU) %>%
  summarise_all(list(func_paste)) %>%
  column_to_rownames("OTU")

OTU <- otu_table(bracken_norm$feature_table, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(bracken_clean_tax)) 

bracken_clean <- phyloseq(OTU, TAX)

save(bracken_clean,
     file = here(
       "results",
       "RData",
       "bracken_norm_clean.rdata"
    )
)

## Calculate bacterial alpha diversity
biom_bacteria <- subset_taxa(bracken_clean, Rank1 == "k__Bacteria")

alpha_bacteria <- estimate_richness(biom_bacteria) %>%
  rownames_to_column("ref") %>%
  mutate(ref = sub("X.+(2020.01.....)", "\\1", ref),
         ref = gsub("\\.", "-", ref)) %>%
  left_join(group_info, by = c("ref" = "saksnr")) %>%
  select(-c(se.chao1, se.ACE, InvSimpson, Fisher, ACE)) %>%
  pivot_longer(cols = c("Observed",
                      "Chao1",
                      "Shannon",
                      "Simpson"),
             names_to = "measure",
             values_to = "value")

## Run stats
### Factorial ANOVA first to see if there are
### any significant differences overall
anova(lm(value ~ measure * origin, data = alpha_bacteria))

#### Doesn't seem to be any significant differences
#### in the data from only bacteria.

stats_data <- alpha_bacteria %>%
  pivot_wider(names_from = "measure",
              values_from = "value")

wilcox.test(Observed ~ origin,
            data = stats_data,
            alternative = "two.sided")

wilcox.test(Chao1 ~ origin,
            data = stats_data,
            alternative = "two.sided")

wilcox.test(Shannon ~ origin,
            data = stats_data,
            alternative = "two.sided")


## Plot figures
palette <- c("Broiler" = colorspace::lighten("#006c89", 
                                             amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", 
                                            amount = 0.5))

p_alpha <- alpha_bacteria %>%
  filter(measure %in% c("Shannon","Observed")) %>%
  ggplot(aes(origin, value)) +
  geom_violin(aes(fill = origin),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  scale_fill_manual(values = palette) +
  labs(y = "Distance measure value",
       fill = NULL,
       title = "A") +
  facet_wrap(~measure, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10))

save(p_alpha,
     file = here(
       "results",
       "RData",
       "alpha_div_plot.rdata"
     )
)

ggsave(
  here(
    "results",
    "figures",
    "05_alpha_diversity.png"
  ),
  p_alpha,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 15
)

# 02. Beta diversity ----
## Ordinate NMDS
ord_bact <- ordinate(
  physeq = biom_bacteria,
  method = "NMDS",
  distance = "bray"
)

## Create plot data
plot_data <- as.data.frame(ord_bact$points) %>%
  rownames_to_column("ref") %>%
  filter(!grepl("zymo", ref)) %>%
  mutate(ref = sub(".+-(2020.+)", "\\1", ref)) %>%
  left_join(group_info, by = c("ref" = "saksnr"))

## Run stats
otu_table <- psmelt(biom_bacteria) %>%
  mutate(Sample = sub(".+(2020.+)", "\\1", Sample)) %>%
  select(-contains("Rank")) %>%
  pivot_wider(names_from = "OTU",
              values_from = "Abundance") %>%
  mutate(Sample = sub(".+-(2020.+)", "\\1", Sample)) %>%
  column_to_rownames("Sample")

stats_metadata <- group_info %>%
  select(-cfu_g_total) %>%
  filter(saksnr %in% rownames(otu_table)) %>%
  column_to_rownames("saksnr")

bray_dist <- vegdist(otu_table, method = "bray")

adonis2(
  bray_dist ~ origin,
  data = stats_metadata,
  permutations = 9999,
  method = "bray"
  )

## Plot figure
p_beta <- ggplot(plot_data, aes(MDS1, MDS2, color = origin)) +
  stat_ellipse(linetype = "dotted",
               linewidth = 0.5) +
  geom_point(aes(shape = result),
             size = 1.5) +
  labs(shape = NULL,
       color = NULL,
       title = "D") +
  annotate(geom = "text",
           label = paste0("Stress: ", round(ord_bact$stress, 4)),
           x = -0.085,
           y = -0.09,
           size = 2) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))

save(p_beta,
     file = here(
       "results",
       "RData",
       "beta_div_plot.rdata"
    )
)

ggsave(
  here(
    "results",
    "figures",
    "05_beta_diversity.png"
  ),
  p_beta,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 16
)



