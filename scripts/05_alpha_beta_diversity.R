# ABSTRACT
# Script for calculating alpha- and beta
# diversity from the kraken data

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(ggsignif)
library(impoRt)
library(readxl)
library(tibble)
library(stringr)
library(vegan)
library(phyloseq)
library(here)
library(rstudioapi)
library(patchwork)

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
bracken_clean_tax <- psmelt(bracken_norm) %>%
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

load(
  here(
    "results",
    "RData",
    "bracken_norm_clean.rdata"
  )
)

## Calculate bacterial alpha diversity
### Define dataset for bacteria only
biom_bacteria <- subset_taxa(bracken_clean, Rank1 == "k__Bacteria")

### Agglomerate for each taxrank
ranks <- c("Rank2",
           "Rank3",
           "Rank4",
           "Rank5",
           "Rank6",
           "Rank7")

tax_names <- c("Phylum",
               "Class",
               "Order",
               "Family",
               "Genus",
               "Species")

biom_subset_list <- lapply(ranks, function(x) {
  tax_glom(biom_bacteria, taxrank = x)
})

names(biom_subset_list) <- tax_names

### Calculate alpha diversity
alpha_bacteria <- lapply(biom_subset_list, function(x) {
  estimate_richness(x) %>%
    rownames_to_column("ref") %>%
    mutate(ref = sub("X.+(2020.01.....)", "\\1", ref),
           ref = gsub("\\.", "-", ref)) %>%
    select(ref, Shannon) %>%
    left_join(group_info, by = c("ref" = "saksnr"))
})

alpha_plot_data <- bind_rows(alpha_bacteria, .id = "taxrank") %>%
  mutate(taxrank = factor(taxrank, levels = tax_names))

### Run stats
aov_stats_origin <- lapply(alpha_bacteria, function(x) {
  anova(lm(Shannon ~ origin, data = x))
}) %>%
  bind_rows(.id = "taxrank") %>%
  mutate(taxrank = factor(taxrank, levels = c("Species",
                                              "Genus",
                                              "Family",
                                              "Order",
                                              "Class",
                                              "Phylum"))) %>%
  rownames_to_column("type") %>%
  filter(grepl("origin",type),
         `Pr(>F)` < 0.05)

aov_stats_group <- lapply(alpha_bacteria, function(x) {
  anova(lm(Shannon ~ group, data = x))
}) %>%
  bind_rows(.id = "taxrank") %>%
  mutate(taxrank = factor(taxrank, levels = c("Species",
                                              "Genus",
                                              "Family",
                                              "Order",
                                              "Class",
                                              "Phylum"))) %>%
  rownames_to_column("type") %>%
  filter(grepl("group",type),
       `Pr(>F)` < 0.05)

#### AOV indicate a difference within the
#### Order and Species taxranks for the groups,
#### and Class and Order taxranks for the host species.
#### Will therefore run pairwise wilcox.tests
#### for each combination within each taxrank,
#### with p-value correction

stats_data <- alpha_plot_data %>%
  split(f = .$taxrank)

tests_origin <- lapply(stats_data, function(x) {
  test <- wilcox.test(
    x$Shannon ~ x$origin,
    alternative = "two.sided",
    correct = FALSE
  )$p.value
}) %>%
  bind_rows() %>%
  pivot_longer(
    cols = everything(),
    names_to = "taxid",
    values_to = "p_val"
  ) %>%
  mutate(p_adj = p.adjust(p_val, 
                          method = "bonferroni"))

#### Significant difference observed on the 
#### order taxon. Broiler samples have less 
#### diversity in the order taxon compared
#### to turkeys

### Generate plots
palette <- c("Broiler" = colorspace::lighten("#006c89", amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", amount = 0.5))

group_palette <- c(
  "broiler_negative" = colorspace::lighten("#006c89", amount = 0.9),
  "broiler_positive" = colorspace::lighten("#006c89", amount = 0.5),
  "turkey_negative" = colorspace::lighten("#3d6721", amount = 0.9),
  "turkey_positive" = colorspace::lighten("#3d6721", amount = 0.5)
)

p_alpha <- alpha_plot_data %>%
  filter(taxrank == "Order") %>%
  ggplot(aes(group, Shannon)) +
  geom_violin(aes(fill = group),
              trim = FALSE,
              linewidth = 0.3) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(width = 0.1,
               outlier.size = 2,
               size = 0.3) +
  labs(y = "Shannon Index",
       fill = NULL,
       title = "A") +
  scale_fill_manual(values = group_palette,
                    labels = c("Broiler negative",
                               "Broiler positive",
                               "Turkey negative",
                               "Turkey positive")) +
  scale_x_discrete(labels = c("Negative","Positive","Negative","Positive")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title.y = element_text(size = 8))


# 02. Beta diversity ----
## Ordinate NMDS
ord_list <- lapply(biom_subset_list, function(x) {
  ordinate(
    physeq = x,
    method = "NMDS",
    distance = "bray"
  )
})


## Clean data
stress_info <- lapply(tax_names, function(x) {
  data.frame("stress" = ord_list[[x]]$stress,
             "taxrank" = x)
}) %>%
  bind_rows() %>%
  mutate(taxrank = factor(taxrank, levels = tax_names),
         stress = round(stress, 4))

clean_ord_list <- lapply(ord_list, function(x) {
  as.data.frame(x$points) %>%
    rownames_to_column("ref") %>%
    filter(!grepl("zymo", ref)) %>%
    mutate(ref = sub(".+-(2020.+)", "\\1", ref)) %>%
    left_join(group_info, by = c("ref" = "saksnr"))
}) %>%
  bind_rows(.id = "taxrank") %>%
  mutate(taxrank = factor(taxrank, levels = c("Species",
                                              "Genus",
                                              "Family",
                                              "Order",
                                              "Class",
                                              "Phylum")))

## Run adonis test on each taxrank
tax_rank_list <- lapply(ranks, function(x) {
  tax_glom(biom_bacteria, taxrank = x) %>%
    psmelt() %>%
    mutate(Sample = sub(".+(2020.+)", "\\1", Sample)) %>%
    select(-contains("Rank")) %>%
    pivot_wider(names_from = "OTU",
                values_from = "Abundance") %>%
    mutate(Sample = sub(".+-(2020.+)", "\\1", Sample)) %>%
    column_to_rownames("Sample")
})

names(tax_rank_list) <- tax_names

bray_dist_list <- lapply(tax_rank_list, function(x) {
  vegdist(x, method = "bray")
})

stats_metadata <- group_info %>%
  select(-cfu_g_total) %>%
  filter(saksnr %in% rownames(tax_rank_list$Phylum)) %>%
  column_to_rownames("saksnr")

adonis_bray <- lapply(bray_dist_list, function(x) {
  as.data.frame(
    adonis2(x ~ origin + group,
            data = stats_metadata,
            permutations = 9999,
            method = "bray")
  )
}) %>%
  bind_rows(.id = "Tax_rank") %>%
  rownames_to_column("group")

## Plot beta diversity
p_beta <- clean_ord_list %>%
  filter(taxrank == "Order") %>%
  ggplot(aes(MDS1, MDS2)) +
  stat_ellipse(aes(color = group),
               linetype = "dotted",
               linewidth = 1) +
  geom_point(aes(color = group),
             size = 2) +
  scale_color_manual(values = group_palette) +
  labs(shape = NULL,
       title = "B") +
  guides(color = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6))


save(
  p_alpha,
  file = here(
    "results",
    "RData",
    "alpha_plot.rdata"
  )
)

save(
  p_beta,
  file = here(
    "results",
    "RData",
    "beta_plot.rdata"
  )
)
