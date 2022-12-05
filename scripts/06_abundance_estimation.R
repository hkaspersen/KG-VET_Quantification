# ABSTRACT
# Script for summarising and running statistics
# on the metagenomics quantification data

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(ggsignif)
library(readxl)
library(here)
library(tibble)
library(impoRt)
library(vegan)
library(janitor)
library(phyloseq)
library(patchwork)
library(gridExtra)
library(ANCOMBC)

# 01. Kraken/Bracken abundance ----
## Abundance estimation Bracken ----
group_info <- read_delim(
  here(
    "results",
    "tables",
    "02_group_info_metagenomics_detailed.txt"
  ),
  delim = "\t"
)

load(
  file = here(
    "results",
    "RData",
    "bracken_norm_clean.rdata"
  )
)

## Filter to only bacteria
bracken_bact <- subset_taxa(bracken_clean, Rank1 == "k__Bacteria")

## Convert to relative abundance
bracken_rel <- transform_sample_counts(
  bracken_bact, function(x) x*100 / sum(x))

## Calculate relative abundance of orders
order_data <- tax_glom(bracken_rel,
                        taxrank = "Rank4") %>%
  psmelt() %>%
  mutate(Sample = sub(".+(2020.+)", "\\1", Sample)) %>%
  left_join(group_info, by = c("Sample" = "saksnr")) %>%
  mutate(Rank4 = sub("o__", "", Rank4),
         plot_group = ifelse(Abundance >= 5, Rank4, "Other")) %>%
  mutate(plot_group = factor(plot_group,
                             levels = c(
                               "Bacteroidales",
                               "Bifidobacteriales",
                               "Coriobacteriales",
                               "Enterobacterales",
                               "Eubacteriales",
                               "Lactobacillales",
                               "Selenomonadales",
                               "Other"
                             )))

## Generate abundance figure ----
palette <- c("Broiler" = colorspace::lighten("#006c89", amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", amount = 0.5))

group_names <- c(
  "broiler_positive" = "Broiler positive",
  "broiler_negative" = "Broiler negative",
  "turkey_positive" = "Turkey positive",
  "turkey_negative" = "Turkey negative"
)

p_order <- ggplot(order_data, aes(Sample, Abundance, fill = plot_group)) +
  geom_col(position = position_fill()) +
  labs(x = "Samples",
       y = "Relative Abundance (%)",
       fill = NULL,
       title = "A") +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 6)) +
  facet_wrap(~group, scales = "free_x", 
             labeller = as_labeller(group_names))


p_order_nolegend <- p_order +
  theme(legend.position = "none")


## Abundance of Klebsiella ----
### Klebsiella genus ----
genera_data <- tax_glom(bracken_bact,
                        taxrank = "Rank6")

klebsiella_results <- psmelt(genera_data) %>%
  filter(Rank6 == "g__Klebsiella") %>%
  mutate(Sample = sub(".+-(2020-.+)", "\\1", Sample)) %>%
  select(Sample, Abundance, Rank6) %>%
  left_join(group_info, by = c("Sample" = "saksnr")) %>%
  mutate(cfu_group = ifelse(is.na(cfu_group), "Negative", cfu_group),
         cfu_group = factor(cfu_group, levels = c("Negative",
                                                  "Low",
                                                  "Medium",
                                                  "High")),
         Abundance = round(Abundance, 4))

klebsiella_summary <- klebsiella_results %>%
  group_by(origin) %>%
  summarise(mean = round(mean(Abundance), 3),
            median = round(median(Abundance), 3),
            sd = round(sd(Abundance), 3),
            range = paste0(min(Abundance), " - ", max(Abundance))) %>%
  mutate(method = "Bracken")

write_delim(
  klebsiella_summary,
  here(
    "results",
    "tables",
    "06_klebsiella_summary_bracken.txt"
  ),
  delim = "\t"
)

test_data <- klebsiella_results %>%
  filter(result == "Positive")

summary(lm(log10(Abundance) ~ log10(cfu_g_total), data = test_data))

p_stats <- ggplot(test_data, aes(log10(cfu_g_total), 
                            log10(Abundance),
                      color = origin)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = palette) +
  labs(x = "Log10 CFU/g",
       y = "Log10 Absolute abundance",
       title = "B") +
  theme_bw() +
  theme(panel.grid = element_blank())


group_palette <- c(
  "broiler_negative" = colorspace::lighten("#006c89", amount = 0.9),
  "broiler_positive" = colorspace::lighten("#006c89", amount = 0.5),
  "turkey_negative" = colorspace::lighten("#3d6721", amount = 0.9),
  "turkey_positive" = colorspace::lighten("#3d6721", amount = 0.5)
)

p_klebs <- ggplot(klebsiella_results, aes(group, Abundance)) +
  geom_violin(aes(fill = group),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  scale_fill_manual(values = group_palette,
                    labels = c("Broiler negative",
                               "Broiler positive",
                               "Turkey negative",
                               "Turkey positive")) +
  labs(y = "Absolute abundance of *Klebsiella* spp.",
       title = "C") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = ggtext::element_markdown())


p_all_nolegend <- p_order_nolegend | (p_stats / p_klebs) +
  plot_layout(guides = "collect")

p_all <- p_order | (p_stats / p_klebs) +
  plot_layout(guides = "collect")

ggsave(
  here(
    "results",
    "figures",
    "06_relative_abundance.png"
  ),
  p_all_nolegend,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 22
)

ggsave(
  here(
    "results",
    "figures",
    "06_relative_abundance_legend.png"
  ),
  p_all,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 22
)


### Klebsiella species ----
klebsiella_species <- psmelt(bracken_rel) %>%
  filter(Rank6 == "g__Klebsiella") %>%
  mutate(Sample = sub(".+-(2020-.+)", "\\1", Sample),
         Rank7 = sub("s__", "", Rank7)) %>%
  left_join(group_info[,c("saksnr","origin","result","group","cfu_group")],
            by = c("Sample" = "saksnr"))

kpsc_species <- psmelt(bracken_rel) %>%
  filter(Rank6 == "g__Klebsiella") %>%
  mutate(Sample = sub(".+-(2020-.+)", "\\1", Sample),
         Rank7 = sub("s__", "", Rank7)) %>%
  left_join(group_info[,c("saksnr","origin","result","group","cfu_group")],
            by = c("Sample" = "saksnr")) %>%
  filter(Rank7 %in% c(
    "pneumoniae",
    "quasipneumoniae",
    "variicola",
    "africana",
    "quasivariicola"
    ))

p_klebs_species <- ggplot(klebsiella_species, aes(reorder(Rank7, -Abundance), Abundance)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(fill = "grey80") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, face = "italic"),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(~result)


ggsave(
  here(
    "results",
    "figures",
    "06_klebs_species.png"
  ),
  p_klebs_species,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 20
)


# 02. Differential abundance analysis ----
## Import data
bracken_data <- import_biom(
  here(
    "data",
    "kraken_biom",
    "bracken_filtered.biom"
  )
)

colnames(tax_table(bracken_data)) <-
  c("kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species")

bracken_bact <- subset_taxa(bracken_data, kingdom == "k__Bacteria")

sample_data <- as.data.frame(bracken_bact@otu_table) %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column("id") %>%
  select(id)

seq_ids <- sample_data$id
metadata_ids <- group_info$saksnr

matches_seq_ids <- c()
matches_metadata_ids <- c()

for (i in seq_ids) {
  for (j in metadata_ids) {
    if (grepl(j, i)) {
      matches_seq_ids <- c(matches_seq_ids, i)
      matches_metadata_ids <- c(matches_metadata_ids, j)
    }
  }
}

sample_df <- data.frame(id = matches_seq_ids,
                        saksnr = matches_metadata_ids) %>%
  left_join(group_info[,c("saksnr","origin","result","group")]) %>%
  select(-saksnr) %>%
  column_to_rownames("id")

bracken_bact@sam_data <- sample_data(sample_df)

## Run differential abundance analysis 
bracken_norm_order <- ancombc2(
  bracken_bact,
  assay_name = "counts",
  fix_formula = "origin",
  tax_level = "order"
)

bracken_clean_tax <- psmelt(bracken_bact) %>%
  select(OTU, contains("Rank")) %>%
  group_by(OTU) %>%
  summarise_all(list(func_paste))

ancom_de_results <- bracken_norm_order$res %>%
  filter(`diff_(Intercept)` == TRUE | 
           diff_originTurkey == TRUE) %>%
  mutate(taxon = sub("o__", "", taxon)) %>%
  t() %>%
  as.data.frame %>%
  rownames_to_column("type") %>%
  mutate(host = case_when(
    grepl("intercept", ignore.case = T, type) ~ "broiler",
    grepl("turkey", ignore.case = T, type) ~ "turkey"
  ),
    type = sub("_.+", "", type),
  host = ifelse(is.na(host), "host", host)) %>%
  row_to_names(row_number = 1,
               remove_rows_above = FALSE) %>%
  rename("type" = taxon) %>%
  pivot_longer(cols = -c(type,host),
               names_to = "species",
               values_to = "value") %>%
  pivot_wider(names_from = "type",
              values_from = "value") %>%
  mutate_at(vars(-c(host, species, diff)),
            ~as.numeric(.)) %>%
  filter(!species %in% c("_4","_5")) %>%
  mutate(sign = case_when(
    q < 0.01 ~ "**",
    q < 0.05 ~ "*"
  ),
  pos = ifelse(
    lfc > 0, lfc + se + 0.08, lfc - se - 0.08
  ),
  host = ifelse(host == "broiler", "Broiler", "Turkey"))

p_diff_abund <- ggplot(ancom_de_results, aes(lfc, reorder(species,lfc), fill = host)) +
  geom_bar(color = "black",
           stat = "identity",
           linewidth = 0.3) +
  geom_errorbar(aes(xmin = lfc - se, xmax = lfc + se),
                width = 0.6,
                linewidth = 0.3) +
  geom_text(aes(label = sign,
                x = pos),
            size = 3,
            color = "#e31a1c",
            vjust = 0.75) +
  scale_fill_manual(values = palette) +
  labs(x = "Log fold change",
       y = "Orders",
       title = "C") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 6),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6))

load(
  here(
    "results",
    "RData",
    "alpha_plot.rdata"
  )
)

load(
  here(
    "results",
    "RData",
    "beta_plot.rdata"
  )
)

p_alpha_nolegend <- p_alpha +
  theme(legend.position = "none")

p_beta_nolegend <- p_beta +
  theme(legend.position = "none")


p_complete <- (p_alpha / p_beta) | p_diff_abund +
  plot_layout(guides = "collect")

p_complete_nolegend <- (p_alpha_nolegend / p_beta_nolegend) | p_diff_abund +
  plot_layout(guides = "collect")


ggsave(
  here(
    "results",
    "figures",
    "06_alpha_beta_diff_order.png"
  ),
  p_complete,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 12,
  width = 15
)

ggsave(
  here(
    "results",
    "figures",
    "06_alpha_beta_diff_order_nolegend.png"
  ),
  p_complete_nolegend,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 12,
  width = 15
)
