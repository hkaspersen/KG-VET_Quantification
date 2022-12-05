# ABSTRACT
# Resistome analysis

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(RColorBrewer)
library(impoRt)
library(readxl)
library(here)
library(patchwork)
library(janitor)
library(vegan)
library(phyloseq)
library(ANCOMBC)

# 01. Import data ----
group_info <- read_delim(
  here(
    "results",
    "tables",
    "02_group_info_metagenomics_detailed.txt"
  ),
  delim = "\t"
)

resistome_reports <- get_data(
  filepath = here(
    "data",
    "groot"
  ),
  pattern = "_report.txt",
  col_names = FALSE
)

# 02. Create phyloseq object ----
resistome_otu <- resistome_reports %>%
  select(X1, ref, X2) %>%
  filter(ref != "45-zymo-std-D6311") %>%
  pivot_wider(names_from = "ref",
              values_from = X2,
              values_fill = 0) %>%
  filter(!is.na(X1)) %>%
  mutate(OTU = paste0("ARG", 1:n())) %>%
  column_to_rownames("OTU") %>%
  select(-X1)

resistome_tax <- resistome_reports %>%
  select(X1) %>%
  group_by(X1) %>%
  summarise_all(list(func_paste)) %>%
  ungroup() %>%
  filter(!is.na(X1)) %>%
  mutate(OTU = paste0("ARG", 1:n())) %>%
  column_to_rownames("OTU") %>%
  rename("gene" = X1) %>%
  mutate(family = case_when(
    grepl("tet", gene) ~ "Tetracycline",
    grepl("sul", gene) ~ "Sulfonamide",
    grepl("erm", gene) ~ "Macrolide",
    grepl("aph|ant|aac|str", gene) ~ "Aminoglycoside",
    grepl("cat", gene) ~ "Phenicol",
    grepl("cep", gene) ~ "Cephalosporin",
    grepl("Van", gene) ~ "Glycopeptide"
  )) %>%
  select(family, gene) %>%
  as.matrix

otu_table_res <- otu_table(resistome_otu, taxa_are_rows = T)
tax_table_res <- tax_table(resistome_tax)

resistome_phyloseq <- phyloseq(otu_table_res,
                               tax_table_res)

# 03. Attach metadata to phyloseq object
sample_data <- as.data.frame(resistome_phyloseq@otu_table) %>%
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
  left_join(group_info[,c("saksnr","origin","result")]) %>%
  select(-saksnr) %>%
  column_to_rownames("id")

resistome_phyloseq@sam_data <- sample_data(sample_df)

# 04. Run normalization and diff abundance ----
resistome_norm <- ancombc2(
  resistome_phyloseq,
  assay_name = "counts",
  fix_formula = "origin"
)

otu_table_norm <- as.data.frame(resistome_norm$feature_table) %>%
  otu_table(taxa_are_rows = T)

tax_table_norm <- as.data.frame(resistome_tax) %>%
  filter(rownames(.) %in% rownames(otu_table_norm)) %>%
  as.matrix %>%
  tax_table()

norm_res <- phyloseq(otu_table_norm,
                     tax_table_norm)

# 05. Beta diversity ----
ord_resistome <- ordinate(
  physeq = norm_res,
  method = "NMDS",
  distance = "bray"
)

plot_data <- as.data.frame(ord_resistome$points) %>%
  rownames_to_column("ref") %>%
  filter(!grepl("zymo", ref)) %>%
  mutate(ref = sub(".+-(2020.+)", "\\1", ref)) %>%
  left_join(group_info, by = c("ref" = "saksnr"))

otu_table <- psmelt(norm_res) %>%
  mutate(Sample = sub(".+(2020.+)", "\\1", Sample)) %>%
  select(-c(family, gene)) %>%
  pivot_wider(names_from = "OTU",
              values_from = "Abundance") %>%
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

palette <- c("Broiler" = colorspace::lighten("#006c89", 
                                             amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", 
                                            amount = 0.5))

p_beta <- ggplot(plot_data, aes(MDS1, MDS2, color = origin)) +
  stat_ellipse(linetype = "dotted",
               linewidth = 0.3) +
  geom_point(aes(shape = result),
             size = 1.2) +
  labs(shape = NULL,
       color = NULL,
       title = "B") +
  annotate(geom = "text",
           label = paste0("Stress: ", round(ord_resistome$stress, 4)),
           x = -0.8,
           y = -0.6,
           size = 2) +
  scale_color_manual(values = palette) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size = 5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6))

p_beta_nolegend <- p_beta + theme(legend.position = "none")

# 06. Differential Abundance ----
merge_data <- as.data.frame(tax_table_norm) %>%
  rownames_to_column("OTU") %>%
  mutate(gene_clean = sub("(.+[0-9])_[A-Z].+", "\\1", gene),
         gene_aggregated = sub("_.+", "", gene))

normalized_counts <- as.data.frame(resistome_norm$feature_table) %>%
  rownames_to_column("OTU") %>%
  left_join(merge_data) %>%
  pivot_longer(cols = -c(OTU, family, gene, gene_clean, gene_aggregated),
               names_to = "Sample",
               values_to = "Abundance") %>%
  select(OTU, Sample, Abundance, family, gene_clean, gene_aggregated) %>%
  group_by(Sample, gene_aggregated) %>%
  summarise(Abundance = sum(Abundance)) %>%
  group_by(Sample) %>%
  mutate(rel_abundance = Abundance/sum(Abundance)*100,
         Sample = sub(".+(2020-)", "\\1", Sample)) %>%
  ungroup() %>%
  left_join(group_info, by = c("Sample" = "saksnr")) %>%
  mutate(gene_aggregated = sub("Van", "van", gene_aggregated),
         plot_group = paste0(origin, " ", result))

palette_genes <- c(
  "aac(6')-aph(2'')" = "#8DD3C7",
  "ant(6)-Ia" = "#FFFFB3",
  "erm(B)" = "#BEBADA",
  "cat" = "#FB8072",
  "ant(6)-Ib" = "#80B1D3",
  "cepA-29" = "#FDB462",
  "aph(3'')-Ib" = "#B3DE69",
  "vanY-B" = "#FCCDE5",
  "sul2" = "#D9D9D9",
  "aph(3')-III" = "#BC80BD",
  "aph(6)-Id" = "#CCEBC5",
  "erm(F)" = "#FFED6F"
)


p_rel_abund_res <- ggplot(normalized_counts, 
                          aes(Sample,
                              rel_abundance,
                              fill = gene_aggregated)) +
  geom_col() +
  labs(x = "Samples",
       y = "Relative abundance (%)",
       fill = NULL,
       title = "A") +
  scale_fill_manual(values = palette_genes) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.2, 'cm'),
        legend.text = element_text(size = 5),
        strip.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text.y = element_text(size = 5)) +
  facet_wrap(~plot_group, scales = "free_x")

p_rel_abund_res_nolegend <- p_rel_abund_res + theme(legend.position = "none")

p_resistome_nolegend <- p_rel_abund_res_nolegend | (p_beta / guide_area()) +
  plot_layout(guides = "collect")

p_resistome_legend <- p_rel_abund_res | (p_beta / guide_area()) +
  plot_layout(guides = "collect")

ggsave(
  here(
    "results",
    "figures",
    "07_diff_abundance_res_nolegend.png"
  ),
  p_resistome_nolegend,
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
    "07_diff_abundance_res.png"
  ),
  p_resistome_legend,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 12,
  width = 15
)


# 07. Counts per host ----
host_counts <- resistome_reports %>%
  select(ref, X1) %>%
  filter(!grepl("zymo", ref)) %>%
  rename("gene" = X1) %>%
  mutate(gene_clean = sub("(.+[0-9])_[A-Z].+", "\\1", gene)) %>%
  mutate(ref = sub(".+(2020.+)", "\\1", ref)) %>%
  left_join(group_info, by = c("ref" = "saksnr")) %>%
  group_by(gene_clean, origin) %>%
  count() %>%
  ungroup()

p_counts <- ggplot(host_counts, aes(reorder(gene_clean, -n), n)) +
  geom_col() +
  geom_text(aes(label = n),
            vjust = -0.3,
            size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.3, 
                                   size = 6),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(~origin)

ggsave(
  here(
    "results",
    "figures",
    "07_res_gene_counts.png"
  ),
  p_counts,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 20,
  width = 20
)
