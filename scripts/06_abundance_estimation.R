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
library(impoRt)
library(vegan)
library(phyloseq)
library(patchwork)
library(gridExtra)

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

## Calculate relative abundance of phyla
phylum_data <- tax_glom(bracken_rel,
                        taxrank = "Rank2") %>%
  psmelt() %>%
  mutate(Sample = sub(".+(2020.+)", "\\1", Sample)) %>%
  left_join(group_info, by = c("Sample" = "saksnr")) %>%
  mutate(Rank2 = sub("p__", "", Rank2)) %>%
  filter(Abundance > 0.1)

## Abundance of Enterobacteriaceae ----
ent_data <- tax_glom(bracken_rel,
                     taxrank = "Rank5") %>%
  psmelt() %>%
  mutate(Sample = sub(".+(2020.+)", "\\1", Sample)) %>%
  left_join(group_info, by = c("Sample" = "saksnr")) %>%
  mutate(Rank5 = sub("f__", "", Rank5)) %>% 
  filter(Rank5 == "Enterobacteriaceae")

wilcox.test(Abundance ~ origin, data = ent_data, exact = FALSE)

## Generate abundance figure ----
palette <- c("Broiler" = colorspace::lighten("#006c89", 
                                             amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", 
                                            amount = 0.5))

p_phylum <- ggplot(phylum_data, aes(Sample, Abundance, fill = Rank2)) +
  geom_col(position = position_fill()) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Samples",
       y = "Relative Abundance (%)",
       fill = NULL,
       title = "B") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 6)) +
  facet_wrap(~origin, scales = "free_x")


p_ent <- ggplot(ent_data, aes(origin, Abundance)) +
  geom_violin(aes(fill = origin),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  geom_signif(comparisons = list(c("Broiler","Turkey")),
              annotations = "~italic(p) < 0.05",
              parse = TRUE) +
  scale_fill_manual(values = palette) +
  labs(y = "Relative Abundance (%)",
       title = "C") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10))


## Import alpha and beta diversity figures
load(
  here(
    "results",
    "RData",
    "alpha_div_plot.rdata"
  )
)

load(
  here(
    "results",
    "RData",
    "beta_div_plot.rdata"
  )
)

## Create full plot
p_all <- p_alpha + 
  plot_spacer() + 
  p_phylum + 
  p_ent + 
  plot_spacer() + 
  p_beta +
  plot_layout(
    ncol = 3, 
    nrow = 2,
    widths = c(
      0.7, 0.05, 1, 0.7, 0.05, 1
      )
    )

ggsave(
  here(
    "results",
    "figures",
    "06_microbiome_results.png"
  ),
  p_all,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 15,
  width = 19
)


## Abundance of Klebsiella ----
### Klebsiella genus ----
genera_data <- tax_glom(bracken_rel,
                        taxrank = "Rank6")

klebsiella_results <- psmelt(genera_data) %>%
  filter(Rank6 == "g__Klebsiella") %>%
  mutate(Sample = sub(".+-(2020-.+)", "\\1", Sample)) %>%
  select(Sample, Abundance, Rank6) %>%
  left_join(group_info[,c("saksnr","origin","result","group","cfu_group")],
            by = c("Sample" = "saksnr")) %>%
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



# 02. Differential abundance analysis ----
## Import data
load(
  here(
    "results",
    "RData",
    
  )
)









# 02. Abundance estimation StrainGE ----
strainge_results <- get_data(
  here("data",
       "strainge"),
  pattern = "stats.tsv",
  delim = "\t"
  ) %>%
  mutate(
    sample = sub(
      ".+(2020.01.....)_S.+",
      "\\1",
      sample
      )) %>%
  left_join(
    group_info[,c("saksnr",
                  "origin",
                  "result",
                  "group",
                  "cfu_g_total")],
            by = c("sample" = "saksnr")) %>%
  filter(sample != "45-zymo-std-D6311_S73")

strainge_summary <- strainge_data %>%
  mutate(`pan%` = round(`pan%`, 4)) %>%
  group_by(origin) %>%
  summarise(mean = round(mean(`pan%`), 3),
            median = round(median(`pan%`), 3),
            sd = round(sd(`pan%`), 3),
            range = paste0(min(`pan%`), " - ", max(`pan%`))) %>%
  mutate(method = "StrainGE")

write_delim(
  strainge_summary,
  here(
    "results",
    "tables",
    "05_strainGE_summary.txt"
  ),
  delim = "\t"
)

# 03. Mapping with bbsplit ----
seq_info <- read_delim(
  here(
    "data",
    "seqkit",
    "filtered_seqkit_stats_report.txt"
    ),
  delim = "\t"
  ) %>%
  filter(!grepl("zymo", file),
         !grepl("R2", file)) %>%
  mutate(file = basename(file),
         file = sub(".+-(2020-01-.+)_S.+", "\\1", file)) %>%
  select(file, num_seqs, sum_len)

mapping_results <- get_data(
  filepath = here(
    "data",
    "bbsplit",
    "filtered"
    ),
  pattern = "ref_report.txt",
  delim = "\t"
) %>%
  mutate(
    ref = sub("_bbsplit_ref_report.txt", "", ref),
    ref = sub(".+-(2020-01-.+)_.+", "\\1", ref),
    reference = case_when(
      `#name` == "GCA_000750555.1_ASM75055v1_genomic_masked" ~ "E. coli",
      `#name` == "GCA_003812345.1_ASM381234v1_genomic_masked" ~ "C. freundii",
      `#name` == "GCA_001518835.1_ASM151883v1_genomic_masked" ~ "L. adecarboxylata",
      `#name` == "GCA_000770155.1_ASM77015v1_genomic_masked" ~ "E. cloacae",
      `#name` == "GCA_013892435.1_ASM1389243v1_genomic_masked" ~ "E. fergusonii",
      `#name` == "GCA_016767755.1_ASM1676775v1_genomic_masked" ~ "E. hormaechei",
      `#name` == "GCA_015137465.1_ASM1513746v1_genomic_masked" ~ "E. kobei",
      `#name` == "GCA_000006945.2_ASM694v2_genomic_masked" ~ "S. enterica",
      `#name` == "GCA_016726285.1_ASM1672628v1_genomic_masked" ~ "S. boydii",
      `#name` == "GCF_000009885.1_ASM988v1_genomic_masked" ~ "K. pneumoniae",
      `#name` == "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic_masked" ~ "Broiler",
      `#name` == "GCA_002984395.1_ASM298439v1_genomic_masked" ~ "K. oxytoca",
      `#name` == "GCA_007632255.1_ASM763225v1_genomic_masked" ~ "K. aerogenes",
      `#name` == "GCA_013305245.1_ASM1330524v1_genomic_masked" ~ "K. variicola",
      `#name` == "GCA_016415705.1_ASM1641570v1_genomic_masked" ~ "K. quasipneumoniae",
      `#name` == "GCF_000146605.3_Turkey_5.1_genomic_masked" ~ "Turkey",
      `#name` == "hg19_main_mask_ribo_animal_allplant_allfungus_masked" ~ "Human"
    )
  ) %>%
  left_join(seq_info, by = c("ref" = "file")) %>%
  filter(!grepl("zymo", ref)) %>%
  select(
    ref,
    reference,
    unambiguousReads,
    ambiguousReads,
    assignedReads,
    num_seqs,
    sum_len
  ) %>%
  mutate(perc_reads_assigned = round(assignedReads/num_seqs* 100, 5)) %>%
  left_join(group_info[,c("saksnr","origin")], by = c("ref" = "saksnr"))


mapping_summary <- mapping_results %>%
  filter(reference %in% c("K. pneumoniae",
                          "K. aerogenes",
                          "K. variicola",
                          "K. oxytoca",
                          "K. quasipneumoniae")) %>%
  group_by(ref) %>%
  mutate(mapped_reads = sum(assignedReads)) %>%
  summarise_all(list(func_paste)) %>%
  ungroup() %>%
  select(ref, origin, num_seqs, sum_len, mapped_reads) %>%
  mutate(mapped_reads = as.numeric(mapped_reads),
         num_seqs = as.numeric(num_seqs)) %>%
  mutate(perc_reads_assigned = round(mapped_reads/num_seqs* 100, 4)) %>%
  group_by(origin) %>%
  summarise(mean = round(mean(perc_reads_assigned), 3),
            median = round(median(perc_reads_assigned), 3),
            sd = round(sd(perc_reads_assigned), 3),
            range = paste0(min(perc_reads_assigned), 
                           " - ",
                           max(perc_reads_assigned))) %>%
  mutate(method = "Mapping")

write_delim(
  mapping_summary,
  here(
    "results",
    "tables",
    "05_mapping_summary.txt"
  ),
  delim = "\t"
)


# 04. Statistics ----

wilcox.test(Abundance ~ origin,
            data = kraken_results,
            alternative = "two.sided")

wilcox.test(`pan%` ~ origin,
            data = strainge_results,
            alternative = "two.sided")

wilcox.test(perc_reads_assigned ~ origin,
            data = mapping_results,
            alternative = "two.sided")

# 05. Figure ----
## Define palette
palette <- c("Broiler" = colorspace::lighten("#006c89", 
                                             amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", 
                                            amount = 0.5))

## Kraken results
p_kraken <- ggplot(kraken_results, 
                   aes(origin, Abundance)) +
  geom_violin(aes(fill = origin),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  labs(y = "Relative abundance (%)",
       title = "Kraken") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

## StrainGE results
p_strainge <- ggplot(strainge_data, 
                     aes(origin, `pan%`)) +
  geom_violin(aes(fill = origin),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  labs(y = "Pan%",
       title = "StrainGE") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())

## Mapping results
p_mapping <- ggplot(mapping_results, 
                    aes(origin, perc_reads_assigned)) +
  geom_violin(aes(fill = origin),
              trim = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  labs(y = "Percent (%) mapped reads",
       title = "Mapping") +
  scale_fill_manual(values = palette) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank())


## Combine plots

stats_table <- rbind(kraken_summary,
                     strainge_summary,
                     mapping_summary) %>%
  select(method, everything(), -range) %>%
  rename("Method" = method,
         "Host" = origin,
         "Mean" = mean,
         "Median" = median,
         "SD" = sd)

summary_table <- ggtexttable(stats_table,
            rows = NULL,
            theme = ttheme("blank")) %>%
  tab_add_hline(at.row = c(1, 2), 
                row.side = "top", 
                linewidth = 3, 
                linetype = 1) %>%
  tab_add_hline(at.row = c(7),
                row.side = "bottom",
                linewidth = 3,
                linetype = 1)

p_all <- p_kraken + p_strainge + p_mapping + summary_table +
  plot_layout(
    nrow = 2,
    ncol = 2,
    guides = "collect")

ggsave(
  here(
    "results",
    "figures",
    "05_abundance_estimations.png"
  ),
  p_all,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 18,
  width = 18
)
