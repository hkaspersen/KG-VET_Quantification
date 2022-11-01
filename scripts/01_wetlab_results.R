# ABSTRACT
# Script for summarising and running statistics
# on the culture quantification data and ZKIR
# qPCR results

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(impoRt)
library(readxl)
library(here)
library(tibble)
library(readr)

# 01. Calculate occurrence ----
quant_results <- read_xlsx(
  here("data",
       "list_quant_results.xlsx"
       )
  ) %>%
  ## These samples are removed, as they represent duplicated
  ## samples from the campylobacter samples
  filter(
    !saksnr %in% c(
      "2020-01-3390",
      "2020-01-3391",
      "2020-01-3985",
      "2020-01-4081",
      "2020-01-4082"
    )
  ) %>%
  mutate(
    result = ifelse(is.na(species_frozen),
                    "Negative",
                    "Positive"),
    result_kp = ifelse(
      grepl("K. pneumoniae|K. variicola",
            species_frozen),
      "Positive",
      "Negative"
    )
  )

## Function for calculating 95% confidence intervals
get_binCI <- function(x, n) as.numeric(
  setNames(
    binom.test(
      x, n, p = 0.5, "two.sided", conf.level = 0.95
      )$conf.int*100,
    c("lwr", "upr")
    )
  )

## Calculate occurrence
### Occurrence of Klebsiella spp.
occurrence_report <- quant_results %>%
  select(origin, result) %>%
  group_by_all() %>%
  count() %>%
  pivot_wider(names_from = "result",
              values_from = "n") %>%
  mutate(Total = Positive + Negative,
         Percent = round(Positive/Total*100, 2)) %>%
  rowwise() %>%
  mutate(lwr = get_binCI(Positive, Total)[1],
         upr = get_binCI(Positive, Total)[2],
         type = "Klebsiella spp.")

### Occurrence report for Kpn
occurrence_report_kp <- quant_results %>%
  select(origin, result_kp) %>%
  group_by_all() %>%
  count() %>%
  pivot_wider(names_from = "result_kp",
              values_from = "n") %>%
  mutate(Total = Positive + Negative,
         Percent = round(Positive/Total*100, 2)) %>%
  rowwise() %>%
  mutate(lwr = get_binCI(Positive, Total)[1],
         upr = get_binCI(Positive, Total)[2],
         type = "K. pneumoniae")

## Combine reports and save to file
total_occurrence <- rbind(occurrence_report,
                          occurrence_report_kp)

write_delim(
  total_occurrence,
  here(
    "results",
    "tables",
    "01_occurrence_data.txt"
  ),
  delim = "\t"
)

# 02. Calculate abundance ----
## Calculate CFU/g material
cfu_calc <- quant_results %>%
  select(
    saksnr,
    origin,
    result,
    result_kp,
    undiluted,
    `10-1`,
    `10-2`,
    `10-3`,
    `10-4`) %>%
  # convert the dilution columns to numeric
  mutate_at(vars(contains("10")),
            ~ as.numeric(.)) %>%
  # Calculate CFU for each dilution
  # Initial sample was diluted 1:10 in peptone water, therefore already
  # diluted
  mutate(cfu_10_1 = ifelse(
           is.na(`10-1`) == TRUE, NA, `10-1` / 0.01
           ),
         cfu_10_2 = ifelse(
           is.na(`10-2`) == TRUE, NA, `10-2` / 0.001
           ),
         cfu_10_3 = ifelse(
           is.na(`10-3`) == TRUE, NA, `10-3` / 0.0001
           ),
         cfu_10_4 = ifelse(
           is.na(`10-4`) == TRUE, NA, `10-4` / 0.00001
           )
         ) %>%
  # Remove all 0 values, to get a correct mean
  mutate_at(vars(contains("cfu_10")),
            ~ ifelse(. == 0, NA, .)) %>%
  # Calculate the mean values
  mutate(cfu_g_total = ifelse(
    rowSums(.[,10:13], na.rm = TRUE) > 0,
      apply(.[,10:13], 1, FUN = mean, na.rm = TRUE),
      NA
    ),
    log_cfu = log10(cfu_g_total)
  )

quant_stats <- cfu_calc %>%
  filter(cfu_g_total > 0) %>%
  # Calculate mean and median for each group
  group_by(origin) %>%
  mutate(median = median(cfu_g_total),
         mean = round(mean(cfu_g_total), 2),
         samples = n()) %>%
  ungroup() %>%
  mutate(cfu_g_total = round(cfu_g_total, 0),
         mean_total = round(mean(cfu_g_total), 2),
         median_total = round(median(cfu_g_total), 2)) %>%
  select(origin, samples, median, mean, mean_total, median_total) %>%
  group_by(origin) %>%
  summarise_all(list(func_paste))

# 03. Summarise ZKIR qPCR results ----
## Import ZKIR qPCR results
zkir_results <- read_xlsx(
  here(
    "data",
    "ZKIR_results.xlsx"
    )) %>%
  filter(!saksnr %in% c(
    "2020-01-3390",
    "2020-01-3391",
    "2020-01-3985",
    "2020-01-4081",
    "2020-01-4082"
    )) %>%
  mutate(
    ZKIR = ifelse(
      ZKIR == "Positive/negative?", 
      "Negative",
      ZKIR
      )
    )

zkir_stats <- zkir_results %>%
  count(origin, ZKIR) %>%
  pivot_wider(names_from = "ZKIR",
              values_from = "n") %>%
  mutate(Total = Negative + Positive,
         Percent = round(Positive/Total * 100, 2)) %>%
  rename("Host" = origin)

zkir_quant_results <- cfu_calc %>%
  select(saksnr, origin, result, cfu_g_total) %>%
  left_join(zkir_results[,c("saksnr","ZKIR")],
            by = "saksnr")
  
write_delim(
  zkir_stats,
  here(
    "results",
    "tables",
    "01_zkir_qpcr_results.txt"
  ),
  delim = "\t"
)

write_delim(
  zkir_quant_results,
  here(
    "results",
    "tables",
    "01_zkir_quant_results.txt"
  ),
  delim = "\t"
)

# 04. Run statistics ----
## Chi Squared test for difference in occurrence
## for Klebsiella spp. and Klebsiella pneumoniae
occurrence_report %>%
  select(origin, Positive, Total) %>%
  column_to_rownames("origin") %>%
  as.matrix %>%
  chisq.test(correct = FALSE)

occurrence_report_kp %>%
  select(origin, Positive, Total) %>%
  column_to_rownames("origin") %>%
  as.matrix %>%
  chisq.test(correct = FALSE)

## Wilcox test for overall difference in abundance
wilcox.test(cfu_g_total ~ origin, data = cfu_calc, alternative = "two.sided")



# 05. Create figure ----
## Define palette
palette <- c("Broiler" = colorspace::lighten("#006c89", amount = 0.5),
             "Turkey" = colorspace::lighten("#3d6721", amount = 0.5))

## Occurrence plot
p_occurrence <- ggplot(occurrence_report, 
                       aes(
                         origin,
                         Percent,
                         fill = origin
                         )
                       ) +
  geom_col(color = "black",
           width = 0.5) +
  geom_errorbar(aes(ymin = lwr,
                    ymax = upr),
                width = 0.3) +
  geom_signif(comparisons = list(c("Broiler","Turkey")),
              annotations = "~italic(p) < 0.01",
              parse = TRUE) +
  scale_fill_manual(values = palette) +
  scale_y_continuous(limits = c(0, 90)) +
  labs(y = "Percent (%) occurrence of *Klebsiella*") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = 7))


## Quantification plot
p_abundance <- cfu_calc %>%
  filter(cfu_g_total > 0) %>%
  ggplot(aes(origin, log_cfu)) +
  geom_violin(trim = FALSE,
              scale = "count",
              adjust = 0.5,
              aes(fill = origin)) +
  stat_boxplot(geom = "errorbar", width = 0.05) +
  geom_boxplot(width = 0.05) +
  geom_jitter(width = 0.05,
              alpha = 0.3) +
  geom_signif(comparisons = list(c("Broiler","Turkey")),
              annotations = "~italic(p) < 0.05",
              parse = TRUE) +
  scale_fill_manual(values = palette) +
  labs(x = NULL,
       y = "Log10 CFU/g *Klebsiella*") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = ggtext::element_markdown(size = 7),
        axis.ticks.x = element_blank())

## Arrange plots
p_all <- p_occurrence / p_abundance +
  plot_layout(ncol = 1, nrow = 2, heights = c(0.4, 0.6))


## Save plot
ggsave(
  here(
    "results",
    "figures",
    "01_wetlab_results.png"
    ),
  p_all,
  device = "png",
  units = "cm",
  dpi = 600,
  height = 14,
  width = 8
)

