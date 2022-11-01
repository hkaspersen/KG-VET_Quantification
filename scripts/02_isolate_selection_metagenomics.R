# ABSTRACT
# Script for isolate selection for metagenomics
# sequencing

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(here)
library(readxl)

# Import data
data <- read_delim(
  here(
    "results",
    "tables",
    "01_zkir_quant_results.txt"
    ),
  delim = "\t") %>%
  # Define groups
  mutate(group = case_when(origin == "Turkey" & 
                             result == "Negative" & 
                             ZKIR == "Negative" ~ "turkey_negative",
                           origin == "Turkey" & 
                             result == "Positive" & 
                             ZKIR == "Positive" ~ "turkey_positive",
                           origin == "Broiler" & 
                             result == "Negative" & 
                             ZKIR == "Negative" ~ "broiler_negative",
                           origin == "Broiler" & 
                             result == "Positive" & 
                             ZKIR == "Positive" ~ "broiler_positive",
                           TRUE ~ "exclude")) %>%
  filter(group != "exclude",
         !saksnr %in% c("2020-01-3390",
                        "2020-01-3391",
                        "2020-01-3985",
                        "2020-01-4081",
                        "2020-01-4082"))

sample_groups <- data %>%
  filter(!is.na(cfu_g_total)) %>%
  group_by(group) %>%
  mutate(mean_val = round(mean(cfu_g_total), 2),
         median_val = median(cfu_g_total),
         range = paste0(min(cfu_g_total), " - ", max(cfu_g_total)),
         n_samples = n()) %>%
  summarise_all(list(func_paste)) %>%
  select(group, origin, n_samples, mean_val, median_val, range)

write_delim(data,
            here("results",
                 "tables",
                 "02_sample_groups_metagenomics.txt"),
            delim = "\t")


groups <- data %>%
  filter(cfu_g_total > 0) %>%
  mutate(cfu_group = case_when(cfu_g_total <= 500 ~ "Low",
                               cfu_g_total > 500 & cfu_g_total <= 5000 ~ "Medium",
                               cfu_g_total > 5000 ~ "High")) %>%
  group_by(group, cfu_group) %>%
  mutate(mean_val = round(mean(cfu_g_total), 2),
         median_val = median(cfu_g_total),
         range = paste0(min(cfu_g_total), " - ", max(cfu_g_total)),
         n_samples = n()) %>%
  summarise_all(list(func_paste)) %>%
  select(group, cfu_group, origin, n_samples, mean_val, median_val, range)

write_delim(groups,
            here("results",
                 "tables",
                 "02_subgroup_info_metagenomics.txt"),
            delim = "\t")


group_info <- data %>%
  mutate(cfu_group = case_when(cfu_g_total <= 500 ~ "Low",
                               cfu_g_total > 500 & cfu_g_total <= 5000 ~ "Medium",
                               cfu_g_total > 5000 ~ "High",
                               TRUE ~ "Negative"))

write_delim(group_info,
            here("results",
                 "tables",
                 "02_group_info_metagenomics_detailed.txt"),
            delim = "\t")


set.seed(12345)

random_sampling_positive <- data %>%
  filter(cfu_g_total > 0) %>%
  mutate(cfu_group = case_when(cfu_g_total <= 500 ~ "Low",
                               cfu_g_total > 500 & cfu_g_total <= 5000 ~ "Medium",
                               cfu_g_total > 5000 ~ "High")) %>%
  group_by(group, cfu_group) %>%
  slice_sample(n = 4) 


random_sampling_negative <- data %>%
  filter(cfu_g_total == 0) %>%
  group_by(group) %>%
  slice_sample(n = 10) 


total_sampled_data <- rbind(random_sampling_negative,
                            random_sampling_positive) %>%
  select(saksnr,
         origin,
         result,
         ZKIR,
         cfu_g_total,
         qubit,
         group,
         cfu_group)

sampled_data_stats <- total_sampled_data %>%
  group_by(group, cfu_group) %>%
  mutate(mean_val = round(mean(cfu_g_total), 2),
         median_val = median(cfu_g_total),
         range = paste0(min(cfu_g_total), " - ", max(cfu_g_total)),
         n_samples = n()) %>%
  summarise_all(list(func_paste)) %>%
  select(group, cfu_group, origin, n_samples, mean_val, median_val, range)

# Find position in freezer
extract_list <- read_xlsx(here("data","list_isol_extraction.xlsx")) %>%
  select(saksnr, source, nd_280, nd_230, Box, Pos, Note)

samples <- total_sampled_data %>%
  left_join(extract_list, by = "saksnr")


# Write to file
write_delim(samples,
            here("results","samples_for_metagenomics.txt"),
            delim = "\t")

write_delim(sampled_data_stats,
            here("results","stats_selected_samples.txt"),
            delim = "\r")
