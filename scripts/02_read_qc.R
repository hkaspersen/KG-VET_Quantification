# ABSTRACT
# Script for general QC analysis of the
# metagenomic shotgun reads

# Libraries ----
library(dplyr)
library(tidyr)
library(readr)
library(funtools)
library(ggplot2)
library(readr)
library(readxl)
library(here)

# 01. Seqkit stats of the reads ----
seqkit_data_filtered <- read_delim(
  here("data",
       "seqkit",
       "filtered_seqkit_stats_report.txt"),
  delim = "\t"
) %>%
  mutate(file = sub(".+\\/.+-(2020.+)_S.+", "\\1", file),
         file = sub(".+\\/45-(zymo.+)_S.+", "\\1", file)) %>%
  group_by(file) %>%
  summarise_all(list(func_paste)) %>%
  select(file,
         num_seqs) %>%
  rename("number_of_reads_filtered" = num_seqs)


seqkit_data <- read_delim(
  here("data",
       "seqkit",
       "raw_seqkit_stats_report.txt"),
  delim = "\t"
) %>%
  mutate(file = sub(".+\\/.+-(2020.+)_S.+", "\\1", file),
         file = sub(".+\\/45-(zymo.+)_S.+", "\\1", file)) %>%
  group_by(file) %>%
  summarise_all(list(func_paste)) %>%
  select(file,
         min_len,
         sum_len,
         num_seqs) %>%
  rename("number_of_reads_raw" = num_seqs,
         "total_bases_raw" = sum_len,
         "read_length" = min_len) %>%
  left_join(seqkit_data_filtered) %>%
  mutate_at(vars(contains("number")),
            ~as.numeric(.)) %>%
  mutate(diff = number_of_reads_raw - number_of_reads_filtered,
         diff_percent = round(diff/number_of_reads_raw*100, 2))

## Write data to file
write_delim(
  seqkit_data,
  here(
    "results",
    "tables",
    "02_read_stats.txt"
  ),
  delim = "\t"
)
