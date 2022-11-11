# ABSTRACT
# Script for normalizing bracken data and 
# running differential analysis 

# Libraries ----
library(readr)
library(phyloseq)
library(here)
library(tibble)
library(dplyr)
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

bracken_data <- import_biom(
  here(
    "data",
    "kraken_biom",
    "bracken_filtered.biom"
  )
)

# 02. Attach sample info to biom file ----
sample_data <- as.data.frame(bracken_data@otu_table) %>%
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

bracken_data@sam_data <- sample_data(sample_df)

# 03. Normalize and run differential analysis ----
# This is for taxa that are different between
# host species
bracken_norm <- ancombc2(
  bracken_data,
  assay_name = "counts",
  fix_formula = "origin"
)

save(bracken_norm,
     file = here(
       "results",
       "RData",
       "bracken_norm.rdata"
  )
)
