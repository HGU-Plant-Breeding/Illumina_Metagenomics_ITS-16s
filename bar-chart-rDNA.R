#!/usr/bin/env Rscript

library(tidyverse)

setwd("/home/max/Metagenomics_Mikrobio")

# input / output
infile <- "species_matrix_rDNA.tsv"
out_prefix <- sub("\\.tsv$", "", infile)

tax_cols <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# load table
df <- read.delim(
  infile,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_cols <- setdiff(colnames(df), tax_cols)

# total reads per sample
sample_totals <- colSums(df[, sample_cols], na.rm = TRUE)

# keep rows that reach >=1% in at least one sample
keep <- apply(df[, sample_cols], 1, function(x) {
  any(x / sample_totals >= 0.005, na.rm = TRUE)
})

df_filt <- df[keep, ]

# make label for plotting
df_filt <- df_filt %>%
  mutate(
    taxon_label = case_when(
      genus != "" & species != "" ~ paste(genus, species),
      genus != "" ~ genus,
      family != "" ~ family,
      TRUE ~ "Unclassified"
    )
  )

# long format
long_abs <- df_filt %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "reads"
  )

long_rel <- long_abs %>%
  group_by(sample) %>%
  mutate(relative_reads = reads / sum(reads) * 100) %>%
  ungroup()

# absolute stacked barplot
p_abs <- ggplot(long_abs, aes(x = sample, y = reads, fill = taxon_label)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = "Sample",
    y = "Read count",
    title = "Absolute read counts"
  )

print(p_abs)

# relative stacked barplot
p_rel <- ggplot(long_rel, aes(x = sample, y = relative_reads, fill = taxon_label)) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = "Sample",
    y = "Relative abundance (%)",
    title = "Relative read abundance"
  )

print(p_rel)

# also save filtered table
write.table(
  df_filt,
  paste0(out_prefix, "_filtered_min0-5percent.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
