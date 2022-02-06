# Generate FASTA File for the UCE

library(tidyverse)
library(Biostrings)

fto <- readDNAStringSet("data/NC_000016.9.fasta")[[1]]

extract_sequences <- function(sequence, start, end) {
  subsequence <- subseq(sequence, start, end)
  toString(subsequence)
}


uce_coordinates <- read_csv("data/uce_coordinates.csv") |>
  select(uce, contains("37"), length) |>
  rename_with(~ str_remove(.x, "_hg37"))

uce_data <- uce_coordinates |>
  mutate(sequence = map2_chr(start_pos, end_pos, ~ extract_sequences(fto, .x, .y)),
         name = str_c("NC_000016.9", "FTO", uce, start_pos, end_pos, length, sep = " | ")) |>
  select(name, everything()) |>
  write_csv("data/uce_data_with_sequences.csv")

uce_data |>
  select(name, sequence) |>
  deframe() |>
  DNAStringSet() |>
  writeXStringSet("data/uce_seqs.fasta")
