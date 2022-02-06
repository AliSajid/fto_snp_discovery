# Random Samples and their FASTA File

library(tidyverse)
library(Biostrings)

files <- list.files("data", "UCE")

filepaths <- file.path("data", files)

fto <- readDNAStringSet("data/NC_000016.9.fasta")[[1]]

extract_sequences <- function(sequence, start, end) {
  subsequence <- subseq(sequence, start, end)
  toString(subsequence)
}

create_fasta <- function(dataset) {

  dataset_name <- str_to_lower(str_extract(dataset, "UCE\\d{4}"))

  sequences <- read_csv(dataset) |>
    separate(uce, into = c("uce", "xce")) |>
    mutate(xce = str_c("XCE", str_pad(row_number(xce), width = 4, pad = "0"))) |>
    mutate(sequence = map2_chr(start_pos, end_pos, ~ extract_sequences(fto, .x, .y)),
           name = str_c("NC_000016.9", "FTO", uce, xce, start_pos, end_pos, length, sep = " | ")) |>
    select(name, everything()) |>
    write_csv(str_glue("data/{dataset_name}_data_with_sequences.csv")) |>
    select(name, sequence) |>
    deframe() |>
    DNAStringSet() |>
    writeXStringSet(str_glue("data/{dataset_name}_xce_seqs.fasta"))
}

filepaths |>
  walk(create_fasta)
