# Retrieve all FTO SNPs

library(tidyverse)
library(reutils)
library(httr)
library(jsonlite)

options(reutils.email = "Ali.Imami@rockets.utoledo.edu",
        reutils.api.key = "a20864ad0756cdd3c694893bedc036a99708")

retmode <- "json"
retmax <- 150000

query <-
  "((FTO[Gene Name]) AND Homo sapiens[Organism]) AND SNV[SNP Class]"

if (!file.exists("data/snp_uids")) {
  uids <-
    esearch(query,
            db = "snp",
            retmode = retmode,
            retmax = retmax) |>
    uid() |>
    str_c(collapse = "\n") |>
    write_file("data/snp_uids") |>
    str_split("\n", simplify = TRUE) |>
    t()
} else {
  uids <- read_file("data/snp_uids") |>
    str_split("\n", simplify = TRUE) |>
    t()
}

make_summary_request <- function(uids) {
  url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
  query <- list(
    db = "snp",
    api_key = options()$reutils.api.key,
    email = options()$reutils.email,
    tool = "SNPInformant",
    id = str_c(uids, collapse = ","),
    retmode = "json",
    version = "2.0"
  )

  r <- POST(url, query = query)

  r
}


get_summary_data <- function(uids, start, step = 10000) {
  slice <- uids[start:(start + step - 1)]
  out <- make_summary_request(slice)

  out
}

parse_results <- function(json) {
  data <- fromJSON(json)
  ids <- data |> pluck("result", 1)

  parsed_data <- ids |>
    map( ~ pluck(data, "result", .x)) |>
    map(enframe) |>
    map( ~ column_to_rownames(.x, "name")) |>
    map( ~ t(.x)) |>
    map(as_tibble) |>
    map_dfr( ~ select(.x, uid, chrpos, chrpos_prev_assm)) |>
    mutate(across(everything(), unlist))

  parsed_data
}


if (file.exists("data/snp_locations.csv")) {
  dataset <- read_csv("data/snp_locations.csv", col_types = cols(.default = col_character()))
} else {
  dataset <-
    map(seq(1, nrow(uids), 112), ~ get_summary_data(uids, .x, 112), ) |>
    map( ~ content(.x, as = "text", encoding = "UTF-8")) |>
    map_dfr( ~ parse_results(.x)) |>
    write_csv("data/snp_locations.csv")
}
