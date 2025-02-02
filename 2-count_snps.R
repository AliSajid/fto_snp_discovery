# Identify how many SNPs are between each UCE

library(tidyverse)
library(broom)

set.seed(1989)

snp_locations <-
  read_csv("data/snp_locations.csv",
           col_types = cols(.default = col_character())) |>
  transmute(
    uid = as.character(uid),
    chrpos_hg38 = as.numeric(str_remove(chrpos, "16:")),
    chrpos_hg37 = as.numeric(str_remove(chrpos_prev_assm, "16:"))
  ) |>
  select(-chrpos_hg38) |>
  mutate(chrpos = chrpos_hg37)

uce_spans <- read_csv("data/uce_coordinates.csv") |>
  select(uce, contains("37"), length) |>
  rename_with(~ str_remove(.x, "_hg37"))



count_snps <- function(uce, start_pos, end_pos, length) {
  res <- snp_locations |>
    filter(chrpos_hg37 <= end_pos & chrpos_hg37 >= start_pos) |>
    nrow()

  unlist(res)
}


generate_random_interval <-
  function(uce,
           length,
           start = 53701692,
           end = 54158512) {
    start <- sample(start:(end - length), 1)
    out <- list(
      uce = str_c(uce, str_replace(uce, "U", "X"), sep = "-"),
      start_pos = start,
      end_pos = start + length - 1,
      length = length
    )
    out
  }


generate_samples <-
  function(n,
           uce,
           length,
           start = 53701692,
           end = 54158512) {
    res <-
      map_dfr(1:n, ~ generate_random_interval(uce, length, start, end))
    res
  }

t_test_comparison <- function(dataset, mean, alternative = "two.sided") {
x <- dataset$snp_count

res <- t.test(x, mu = mean, alternative = alternative)
res
}

uce_randoms <- uce_spans |>
  mutate(randoms = map2(uce, length, ~ generate_samples(1000, .x, .y)),
         filename = str_glue("data/{uce}-samples_snp_count_data.csv"),
         counted_randoms = map(randoms, \(x) {
           x |> mutate(snp_count = pmap_dbl(
             list(uce, start_pos, end_pos, length),
             ~ count_snps(..1, ..2, ..3, ..4)
           ))
         })
         ) |>
  mutate(x = map2(counted_randoms, filename, ~ write_csv(.x, .y)),) |>
  select(-x, -filename, -randoms)



uce_data <- uce_spans |>
  mutate(
    snp_count = pmap_dbl(
      list(uce, start_pos, end_pos, length),
      ~ count_snps(..1, ..2, ..3, ..4)
    )) |>
  write_csv("data/uce_snp_count_data.csv") |>
  inner_join(uce_randoms, by = c("uce", "start_pos", "end_pos", "length")) |>
  expand_grid(alternative = c("two.sided", "less", "greater")) |>
  mutate(comparison_results = pmap(list(counted_randoms, snp_count, alternative), ~ t_test_comparison(..1, ..2, ..3)),
         results = map(comparison_results, glance)) |>
  select(-alternative, -counted_randoms, -comparison_results) |>
  unnest_wider(results) |>
  write_csv("snp_count_comparison.csv")

