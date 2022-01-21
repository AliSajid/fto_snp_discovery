# Create the UCE region marker data

library(tidyverse)

uce <- str_c("UCE", str_pad(1:10, width = 2, pad = "0"))
start_pos <- c(53756950, 53929295, 53957123, 53957123, 54020764, 54056490, 54093096, 54093563, 54101268, 54144134)

end_pos <- c(53757285, 53929573, 53957516, 53957516, 54020989, 54056818, 54093465, 54093797, 54101623, 54144343)

uce_data <- tibble(uce, start_pos, end_pos) |>
  mutate(length = end_pos - start_pos + 1) |>
  write_csv("data/uce_coordinates.csv") |>
  unique()
