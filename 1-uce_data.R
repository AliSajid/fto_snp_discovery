# Create the UCE region marker data

library(tidyverse)

uce <- str_c("UCE", str_pad(1:10, width = 4, pad = "0"))
start_pos_hg37 <- c(53756950, 53929295, 53957123, 54020764, 54056490, 54093096, 54093563, 54101268, 54127008, 54144134)

end_pos_hg37 <- c(53757285, 53929573, 53957516, 54020989, 54056818, 54093465, 54093797, 54101623, 54127377, 54144343)


start_pos_hg38 <- c(53723038, 53895383, 53923211, 53986852, 54022578, 54059184, 54059651, 54067356, 54093096, 54110222)

end_pos_hg38 <- c(53723373, 53895661, 53923604, 53987077, 54022906, 54059553, 54059885, 54067711, 54093465, 54110431)

uce_data <- tibble(uce, start_pos_hg37, end_pos_hg37, start_pos_hg38, end_pos_hg38) |>
  mutate(length = end_pos_hg37 - start_pos_hg37 + 1) |>
  unique() |>
  write_csv("data/uce_coordinates.csv")
