library("readr")
library("janitor")
library("tibble")
library("dplyr")

small <- readr::read_rds("small-water-isotope-data.rds")
small <- small |>
  janitor::clean_names() |>
  tibble::as_tibble()

write_csv(small, "small-water-isotopes.csv")

