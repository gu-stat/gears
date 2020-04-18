## code to prepare `commodities_prices` dataset goes here

commodities_prices <- read.csv(
  file = "data-raw/commoditiesprices.csv",
  stringsAsFactors = FALSE
)

commodities_prices_data_dictionary <- read.csv(
  file = "data-raw/commoditiesprices_data_dictionary.csv",
  stringsAsFactors = FALSE
)

usethis::use_data(commodities_prices, overwrite = TRUE)

usethis::use_data(commodities_prices_data_dictionary, overwrite = TRUE)
