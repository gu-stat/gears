#' Primary Commodities Price Data
#'
#' The most recent data from the International Monetary Fund on Primary
#' Commodities Prices. In particular, data on Beef, Swine (pork) meat,
#' Poultry (chicken) meat, Maize (corn), and Wheat.
#' Data is collected monthly, starting on January, 1980.
#'
#' @format A dataframe with 483 observations, and 6 variables.
#'    The variables are described in the
#'    \code{\link{commodities_prices_data_dictionary}}.
#' @source \url{https://www.imf.org/~/media/Files/Research/CommodityPrices/Monthly/Table3.ashx}
"commodities_prices"

#' Descriptions of the commoditiesprices data.
#'
#' A data dictonary describing the commoditiesprices data.
#'
#' @format A dataframe with 6 observations of fields, and 8 variables.
#'    \describe{\item{variable_name}{The variable name as it appears in the
#'    \code{\link{commodities_prices}} data frame.}
#'    \item{variable_description}{Detailed description of the contents of each
#'    variable.}
#'    \item{data_type}{Variable type: character (char) or numeric (num).}
#'    \item{units}{Monetary unit.}
#'    \item{frequency}{Time frequency of the data.}
#'    \item{source}{Name of the data source.}
#'    \item{source_link}{URL Link to the data source.}
#'    \item{access_on}{Date in which the data used to create the
#'    \code{\link{commodities_prices}} data frame was obtained.}}
"commodities_prices_data_dictionary"
