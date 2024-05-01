#' Distribution patterns of the crab \emph{Ucides cordatus} at different spatial scales
#'
#' Distribution patterns of the mangrove crab \emph{Ucides cordatus} were assessed
#' at four levels of spatial hierarchy (from meters to tens of km) using a nested
#' ANOVA and variance components measures. The design included four spatial scales
#' of variation: regions (10s km), locations (km), sites (10s m), and quadrats (m).
#' Three regions tens of kilometers apart were sampled across the Paranagua Bay
#' salinity-energy gradient. Within each region, three locations
#' (1.5 - 3.5 km apart from each other) were randomly chosen. In each location,
#' three sites of 10 x 10 m were randomized, and within these, five quadrats
#' of 2 x 2 m were sampled (Sandrini-Neto and Lana, 2012).
#'
#' \itemize{
#'   \item Region: a random factor with three levels (\code{R1, R2, R3})
#'   \item Location: a random factor with three levels (\code{L1, L2, L3}) nested in Region
#'   \item Site: a random factor with three levels (\code{S1, S2, S3}) nested in Location
#'   \item Quadrat: sample size
#'   \item Crabs: response variable; number of crabs per quadrat
#'   \item Density: response variable; number of crabs per square meter
#' }
#'
#' @docType data
#' @keywords datasets
#' @name crabs
#' @usage data(crabs)
#' @format A data frame with 135 rows and 6 variables
#' @references Sandrini-Neto, L., Lana, P.C. 2012. Distribution patterns of the crab \emph{Ucides cordatus} (Brachyura, Ucididae) at different spatial scales in subtropical mangroves of Paranagua Bay (southern Brazil). Helgoland Marine Research 66, 167-174.
#' @examples
#' library(GAD)
#' data(crabs)
#' crabs
NULL
