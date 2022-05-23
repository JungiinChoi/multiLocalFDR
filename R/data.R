#' @importFrom tibble tibble
NULL

#' P-values for Pathways.
#'
#' The list of significantly upregulated genes in three published studies:
#' (1) peripheral leukocytes(L), (2) orbital inflammatory disease(O), and
#' (3) sinus brushings(S) compared to healthy controls.
#'
#' @format A data frame with four variables: \code{Gene.Set}, \code{L}
#' (peripheral leukocytes), \code{O} (orbital inflammatory), and
#' \code{S} (sinus brushings).
"pathways"

#' Velocity of Galaxies
#'
#' Radial velocities of globular clusters of M104 with average v=1121 km/s,
#' and Milky Way stars which have much lower radial velocities.
#'
#'#' @format A data frame with two variables:
#' \describe{
#' \item{\code{radius}}{Angular measurement of each globular cluster and star (arcmin)}
#' \item{\code{velocity}}{Radial velocities (km/s)}
#' }
"galaxy"

#' P-values for Microarray Data
#'
#' List of p-values for three different experiments.
#'
#' @format A data frame with three variables: \code{expm1}, \code{expm2}, and
#' \code{expm3}.
#'
"microarrays"

