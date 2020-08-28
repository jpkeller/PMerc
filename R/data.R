#' Posterior parameters for PM-ALRI Curve
#'
#' A list containing the posterior samples and spline information needed
#' for computing values from the exposure-response curve for PM and ALRIs
#' in children.
#'
#' The parameters are output from a model fit via the `bercs` package.
#'
#' @format A list with two elements
#' \describe{
#'   \item{posterior_params}{List containing posterior samples for `beta` and `bS` parameters}
#'   \item{model_data}{List with minimal model information needed for computing ERC.}
#'   ...
#' }
#' @source (To be added)
"nepal_pm_alri"
