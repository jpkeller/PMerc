##' @title Extract fitted Exposure-Response Curve (ERC)
##' @description Computes the exposure-response curve from a fitted outcome model
##' @param exprange range of exposure values over which to compute the curve
##' @param expsequence sequence of exposure values at which the curve should be evaluated. For plotting, it is preferable to use \code{exprange}, but for specific known exposure values use \code{expsequence}.
##' @param ciband width of credible interval band
##' @param inclInterceptUncertainty Should intercept uncertainty be included uncertainty estimates? See details for more information.
##' @param inclIntercept Should the intercept term be included in the curve?
##' @param intercept_prop Proportions used in calculating the "average" intercept. Defaults to equal proportions for each study (\code{"equal"}) and can be set to be proportional to the number of observations in each study (\code{"obs"}). A vector of proportions can also be given.
##' @param study Return the curve for a specific study, or the
##' @param beta_post vector or matrix of posterior samples for 'beta' parameter. Only needed if \code{stanfit} not provided.
##' @param bs_post vector or matrix of posterior samples for 'bS' parameter. Only needed if \code{stanfit} not provided.
##' @param nS Number of studies. Only needed if \code{standata} not provided.
##' @param Mx Spline matrix for exposure. Only needed if \code{standata} not provided.
##' @param Mx_attributes Attributes for the spline matrix. Only needed if \code{standata} not provided.
##' @param xdf Degrees of freedom in exposure splines
##' @param ... Passed to \code{\link{seq}} to control sequence of exposure values.
##'
##' @details These functions are derived from the functions of similar name in the `bercs` package (e.g. `compute_ERC2()`). They are included here in `PMerc` to facilitate calculating ERC values without requiring the full `bercs` package (which requires compiled C code for the STAN model objects). This package is recommended only for calculating specific risk measures from the accompanying data files. If developing your own models, it is recommended to use the functions in the `bercs` package and not these. To avoid namespace issues, these functions have `2` appended in their name.
##'
##' This function creates a data frame containing the values of the exposure-response curve over a given range of values. The output is designed for easy plotting.
##' Currently, the fitted curve is plotted based upon the posterior means of the parameters, with the uncertainty bands based upon quantiles.
##'
##' Uncertainty from the intercepts is included in the confidence bands by default, since this corresponds to a common interpretation of such intervals. This requires picking a single value for the intercept.
##' For models fit to multiple studies, a separate set of uncertainty will be created for each study. Additionally, a set of results corresponding to "average" intercept is created. The \code{intercept_prop} argument controls the relative contribution of the intercepts from each model.
##'
##' @seealso \code{\link[bercs]{compute_ERC}}
##' @export
#' @importFrom stats quantile
#' @importFrom utils getS3method
compute_ERC2 <- function (  exprange = c(0, 100),
                            expsequence=NULL,
                            ref_exposure=NULL,
                            ciband = 0.95,
                            inclInterceptUncertainty=TRUE,
                            inclIntercept=FALSE,
                            intercept_prop=c("equal"),
                            study=NULL,
                            beta_post,
                            bs_post,
                            nS,
                            Mx,
                            Mx_attributes=attributes(Mx),
                            xdf=ncol(Mx),
                            ...)
{
    if (ciband < 0 || ciband > 1)
        stop("'ciband' must be between 0 and 1.")
    # If multiple studies, will add column with "average" intercept
    nSout <- ifelse(nS>1 & length(dim(beta_post))==2, nS+1, nS)
    if (is.null(expsequence)){
      expsequence <- seq(exprange[1], exprange[2],...)
    } else {
      if (!is.null(ref_exposure)){
        if (!ref_exposure %in% expsequence) {
          expsequence <- c(ref_exposure, expsequence)
        }
      }
    }

    if (nSout > nS){
        # Add column that has "average" intercept, if there are multiple studies
        if (intercept_prop=="equal"){
            intercept_prop <- rep(1/nS, length=nS)
        } else if (!is.numeric(intercept_prop)){
            stop("'intercept_prop' must be 'equal', 'open' or a numeric vector")
        }
        if (sum(intercept_prop)!=1) stop("The values in 'intercept_prop' should equal 1.")
        if (any(intercept_prop < 0)) stop("The values in 'intercept_prop' should be non-negative.")


        bs_post <- cbind(bs_post, bs_post %*% intercept_prop)
    }
    if (length(dim(beta_post))==2) {
        beta_post2 <- array(dim=c(nrow(beta_post),
                                  nSout,
                                  ncol(beta_post)))
        for (j in 1:dim(beta_post2)[2]){
            beta_post2[, j, ] <-beta_post
        }
        beta_post <- beta_post2
    }

    if (xdf == 1) {
        if (!is.null(Mx_attributes$`scaled:center`)) {
            exposure_seq_scaled <- (expsequence - Mx_attributes$`scaled:center`)/Mx_attributes$`scaled:scale`
        }
        else {
            exposure_seq_scaled <- expsequence
        }
    }
    else {
        Mxtemp <- Mx
        attributes(Mxtemp) <- Mx_attributes
        predfn <- utils::getS3method(f = "predict",
                                     class(Mxtemp)[class(Mxtemp) !="matrix"][1])
        exposure_seq_scaled <- predfn(Mxtemp, newx = expsequence)
    }
    fitted_seq <- array(dim=c(length(expsequence),
                              dim(beta_post)[1], # B from stan fit
                              nSout))
    for (i in 1:nSout){
        fitted_seq[, , i] <- exposure_seq_scaled %*% t(beta_post[, i, ])
    }

    if (inclIntercept & !inclInterceptUncertainty) {
        inclInterceptUncertainty <- TRUE
        warning("inclIntercept set to TRUE. This requires inclInterceptUncertainty to be TRUE.")
    }
    if (inclInterceptUncertainty) {
        if (!inclIntercept){
            bs_post <- sweep(bs_post, 2, colMeans(bs_post), FUN = "-")
        }
        bs_post <- rep(1, length(expsequence)) %o% bs_post

        fitted_seq <- fitted_seq + bs_post
    }
    fitted_seq_mean <- apply(fitted_seq, c(1, 3), mean)
    fitted_seq_low <- apply(fitted_seq, c(1,3), stats::quantile,
                                  probs = (1 - ciband)/2)
    fitted_seq_high <- apply(fitted_seq, c(1,3), stats::quantile,
                                   probs = 0.5 + ciband/2)

    obj <- list(mean = fitted_seq_mean,
                low = fitted_seq_low,
                high = fitted_seq_high,
                exposure = expsequence)

    dflist <- vector("list", nSout)
    for (j in 1:nSout){
      dflist[[j]] <- data.frame(mean=obj$mean[, j],
                                low=obj$low[, j],
                                high=obj$high[, j],
                                exposure=obj$exposure)
    }
    dflist <- dflist[[length(dflist)]] # only the combined intercept curve
    if (!is.null(ref_exposure)){
      dflist <- center_ERC2(dflist,
                         ref_exposure=ref_exposure)
    }
    dflist
}

##' @rdname compute_ERC2
##' @param obj Data frame, or list of data frames, containing \code{exposure}, \code{mean}, \code{low}, and \code{high}. Typically generated from \code{\link{compute_ERC2}}.
##' @param incS If model has a single curve but with different selections of intercept uncertainty, this indicates which choice of uncertainty to use. Defaults to the last column of the curve in \code{obj}, which is typically the averaged intercept uncertainty. If multiple curves are fit, this selects which curve(s) is plotted.
##' @param expERC Should the fitted curve be exponentiated (TRUE) or not (FALSE).
##' @param ylab String providing y-axis label.
##' @param xlab String providing x-axis label.
##' @param ribbon Should the uncertainty be represented as a filled ribbon (TRUE) or lines without fill (FALSE).
##' @export
##' @import ggplot2
plot_ERC2 <- function (obj,
                             incS=NULL,
                             expERC=TRUE,
                             ylab = "Relative Risk",
                             xlab = "Exposure",
                             ribbon=FALSE)
{
    if (!is.data.frame(obj)){
      nS <- length(obj)
    } else {
      nS <- 1
    }
    fulldf <- do.call(rbind, obj)
    if (expERC) {
        fulldf$mean <- exp(fulldf$mean)
        fulldf$low <- exp(fulldf$low)
        fulldf$high <- exp(fulldf$high)
    }

    if(is.null(incS)){
        incS <- nS
    }
    fulldf <- subset(fulldf, study %in% incS)
    fulldf$study <- factor(fulldf$study)
    g <- ggplot(fulldf) + theme_bw()
    if (ribbon){
        g <- g + geom_ribbon(aes(x = .data$exposure,
                                 ymin = .data$low,
                                 ymax = .data$high,
                                 group=.data$study),
                             fill = "grey80")
    } else {
        g <- g+   geom_line(aes(x=.data$exposure,
                                y=.data$low,
                                group = .data$study,
                                col=.data$study),
                            data=fulldf) +
            geom_line(aes(x=.data$exposure,
                          y=.data$high,
                          group = .data$study,
                          col=.data$study),
                      data=fulldf)
    }
      g   + geom_line(aes(x = .data$exposure,
                      y = .data$mean,
                      group=.data$study,
                      col=.data$study),
                  lwd = 1.5) +
        geom_hline(yintercept = 1,
                   lty = 2) + xlab(xlab) + ylab(ylab)

}


##' @rdname compute_ERC2
##' @param ref_exposure Exposure concentration at which the spline should be shifted to have value 0 on log scale (value 1 on exponentiated scale). The largest exposure value less than or equal to \code{ref_exposure} is used as the new reference concentration.
##' @export
center_ERC2 <- function(obj, ref_exposure=min(obj$exposure)){
  obj_new <- obj
  if (!is.data.frame(obj_new)){
    for (j in 1:length(obj_new)){
      ind_cut <- max(which(obj_new[[j]]$exposure <= ref_exposure))
      obj_new[[j]]$low <- obj_new[[j]]$low - obj_new[[j]]$mean[ind_cut]
      obj_new[[j]]$high <- obj_new[[j]]$high - obj_new[[j]]$mean[ind_cut]
      obj_new[[j]]$mean <- obj_new[[j]]$mean - obj_new[[j]]$mean[ind_cut] # do last
    }} else {
      ind_cut <- max(which(obj_new$exposure <= ref_exposure))
      obj_new$low <- obj_new$low - obj_new$mean[ind_cut]
      obj_new$high <- obj_new$high - obj_new$mean[ind_cut]
      obj_new$mean <- obj_new$mean - obj_new$mean[ind_cut] # do last
    }
  obj_new
}
