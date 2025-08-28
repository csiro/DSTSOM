#' Convert ordinal data frame for binary data analysis
#'
#' @param x dataframe with each row corresponding to a numeric ordinal valued
#' observation with name \code{yname} with labels 1, 2, ...
#' @param yname character, name of ordinal response variable on positive integers
#' @return an expanded dataframe with additional rows for implied zeroes in
#' binary likelihood. The new response variable will be named \code{y_binary}.
#' The cut levels are appended with names \code{cut_1} ... \code{cut_k},
#' where \eqn{k} is the maximum observed category.
#' Additionally, a column \code{index} is appended that records the
#' corresponding index of the original ordinal observation vector.
#' @export
dfBinary <- function(x, yname) {

  if (!is.data.frame(x)) stop("x should be a data.frame")
  if (!yname%in%names(x)) stop("yname should be name of variable in x")

  y <- x[[yname]]
  if (!is.numeric(y) && !is.integer(y))
    stop("response variable should be numeric")
  if (any(is.na(y))) stop("no missing values allowed")
  if (any(y < 1)) stop("response variable should be positive")
  if (any(!is.wholenumber(y))) stop("response variable should be integer valued")

  # max cutpoint level
  q <- max(y)

  # Loop by the ordinal levels
  indx <- NA
  y_binary <- NA
  cuts <- matrix(NA, nrow=1, ncol=q+1)

  for (ii in 1:q) {
    aa <- which(y==ii)
    indx <- c(indx, rep(aa, each=ii))
    y_binary <- c(y_binary, rep(c(rep(0,ii-1),1),length(aa)))

    gg <- matrix(0, nrow=length(aa)*ii, ncol=q+1)
    for (jj in 1:ii) {
      gg[ jj + (0:(length(aa)-1))*ii,jj] <- 1
    }
    gg[,q+1] <- rep(aa, each=ii)
    cuts <- rbind(cuts,gg)

    rm(aa, gg)
  }
  indx <- indx[-1]
  y_binary <- y_binary[-1]
  cuts <- cuts[-1,]

  df <- x[indx,]
  df$y_binary <- y_binary

  # cut levels
  colnames(cuts) <- c(paste("cut", 1:q, sep = "_"), "index")
  df <- data.frame(
    df,
    cuts,
    stringsAsFactors = FALSE
  )

  df <- df[df[ , paste("cut", q, sep ="_")] != 1, ]
  df <- df[ , -which(names(df) ==  paste("cut", q, sep ="_"))]
  df <- df[order(df$index, as.numeric(row.names(df))), ]

  df

}


#' Verify if the predicted spatial random effects for the input ordinal data
#' points at the same location and time point are the same.
#'
#' @param inlaObj the inla output object of a spatial or spatial-temporal
#' ordinal model.
#'
#' @param Aproj the projection matrix A used to fit the spatial or
#' spatial-temporal ordinal model which produces the inlaObj.
#'
#' @param DataLocsTime the matrix of locations and/or time in the input ordinal
#' data fitted to a spatial or spatial-temporal model. DataLocsTime should
#' be a matrix of two columns representing locations of the data fitted to the
#' spatial model. Otherwise, DataLocsTime should be a matrix of three columns
#' representing locations and time of the data fitted to the spatial-temporal
#' model.
#'
#' @param spdeName the name of spatial random effect
#'
#' @return a logic value. If TRUE, the predicted spatial random effects are the
#' same; if FALSE, they are different.
#' @export
is_SpatialRandom_consistent <- function(inlaObj=NULL, Aproj=NULL,
                                   DataLocsTime=NULL, spdeName=NULL) {

  if (any(c(is.null(inlaObj),is.null(Aproj),is.null(DataLocsTime),is.null(spdeName)))) {
    stop("Cannot find inlaObj, Aproj, DataLocsTime, or spdeName!")
  }

  # check if nrow(Aproj) = nrow(DataLocsTime)
  if (nrow(Aproj) != nrow(DataLocsTime)) {
    stop("The number of locations and the number of the Aproj rows are different!")
  }

  # check if the given spdeName is correct
  if (!(spdeName %in% names(inlaObj$summary.random))) {
    stop("Cannot find the given spdeName!")
  }

  # random effect of knot
  rand_knot <- matrix(inlaObj$summary.random[[spdeName]]$mean, ncol=1)

  # check if ncol(Aproj) = the number of knots
  if (ncol(Aproj) != nrow(rand_knot)) {
    stop("The number of knots and the number of the Aproj columns are different!")
  }

  # Spatial effects at the data locations
  aa <- Aproj %*% rand_knot

  temp00 <- data.frame(xyt=apply(DataLocsTime,1, function(x) paste(x,collapse=",")),
                    spEffect=as.character(aa@x))

  bb <- with(temp00, tapply(spEffect, list(xyt), function(x) length(unique(x))))

  return(ifelse(sum(bb>1)==0, TRUE, FALSE))

}


#' Predict marginal probabilities of ordinal levels in space-time.
#'
#' @param fm_inla spatio-temporal ordinal model fit by INLA. The model fixed
#' effects include cut levels labelled by cut_1, cut_2, ... cut_k - 1 for k
#' ordinal levels that correspond to first k - 1 columns in model matrix. The
#' model formula should include no intercept, e.g., \code{y ~ -1 + cut_1 +
#' cut_2 + ...}. The name of first argument in inla.spde.make.index is "s",
#' e.g., \code{inla.spde.make.index('s', n.spde, n.group)}. Please refer to
#' paper and worked examples in package vignette for further details on INLA
#' model formulation.
#' @param nSample number of Monte Carlo realisations sampled from ordinal model
#' by INLA
#' @param arrPred array of sites (first dimension) by fixed covariates less cut
#' levels (second dimension) by time levels (e.g., years, third dimension)
#' conformal to model matrix so that covariate names  match, in order, "fixed"
#' covariates from INLA model fit excluding cut levels.
#' @param A_s_Pred prediction weight matrix output from
#' \code{INLA::inla.spde.make.A} based on spatial mesh for prediction
#' @param INLAseed seed passed to INLA posterior sampling function
#' @return List of 1) MCsummary object that summarises predictive posterior
#' quantiles by site, ordinal level and year, and 2) a summary matrix maxDev.
#' The latter is a numerical check to show that max absolute deviation from one
#' after summation of predicted marginal probabilities by ordinal level across
#' all realisations of each site and year combination.
#' @export
predict_STordinal <- function(fm_inla = NULL,
                              nSample = 50,
                              arrPred = NULL,
                              A_s_Pred = NULL,
                              INLAseed = 0L) {

  # argument checking ----------------------------------------------------------

  stopifnot(length(nSample) == 1 && is.numeric(nSample) && nSample > 0)
  if (is.null(fm_inla)) stop("Cannot find fm_inla!")
  if (is.null(arrPred)) stop("Cannot find arrPred!")
  if (is.null(A_s_Pred)) stop("Cannot find A_s_Pred!")

  # collect crucial info
  FixedName <- fm_inla$names.fixed
  if (FixedName[1] != "cut_1")
    stop("expecting first covariates in INLA model to be cut_1, cut_2 etc ...")

  cutnames <- fm_inla$names.fixed[grepl("cut", fm_inla$names.fixed)]
  if (!all(cutnames == paste("cut", 1:length(cutnames), sep = "_")))
    stop("expecting first covariates in INLA model to be cut_1, cut_2 etc ...")
  if (!all(fm_inla$names.fixed[1:length(cutnames)] == cutnames))
    stop("expecting first covariates in INLA model to be cut_1, cut_2 etc ...")

  # check any fixed covariate names
  noncutnames <- fm_inla$names.fixed[(length(cutnames) + 1):length(fm_inla$names.fixed)]
  if (!all(dimnames(arrPred)[[2]] == noncutnames))
    stop("discrepancy between array of predictors and INLA model object covariates")

  # prediction function is optimised for homogeneous covariate structure
  if ("cut"%in%noncutnames)
    stop("prediction function optimised for no interactions with cut levels")

  ngroup <- fm_inla[["size.random"]][[1]][["ngroup"]]
  if (dim(arrPred)[3] != ngroup)
    stop("temporal dimension in arrPred should match n.group in INLA fit")

  # INLA prediction ------------------------------------------------------------

  STindex <- as.data.frame(
    INLA::inla.spde.make.index('s', n.spde = ncol(A_s_Pred), n.group = dim(arrPred)[3])
  )

  # need config = TRUE in control.compute of inla and require the sn package installed
  Samples <- INLA::inla.posterior.sample(n = nSample,
                                         result = fm_inla,
                                         intern = FALSE,
                                         use.improved.mean = TRUE,
                                         skew.corr = TRUE,
                                         seed = INLAseed)

  # grab MC samples
  ## All fixed (row: fixed, column: MC samples)
  FixSamples <- INLA::inla.posterior.sample.eval(FixedName,
                                                 Samples,
                                                 return.matrix = TRUE)
  row.names(FixSamples) <- FixedName
  tCutSamples <- t(FixSamples[cutnames, ])
  tNonCutSamples <- t(FixSamples[noncutnames, ])

  rm(FixSamples)
  invisible(gc())

  ## spatial-temporal effect in the feature space
  # knots * years
  STeffSamples <- INLA::inla.posterior.sample.eval("s", Samples,
                                                   return.matrix = TRUE)
  rm(Samples)
  invisible(gc())

  # collect spatial temporal effects by step
  tSTeff.ls <- vector("list", length = ngroup)
  for (tt in 1:ngroup) {
    tSTeff.ls[[tt]] <- t(STeffSamples[STindex$s.group == tt, ])
  }
  rm(STeffSamples)
  invisible(gc())

  # preallocate storage for results --------------------------------------------

  ## Collect the MC outputs (year and spatial): tail quantiles and median
  MCsummary  <- array(data = NA,
                      dim = c(dim(arrPred)[1],
                              length(cutnames) + 1,
                              3,
                              ngroup),
                      dimnames = list(1:dim(arrPred)[1],
                                      c(cutnames, "k"), # append last category (highest ordinal level)
                                      c("quant05", "median", "quant95"),
                                      paste0("year", 1:ngroup)))

  # maximum abs deviation from one by sum of marginal probs over categories
  # these should always equal one
  # max abs deviation recorded in each year for each site
  maxDev <- matrix(NA, nrow = dim(arrPred)[1], ncol = ngroup,
                   dimnames = list(
                     1:dim(arrPred)[1],
                     paste0("year", 1:ngroup)
                   ))

  # Posterior predictive marginal probabilities of ordinal levels --------------

  # by time level (year)
  for (tt in 1:ngroup) {

    # predictions for fixed effects for all time levels
    # sims by sites
    fpreds <- tcrossprod(tNonCutSamples, arrPred[ , , tt])

    # random effects
    # sims by sites
    w <- as.matrix(tcrossprod(tSTeff.ls[[tt]], A_s_Pred))

    # linear predictor for all cut levels
    # sims by cut levels
    Eta <- fpreds + w

    for (l in 1:length(cutnames)) {

      # inverse cloglog transform to sequential probabilities
      # after adding cut level to linear predictor
      # sims by site for each cut level in time level
      seqProbs <- -expm1(-exp(Eta + tCutSamples[ , l]))

      # across Monte Carlo realisations:
      # marginal probabilities that correspond to
      # successive cut levels depend on products
      # of conditional and cumulative probabilities
      # store cumulative probs for use at next cut level

      if (l == 1) {
        MarginProbs <- seqProbs
        CumulProbs <- 1 - seqProbs
        sumMP <- MarginProbs
      } else {
        MarginProbs <- seqProbs*CumulProbs
        CumulProbs <- CumulProbs*(1 - seqProbs)
        sumMP <- sumMP + MarginProbs
      }

      MCsummary[ ,
                 cutnames[l],
                 c("quant05", "median", "quant95"),
                 paste0("year", tt)] <- matrixStats::colQuantiles(
                   MarginProbs, probs = c(0.05, 0.50, 0.95)
                 )

    }

    # last category
    MCsummary[ ,
               "k",
               c("quant05", "median", "quant95"),
               paste0("year", tt)] <- matrixStats::colQuantiles(
                 CumulProbs, probs = c(0.05, 0.50, 0.95)
               )

    # should sum to one for each realised site x time combo
    # record max absolute deviation from one
    maxDev[ , tt] <- matrixStats::colMaxs(abs(sumMP + CumulProbs - 1))

  }

  # return predictions
  list(MCsummary = MCsummary, maxDev = maxDev)

}

