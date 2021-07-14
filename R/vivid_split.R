
#' VIVID method when the data must be split into multiple groups.
#'
#' @param x input matrix, rows represent observations and columns represent
#' features. column names and row names are recommended.
#' @param y response variable, restricted to two classes.
#' @param bootstraps integer, the number of bootstraps to run.
#' @param cores integer, number of cores to use for parallel processing. VIVID
#' makes uses of the "future" and "furrr" packages and will set the plan to
#' "multiprocess". The cores parameter is directly fed to workers arguement of
#' future::plan.
#' @param seed integer, passed to set.seed function, a full L'Ecuyer-CMRG RNG seed (vector of 7 integers).
#' @param minSize minimum cluster size to consider when clustering.
#' @param lambda Used in cv.glmnet
#' @param compareMethod AIC,  AICc, BIC or EBIC.
#' @param gamma gamma value for EBIC.
#' @param groups Integer value for the total number of groups.
#' @param disjoint TRUE or FALSE value on whether the groups are disjoint.
#' @param repFeatures Vector of features that are repeated in each group.
#' @importFrom furrr future_options
#' @importFrom parallel detectCores
#'
#' @return
#' @export
#'
#' @examples
#' library('ropls')
#' data("sacurine") #Load sacurine dataset from the 'ropls' package
#'
#' dat = sacurine$dataMatrix
#' outcomes = sacurine$sampleMetadata$gender
#' 
#' vividResults = vivid_split(x = dat,
#'       y = outcomes,
#'       groups = 2)
#'       
#' vividResults$optFeatures

vivid_split = function(x,
                       y,
                       bootstraps = 100,
                       cores = parallel::detectCores() - 1,
                       seed = 1234567,
                       minSize = 2,
                       lambda = 'lambda.1se',
                       compareMethod = 'BIC',
                       gamma = 1,
                       groups,
                       disjoint = TRUE,
                       repFeatures) {
  # string to factor function
  if (base::is.vector(y)) {
    y = base::factor(y)
  }

  nClass = length(levels(y))

  if (nClass != 2) {
    stop(
      base::paste(
        "The response vector contains",
        nClass,
        "classes and should only contain 2"
      )
    )
  }

  # params
  p = base::NCOL(x)

  base::RNGkind("L'Ecuyer-CMRG")
  base::set.seed(seed)



  output = base::rep(x = base::list(NA_real_),
                     times = groups + 1)
  featureTrack = base::rep(base::list(NA_real_),
                           times = groups + 1)

  if (disjoint == TRUE) {
    groupSize = base::floor(p / groups)
    remainder = p - groupSize * groups
    if (remainder == 0) {
      allocation = base::c(base::rep(x = 1:groups,
                                     times = groupSize))
    }

    if (remainder != 0) {
      allocation = base::c(base::rep(x = 1:groups, 
                                     times = groupSize), 1:remainder)
    }

    allocation = base::sample(x = allocation,
                               size = p,
                               replace = FALSE)

    for (i in 1:groups) {
      xNew = x[, allocation == i]
      nVar = base::ncol(xNew)
      output[[i]] = vivid(
        x = xNew,
        y = y,
        bootstraps = bootstraps,
        cores = cores,
        seed = seed,
        minSize = minSize,
        lambda = lambda,
        compareMethod = 'AIC'
      )[-2]
      featureTrack[[i]] = base::which(allocation == i)
    }

  }

  if (disjoint == FALSE) {
    ncol = p - base::length(repFeatures)
    groupSize = base::floor(ncol / groups)
    remainder = ncol - groupSize * groups
    allocation = base::c(base::rep(x = 1:groups,
                                   times = groupSize),
                         1:remainder)
    allocation = base::sample(x = allocation,
                              size = ncol,
                              replace = FALSE)
    xFixed = x[, repFeatures]
    xChange = x[,-repFeatures]

    for (i in 1:groups) {
      xNew = base::cbind(xFixed, xChange[, allocation == i])
      nVar = base::ncol(xNew)
      output[[i]] = vivid(
        x = xNew,
        y = y,
        bootstraps = bootstraps,
        cores = cores,
        seed = seed,
        minSize = minSize,
        lambda = lambda,
        compareMethod = 'AIC'
      )[-2]
      featureTrack[[i]] = base::c(repFeatures, base::which(allocation == i))
    }

  }

  finalPool = base::matrix(data = NA_real_,
                           nrow = groups,
                           ncol = p)

  for (j in 1:groups) {
    whichVar = base::unlist(output[[j]]$optModel)
    finalPool = base::c(finalPool, featureTrack[[j]][whichVar])
  }

  finalPool = base::sort(base::unique(finalPool))
  xFinal = x[, finalPool]
  featureTrack[[(groups + 1)]] = finalPool

  output[[(groups + 1)]] = vivid(
    x = xFinal,
    y = y,
    bootstraps = bootstraps,
    cores = cores,
    seed = seed,
    minSize = minSize,
    lambda = lambda,
    compareMethod = compareMethod,
    gamma = gamma
  )

  output$vividSplit = TRUE

  return(output)

}
