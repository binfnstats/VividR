#' Adjust variables in a VIVID object
#'
#' @param vividObj An object passed from the vivid() funciton.
#' @param minFinalFeatures Integer value for the minimum number of final features selected.
#'
#' @return A list of results:
#' \itemize{
#' \item{optFeatures}: A character vector containing all features in the final vivid model, adjusted for minFinalFeatures.
#' \item{n}: A count of the new number of optimum features.
#' }
#' @export
#'
#' @examples
#' library('ropls')
#' data("sacurine") #Load sacurine dataset from the 'ropls' package
#'
#' dat = sacurine$dataMatrix
#' outcomes = sacurine$sampleMetadata$gender
#'
#' vividResults = vivid(x = dat,
#'                      y = outcomes)
#'
#' vivid_adj(vividResults,
#' minFinalFeatures = 10)

vivid_adj = function(vividObj, minFinalFeatures) {
  features = vividObj$selection
  sizes = vividObj$sizes
  compareValues = vividObj$compareValues

  possModels = base::which(sizes >= minFinalFeatures)

  if (length(possModels) == 0) {
    stop("Size on minimum features exceeds that of possible feature groups found.")
  }

  optModel = features[base::which.min(compareValues[possModels]), ]

  optFeatures = base::names(optModel)[base::unlist(optModel)]
  newLength = base::length(optFeatures)

  return(list(
    optFeatures = optFeatures,
    n = newLength))
}
