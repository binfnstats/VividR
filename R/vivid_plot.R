
#' A funciton that creates graphics for the VIVID method
#'
#' @param vividObj An object passed from the vivid() funciton.
#' @param log A TRUE or FALSE variable that defines whether the values plotted are log-transformed.
#' @importFrom dendsort dendsort
#' @importFrom stringr str_trunc
#' @importFrom grDevices colorRampPalette
#'
#' @return A VIVID plot.
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
#' vivid_plot(vividObj = vividResults)

vivid_plot = function(vividObj, log = FALSE, topN = 0) {
  vividSplit  = vividObj$vividSplit
  if (vividSplit == TRUE) {
    vividObj = vividObj[[length(vividObj)-1]]
  }

  if (topN == 0){
    topN = length(vividObj$varClust$labels)
  }

  vividObj$varClust = dendsort::dendsort(d = vividObj$varClust)

  mat = vividObj$varMat

  vividObj$varClust$labels = stringr::str_trunc(
    string = vividObj$varClust$labels,
    width = 12,
    side = "right"
  )

  ddRow = stats::as.dendrogram(vividObj$varClust)
  rowOrd = stats::order.dendrogram(ddRow)
  ddCol = ddRow
  colOrd = rowOrd

  colnames(vividObj$varMat) = stringr::str_trunc(
    string = colnames(vividObj$varMat),
    width = 12,
    side = "right"
  )
  rownames(vividObj$varMat) = stringr::str_trunc(
    string = rownames(vividObj$varMat),
    width = 12,
    side = "right"
  )

  cc = grDevices::colorRampPalette(c("darkgreen", "green", "yellow", "red", "darkred"))
  lattice::trellis.par.set(regions = list(col = cc(100)))


  base::switch(
    2-log*1,
    "TRUE" = lattice::levelplot(xlab = "Features",
                                ylab = "Features",
                                main = "Mosaic Pairwise Variance Plot",
                                log(vividObj$varMat[rowOrd[1:topN], colOrd[1:topN]]+1),
                                col.regions = cc,
                                scales = list(cex = 0.5, x = list(rot =
                                                                    90)),
                                aspect = "fill",
                                colorkey = list(space = "left"),
    ),
    "FALSE" = lattice::levelplot(xlab = "Features",
                                 ylab = "Features",
                                 main = "Mosaic Pairwise Variance Plot",
                                 vividObj$varMat[rowOrd[1:topN], colOrd[1:topN]],
                                 col.regions = cc,
                                 scales = list(cex = 0.5, x = list(rot =
                                                                     90)),
                                 aspect = "fill",
                                 colorkey = list(space = "left")
    )
  )

}
