#' Subset S4? object according to expression
#' 
#' @encoding UTF-8
#' 
#' @param object.expr Expression object
#' @param inputAssay Options: c("counts", "tpm"). Default sets to "counts".
#' @param minExpr Minimum expression that would like to include in the subset condition. Default sets to 1.
#' @param minSamples Minimum library that would like to include in the subset condition. Default sets to 2.
#' @param across Condition that would like to subset including 
#' c("allSamples", "eachClass", "sumAllSamples", "sumInEachClass"). 
#' Default sets to "allSamples".
#' @param groupBy Column name of metadata. this will be set only for "eachClass", "sumInEachClass". 
#' Default sets to FALSE.
#' 
#' @return new expression object according to the preferred condition
#' 
#' @examples
#' FANTOM setting (Is it the FANTOM setting?, Robins setting) 
#' subset.expr <- subsetter(cage.expr, minExpr=2, minSamples=2, across="allSamples")
#' More than 1 tags in at least 2 libraies for at least 1 condition
#' subset.expr <- subsetter(cage.expr, minExpr=2, minSamples=2, across="eachClass", groupBy="Class")
#' More than 1 tags in at least 1 condition
#' subset.expr <- subsetter(cage.expr, minExpr=2, across="sumInEachClass", groupBy="Class")
#' 
#' @export

subsetter <- function(object.expr, inputAssay="counts",
                      minExpr=1, minSamples=2,
                      across="allSamples", groupBy=FALSE) {

  # across = c("allSamples", "eachClass", "sumAllSamples", "sumInEachClass")

  # "sumInEachClass" and "eachClass" : one of these classes passes the minimum setting

  # expressed:         expr >= minExpr, TRUE/FALSE assigned
  # expressed.byClass: rowSum data according to class indicated in design.matrix
  # expressed.all:     rowSum across all samples
  # keep:              according to the each << across >>

  # no support column as in subsetBySupport (add?)

  ch = FALSE

  # Pre-checks
  assertthat::assert_that(is.numeric(minExpr))
  assertthat::assert_that(is.numeric(minSamples))
  

  if (across=="allSamples")  {
    message("Subsetting data across all samples ...")
    expressed <- assay(object.expr, inputAssay) >= minExpr
    keep <- rowSums(expressed) >= minSamples
    ch = TRUE
  }

  else if (across=="eachClass") {
    message("Subsetting data across each class ...")
    expressed <- assay(object.expr, inputAssay) >= minExpr
    expressed.byClass <- sapply(unique(colData(object.expr)[[groupBy]]),
                                function(n) rowSums(expressed[, colData(object.expr)[[groupBy]]==n]))
    keep <- rowMax(expressed.byClass) >= minSamples
    ch = TRUE
  }

  else if (across=="sumAllSamples") {
    message("Subsetting data on summation of all samples ...")
    message("try out")
    expressed.all <- rowSums(assay(object.expr, inputAssay))
    keep <- expressed.all >= minExpr
    ch = TRUE
  }

  else if (across=="sumInEachClass") {
    message("Subsetting data on summation of samples of each class ...")
    expressed.byClass <- sapply(unique(colData(object.expr)[[groupBy]]),
                                function(n) rowSums(assay(object.expr, inputAssay)[, colData(object.expr)[[groupBy]]==n]))
    keep <- rowMax(expressed.byClass) >= minExpr
    ch = TRUE
  }

  else {
    message("Wrong 'across' is given. Options: 'allSamples', 'eachClass', 'sumAllSamples', or 'sumInEachClass'")
  }

  if (ch == TRUE) {
    new.object.expr <- object.expr[keep, ]
    before <- nrow(object.expr)
    after <- nrow(new.object.expr)
    message(after, " regions (", round(after/before * 100, digits = 1), "%) out of ", before, " regions were kept.")
  }

  return(new.object.expr)

}
