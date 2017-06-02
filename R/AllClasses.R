#' @title svplsSurr
#' @exportClass svplsSurr
setClass("svplsSurr", representation(surr = "matrix", prop.vars = "numeric"))

#' @title svplsTest
#' @exportClass svplsTest
setClass("svplsTest", representation(pvs.unadj = "numeric", pvs.adj = "numeric", sig.features = "character"))
