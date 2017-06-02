#' Accessor for the 'surr' slot of a 'svplsSurr' object
#'
#' @usage
#' \S4method{surr}{svplsSurr}(object)
#'
#' @docType methods
#' @name surr
#' @rdname surr
#' @aliases surr surr,svplsSurr-method
#' @param object a \code{svplsSurr} object
#'
#' @examples
#' data(sim.dat)
#' group = as.factor(c(rep(1, 10), rep(-1, 10)))
#' sv = svplsSurr(sim.dat, group, surr.select = "automatic")
#' surr(sv)
surr.svplsSurr <- function(object) object@surr

#' @export
setMethod("surr", signature(object = "svplsSurr"), surr.svplsSurr)

#' Accessor for the 'prop.vars' slot of a 'svplsSurr' object
#'
#' @usage
#' \S4method{prop.vars}{svplsSurr}(object)
#'
#' @docType methods
#' @name prop.vars
#' @rdname prop.vars
#' @aliases prop.vars prop.vars,svplsSurr-method
#' @param object a \code{svplsSurr} object
#'
#' @examples
#' data(sim.dat)
#' group = as.factor(c(rep(1, 10), rep(-1, 10)))
#' sv = svplsSurr(sim.dat, group, surr.select = "automatic")
#' prop.vars(sv)
pvars.svplsSurr <- function(object) object@prop.vars

#' @export
setMethod("prop.vars", signature(object = "svplsSurr"), pvars.svplsSurr)

## new S4 summary function for 'svplsSurr' objects
setMethod("summary", "svplsSurr", function(object){
  cat("Significant surrogate variables: \n")
  print(object@surr)
  cat("Proportion of total variance explained: \n")
  print(object@prop.vars)
})

## new S4 print function for 'svplsSurr' objects
setMethod("print", "svplsSurr",  function(x){
  cat("Significant surrogate variables: \n")
  print(x@surr)
})

#' Accessor for the 'pvs.unadj' slot of a 'svplsTest' object
#'
#' @usage
#' \S4method{pvs.unadj}{svplsTest}(object)
#'
#' @docType methods
#' @name pvs.unadj
#' @rdname pvs.unadj
#' @aliases pvs.unadj pvs.unadj,svplsTest-method
#' @param object a \code{svplsTest} object
#'
#' @examples
#' data(sim.dat)
#' group = as.factor(c(rep(1, 10), rep(-1, 10)))
#' sv = svplsSurr(sim.dat, group, surr.select = "automatic")
#' surr = surr(sv)
#' fit = svplsTest(dat = sim.dat, group = group, surr = surr, normalization = "TMM", test = "t-test")
#' pvs.unadj(fit)
pvs.unadj.svplsTest <- function(object) object@pvs.unadj

#' @export
setMethod("pvs.unadj", signature(object = "svplsTest"), pvs.unadj.svplsTest)

#' Accessor for the 'pvs.adj' slot of a 'svplsTest' object
#'
#' @usage
#' \S4method{pvs.adj}{svplsTest}(object)
#'
#' @docType methods
#' @name pvs.adj
#' @rdname pvs.adj
#' @aliases pvs.adj pvs.adj,svplsTest-method
#' @param object a \code{svplsTest} object
#'
#' @examples
#' data(sim.dat)
#' group = as.factor(c(rep(1, 10), rep(-1, 10)))
#' sv = svplsSurr(sim.dat, group, surr.select = "automatic")
#' surr = surr(sv)
#' fit = svplsTest(dat = sim.dat, group = group, surr = surr, normalization = "TMM", test = "t-test")
#' pvs.adj(fit)
pvs.adj.svplsTest <- function(object) object@pvs.adj

#' @export
setMethod("pvs.adj", signature(object = "svplsTest"), pvs.adj.svplsTest)

#' Accessor for the 'sig.features' slot of a 'svplsTest' object
#'
#' @usage
#' \S4method{sig.features}{svplsTest}(object)
#'
#' @docType methods
#' @name sig.features
#' @rdname sig.features
#' @aliases sig.features sig.features,svplsTest-method
#' @param object a \code{svplsTest} object
#'
#' @examples
#' data(sim.dat)
#' group = as.factor(c(rep(1, 10), rep(-1, 10)))
#' sv = svplsSurr(sim.dat, group, surr.select = "automatic")
#' surr = surr(sv)
#' fit = svplsTest(dat = sim.dat, group = group, surr = surr, normalization = "TMM", test = "t-test")
#' sig.features(fit)
sig.features.svplsTest <- function(object) object@sig.features

#' @export
setMethod("sig.features", signature(object = "svplsTest"), sig.features.svplsTest)

## new S4 summary function for 'svplsTest' objects
setMethod("summary", "svplsTest", function(object){
  cat("Unadjusted pvalues: \n")
  print(object@pvs.unadj)
  cat("FDR corrected pvalues: \n")
  print(object@pvs.adj)
})

## new S4 print function for 'svplsTest' objects
setMethod("print", "svplsTest", function(x){
  cat("The significantly differentially expressed features are: \n")
  print(x@sig.features)
})
