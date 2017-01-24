
#' @title svplsSurr
#'
#' @description This function extracts the surrogated estimates of the hidden variables in the 
#' data by using the partial least squares (PLS) algorithm 
#' on two multivariate random matrices. It provides the user with two options: 
#'
#' (1) \bold{Unsupervised SVAPLS}: Here a standard linear regression model is first used on
#' a transformed version of the expression count matrix to estimate the primary signals 
#' of differential expression for all the genes. The fitted model residuals and the 
#' transformed count matrix are then organized respectively into two multivariate matrices 
#' \code{E} and \code{Y}, in such a way that each column corresponds to a certain gene.  
#' \code{E} is then regressed on \code{Y} using a Non-linear partial least squares (NPLS) 
#' algorithm and the extracted scores in the column-space of \code{Y} are deemed as the surrogate 
#' variables.
#'
#' (2) \bold{Supervised SVAPLS}: In case information on a set of control genes (probes) is provided,
#' this function uses a Non-linear partial least squares (NPLS) algorithm to regress \code{Y} on a submatrix 
#' of \code{Y} (\code{Y.sub}) corresponding to the set of controls and scores in the column-space 
#' of \code{Y.sub} are considered as the surrogate variables.
#' 
#' The function then regresses the first eigenvector of the residual matrix \code{E} (for Unsupervised
#' SVAPLS or the control matrix \code{Y.sub} for Supervised SVAPLS) on these surrogate variables and 
#' tests them for statistical significance with a certain user-specified pvalue cutoff. The variables 
#' yielding a pvalue below the cutoff are returned.
#'
#' @param dat The original gene expression count matrix.
#' @param group a factor representing the sample indices belonging to the two 
#' different groups.
#' @param controls The set of control genes with no differential expression 
#' between the two groups (set to NULL by default).
#' @param phi The transforming function to be applied on the original gene 
#' expression count data (set to be log function with an offset \code{const}).
#' @param const The offset parameter for the transforming function \code{phi} 
#' (set to 1 by default).
#' @param pls.method The non-linear partial least squares method to be used. 
#' The different options available are: the classical orthogonal scores 
#' algorithm ("oscorespls, default), the kernel algorithm ("kernelpls") and 
#' wide kernel algorithm ("widekernelpls"). Using the "oscorespls" option is 
#' recommended for producing mutually orthogonal surrogate variables.
#' @param max.surrs The maximum number if surrogate variables to be extracted 
#' from the NPLS algorithm (set to 3 by default).
#' @param cutoff The user-specified pvalue cutoff for testing the significance
#' of the extracted surrogate variables (set to 1e-07 by default).
#' @param parallel Logical, indicating if the computations should be 
#' parallelized or not (set to \code{FALSE} by default).
#' @param num.cores The requested number of cores to be used in the parallel 
#' computations inside the function (used only when \code{parallel} is 
#' \code{TRUE}, \code{NULL} by default).
#' @param plot Logical, if \code{TRUE} a barplot of the variance proportions 
#' explained by the significant surrogate variables is returned (set to 
#' \code{FALSE} by default).
#'
#' @return surr A \code{data.frame} of the significant surrogate variables.
#' @return prop.vars A vector of the variance proportions explained by the
#' variables in \code{surr}.
#'
#' @examples
#' ##Loading the simulated dataset
#' data(sim.dat)
#'
#' ##Extracting the significant surrogate variables
#' group = as.factor(c(rep(1, 10), rep(-1, 10)))
#' sv <- svplsSurr(dat = sim.dat, group = group)
#' slotNames(sv)
#' head(surr(sv))             
#' head(prop.vars(sv)) 
#'
#' @rdname svplsSurr
#' @export
svplsSurr <-
function(dat, group, controls = NULL, phi = function(x) log(x + const), const = 1, 
            pls.method = "oscorespls", max.surrs = 3, cutoff = 10^-7, parallel = FALSE, 
            num.cores = NULL, plot = FALSE){

if (class(dat) == "matrix") data = dat
if (class(dat) == "SummarizedExperiment") data = assay(dat)
if (class(dat) == "DGEList") data = dat$counts
Y = phi(data)

###Fit an Initial Model to Y and organize the fitted model residuals into a matrix
if (parallel){
   e = unlist(mclapply(1:nrow(Y), function(u) lm(Y[u, ] ~ group)$resid, mc.cores = num.cores))
}
if (!parallel) e = unlist(lapply(1:nrow(Y), function(u) lm(Y[u, ] ~ group)$resid))

E = matrix(e, nrow(Y), ncol(Y), byrow = TRUE)

###Extract the signatures of hidden variation using Partial Least Squares 
if (!is.null(controls)) pls.fit = mvr(t(Y) ~ t(Y[controls, ]), ncomp = max.surrs, method = pls.method)   
if (is.null(controls)) pls.fit = mvr(t(E) ~ t(Y), ncomp = max.surrs, method = pls.method)  

sc = scores(pls.fit)
surr.effect = paste("+", paste("sc[, ", paste(as.character(1:ncol(sc)), "]", sep = ""), sep = "", collapse = "+"), sep = "")

if (!is.null(controls)) pvals = anova(lm(as.formula(paste("svd(Y[controls, ])$v[, 1] ~", surr.effect, sep = " "))))[, 5]
if (is.null(controls)) pvals = anova(lm(as.formula(paste("svd(E)$v[, 1] ~", surr.effect, sep = " "))))[, 5]

index = which(pvals < cutoff)      
if (length(index) == 0) stop("No Significant Surrogate Variables")      
if (length(index) > 0){
   	surr = sc[, index]

	if (length(index) == 1) vars = var(surr)
	if (length(index) > 1)  vars = apply(surr, 2, var)
	prop.vars = vars/sum(vars)

	if (plot == TRUE){
  	   surr.df = as.data.frame(surr[, 1:max.surrs])
	   colnames(surr.df) = paste("surr", as.character(1:max.surrs), sep = "")
	   grp = factor(colnames(surr.df), levels = colnames(surr.df))

	   bp.df = data.frame(grp, prop.vars)
	   bp <- qplot(grp, data = bp.df, geom="bar", weight = prop.vars) + scale_x_discrete("Surrogate Variable") + scale_y_continuous("Explained Variance Proportion")
	   bp + theme(axis.text.x = element_text(angle = 0, hjust = 1))
	}

	res = new("svplsSurr", surr = surr, prop.vars = prop.vars)
}
return(res)
}