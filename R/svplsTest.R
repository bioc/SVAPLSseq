#' @title svplsTest
#'
#' @description This function incorporates the significant surrogate variables
#' returned by the function \code{svplsSurr} in a linear model along with the 
#' group variable in order to estimate the group effect more accurately. The 
#' reestimated primary signal (group) effects are then used to test the genes 
#' for differential expression. The resulting pvalues are further corrected 
#' for multiple hypothesis testing at a prespecified FDR level. The 
#' significantly differentially expressed genes are finally returned along 
#' with their uncorrected and corrected pvalues.
#'
#' @param dat The original gene expression count matrix.
#' @param phi The transforming function to be applied on the original gene 
#' expression count data (set to be log function with an offset \code{const}).
#' @param const The offset parameter for the transforming function \code{phi}
#' (set to 1 by default).
#' @param group a factor representing the sample indices belonging to the two 
#' different groups.
#' @param surr A \code{data.frame} of the significant surrogate variables.
#' @param test The test to be used for detecting the differentially expressed 
#' genes. Options are "Wald" (Wald test with the gene-specific estimated group
#' effects after asjusting for the surrogate variables) and "LRT" (Likelihood 
#' Ratio Test).
#' @param mht.method The method to be used for the multiple hypothesis 
#' correction (set to the Benjamini-Hochberg procedure ("BH") by default).
#' @param fdr.level The specified level of the False Discovery Rate (FDR) for 
#' the multiple hypothesis testing (set to 0.05 by default).
#' @param parallel Logical, indicating if the computations should be 
#' parallelized or not (set to \code{FALSE} by default).
#' @param num.cores The requested number of cores to be used in the parallel 
#' computations inside the function (used only when \code{parallel} is 
#' \code{TRUE}, \code{NULL} by default).
#'
#' @return pvs.unadj The uncorrected pvalues corresponding to the genes after 
#' adjusting for the signatures of hidden variability.
#' @return pvs.adj The multiple hypothesis corrected pvalues after adjusting
#' for the signatures of hidden variability.
#' @return sig.genes The genes detected to be significantly differentially 
#' expressed between the two groups.
#'
#' @examples
#' ##Loading the simulated dataset
#' data(sim.dat)
#'
#' ##Fitting a linear model with the surrogate variables and detecting the 
#' ##differentially expressed genes
#' group = as.factor(c(rep(1, 10), rep(-1, 10)))
#' sv <- svplsSurr(dat = sim.dat, group = group)$surr
#' surr = surr(sv)
#' fit = svplsTest(dat = sim.dat, group = group, surr = surr, test = "Wald")
#'
#' ##The detected genes, hidden effect adjusted pvalues, FDR-corrected pvalues and the positive genes detected from the fitted model are given by:
#' head(sig.genes(fit))
#'
#' head(pvs.unadj(fit))
#'
#' head(pvs.adj(fit))
#'
#'
#' @rdname svplsTest 
#' @export
svplsTest <-
function(dat, phi = function(x) log(x + const), const = 1, group, surr, test = c("Wald", "LRT"), 
             mht.method = "BH", fdr.level = 0.05, parallel = FALSE, num.cores = NULL){

   if (class(dat) == "matrix") data = dat
   if (class(dat) == "SummarizedExperiment") data = assay(dat)
   if (class(dat) == "DGEList") data = dat$counts
   Y = phi(data)

   sv = paste("+", paste("surr[, ", paste(as.character(1:ncol(surr)), "]", sep = ""), sep = "", collapse = "+"), sep = "")

   if (test == "Wald"){
     	design = model.matrix(as.formula(paste0("~ group", sv)))
        dge = DGEList(data)
	dge = calcNormFactors(dge)
  	v <- voom(dge,design,plot=FALSE)
  	fit <- lmFit(v,design)
  	fit <- eBayes(fit)
      pvs = fit$p.value[,2]
   }

   if (test == "LRT"){
      lrt.test = function(row){
            fit1 = lm(as.formula(paste("Y[row, ] ~", sv, sep = " ")))
            fit2 = lm(as.formula(paste(paste("Y[row, ] ~ group", sv, sep = " "), sep = "")))
            pval = lrtest(fit1, fit2)[, 5][2]
            return(pval)
      }
    
      if (parallel){
          pvs = unlist(mclapply(1:nrow(dat), lrt.test, mc.cores = num.cores))
      }
      if (!parallel) pvs = unlist(lapply(1:nrow(dat), lrt.test))
   }

   pvs.adj = p.adjust(pvs, method = mht.method)
   if (is.null(row.names(dat))) sig.genes = as.character(which(pvs.adj < fdr.level))
   if (!is.null(row.names(dat))) sig.genes = row.names(dat)[which(pvs.adj < fdr.level)]   

   res = new("svplsTest", sig.genes = sig.genes, pvs.unadj = pvs, pvs.adj = pvs.adj)
   return(res)
}