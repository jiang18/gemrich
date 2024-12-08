#' Compute Empirical Distributions of Enrichments with Bootstrapping
#'
#' This function performs bootstrapping to construct empirical distributions of genetic effect enrichments in functional annotation categories.
#' 
#' @param bfmap A `data.frame` or `data.table` containing fine-mapping summary statistics generated by BFMAP (forward selection).
#'        Required columns are signal, SNPindex, Pval, normedProb, and SNPname.
#' @param snpinfo A `data.frame` or `data.table` containing SNP functional annotation data, which can be generated using `map_snp_annotation()`.
#' @param cat_prop A `data.frame` or `data.table` containing the proportion of SNPs in each annotation category. 
#'        This can be generated using either `calc_snp_category_prop()` or `calc_category_coverage()`.
#' @param annot A character string specifying the functional annotation column to analyze (default: "multi_cat"). 
#'        The specified name must match a column name in the snpinfo input.
#' @param pvalue_threshold A numeric threshold for filtering fine-mapped signals based on p-values (default: 5e-5). 
#'        Signals with lead SNP p-values below this threshold are included in the analysis.
#' @param n_bootstraps A numeric value specifying the number of bootstrap replicates (default: 1000).
#' @param seed A numeric value specifying the random seed for resampling signals (default: 1).
#' 
#' @return A list containing:
#'   \item{prob_mle}{Bootstrap distributions of probabilities of causal variants being in annotation categories}
#'   \item{enrichment_mle}{Bootstrap distributions of genetic effect enrichments in annotation categories}
#'   \item{loglik}{Bootstrap distribution of log-likelihood value at the maximum}
#' 
#' @examples
#' data("dairy_example")
#' # Use 5 bootstraps for quick testing
#' bootstrap_result <- bootstrap_category_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot,
#'   dairy_example$cat_prop,
#'   pvalue_threshold = 5e-5,
#'   n_bootstraps = 5
#' )
#' 
#' @export
#' 
bootstrap_category_enrichment <- function(bfmap, snpinfo, cat_prop, annot = "multi_cat", pvalue_threshold = 5e-5, n_bootstraps = 1000, seed = 1) {
  # Validate inputs
  pvalue_threshold = as.numeric(pvalue_threshold)
  if(is.na(pvalue_threshold)) {
    stop("'pvalue_threshold' must be numeric.\n", call.=FALSE)
  }
  # Load data
  kept_signals <- bfmap$signal[bfmap$SNPindex == 0 & bfmap$Pval < pvalue_threshold]
  if(length(kept_signals) != length(unique(kept_signals))) {
    stop("Multiple signals have identical signal index.", call.=FALSE)
  }
  bfmap <- bfmap[bfmap$signal %in% kept_signals, ]
  setDT(bfmap)
  nloci = length(kept_signals)

  if (!(annot %in% colnames(snpinfo))) stop(paste(annot, "is not one of the column names in 'snpinfo'"))

  cat_prop[[1]] = factor(cat_prop[[1]], levels=cat_prop[[1]])
  cat_names = levels(cat_prop[[1]])

  snpinfo <- copy(as.data.table(snpinfo))
  snp_colname = colnames(snpinfo)[1]
  snpinfo = copy(snpinfo)[, .SD, .SDcols = c(snp_colname, annot)]
  snpinfo = snpinfo[snpinfo[[snp_colname]] %in% unique(bfmap$SNPname), ]
  snpinfo[is.na(snpinfo[[annot]]),2] = "remaining"
  if(! all(unique(snpinfo[[annot]]) %in% cat_names) ) {
    stop(paste("Some categories in 'snpinfo' are missing in 'cat_prop'.\n"), call.=FALSE)
  }
  snpinfo[[annot]] = factor(snpinfo[[annot]], levels=cat_names)
  cat_cnt = table(snpinfo[[annot]])
  
  levels(snpinfo[[annot]])= c(1:length(cat_names))
  var2cat = new.env(hash=TRUE)
  cats = apply(snpinfo, 1, function(x) var2cat[[x[1]]] = x[annot])
  bfmap[, cat_idx := as.numeric(unlist(sapply(SNPname, function(x) var2cat[[x]])))]

  cat("Number of unique model SNPs in each category:")
  print(cat_cnt)
  if(any(cat_cnt < 50)) {
    cat("\nWarning: Some categories involve too few signal SNPs, which may cause ML estimation problems.\n")
  }

  if(any(cat_prop[[2]]<=0)) {
    stop(paste("Category proportions in cat_prop must be positive.\n"), call.=FALSE)
  }
  freq = cat_prop[[2]]
  if(sum(freq) != 1) {
    cat("\nSum of proportions in cat_prop is not equal to 1. Scaled to 1.\n")
    freq = freq / sum(freq)
  }

  cat(paste("\nCompleted reading all data files.\n",
      "The ML estimation uses ", nloci, " loci and ", nrow(bfmap)," variants.\n",sep=""))
  flush.console()

  logLL = function(par) {
    par <- c(par, 1-sum(par))   # sum(par) = 1
    if(!any(par<0)) {
      ret = bfmap[, .(li = sum(normedProb * par[cat_idx] / freq[cat_idx])), by = signal][, sum(log(li))]
      return(ret)
    } else return(-100000)
  }

  ###############################################################################
  ## Bootstrapping
  ###############################################################################

  # Function to perform one-level resampling of signals
  resample_data <- function(dat) {
    # Get unique signals
    unique_signals <- unique(dat$signal)

    # Sample signals with replacement
    sampled_signals <- sample(unique_signals, replace = TRUE)

    # Re-index signals
    resampled_dat <- do.call(rbind, lapply(seq_along(sampled_signals), function(i) {
      # Get original signal ID
      orig_signal <- sampled_signals[i]

      # Get data for this signal
      signal_data <- dat[dat$signal == orig_signal, ]

      # Assign new unique signal index
      signal_data$signal <- i - 1  # To maintain 0-based indexing

      return(signal_data)
    }))

    return(resampled_dat[order(resampled_dat$signal), ])
  }

  set.seed(seed)

  cat("Started bootstrapping.\n")
  flush.console()

  dat0 = bfmap
  bootstrap = as.data.frame(matrix(NA, ncol=length(freq)+1, nrow=n_bootstraps))
  npars = length(freq) - 1
  for(i in 1:n_bootstraps) {
    bfmap = resample_data(dat0)
    kept_signals = c(0:(nloci-1))
    par = freq[1:(length(freq)-1)]
    if(npars == 1) {
      result <- optimize(logLL, c(0,1), maximum=TRUE)
      bootstrap[i,] = c(result$maximum, 1-result$maximum, result$objective)
    } else {
      result <- optim(par, logLL, method="L-BFGS-B", lower=rep(1e-6, npars), upper=rep(0.999999, npars), control=list(fnscale=-1))
      # ui = rbind(diag(npars), rep(-1, npars))
      # ci = c(rep(1e-6, npars), -1)
      # result <- constrOptim(theta=par, f=logLL, grad=NULL, ui=ui, ci=ci, control=list(fnscale=-1))
      bootstrap[i,] = c(result$par, 1-sum(result$par), result$value)
    }
    cat(sprintf("\rCompleted replicate %d", i))
    flush.console()
  }
  cat("\nCompleted bootstrapping.\n")
  flush.console()

  colnames(bootstrap) = c(cat_names, "loglik")
  output <- list()
  output[["prob_mle"]] = bootstrap[, 1:length(freq)]
  output[["enrichment_mle"]] = as.data.frame( t( t(bootstrap[, 1:length(freq)]) / freq ) )
  output[["loglik"]] = bootstrap$loglik

  return(output)
}

