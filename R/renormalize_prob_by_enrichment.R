#' Renormalize Fine-Mapping Probabilities Using Functional Enrichment Estimates
#'
#' This function takes SNP probabilities from fine-mapping and renormalizes them within each
#' signal using functional enrichment estimates. For each signal, SNP probabilities are 
#' multiplied by enrichment factors and renormalized to their original sum.
#'
#' @param bfmap A `data.frame` containing fine-mapping summary data generated by BFMAP (forward selection).
#'        Required columns are signal, normedProb, and SNPname.
#' @param snpinfo A `data.frame` or `data.table` containing SNP functional annotation data, which can be generated using `map_snp_annotation()`.
#' @param enrichment_mle A `data.frame` of enrichment estimates from estimate_category_enrichment(),
#'        containing fold enrichments of genetic effects across functional annotation categories.
#' @param annot A character string specifying a functional annotation column in snpinfo (default: "multi_cat").
#' 
#' @return A `data.table` with renormalized probabilities (renormedProb). Other columns from the input bfmap are preserved.
#' 
#' @examples
#' data("dairy_example")
#' mle_result <- estimate_category_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot,
#'   dairy_example$cat_prop
#' )
#' renormed_bfmap <- renormalize_prob_by_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot,
#'   mle_result$enrichment_mle
#' )
#' 
#' @export
#' 
renormalize_prob_by_enrichment <- function(bfmap, snpinfo, enrichment_mle, annot = "multi_cat") {
  # Column checks
  required_bfmap <- c("signal", "SNPname", "normedProb")
  
  if (!all(required_bfmap %in% colnames(bfmap))) {
    stop("Missing bfmap columns: ", paste(setdiff(required_bfmap, colnames(bfmap)), collapse=", "))
  }
  if (!(annot %in% colnames(snpinfo))) stop(paste(annot, "is not one of the column names in 'snpinfo'"))

  bfmap <- copy(bfmap)
  setDT(bfmap)
  signals = unique(bfmap$signal)
  nloci = length(signals)

  enrichment_mle[[1]] = factor(enrichment_mle[[1]], levels=enrichment_mle[[1]])
  cat_names = levels(enrichment_mle[[1]])

  snp_colname = colnames(snpinfo)[1]
  snpinfo = snpinfo[,mget(c(snp_colname, annot))]
  snpinfo = snpinfo[snpinfo[[snp_colname]] %in% unique(bfmap$SNPname), ]
  snpinfo[is.na(snpinfo[[annot]]),2] = "remaining"
  if(! all(unique(snpinfo[[annot]]) %in% cat_names) ) {
    stop(paste("Some categories in 'snpinfo' are missing in 'enrichment_mle'.\n"), call.=FALSE)
  }
  snpinfo[[annot]] = factor(snpinfo[[annot]], levels=cat_names)
  
  levels(snpinfo[[annot]])= c(1:length(cat_names))
  var2cat = new.env(hash=TRUE)
  cats = apply(snpinfo, 1, function(x) var2cat[[x[1]]] = x[annot])
  bfmap[, cat_idx := as.numeric(unlist(sapply(SNPname, function(x) var2cat[[x]])))]

  if(any(enrichment_mle[[2]]<=0)) {
    stop(paste("Fold enrichments in enrichment_mle must be positive.\n"), call.=FALSE)
  }
  fold = enrichment_mle[[2]]

  cat(paste("\nCompleted reading all data files.\n",
      "Renormalization will cover ", nloci, " loci and ", nrow(bfmap)," variants.\n",sep=""))
  flush.console()
  
  bfmap[, category := cat_names[cat_idx]]
  bfmap[, renormedProb := {
    ff <- fold[cat_idx]
    normedProb * ff / sum(normedProb * ff) * sum(normedProb)
  }, by = signal]
  bfmap[, cat_idx := NULL]

  cat("Renormalization completed normally.\n")
  
  bfmap[order(-renormedProb), SNPindex2 := seq_len(.N)-1L, by=signal]
  setorder(bfmap, signal, SNPindex2)
  return(bfmap[])
}
