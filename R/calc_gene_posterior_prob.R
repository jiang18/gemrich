#' Calculate Gene-level Posterior Probabilities from Fine-Mapping
#'
#' This function calculates the posterior probability of genes containing causal variants
#' by integrating BFMAP forward-selection fine-mapping results with gene annotations. 
#' The probability for each gene is computed by summing the normalized probabilities of SNPs within its region.
#'
#' @param bfmap A `data.frame` or `data.table` containing fine-mapping summary statistics generated by BFMAP (forward selection).
#'        Required columns are signal, SNPindex, Chr, Pos, Pval, and normedProb.
#' @param gene_annot A `data.frame` or `data.table` containing gene coordinates and identifiers.
#'        Required columns are chr, start, end, and gene_id.
#' @param extension A numeric value specifying the number of base pairs to extend upstream and downstream of each gene (default: 0).
#' @param pvalue_threshold A numeric threshold for filtering fine-mapped signals based on p-values (default: 5e-5). 
#'        Signals with lead SNP p-values below this threshold are included in the analysis.
#' 
#' @return A `data.table` containing genes and their posterior probabilities of containing causal variants.
#'
#' @examples
#' data("dairy_example")
#' gene_probs = calc_gene_posterior_prob(
#'   dairy_example$bfmap, 
#'   dairy_example$gene_annot, 
#'   extension = 1000, 
#'   pvalue_threshold = 5e-7
#' )
#' 
#' @export
#'
calc_gene_posterior_prob <- function(bfmap, gene_annot, extension = 0, pvalue_threshold = 5e-5) {
  setDT(bfmap)
  setDT(gene_annot)
  
  bfmap[, Chr := as.character(Chr)]
  message("Starting calculations..."); flush.console()
  
  # Column checks
  required_bfmap <- c("signal", "SNPindex", "Chr", "Pos", "Pval", "normedProb")
  required_gene <- c("chr", "start", "end", "gene_id")
  
  if (!all(required_bfmap %in% colnames(bfmap))) {
    stop("Missing bfmap columns: ", paste(setdiff(required_bfmap, colnames(bfmap)), collapse=", "))
  }
  if (!all(required_gene %in% colnames(gene_annot))) {
    stop("Missing gene_annot columns: ", paste(setdiff(required_gene, colnames(gene_annot)), collapse=", "))
  }
  
  # Process reverse strand
  message("Processing gene coordinates..."); flush.console()
  reverse_idx <- gene_annot[, which(end < start)]
  if (length(reverse_idx) > 0) {
    message(sprintf("Found %d reverse strand genes", length(reverse_idx))); flush.console()
    gene_annot[reverse_idx, c("start", "end") := .(end, start)]
  }
  
  # Filter signals and extend boundaries
  message("Filtering signals..."); flush.console()
  sig_signals <- bfmap[SNPindex == 0 & Pval <= pvalue_threshold, unique(signal)]
  message(sprintf("Found %d significant signals", length(sig_signals))); flush.console()
  
  bfmap_filtered <- bfmap[signal %in% sig_signals]
  gene_annot[, `:=`(
    start_ext = start - extension,
    end_ext = end + extension
  )]
  
  # Calculate overlaps and sum probabilities
  message("Calculating overlaps..."); flush.console()
  results <- bfmap_filtered[gene_annot, on = .(Chr = chr), 
    allow.cartesian = TRUE, nomatch = 0][
      Pos >= start_ext & Pos <= end_ext,
      .(summed_prob = sum(normedProb)), 
      by = .(gene_id, signal)
    ]
  
  # Add gene information and lead p-values
  results <- gene_annot[results, on = "gene_id"]
  lead_pvals <- bfmap_filtered[SNPindex == 0, .(signal, lead_pvalue = Pval)]
  results <- lead_pvals[results, on = "signal"]
  
  message(sprintf("Found %d gene-signal pairs", nrow(results))); flush.console()

  setcolorder(results, c("gene_id", "chr", "start", "end", "signal", "summed_prob", "lead_pvalue", "start_ext", "end_ext"))
  setorder(results, -summed_prob, chr, start, signal)
  return(results[])
}
