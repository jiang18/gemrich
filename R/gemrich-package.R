#' gemrich: Genetic Effect Mapping enRICHment
#'
#' @description
#' An R package for estimating genetic effect enrichments in functional annotations 
#' from fine-mapping summary statistics. The package takes BFMAP results as input.
#'
#' @details
#' The package provides tools for genetic effect enrichment analysis and 
#' fine-mapping integration. Key functionalities include enrichment estimation, 
#' bootstrapping for empirical distributions, enrichment-informed probability renormalization, 
#' and genomic-feature posterior probability calculations.
#'
#' @section Main Functions:
#' \strong{Enrichment Analysis}
#' \itemize{
#'   \item \code{estimate_category_enrichment()}: Estimate genetic effect enrichments 
#'         for annotation categories using maximum likelihood
#'   \item \code{bootstrap_category_enrichment()}: Compute empirical distributions 
#'         of enrichments with bootstrapping
#' }
#'
#' \strong{Fine-Mapping Integration}
#' \itemize{
#'   \item \code{renormalize_prob_by_enrichment()}: Update SNP probabilities using 
#'         functional enrichment estimates
#'   \item \code{calc_feature_posterior_prob()}: Calculate genomic-feature posterior probabilities
#' }
#'
#' \strong{Annotation Processing}
#' \itemize{
#'   \item \code{map_snp_annotation()}: Map SNPs to functional annotations
#'   \item \code{calc_category_coverage()}: Calculate genomic coverage of annotation categories
#'   \item \code{calc_snp_category_prop()}: Calculate SNP proportions in annotation categories
#' }
#'
#' @section Examples:
#' \strong{Example Data}
#' 
#' The package includes a dairy cattle example dataset (\code{dairy_example}) containing:
#' \itemize{
#'   \item \code{bfmap}: BFMAP forward-selection summary statistics for five dairy production traits
#'   \item \code{snp2annot}: SNP functional annotations
#'   \item \code{cat_prop}: SNP proportions in annotation categories
#'   \item \code{gene_annot}: Gene annotations
#' }
#'
#' \strong{Enrichment Analysis}
#' \preformatted{
#' # Estimate enrichments using MLE
#' mle_result <- estimate_category_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot, 
#'   dairy_example$cat_prop,
#'   pvalue_threshold = 5e-5
#' )
#'
#' # Get empirical distributions with bootstrapping  
#' bootstrap_result <- bootstrap_category_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot,
#'   dairy_example$cat_prop,
#'   pvalue_threshold = 5e-5,
#'   n_bootstraps = 100
#' )
#' }
#'
#' \strong{Fine-Mapping Integration}
#' \preformatted{
#' # Renormalize probabilities using enrichment estimates
#' renormed_bfmap <- renormalize_prob_by_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot,
#'   mle_result$enrichment_mle
#' )
#'
#' # Calculate genomic-feature posterior probabilities
#' gene_probs <- calc_feature_posterior_prob(
#'   dairy_example$bfmap,
#'   dairy_example$gene_annot,
#'   extension = 1000, # extend gene boundaries by 1kb
#'   pvalue_threshold = 5e-7
#' )
#' }
#'
#' @author Jicai Jiang \email{jjiang26@ncsu.edu}
#'
#' @seealso
#' \href{https://github.com/jiang18/bfmap}{BFMAP} for the fine-mapping method
#' used to generate input data.
#'
#' @docType package
#' @name gemrich
"_PACKAGE"
