% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_category_enrichment.R
\name{estimate_category_enrichment}
\alias{estimate_category_enrichment}
\title{Estimate Genetic Effect Enrichments for Annotation Categories Using Maximum Likelihood}
\usage{
estimate_category_enrichment(
  bfmap,
  snpinfo,
  cat_prop,
  annot = "multi_cat",
  pvalue_threshold = 5e-05
)
}
\arguments{
\item{bfmap}{A \code{data.frame} or \code{data.table} containing fine-mapping summary statistics generated by BFMAP (forward selection).
Required columns are signal, SNPindex, Pval, normedProb, and SNPname.}

\item{snpinfo}{A \code{data.frame} or \code{data.table} containing SNP functional annotation data, which can be generated using \code{map_snp_annotation()}.}

\item{cat_prop}{A \code{data.frame} or \code{data.table} containing the proportion of SNPs in each annotation category.
This can be generated using either \code{calc_snp_category_prop()} or \code{calc_category_coverage()}.}

\item{annot}{A character string specifying the functional annotation column in snpinfo (default: "multi_cat").}

\item{pvalue_threshold}{A numeric threshold for filtering fine-mapped signals based on p-values (default: 5e-5).
Only signals with lead SNP p-values below this threshold are included in the analysis.}
}
\value{
A list containing:
\item{prob_mle}{Maximum likelihood estimates of the probabilities of causal variants being in annotation categories}
\item{prob_cov_matrix}{Covariance matrix of the probability estimates}
\item{loglik}{Log-likelihood value at the maximum}
\item{enrichment_mle}{Maximum likelihood estimates of genetic effect enrichments}
}
\description{
This function performs maximum likelihood estimation (MLE) for the probability of causal variants being in each functional annotation category.
The output includes estimated probabilities, standard errors, enrichment fold values,
and profile-likelihood confidence intervals of probabilities (available only for binary functional annotations).
}
\examples{
data("dairy_example")
mle_result <- estimate_category_enrichment(
  dairy_example$bfmap,
  dairy_example$snp2annot,
  dairy_example$cat_prop,
  pvalue_threshold = 5e-5
)

}
