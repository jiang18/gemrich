% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_snp_category_prop.R
\name{calc_snp_category_prop}
\alias{calc_snp_category_prop}
\title{Calculate SNP Proportions in Annotation Categories}
\usage{
calc_snp_category_prop(snplist, bed, category_list = NULL)
}
\arguments{
\item{snplist}{A \code{data.frame} or \code{data.table} containing genomic positions across all SNPs.
Required columns are chr and pos.}

\item{bed}{A \code{data.frame} or \code{data.table} containing functional annotations from a BED file.
The first four columns are chr, start, end, and category.}

\item{category_list}{An optional character vector specifying the priority order of categories for SNP assignment.
Required if there is more than one category in the bed input.
SNPs will be assigned to mutually exclusive categories, where categories listed earlier take precedence.
A SNP in multiple categories will be assigned to the category with the highest priority.}
}
\value{
A \code{data.table} containing the proportion of SNPs in each annotation category.
}
\description{
This function calculates the proportion of SNPs that fall within each functional annotation category, relative to the total number of SNPs across the genome.
}
\examples{
snplist <- data.frame(
  chr = c("1", "1", "2"),
  pos = c(1000, 2000, 2000)
)
bed <- data.frame(
  chr = c("1", "1", "2"),
  start = c(500, 1800, 1400),
  end = c(1500, 2200, 1600),
  category = c("CDS", "intron", "CDS")
)
cat_prop = calc_snp_category_prop(snplist, bed, category_list = c("CDS", "intron"))

}
