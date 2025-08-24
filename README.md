# gemrich: Genetic Effect Mapping enRICHment

An R package for estimating genetic effect enrichments in functional annotations from fine-mapping summary statistics.

The package takes [BFMAP](https://github.com/jiang18/bfmap) results as input. 

## Installation
```r
# Dependencies: data.table and Rcpp

# Install from GitHub
devtools::install_github("jiang18/gemrich")

# Coming to CRAN soon
# install.packages("gemrich")
```

## Documentation
```r
# Package overview
?gemrich

# Detailed function documentation
?estimate_category_enrichment
?bootstrap_category_enrichment
?renormalize_prob_by_enrichment
?calc_feature_posterior_prob
?map_snp_annotation
?calc_category_coverage
?calc_snp_category_prop
```

## Main Functions

### Enrichment Analysis
- `estimate_category_enrichment()`: Estimate genetic effect enrichments for annotation categories using maximum likelihood
- `bootstrap_category_enrichment()`: Compute empirical distributions of enrichments with bootstrapping (not available for SSS)

### Fine-Mapping Integration 
- `renormalize_prob_by_enrichment()`: Update posterior probabilities using functional enrichment estimates
- `calc_feature_posterior_prob()`: Calculate genomic-feature posterior probabilities

### Annotation Processing
- `map_snp_annotation()`: Map SNPs to functional annotations
- `calc_category_coverage()`: Calculate genomic coverage of annotation categories
- `calc_snp_category_prop()`: Calculate SNP proportions in annotation categories

## Usage Examples

### Example Data
The package includes a dairy cattle example dataset (`dairy_example`) containing:
- `bfmap`: BFMAP forward-selection summary statistics for five dairy production traits
- `snp2annot`: SNP functional annotations
- `cat_prop`: SNP proportions in annotation categories
- `gene_annot`: Gene annotations

### Enrichment Analysis
```r
library(gemrich)

# Load example data
data("dairy_example")

# Estimate enrichments using MLE
mle_result <- estimate_category_enrichment(
  dairy_example$bfmap,
  dairy_example$snp2annot, 
  dairy_example$cat_prop,
  pvalue_threshold = 5e-5
)

# Get empirical distributions with bootstrapping  
bootstrap_result <- bootstrap_category_enrichment(
  dairy_example$bfmap,
  dairy_example$snp2annot,
  dairy_example$cat_prop,
  pvalue_threshold = 5e-5,
  n_bootstraps = 100
)
```

### Fine-Mapping Integration
```r
# Renormalize probabilities using enrichment estimates
renormed_bfmap <- renormalize_prob_by_enrichment(
  dairy_example$bfmap,
  dairy_example$snp2annot,
  mle_result$enrichment_mle
)

# Calculate gene-level probabilities
gene_probs <- calc_feature_posterior_prob(
  dairy_example$bfmap,
  dairy_example$gene_annot,
  extension = 1000, # extend gene boundaries by 1kb
  pvalue_threshold = 5e-7
)
```

## Data Resources
Coming soon

## Citation
Coming soon

## Contact
Jicai Jiang (jjiang26@ncsu.edu)
