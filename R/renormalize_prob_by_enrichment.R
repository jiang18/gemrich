#' Renormalize Fine-Mapping Probabilities Using Functional Enrichment Estimates
#'
#' This function takes probabilities from fine-mapping and renormalizes them using functional 
#' enrichment estimates. For forward selection, SNP probabilities are renormalized within each
#' signal. For SSS, model probabilities are renormalized within each locus and variant PIPs
#' are calculated.
#'
#' @param bfmap A `data.frame` or `data.table` containing fine-mapping summary data.
#'        For forward selection: Required columns are signal, normedProb, and SNPname.
#'        For SSS: columns should be locus ID, variant names, logbf, and model probabilities.
#' @param snpinfo A `data.frame` or `data.table` containing SNP functional annotation data, which can be generated using `map_snp_annotation()`.
#' @param enrichment_mle A `data.frame` of enrichment estimates from estimate_category_enrichment(),
#'        containing fold enrichments of genetic effects across functional annotation categories.
#' @param annot A character string specifying a functional annotation column in snpinfo (default: "multi_cat").
#' @param input_type Character string specifying the input format: "forward_selection" (default) or "sss".
#' 
#' @return For forward selection: A `data.table` with renormalized probabilities (renormedProb).
#'         For SSS: A list containing $model (renormalized model probabilities) and $pip (posterior inclusion probabilities for variants).
#' 
#' @examples
#' data("dairy_example")
#' mle_result <- estimate_category_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot,
#'   dairy_example$cat_prop
#' )
#' # Forward selection
#' renormed_bfmap <- renormalize_prob_by_enrichment(
#'   dairy_example$bfmap,
#'   dairy_example$snp2annot,
#'   mle_result$enrichment_mle
#' )
#' # SSS
#' # renormed_sss <- renormalize_prob_by_enrichment(
#' #   sss_data,
#' #   snp2annot,
#' #   mle_result$enrichment_mle,
#' #   input_type = "sss"
#' # )
#' 
#' @export
#' 
renormalize_prob_by_enrichment <- function(bfmap, snpinfo, enrichment_mle, annot = "multi_cat", input_type = "forward_selection") {
  
  # Handle SSS format
  if (input_type == "sss") {
    return(renormalize_prob_by_enrichment_sss(bfmap, snpinfo, enrichment_mle, annot))
  }
  
  # Forward selection code (unchanged from original)
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
  bfmap[, cat_idx := as.integer(unlist(sapply(SNPname, function(x) var2cat[[x]])))]

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


# SSS renormalization function
renormalize_prob_by_enrichment_sss <- function(bfmap, snpinfo, enrichment_mle, annot = "multi_cat") {
  
  setDT(bfmap)
  setDT(snpinfo)
  
  # Annotation processing
  if(!(annot %in% colnames(snpinfo))) {
    stop(paste("Annotation [", annot, "] is not available in annotation data.\n"), call.=FALSE)
  }
  
  # Filter annotation to only variants present in bfmap
  pcol <- ncol(bfmap)
  bfmap_variants <- unique(unlist(bfmap[, 2:(pcol-2), with=FALSE], use.names=FALSE))
  bfmap_variants <- bfmap_variants[!is.na(bfmap_variants)]

  snp_colname = colnames(snpinfo)[1]
  snpinfo = copy(snpinfo)[, .SD, .SDcols = c(snp_colname, annot)]
  snpinfo <- snpinfo[get(snp_colname) %in% bfmap_variants]
  set(snpinfo, which(is.na(snpinfo[[annot]])), annot, "remaining")
  
  # Process enrichment data
  enrichment_mle[[1]] = factor(enrichment_mle[[1]], levels=enrichment_mle[[1]])
  cat_names = levels(enrichment_mle[[1]])
  n_cats <- length(cat_names)
  
  # Validate that annotation categories match enrichment_mle
  if(!all(unique(snpinfo[[annot]]) %in% cat_names)) {
    stop("Categories in annotation data are not all present in enrichment_mle.\n", call.=FALSE)
  }
  
  # Create environment-based lookup for category indices
  snpinfo[[annot]] = factor(snpinfo[[annot]], levels=cat_names)
  levels(snpinfo[[annot]]) = c(1:length(cat_names))
  
  var2cat = new.env(hash=TRUE)
  cats = apply(snpinfo, 1, function(x) var2cat[[x[1]]] = x[annot])
  
  # Get enrichment factors
  if(any(enrichment_mle[[2]] <= 0)) {
    stop("Fold enrichments in enrichment_mle must be positive.\n", call.=FALSE)
  }
  fold = enrichment_mle[[2]]
  
  # Data preprocessing
  locus_levels <- unique(bfmap$locus)
  nloci <- length(locus_levels)
  prob_col <- names(bfmap)[pcol]
  
  cat("Preprocessing data for renormalization...\n")
  flush.console()
  
  # Process each locus
  model_results <- vector("list", nloci)
  pip_results <- vector("list", nloci)
  
  for(i in seq_len(nloci)) {
    locus_subset <- bfmap[locus == locus_levels[i]]
    
    # Extract variants and model probabilities
    var_cols <- 2:(pcol-2)
    model_data <- as.matrix(locus_subset[, var_cols, with=FALSE])
    probs <- locus_subset[[pcol]]
    original_sum <- sum(probs)
    
    # Calculate enrichment factors for each model
    model_enrichments <- numeric(nrow(model_data))
    
    for(m in seq_len(nrow(model_data))) {
      variants <- model_data[m, ]
      variants <- variants[!is.na(variants)]
      
      if(length(variants) > 0) {
        # Get category indices for variants in this model
        cats <- mget(as.character(variants), envir = var2cat, ifnotfound = list(NA))
        cats <- as.integer(unlist(cats))
        
        # Calculate product of enrichment factors
        if(length(cats) > 0) {
          model_enrichments[m] <- prod(fold[cats])
        } else {
          model_enrichments[m] <- 1.0
        }
      } else {
        model_enrichments[m] <- 1.0
      }
    }
    
    # Apply enrichment factors and renormalize
    enriched_probs <- probs * model_enrichments
    renormed_probs <- enriched_probs / sum(enriched_probs) * original_sum
    
    # Update model data with renormalized probabilities
    locus_result <- copy(locus_subset)
    locus_result[, renormedProb := renormed_probs]
    model_results[[i]] <- locus_result
    
    # Single-pass PIP calculation for this locus
    current_locus <- locus_levels[i]
    variant_pips <- new.env(hash=TRUE)
    
    for(m in seq_len(nrow(model_data))) {
      variants <- model_data[m, ]
      variants <- variants[!is.na(variants)]
      model_prob <- renormed_probs[m]
      
      # Add this model's probability to each variant's PIP
      for(variant in variants) {
        variant_key <- as.character(variant)
        if(exists(variant_key, envir = variant_pips)) {
          variant_pips[[variant_key]] <- variant_pips[[variant_key]] + model_prob
        } else {
          variant_pips[[variant_key]] <- model_prob
        }
      }
    }
    
    # Collect PIP results for this locus
    if(length(ls(variant_pips)) > 0) {
      pip_variants <- ls(variant_pips)
      pip_values <- sapply(pip_variants, function(v) variant_pips[[v]])
      
      pip_results[[i]] <- data.table(
        variant = pip_variants,
        locus = current_locus,
        PIP = pip_values
      )
    }
  }
  
  # Combine all results
  model_data_final <- rbindlist(model_results)
  variant_pips_dt <- rbindlist(pip_results)
  setorder(variant_pips_dt, locus, -PIP)
  
  total_models <- nrow(model_data_final)
  total_variants <- nrow(variant_pips_dt)
  
  cat(paste("Renormalization complete. Processed", nloci, "loci,", 
            total_models, "models, and calculated PIPs for", 
            total_variants, "variants.\n"))
  flush.console()
  
  return(list(
    model = model_data_final,
    pip = variant_pips_dt
  ))
}
