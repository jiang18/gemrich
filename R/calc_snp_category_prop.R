#' Calculate SNP Proportions in Annotation Categories
#'
#' This function calculates the proportion of SNPs that fall within each functional annotation category, relative to the total number of SNPs across the genome.
#'
#' @param snplist A `data.frame` or `data.table` containing genomic positions across all SNPs. 
#'        Required columns are chr and pos.
#' @param bed A `data.frame` or `data.table` containing functional annotations from a BED file. 
#'        The first four columns are chr, start, end, and category.
#' @param category_list An optional character vector specifying the priority order of categories for SNP assignment. 
#'        Required if there is more than one category in the bed input.
#'        SNPs will be assigned to mutually exclusive categories, where categories listed earlier take precedence.
#'        A SNP in multiple categories will be assigned to the category with the highest priority.
#' 
#' @return A `data.table` containing the proportion of SNPs in each annotation category.
#' 
#' @examples
#' snplist <- data.frame(
#'   chr = c("1", "1", "2"),
#'   pos = c(1000, 2000, 2000)
#' )
#' bed <- data.frame(
#'   chr = c("1", "1", "2"),
#'   start = c(500, 1800, 1400),
#'   end = c(1500, 2200, 1600),
#'   category = c("CDS", "intron", "CDS")
#' )
#' cat_prop = calc_snp_category_prop(snplist, bed, category_list = c("CDS", "intron"))
#' 
#' @export
#' 
calc_snp_category_prop <- function(snplist, bed, category_list = NULL) {
  # Ensure the input is a data.table
  setDT(bed)
  setDT(snplist)

  if (ncol(bed) < 4) {
    stop("Input 'bed' should have the following first four columns: chr, start, end, category.")
  }
  if(!all( c("chr", "pos") %in% colnames(snplist) )) {
    stop("Input 'snplist' needs to have the following two columns: chr and pos.")
  }

  cat_list = unique(bed[[4]])

  # Handle category_list
  if (is.null(category_list)) {
    if( length(cat_list) == 1 ) {
      category_list = cat_list
    } else {
      stop(paste(c("'category_list' is required when there are multiple categories in 'bed'"), collapse=" "))
    }
  } else {
    category_list = unique(category_list)
    if(any(!(category_list %in% cat_list))) {
      stop(paste(c("Categories in 'category_list' should be available in 'bed'"), collapse=" "))
    }
  }

  cat("There are", nrow(snplist), "SNPs in 'snplist'.\n")
  flush.console()

  if (!is.numeric(snplist$pos)) {
    stop("Error: column 'pos' should be numeric.")
  }

  snplist <- snplist[order(chr, pos)]
  snplist[, `:=`(pos2 = pos, id = .I)]

  bed = bed[, 1:4]
  setnames(bed, c("chr", "start", "end", "category"))
  bed$chr <- gsub("^chr", "", bed$chr, ignore.case = TRUE)
  bed <- bed[order(chr, start)]
  bed <- bed[category %chin% category_list]

  snplist[, chr := as.character(chr)]
  bed[, chr := as.character(chr)]

  cat("Setting keys for 'snplist' and 'bed'...\n")
  flush.console()

  # Perform a join using data.table
  setkey(snplist, chr, pos, pos2)
  setkey(bed, chr, start, end)

  cat("Finding overlaps of 'snplist' and 'bed'...\n")
  flush.console()

  # Join SNPs to genomic features
  overlaps <- foverlaps(
    snplist,
    bed,
    by.x = c("chr", "pos", "pos2"),
    by.y = c("chr", "start", "end"),
    nomatch = 0
  )

  if (nrow(overlaps) == 0) {
    cat("No SNPs are within genomic features. No SNP info file generated.\n")
    return(NULL)
  }
  cat(length(unique(overlaps$id)), "SNPs are within genomic annotations.\n")
  flush.console()

  categories <- unique(overlaps$category)
  snp_info <- data.table(SNP = snplist$id)

  # Add categories with NA by default
  snp_info[, (categories) := NA_real_]

  # Create unique id-category pairs for efficiency
  unique_pairs <- unique(overlaps[, .(id, category)])
  # Set values using data.table syntax
  for (cat in categories) {
      matched_ids <- unique_pairs[category == cat, id]
      set(snp_info, i = which(snp_info$SNP %in% matched_ids), j = cat, value = 1)
  }

  if (length(categories) == 1){
    snp_info_merged = snp_info
    colnames(snp_info_merged)[2] = "multi_cat"
    snp_info_merged$multi_cat[is.na(snp_info_merged$multi_cat)] = "remaining"
  }else{
    # Create output data.table
    snp_info_merged = data.table(SNP = snp_info$SNP, multi_cat = NA_character_)  # Changed to character
    
    # Process all categories at once
    for (cat in category_list) {
      if(!(cat %in% categories)) next
      # Use set() for more efficient updates
      set(snp_info_merged, 
          i = which(snp_info[[cat]] == 1 & is.na(snp_info_merged$multi_cat)), 
          j = "multi_cat",
          value = cat)
    }
    snp_info_merged$multi_cat[is.na(snp_info_merged$multi_cat)] = "remaining"
  }

  result = table(snp_info_merged$multi_cat) / nrow(snp_info_merged)
  if (length(categories) == 1) names(result)[names(result) == "1"] = categories[1]
  result = as.data.table(result[c(category_list, "remaining")])
  colnames(result) = c("category", "prop")
  result[,1] = c(category_list, "remaining")
  result[is.na(prop), prop := 0]
  return(result[])
}

