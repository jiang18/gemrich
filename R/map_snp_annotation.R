#' Map SNPs to Functional Annotations
#'
#' This function maps SNPs to functional annotations based on their genomic positions. 
#' It takes fine-mapping summary data and functional annotation BED data as input and 
#' outputs a data table of SNP-to-annotation mappings.
#'
#' @param bfmap A `data.frame` or `data.table` containing fine-mapping summary statistics generated by BFMAP (forward selection).
#'        Required columns are Chr, Pos, and SNPname.
#' @param bed A `data.frame` or `data.table` containing functional annotations from a BED file. 
#'        The first four columns are chr, start, end, and category.
#' @param multi_cat A logical value indicating whether to assign SNPs to each category independently or to a set of mutually exclusive categories. 
#'        If FALSE (default), SNPs will be assigned to each category independently, resulting in binary indicators for each category. 
#'        If TRUE, SNPs will be assigned to a set of mutually exclusive categories, and category_list is required to provide the priority order.
#' @param category_list An optional character vector specifying the priority order of categories for SNP assignment. 
#'        Categories listed earlier will take precedence. 
#'        A SNP in multiple categories will be assigned to the category with the highest priority.
#' 
#' @return A `data.table` containing SNP-to-annotation mappings. For multi_cat = FALSE, 
#'         the output includes a binary indicator column for each category. For multi_cat = TRUE, 
#'         it includes a single category assignment column based on the provided priority order.
#' 
#' @examples
#' bfmap <- data.frame(
#'   Chr = c("1", "1", "2"),
#'   Pos = c(1000, 2000, 1500),
#'   SNPname = c("rs1", "rs2", "rs3")
#' )
#' bed <- data.frame(
#'   chr = c("1", "1", "2"),
#'   start = c(500, 1800, 1400),
#'   end = c(2200, 2200, 1600),
#'   category = c("CDS", "intron", "CDS")
#' )
#' snp2annot = map_snp_annotation(bfmap, bed)
#' snp2annot = map_snp_annotation(bfmap, bed, category_list = c("CDS", "intron"), multi_cat = TRUE)
#' 
#' @export
#' 
map_snp_annotation <- function(bfmap, bed, category_list = NULL, multi_cat  = FALSE) {
  
  # Copies are not needed.
  setDT(bfmap)
  setDT(bed)

  required_columns <- c("Chr", "Pos", "SNPname")

  if (!all(required_columns %in% colnames(bfmap))) {
    missing_cols <- required_columns[!required_columns %in% colnames(bfmap)]
    stop(
      paste0(
        "Input 'bfmap' needs to contain the following columns: ",
        paste(required_columns, collapse = ", "),
        ". Missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  if (ncol(bed) < 4) {
    stop("Input 'bed' should have the following first four columns: chr, start, end, category.")
  }

  # Step 1: Process SNP positions
  bfmap <- bfmap[, .(chr = as.character(Chr), pos = as.integer(Pos), id = as.character(SNPname))]
  bfmap <- bfmap[order(chr, pos)]
  bfmap <- unique(bfmap)
  bfmap[, pos2 := pos]  # Duplicate position for foverlaps

  cat("There are", nrow(bfmap), "SNPs in 'bfmap'.\n")
  flush.console()

  # Step 2: Genomic annotation data
  bed = bed[, 1:4]
  setnames(bed, c("chr", "start", "end", "category"))
  bed$chr <- gsub("^chr", "", bed$chr, ignore.case = TRUE)
  bed <- bed[order(chr, start)]


  # Step 3: Map SNPs to genomic annotations
  cat("Mapping SNPs to genomic annotations...\n")
  flush.console()

  # Ensure positions are numeric for efficient joining
  bfmap[, chr := as.character(chr)]
  bed[, chr := as.character(chr)]

  # Perform a join using data.table
  setkey(bfmap, chr, pos, pos2)
  setkey(bed, chr, start, end)

  # Join SNPs to genomic annotations
  overlaps <- foverlaps(
    bfmap,
    bed,
    by.x = c("chr", "pos", "pos2"),
    by.y = c("chr", "start", "end"),
    nomatch = 0
  )

  if (nrow(overlaps) == 0) {
    cat("No SNPs are within genomic annotations. No SNP info file generated.\n")
    return(NULL)
  }
  cat(length(unique(overlaps$id)), "SNPs are within genomic annotations.\n")

  # Step 4: Generate output
  cat("Generating SNP-to-annotation mapping output...\n")
  flush.console()

  categories <- unique(overlaps$category)
  snp_info <- data.table(SNP = bfmap$id)

  # Add categories with NA by default
  snp_info[, (categories) := NA_character_]

  # Create unique id-category pairs for efficiency
  unique_pairs <- unique(overlaps[, .(id, category)])
  # Set values using data.table syntax
  for (cat in categories) {
      matched_ids <- unique_pairs[category == cat, id]
      set(snp_info, i = which(snp_info$SNP %in% matched_ids), j = cat, value = cat)
  }

  if (!multi_cat ){
    return(snp_info)
  }else{
    if (is.null(category_list)) {
      stop("Input 'category_list' is required when multi_cat = TRUE.")
    }
    category_list = unique(category_list)
    if (!all(category_list %in% categories)) {
      missing_categories <- category_list[!category_list %in% categories]
      stop(
        paste0(
          "Categories in 'category_list' need to be included in 'bed': ",
          paste(categories, collapse = ", "),
          ". Missing category: ",
          paste(missing_categories, collapse = ", ")
        )
      )
    } else {
      # Create output data.table
      snp_info_merge = data.table(SNP = bfmap$id, multi_cat = NA_character_)  # Changed to character
      
      # Process all categories at once
      for (cat in category_list) {
        # Use set() for more efficient updates
        set(snp_info_merge, 
            i = which(snp_info[[cat]] == cat & is.na(snp_info_merge$multi_cat)), 
            j = "multi_cat",
            value = cat)
      }
    }
    return(snp_info_merge)
  }
}

