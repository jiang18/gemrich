#' Calculate Genomic Coverage of Annotation Categories
#'
#' This function calculates the proportion of genomic sequence covered by each functional annotation category, relative to the total genome length.
#'
#' @param bed A `data.frame` or `data.table` containing functional annotations from a BED file. 
#'        The first four columns are chr, start, end, and category.
#' @param category_list An optional character vector specifying the priority order of categories. 
#'        Genomic regions will be assigned to mutually exclusive categories, where categories listed earlier take precedence.
#'        A region in multiple categories will be assigned to the category with the highest priority.
#' @param bp A numeric value specifying the total genome length in base pairs (default is 2489385779, for cattle autosomes).
#' 
#' @return A `data.table` containing the proportion of genomic sequence covered by each functional annotation category.
#' 
#' @examples
#' bed <- data.frame(
#'   chr = c("1", "1", "2", "1", "2"),
#'   start = c(1000, 1500, 2000, 1800, 2200),
#'   end = c(2000, 3000, 2500, 2200, 2500),
#'   category = c("CDS", "intron", "CDS", "CDS", "intron")
#' )
#' # Calculate the proportion of genome covered by each category
#' cat_prop = calc_category_coverage(bed, category_list = c("CDS", "intron"))
#' 
#' @export
#' 
calc_category_coverage <- function(bed, category_list = NULL, bp = 2489385779) {
  bed <- copy(bed)
  setDT(bed)
  if (ncol(bed) < 4) {
    stop("Input 'bed' should have the following first four columns: chr, start, end, category.")
  }
  
  # Set correct column names
  setnames(bed, names(bed), c("chr", "start", "end", "category"))
  
  cat_list = unique(bed$category)
  
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
  
  cat(paste(c("Input 'bed' contains", length(cat_list), "annotation categories:", cat_list, "\n"), collapse=" "))
  flush.console()

  # Initialize result table
  result <- data.table(category = category_list, length = 0)
  
  # First merge overlapping regions within each category
  merged_regions <- data.table()
  for(cat in category_list) {
    cat_regions <- copy(bed[category == cat])
    if(nrow(cat_regions) > 0) {
      # Prepare for foverlaps
      setorder(cat_regions, chr, start)
      cat_regions[, end2 := end]
      setkey(cat_regions, chr, start, end2)
      
      # Find overlaps within category
      overlaps <- foverlaps(cat_regions, cat_regions, type="any")
      overlaps <- overlaps[i.start <= end & start <= i.end]
      
      # Merge overlapping regions
      overlaps[, merge_group := .GRP, by = .(chr, i.start)]
      merged <- overlaps[, .(start = min(pmin(start, i.start)),
                           end = max(pmax(end, i.end))),
                       by = .(chr, merge_group, category)]
      merged[, merge_group := NULL]
      
      merged_regions <- rbind(merged_regions, unique(merged))
    }
  }
  
  # Prepare merged regions for overlap detection
  merged_regions[, end2 := end]
  setkey(merged_regions, chr, start, end2)
  
  # Process categories in priority order
  for(cat in category_list) {
    print(paste("Processing category:", cat))
    cat_regions <- merged_regions[category == cat]
    print(paste("Number of regions:", nrow(cat_regions)))
    
    if(nrow(cat_regions) > 0) {
      cat_length <- sum(cat_regions$end - cat_regions$start)
      print(paste("Initial length:", cat_length))
      
      higher_cats <- category_list[1:which(category_list == cat) - 1]
      print(paste("Higher priority categories:", paste(higher_cats, collapse=", ")))
      
      if(length(higher_cats) > 0) {
        higher_regions <- merged_regions[category %in% higher_cats]
        print(paste("Number of higher priority regions:", nrow(higher_regions)))
        
        if(nrow(higher_regions) > 0) {
          overlaps <- foverlaps(cat_regions, higher_regions, type="any")
          print(paste("Number of overlaps found:", nrow(overlaps)))
          
          if(nrow(overlaps) > 0) {
            # Calculate overlap lengths
            overlaps[, overlap_length := pmax(0, pmin(i.end, end) - pmax(i.start, start))]
            total_overlap <- sum(overlaps[overlap_length > 0, overlap_length])
            print(paste("Total overlap length:", total_overlap))
            cat_length <- cat_length - total_overlap
            print(paste("Final length after removing overlaps:", cat_length))
          }
        }
      }
      result[category == cat, length := cat_length]
    }
  }
  
  result[, prop := length / bp]
  result[, length := NULL]
  if(nrow(result) > 0) {
    result <- rbind(result, data.table(category = "remaining", prop = 1-sum(result[!is.na(prop), prop])))
  }
  return(result[])
}
