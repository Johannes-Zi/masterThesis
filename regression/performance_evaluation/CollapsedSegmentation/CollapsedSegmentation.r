if (!require(pbapply)) {
  install.packages("pbapply")
}
library(pbapply)

if (!require(crayon)) {
  install.packages("crayon")
}
library(crayon)

if (!require(glmnet)) {
  install.packages("glmnet")
}
library(glmnet)

if (!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!require(tidyverse)) {
  install.packages("tidyverse")
}
library(tidyverse)

if (!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)

if (!require(tidyr)) {
  install.packages("tidyr")
}
library(tidyr)

if (!require(conflicted)) {
  install.packages("conflicted")
}
library(conflicted)


extract_elasticnet_model_segments <- function(elasticnet_model_file_path){
  load(elasticnet_model_file_path) # Load the .RData file - accessible as 'elasticnet_model'

  filename <- basename(elasticnet_model_file_path)
  
  gene_id <- strsplit(filename, "_")[[1]][5]

  coefficients <- coef(elasticnet_model$model, s = elasticnet_model$model$lambda.min)   #nolint

  coefficients_matrix <- as.matrix(coefficients)

  coefficients_df <- as.data.frame(coefficients_matrix)

  row_names <- rownames(coefficients_df)

  chr_indices <- grep("^chr", row_names)

  non_zero_s1_indices <- which(coefficients_df$s1 != 0)

  if (length(non_zero_s1_indices) == 1) {
    return(NULL)
  }

  #indices that are both in chr_indices and non_zero_s1_indices
  indices_to_keep <- intersect(chr_indices, non_zero_s1_indices)

  # only the rows of the coefficients object that start with "chr" and have != 0 coefficient
  coefficients_subset_df <- coefficients_df[indices_to_keep, , drop = FALSE]
  colnames(coefficients_subset_df)[colnames(coefficients_subset_df) == 's1'] <- 'coefficient'

  # row names to a column 'segment'
  coefficients_subset_df$segment <- rownames(coefficients_subset_df)

  # the row names
  rownames(coefficients_subset_df) <- NULL

  # gene ID
  coefficients_subset_df$gene_id <- gene_id

  # order of row names
  coefficients_subset_df <- coefficients_subset_df %>% select(gene_id, segment,)

  return(coefficients_subset_df)
}


add_atac_to_segment_df <- function(segment_df, segmentation_outptut_dir){
  # gene id
  gene_id <- segment_df$gene_id[1]

  # the gene_id specific segmentation output
  matching_segmentation_output_file <- list.files(segmentation_outptut_dir, pattern = paste0(".*", gene_id, ".*_Pearson.txt"), full.names = TRUE)[1]

  # segmentation output to df
  segmentation_output_df <- read.table(matching_segmentation_output_file, header = TRUE, sep = "\t")

  # row names to be the X column
  rownames(segmentation_output_df) <- segmentation_output_df$X

  # X column
  segmentation_output_df$X <- NULL

  # expression column
  segmentation_output_df$Expression <- NULL

  # data frame
  segmentation_output_df <- t(segmentation_output_df)

  segmentation_output_df <- as.data.frame(segmentation_output_df)

  # row names to column segment
  segmentation_output_df$segment <- rownames(segmentation_output_df)

  # row names to standard
  rownames(segmentation_output_df) <- NULL

  # segment column to be the first column
  segmentation_output_df <- segmentation_output_df[, c("segment", setdiff(names(segmentation_output_df), "segment"))]

  # segment_df and segmentation_output_df based on the segment column
  combined_df <- merge(segment_df, segmentation_output_df, by = "segment", all = FALSE)

  return(combined_df)
}

collapse_segment_df <- function(combined_df){
  col_names <- colnames(combined_df)

  cols_to_rename <- setdiff(col_names, c("segment", "gene_id"))

  new_col_names <- ifelse(col_names %in% cols_to_rename, paste0(col_names, "-ATAC"), col_names)

  colnames(combined_df) <- new_col_names

  segments <- strsplit(combined_df$segment, "\\.")
  combined_df$chromosome <- sapply(segments, "[", 1)
  combined_df$start <- as.integer(sapply(segments, "[", 2))
  combined_df$end <- as.integer(sapply(segments, "[", 3))

  combined_df$segment <- NULL

  combined_df <- combined_df[order(combined_df$chromosome, combined_df$start), ]

  collapsed_df <- combined_df[1, ]
  
  if (nrow(combined_df) >= 2) {
    for (i in 2:nrow(combined_df)) {

      # If the current segment is adjacent to the previous segment
      if (combined_df$chromosome[i] == collapsed_df$chromosome[nrow(collapsed_df)] && 
        combined_df$start[i] == (collapsed_df$end[nrow(collapsed_df)] + 1)) {

        # lengths of the current and previous segments
        curr_length <- combined_df$end[i] - combined_df$start[i] + 1

        # length of the previous segment (-curr_length, because end position is already updated)
        prev_length <- collapsed_df$end[nrow(collapsed_df)] - collapsed_df$start[nrow(collapsed_df)] + 1

        # end of the last row of the combined data frame
        collapsed_df$end[nrow(collapsed_df)] <- combined_df$end[i]

        # column names that end with "-ATAC"
        atac_cols <- grep("-ATAC$", names(collapsed_df), value = TRUE)
        # with the segment length weighted average for each ATAC column
        collapsed_df[nrow(collapsed_df), atac_cols] <- unlist(lapply(atac_cols, function(col) {
          # patient specific weighted average
          (combined_df[i, col] * curr_length + collapsed_df[nrow(collapsed_df), col] * prev_length) / 
            (curr_length + prev_length)
        }))
      }

      # Not adjacent case
      else {
        # current row to the combined data frame
        collapsed_df <- rbind(collapsed_df, combined_df[i, ])
      }
    }
  }
  # segment column based on chromosome, start, and end
  collapsed_df$segment <- paste(collapsed_df$chromosome, collapsed_df$start, collapsed_df$end, sep = ".")
  combined_df$segment <- paste(combined_df$chromosome, combined_df$start, combined_df$end, sep = ".")

  # chromosome, start, and end columns
  collapsed_df$chromosome <- NULL
  collapsed_df$start <- NULL
  collapsed_df$end <- NULL

  # row names
  rownames(collapsed_df) <- NULL

  # cat("\nFinal df: ", "###########################\n")

  #   cat("\nCombined df\n")
  #   print(combined_df[, c("segment", "healthy-6818-ATAC", "pah-6747-ATAC", "pah-6744-ATAC")])
  #   cat("\n")
    
  #   cat("Collapsed df\n")
  #   print(collapsed_df[, c("segment", "healthy-6818-ATAC", "pah-6747-ATAC", "pah-6744-ATAC")])
  #   cat("\n")

  return(collapsed_df)

}

create_elnet_segment_df <- function(elnet_model_dir, file_limit = Inf, segmentation_outptut_dir){
  elnet_model_files <- list.files(elnet_model_dir, pattern = "_Pearson.RData", full.names = TRUE)
  
  # number of files processed
  if(length(elnet_model_files) > file_limit) {
    elnet_model_files <- elnet_model_files[1:file_limit]
  }
  
  # progress bar options
  pboptions(char = "=")

  # elastic net model segments
  message(blue("  |-- Load and collapse elastic net model segments --|"))
  elnet_model_segment_df <- pblapply(elnet_model_files, function(file) {
    # elastic net model segments of the current file
    segment_df <- extract_elasticnet_model_segments(file)

    # no non-zero s1 values
    if (is.null(segment_df)) {
      return(data.frame())
    }

    # ATAC signal to the segment dataframe
    combined_df <- add_atac_to_segment_df(segment_df, segmentation_outptut_dir)

    # clos segments
    collapsed_df <-  collapse_segment_df(combined_df)
    return(collapsed_df)
  })
  
  message(blue("  |---- Combining all elastic net model segments ----|"))

  elnet_model_segments_df <- do.call(rbind, pblapply(elnet_model_segment_df, function(x) x))
  return(elnet_model_segments_df)
}




#'
#' Processing workflow
#' 

#' all segements that were considered in the regression models via cofficients into a single dataframe
elnet_model_dir <- "C:/Users/johan/Desktop/local_master_thesis_data/regression/LOneOCV_regression/regression_output/"
segmentation_outptut_dir <- "C:/Users/johan/Desktop/local_master_thesis_data/segmentation/combined_segmentation_output"
elnet_model_segments_df <- create_elnet_segment_df(elnet_model_dir, file_limit = 500, segmentation_outptut_dir = segmentation_outptut_dir)

# segment column to the first column
elnet_model_segments_df <- elnet_model_segments_df[, c("segment", setdiff(names(elnet_model_segments_df), "segment"))]

# column for chromosome
elnet_model_segments_df$chromosome <- sapply(strsplit(elnet_model_segments_df$segment, "\\."), "[", 1)

# column for start
elnet_model_segments_df$start <- as.integer(sapply(strsplit(elnet_model_segments_df$segment, "\\."), "[", 2))

# dataframe by chromosome and start
elnet_model_segments_df <- elnet_model_segments_df[order(elnet_model_segments_df$chromosome, elnet_model_segments_df$start), ]

# chromosome and start columns
elnet_model_segments_df$chromosome <- NULL
elnet_model_segments_df$start <- NULL

# dataframe as .tsv file in the cwd
write.table(elnet_model_segments_df, "elnet_model_segments_df.tsv", sep = "\t", row.names = FALSE)



# single elnet example
#test_path <- "C:/Users/johan/Desktop/local_master_thesis_data/regression/LOneOCV_regression/regression_output/Elasticnet_Regression_Model_Segmentation_ENSG00000035664_10_Pearson.RData"
#extract_elasticnet_model_segments(test_path)