# Load required libraries
library(dplyr)  # For data manipulation
library(MASS)  # For kernel density estimation
library(ggplot2)  # For data visualization
library(viridis)  # For color palettes


preprocess_performance_evaluation_df <- function(input_df) {
  # original filenames for each df entry based on the Sample_Name column
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")  # Split the Sample_Name column by "_"
  #head(filename_parts)  # Print the first few filename parts

  # segmentation preselection correlation type
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))

  # new column in the dataframe for the segmentation preselection correlation type
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # Eentry specific gene identifiers
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Sample_Name column to gene_name (identifiers)
  adapted_df <- rename(adapted_df, gene_name = Sample_Name)

  return(adapted_df)
}


process_file <- function(file_path) {
  load(file_path)

  # filename without extension
  filename <- tools::file_path_sans_ext(basename(file_path))

  # ident of lambda min based on the list with all lambdas
  lambda_min_value <- elasticnet_model$model$lambda.min
  lambda_min_index <- which(c == lambda_min_value)

  # corresponding nzero value
  nzero_value <- elasticnet_model$model$nzero[lambda_min_index]
  number_of_features <- elasticnet_model$model$glmnet.fit$dim[1]

  # data frame with the results
  result <- data.frame(filename = filename,
                       number_of_features = number_of_features,
                       nzero_value = nzero_value)
  return(result)
}


process_directory <- function(directory_path, import_limit, file_pattern) {
  # list of .RData files in the directory
  file_list <- list.files(directory_path, pattern = file_pattern,
                          full.names = TRUE)
  print(paste("Number of detected files:", length(file_list)))
  print(paste("Used filename pattern", file_pattern))

  # number of files to import
  if (length(file_list) > import_limit) {
    file_list <- file_list[1:import_limit]
    print(paste("Number of imported files limited to", import_limit))
  }

  # empty data frame to store the results
  result_df <- data.frame(filename = character(),
                          number_of_features = numeric(),
                          nzero_value = numeric(),
                          stringsAsFactors = FALSE)

  for (i in seq_along(file_list)) {
    file_path <- file_list[i]
    result <- process_file(file_path)
    result_df <- rbind(result_df, result)
    
    if (i %% 50 == 0) {
      print(paste("Processed", i, "files"))
    }
  }

  return(result_df)
}


process_regression_models_df <- function(input_df) {
  # new dataframe using the filename as key - extracting only the filename
  new_df <- data.frame(filename = unique(input_df$filename))

  # new dataframe with the input_df dataframe
  merged_df <- merge(new_df, input_df, by = "filename", all.x = TRUE)

  # filename column
  filename_parts <- strsplit(as.character(merged_df$filename), "_")

  # segmentation preselection correlation type
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(x[7], sep = "_"))

  merged_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # entry specific gene identifiers
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[5], x[6], sep = "_"))

  merged_df$filename <- extracted_gene_names

  # filename column to gene_name (identifiers)
  merged_df <- rename(merged_df, gene_name = filename)

  return(merged_df)
}


import_and_preprocess_performance_evaluation_data <- function(input_file) {
    # input file into a dataframe
    df <- read.table(input_file, header = TRUE, sep = "\t")
    print(head(df))
    
    df_preprocessed_1 <- preprocess_performance_evaluation_df(df)
    print(head(df_preprocessed_1))
    
    # entries based on Pearson based feature selection
    df_preprocessed_2 <- df_preprocessed_1[df_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
    print(head(df_preprocessed_2))
    
    return(df_preprocessed_2)
}


import_and_preprocess_regression_output_data <- function(input_directory_path, import_limit, file_pattern) {
    cat("imported directory path:\n", input_directory_path, "\n")

    # Process the directory and load the results - only for Pearson based feature (segment) preselection input
    output_df <- process_directory(input_directory_path, import_limit, file_pattern)
    print(head(output_df))
    
    # Second preprocessing step before merging the dataframes
    processed_regression_models_df <- process_regression_models_df(output_df)
    print(head(processed_regression_models_df))
    return(processed_regression_models_df)
}


combine_performance_evaluation_and_regression_data <- function(performance_evaluation_df, regression_models_df) {
    # performance evaluation and regression models dataframes
    combined_df <- merge(performance_evaluation_df, regression_models_df, by = c("gene_name", "segmentation_preselection_corellation"), all = TRUE)
    
    # Spearman based entries
    combined_filtered_df <- combined_df[combined_df$segmentation_preselection_corellation == "Pearson", ]
    
    return(combined_filtered_df)
}


process_bed_file <- function(file_path) {
    
    lines <- readLines(file_path)

    # lines that start with "chr"
    chr_lines <- lines[grep("^chr", lines)]

    # headers for the data frame
    headers <- c("chr", "start", "end", "coeffs", "idk")

    # filtered lines back into a data frame
    temp_df <- read.delim(text = chr_lines, header = FALSE, col.names = headers)

    # number of lines in the temp_df dataframe
    num_lines <- nrow(temp_df)

    # nzero value from the temp_df dataframe
    nzero_value <- sum(temp_df$coeffs != 0)

    # gene name from the filename
    filename <- tools::file_path_sans_ext(basename(file_path))
    filename_parts <- strsplit(as.character(filename), "_")  # Split the Sample_Name column by "_"
    extracted_gene_name <- sapply(filename_parts, function(x) paste(x[4], x[5], sep = "_"))

    # new row for the dataframe with the extracted information
    result_row <- data.frame(gene_name = extracted_gene_name,
                            number_of_features_ols = num_lines,
                            nzero_ols = nzero_value,
                            only_zero_coefficients_ols = sum(!(num_lines == nzero_value)))
    return(result_row)
}


process_ols_data <- function(directory_path, import_limit, file_pattern) {
  # list of .bed files in the directory
  file_list <- list.files(directory_path, pattern = file_pattern,
                          full.names = TRUE)
  print(paste("Number of detected files:", length(file_list)))
  print(paste("Used filename pattern", file_pattern))

  # number of files to import
  if (length(file_list) > import_limit) {
    file_list <- file_list[1:import_limit]
    print(paste("Number of imported files limited to", import_limit))
  }
  
  result_df <- data.frame(filename = character(),
                           number_of_features_ols = numeric(),
                           nzero_ols = numeric(),
                           only_zero_coefficients_ols = numeric())

  for (i in seq_along(file_list)) {
    file_path <- file_list[i]
    result <- process_bed_file(file_path)
    result_df <- rbind(result_df, result)
    
    if (i %% 50 == 0) {
      print(paste("Processed", i, "files"))
    }
  }

  return(result_df)
}


combine_ols_and_elnet_data <- function(ols_df, elnet_df) {
  # performance evaluation and regression models dataframes
  combined_df <- merge(elnet_df, ols_df, by = c("gene_name"), all = TRUE)
    
  return(combined_df)
}



# input directory paths for standard and LeaveOneOut regression performance evaluation data
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"

# Import and preprocess performance evaluation data for standard regression
df_standard_preprocessed_2 <- import_and_preprocess_performance_evaluation_data(input_directory_path_standard)
# Import and preprocess performance evaluation data for LeaveOneOut regression
df_LeaveOneOut_preprocessed_2 <- import_and_preprocess_performance_evaluation_data(input_directory_path_LeaveOneOut)


input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/regression_output/regression_output/"
#input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/regression_output_example/"
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/regression_output/"

# Import and preprocess regression output data for standard regression
processed_regression_models_df_standard <- import_and_preprocess_regression_output_data(input_directory_path_standard, 50000, "\\Pearson.RData$")
# Import and preprocess regression output data for LeaveOneOut regression
processed_regression_models_df_LeaveOneOut <- import_and_preprocess_regression_output_data(input_directory_path_LeaveOneOut, 50000, "\\Pearson.RData$")

# Combine performance evaluation and regression data for standard regression
combined_elnet_df_standard <- combine_performance_evaluation_and_regression_data(df_standard_preprocessed_2, processed_regression_models_df_standard)
# Combine performance evaluation and regression data for LeaveOneOut regression
combined_elnet_df_LeaveOneOut <- combine_performance_evaluation_and_regression_data(df_LeaveOneOut_preprocessed_2, processed_regression_models_df_LeaveOneOut)

print(head(combined_elnet_df_standard))
print(head(combined_elnet_df_LeaveOneOut))

ols_df_standard <- process_ols_data(input_directory_path_standard, 50000, "\\Pearson.bed$")
ols_df_LeaveOneOut <- process_ols_data(input_directory_path_LeaveOneOut, 50000, "\\Pearson.bed$")
print(head(ols_df_standard))
print(head(ols_df_LeaveOneOut))

final_df_standard <- combine_ols_and_elnet_data(ols_df_standard, combined_elnet_df_standard)
final_df_LeaveOneOut <- combine_ols_and_elnet_data(ols_df_LeaveOneOut, combined_elnet_df_LeaveOneOut)
print(head(final_df_standard))
print(head(final_df_LeaveOneOut))



#write.table(final_df_standard, file = "file.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Select standrd entries that are potentially problematic
selected_weird_entries_standard <- subset(final_df_standard, only_zero_coefficients_ols == 1 & pVal < 0.05 & qVal < 0.05)
print(nrow(selected_weird_entries_standard))
# Save weird OLS Models as file
write.table(selected_weird_entries_standard, file = "selected_weird_entries_standard.csv", sep = ",", quote = FALSE, row.names = FALSE)

