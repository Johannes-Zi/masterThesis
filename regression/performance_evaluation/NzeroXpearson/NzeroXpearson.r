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

  # entry specific gene identifiers
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Sample_Name column to gene_name (identifiers)
  adapted_df <- rename(adapted_df, gene_name = Sample_Name)

  return(adapted_df)
}

process_file <- function(file_path) {
  # .RData file - accessible as 'elasticnet_model'
  load(file_path)

  # filename without extension
  filename <- tools::file_path_sans_ext(basename(file_path))

  # ident of lambda min based on the list with all lambdas
  lambda_min_value <- elasticnet_model$model$lambda.min
  lambda_min_index <- which(elasticnet_model$model$lambda == lambda_min_value)

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


get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}


create_dotplot <- function(df, output_path) {
  # density for each point in your data
  df$density <- get_density(df$nzero_value, df$Pearson, n = 1000)

  dotplot <- ggplot(df, aes(x = nzero_value, y = Pearson, color = density)) +
    geom_point(size = 2) +
    theme_bw() +
    #labs(title = "Utilized features vs. respective CV Pearson coefficients",
    labs(title = "",
       x = "Number of utilized features",
       y = "Cross validation Pearson coefficient") +
    scale_color_viridis() +
    scale_x_continuous(breaks = seq(min(df$nzero_value), max(df$nzero_value), by = 2)) +
    coord_cartesian(xlim = c(0, 29)) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.title.y = element_text(margin = margin(r = 10,)),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.text = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )

  ggsave(filename = paste(output_path, "/dotplot.png", sep = ""), plot = dotplot, width = 9, height = 6)
  
  ggsave(filename = paste(output_path, "/dotplot.svg", sep = ""), plot = dotplot, width = 9, height = 6)
}


input_directory_path_LeaveOneOut <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
#df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
#head(df_standard)

# Preprocess the standard CV dataframe
#df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
#head(df_standard_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
#df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
#head(df_standard_preprocessed_2)
# Preprocess the LeaveOneOut CV dataframe
df_LeaveOneOut_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_LeaveOneOut_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_LeaveOneOut_preprocessed_2)


#  For regression LeaveOneOut CV
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/temp_inputfile_for_MA_plots/RegFeaturesXnzero/regression_output/regression_output/"
cat("imported directory path:\n", input_directory_path_LeaveOneOut, "\n")
# Process the directory and load the results - only for Pearson based feature (segment) preselection input
output_df_LeaveOneOut <- process_directory(input_directory_path_LeaveOneOut, 20000, "\\Pearson.RData$")
head(output_df_LeaveOneOut) # Crashes sometimes at first iteration and runs only when data is already loaded in cache

# Second preprocessing step before merging the two dataframes
#processed_regression_models_df_standard <- process_regression_models_df(output_df_standard)
#head(processed_regression_models_df_standard)
processed_regression_models_df_LeaveOneOut <- process_regression_models_df(output_df_LeaveOneOut)
head(processed_regression_models_df_LeaveOneOut)


#combined_df_standard <- merge(df_standard_preprocessed_2, processed_regression_models_df_standard , by = c("gene_name", "segmentation_preselection_corellation"), all = TRUE)

#combined_filtered_df_standard <- combined_df_standard[combined_df_standard$segmentation_preselection_corellation == "Pearson", ]
#head(combined_filtered_df_standard)

combined_df_LeaveOneOut <- merge(df_LeaveOneOut_preprocessed_2, processed_regression_models_df_LeaveOneOut , by = c("gene_name", "segmentation_preselection_corellation"), all = TRUE)

combined_filtered_df_LeaveOneOut <- combined_df_LeaveOneOut[combined_df_LeaveOneOut$segmentation_preselection_corellation == "Pearson", ]
head(combined_filtered_df_LeaveOneOut)


# Filters the data frame to only show rows with missing values in the nzero_value or Pearson columns
# Filters out datapoints where the regression model was created due to error and thus resulted in missin entries after merge
missing_values_df_LeaveOneOut <- combined_filtered_df_LeaveOneOut[is.na(combined_filtered_df_LeaveOneOut$nzero_value) | is.na(combined_filtered_df_LeaveOneOut$Pearson), ]
head(missing_values_df_LeaveOneOut)

combined_filtered2_df_LeaveOneOut <- combined_filtered_df_LeaveOneOut[!is.na(combined_filtered_df_LeaveOneOut$nzero_value) & !is.na(combined_filtered_df_LeaveOneOut$Pearson), ]
head(combined_filtered2_df_LeaveOneOut)

create_dotplot(combined_filtered2_df_LeaveOneOut, output_path)
