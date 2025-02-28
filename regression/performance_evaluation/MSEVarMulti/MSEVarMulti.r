library(ggplot2)  # Required for creating boxplots
library(RColorBrewer)  # Required for color palettes in boxplots
library(rstudioapi)  # Required for getting the active document context
library(dplyr)  # Required for data manipulation operations

# Import Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# New location
input_directory_path_LeaveOneOut <- "c:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "c:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
print("LeaveOneOut")
print(head(df_LeaveOneOut))
print("Standard")
print(tail(df_standard))


preprocess_performance_evaluation_df <- function(input_df) {
  # original filenames for each df entry based on the Sample_Name column
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")  # Split the Sample_Name column by "_"
  #head(filename_parts)  # Print the first few filename parts

  # segmentation preselection correlation type for each df entry
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

  # new column in the dataframe for the segmentation preselection correlation type
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # entry specific gene names
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Sample_Name column to gene_name
  adapted_df <- rename(adapted_df, gene_name = Sample_Name)

  return(adapted_df)
}


create_boxplots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  # x axis names
  df_LeaveOneOut_name <- strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]][2]
  df_standard_name <- strsplit(deparse(substitute(df_standard)), "_")[[1]][2]

  # new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])

  # lower and upper whiskers for both data frames
  ylim1 = boxplot.stats(df_LeaveOneOut[[target_column]])$stats[c(1, 5)]
  ylim2 = boxplot.stats(df_standard[[target_column]])$stats[c(1, 5)]

  # minimum lower whisker and maximum upper whisker as the y limits
  ylim = c(min(ylim1[1], ylim2[1]), max(ylim1[2], ylim2[2]))

  # new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)

  print("Combined data frame")
  print(tail(df_combined))
  boxplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste(target_column, "across Cross-Validation Techniques", sep = " "), y = target_column, x = "Cross Validation Type") +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold")) +
    coord_cartesian(ylim = ylim*1.1) +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "none")

  ggsave(filename = paste(output_directory, "/", target_column, "_boxplot.png", sep = ""), plot = boxplot)

  return(boxplot)
}


output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

df_LeaveOneOut_preprocessed_1 <- preprocess_performance_evaluation_df(df_LeaveOneOut)
head(df_LeaveOneOut_preprocessed_1)
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
print("LeaveOneOut")
print(head(df_LeaveOneOut_preprocessed_2))

df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
print("Standard")
print(head(df_standard_preprocessed_2))

create_boxplots(df_LeaveOneOut_preprocessed_2, df_standard_preprocessed_2, target_column = "MSEVar", output_directory = output_path)
