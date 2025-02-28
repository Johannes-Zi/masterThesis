library(ggplot2)
library(RColorBrewer)
library(rstudioapi)
library(dplyr)

input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
head(df_standard)


preprocess_performance_evaluation_df <- function(input_df) {
  # original filenames for each df entry based on the Sample_Name column
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")  # Split the Sample_Name column by "_"
  #head(filename_parts)  # Print the first few filename parts

  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # entry specific gene names
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Sample_Name column to gene_name
  adapted_df <- rename(adapted_df, gene_name = Sample_Name)

  return(adapted_df)
}


create_violin_plots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  
  df_LeaveOneOut_name <- strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]][2]
  df_standard_name <- strsplit(deparse(substitute(df_standard)), "_")[[1]][2]

  # new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])
  print(head(df_LeaveOneOut_new))
  print(head(df_standard_new))

  # new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)

  violinplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) +
    geom_violin(trim = FALSE, width = 0.5) +
    labs(title = paste(target_column, "across Cross-Validation Techniques", sep = " "), y = target_column, x = "Cross Validation Type") +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold")) +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "none")

  ggsave(filename = paste(output_directory, "/", target_column, "_violinplot.png", sep = ""), plot = violinplot)
  return(violinplot)
}


output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# LeaveOneOut CV dataframe
df_LeaveOneOut_preprocessed_1 <- preprocess_performance_evaluation_df(df_LeaveOneOut)
head(df_LeaveOneOut_preprocessed_1)
# df entries which were based on Pearson based feature selection at the end of the segementation
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_LeaveOneOut_preprocessed_2)

# standard CV dataframe
df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
# df entries which were based on Pearson based feature selection at the end of the segementation
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_standard_preprocessed_2)

create_violin_plots(df_LeaveOneOut_preprocessed_2, df_standard_preprocessed_2, target_column = "Pearson", output_directory = output_path)
