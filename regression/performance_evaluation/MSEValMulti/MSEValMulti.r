library(ggplot2)  # Required for creating boxplots
library(RColorBrewer)  # Required for color palettes in boxplots
library(rstudioapi)  # Required for getting the active document context
library(dplyr)  # Required for data manipulation operations

# Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# New location
input_directory_path_LeaveOneOut <- "c:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "c:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/standard_regression/performance_evaluation/Performance_Overview.txt"

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

  print(head(df_combined))

  df_combined$CV_Type <- recode(df_combined$CV_Type, "LeaveOneOut" = "LOOCV")
  df_filtered <- df_combined %>% filter(target < quantile(target, 0.99)) 

  print("Combined data frame")
  print(tail(df_combined))
  boxplot <- ggplot(df_filtered, aes(y = CV_Type, x = target, fill = CV_Type)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = "Cross-Validation Techniques MSE values distribution", x = target_column, y = "Cross Validation Type") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12, angle = 30)) +
    coord_cartesian(xlim = c(0.1, 2)) + 
    scale_fill_manual(values = c("LOOCV" = "#00aeba", "standard" = "#e7b800")) +
    theme(legend.position = "none")

  ggsave(filename = paste(output_directory, "/", target_column, "_boxplot.svg", sep = ""), plot = boxplot, width = 7, height = 4)
  ggsave(filename = paste(output_directory, "/", target_column, "_boxplot.png", sep = ""), plot = boxplot, width = 7, height = 4)

  violin_plot <- ggplot(df_filtered, aes(y = CV_Type, x = target, fill = CV_Type)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 1, fill = "white") +
    labs(title = "Cross-Validation techniques MSE values distribution", x = "Mean squared error", y = "Cross validation type") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 15, color = "black"),
          axis.text = element_text(color = "black", size = 14),
          axis.text.y = element_text(angle = 30)) +
    coord_cartesian(xlim = c(0, 3)) + 
    scale_fill_manual(values = c("LOOCV" = "#00aeba", "standard" = "#e7b800")) +
    theme(legend.position = "none")
  
  ggsave(filename = paste(output_directory, "/", target_column, "_violin_plot.svg", sep = ""), plot = violin_plot, width = 7, height = 4)
  ggsave(filename = paste(output_directory, "/", target_column, "_violin_plot.png", sep = ""), plot = violin_plot, width = 7, height = 4)
      
  return(boxplot)
}


output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# Preprocess the LeaveOneOut CV dataframe
df_LeaveOneOut_preprocessed_1 <- preprocess_performance_evaluation_df(df_LeaveOneOut)
head(df_LeaveOneOut_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
print("LeaveOneOut")
print(head(df_LeaveOneOut_preprocessed_2))

# Preprocess the standard CV dataframe
df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
print("Standard")
print(head(df_standard_preprocessed_2))

create_boxplots(df_LeaveOneOut_preprocessed_2, df_standard_preprocessed_2, target_column = "MSE", output_directory = output_path)
