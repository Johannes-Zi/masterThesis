library(ggplot2)
library(RColorBrewer)
library(rstudioapi)

input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
head(df_standard)

create_boxplots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  
  df_LeaveOneOut_name <- tail(strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]], 1)
  df_standard_name <- tail(strsplit(deparse(substitute(df_standard)), "_")[[1]], 1)
  
  # new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])
  
  # lower and upper whiskers for both dataframes
  ylim1 = boxplot.stats(df_LeaveOneOut[[target_column]])$stats[c(1, 5)]
  ylim2 = boxplot.stats(df_standard[[target_column]])$stats[c(1, 5)]

  # minimum lower whisker and maximum upper whisker as the y limits
  ylim = c(min(ylim1[1], ylim2[1]), max(ylim1[2], ylim2[2]))

  # new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)
  
  boxplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = paste(target_column, "Comparison", sep = " "), y = target_column, x = "CV Type") +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold")) +
    coord_cartesian(ylim = ylim*1.1) +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "none")

  ggsave(filename = paste(output_directory, "/", target_column, "_violinplot.png", sep = ""), plot = boxplot)

  return(boxplot)
}

create_violin_plots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  
  df_LeaveOneOut_name <- tail(strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]], 1)
  df_standard_name <- tail(strsplit(deparse(substitute(df_standard)), "_")[[1]], 1)
  
  # new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])

  # new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)
  
  violinplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) +
    geom_violin(trim = FALSE) +
    labs(title = paste(target_column, "Comparison", sep = " "), y = target_column, x = "CV Type") +
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

create_violin_plots(df_LeaveOneOut, df_standard, target_column = "Pearson", output_directory = output_path)

create_violin_plots(df_LeaveOneOut, df_standard, target_column = "Spearman", output_directory = output_path)

create_boxplots(df_LeaveOneOut, df_standard, target_column = "MSE", output_directory = output_path)

#create_violin_plots(df_LeaveOneOut, df_standard, target_column = "pVal")

#create_violin_plots(df_LeaveOneOut, df_standard, target_column = "qVal")

