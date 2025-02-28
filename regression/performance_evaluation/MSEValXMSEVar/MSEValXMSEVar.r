# Load necessary libraries
library(ggplot2)  # For creating plots and visualizations
library(RColorBrewer)  # For color palettes in plots
library(rstudioapi)  # For interacting with RStudio IDE
library(dplyr)  # For data manipulation and transformation
library(MASS)  # For kernel density estimation
library(viridis)  # For color scales in plots

# Import Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
head(df_standard)


preprocess_performance_evaluation_df <- function(input_df) {
  # original filenames for each df entry based on the Sample_Name column
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")  # Split the Sample_Name column by "_"
  #head(filename_parts)  # Print the first few filename parts

  # segmentation preselection correlation
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

  # new column in the dataframe for the segmentation preselection correlation
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # entry specific gene names
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Sample_Name column to gene_name
  adapted_df <- rename(adapted_df, gene_name = Sample_Name)

  return(adapted_df)
}

get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}


create_dotplot <- function(df_standard, output_directory) {
    # x axis names
    df_standard_name <- strsplit(deparse(substitute(df_standard)), "_")[[1]][2]
    
    # new data frames that only contain the target columns
    df_standard_new <- data.frame(MSE = df_standard$MSE, MSEVar = df_standard$MSEVar)

    # IQR for MSE and MSEVar
    Q1_MSE <- quantile(df_standard_new$MSE, 0.25)
    Q3_MSE <- quantile(df_standard_new$MSE, 0.75)
    IQR_MSE <- IQR(df_standard_new$MSE)

    Q1_MSEVar <- quantile(df_standard_new$MSEVar, 0.25)
    Q3_MSEVar <- quantile(df_standard_new$MSEVar, 0.75)
    IQR_MSEVar <- IQR(df_standard_new$MSEVar)

    # outliers
    df_standard_new <- df_standard_new[!(df_standard_new$MSE < (Q1_MSE - 1.5 * IQR_MSE) | 
                                                                             df_standard_new$MSE > (Q3_MSE + 1.5 * IQR_MSE) |
                                                                             df_standard_new$MSEVar < (Q1_MSEVar - 1.5 * IQR_MSEVar) |
                                                                             df_standard_new$MSEVar > (Q3_MSEVar + 1.5 * IQR_MSEVar)), ]
        
    # density for each point
    df_standard_new$density <- get_density(df_standard_new$MSE, df_standard_new$MSEVar, n = 100)

    dotplot <- ggplot(df_standard_new, aes(x = MSE, y = MSEVar, color = density)) +
        geom_point() +
        scale_color_viridis() +
        labs(title = "MSE vs MSEVar", x = "MSE", y = "MSEVar") +
        theme_gray() +
        theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                    axis.title = element_text(size = 14, face = "bold"),
                    axis.text = element_text(size = 12, face = "bold"))
    
    ggsave(filename = paste(output_directory, "/dotplot.png", sep = ""), plot = dotplot)
    
    return(dotplot)
}


# output path to the directory of the active document
output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# Performance_Evaluation dataframes
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

#create_dotplot(df_standard_preprocessed_2, output_directory = output_path)
create_dotplot(df_LeaveOneOut_preprocessed_2, output_directory = output_path)
 