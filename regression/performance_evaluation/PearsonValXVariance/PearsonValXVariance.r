
library(ggplot2)  # For creating plots and visualizations
library(RColorBrewer)  # For color palettes in plots
library(rstudioapi)  # For interacting with RStudio IDE
library(dplyr)  # For data manipulation and transformation
library(MASS)  # For kernel density estimation
library(viridis)  # For color scales in plots

input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_standard)


preprocess_performance_evaluation_df <- function(input_df) {
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")  # Split the Sample_Name column by "_"
  #head(filename_parts)  # Print the first few filename parts

  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

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


get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}


create_dotplot <- function(df_standard, output_directory) {
    df_standard_name <- strsplit(deparse(substitute(df_standard)), "_")[[1]][2]
    
    # new data frames that only contain the target columns
    df_standard_new <- data.frame(Pearson = df_standard$Pearson, PearsonVar = df_standard$PearsonVar)

    # IQR for Pearson and PearsonVar
    Q1_Pearson <- quantile(df_standard_new$Pearson, 0.25)
    Q3_Pearson <- quantile(df_standard_new$Pearson, 0.75)
    IQR_Pearson <- IQR(df_standard_new$Pearson)

    Q1_PearsonVar <- quantile(df_standard_new$PearsonVar, 0.25)
    Q3_PearsonVar <- quantile(df_standard_new$PearsonVar, 0.75)
    IQR_PearsonVar <- IQR(df_standard_new$PearsonVar)

    # outliers
    df_standard_new <- df_standard_new[!(df_standard_new$Pearson < (Q1_Pearson - 1.5 * IQR_Pearson) | 
                                                                             df_standard_new$Pearson > (Q3_Pearson + 1.5 * IQR_Pearson) |
                                                                             df_standard_new$PearsonVar < (Q1_PearsonVar - 1.5 * IQR_PearsonVar) |
                                                                             df_standard_new$PearsonVar > (Q3_PearsonVar + 1.5 * IQR_PearsonVar)), ]
        
    # density for each point
    # n specifies the number of points to evaluate the density with/ defines number of grid cells
    df_standard_new$density <- get_density(df_standard_new$Pearson, df_standard_new$PearsonVar, n = 1000)

    dotplot <- ggplot(df_standard_new, aes(x = Pearson, y = PearsonVar, color = density)) +
        geom_point() +
        scale_color_viridis() +
        labs(title = "Pearson vs PearsonVar", x = "Pearson", y = "PearsonVar") +
        theme_gray() +
        theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
                    axis.title = element_text(size = 14, face = "bold"),
                    axis.text = element_text(size = 12, face = "bold"))
    
    ggsave(filename = paste(output_directory, "/dotplot.png", sep = ""), plot = dotplot)
    
    return(dotplot)
}

output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

df_standard_preprocessed_1 <- preprocess_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
# Select df entries which were based on Pearson based feature selection at the end of the segementation
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_standard_preprocessed_2)

create_dotplot(df_standard_preprocessed_2, output_directory = output_path)
 