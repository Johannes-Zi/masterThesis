library(ggplot2)
library(RColorBrewer)
library(rstudioapi)
library(dplyr)
library(MASS)
library(viridis)

input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_standard)

preprocess_performance_evaluation_df <- function(input_df) {
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")
  #head(filename_parts)  # Print the first few filename parts

  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  #head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

  # new column in the dataframe for the segmentation preselection correlation type
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # entry specific gene names
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
    
    df_standard_new <- data.frame(Spearman = df_standard$Spearman, SpearmanVar = df_standard$SpearmanVar)

    # IQR for Spearman and SpearmanVar
    Q1_Spearman <- quantile(df_standard_new$Spearman, 0.25)
    Q3_Spearman <- quantile(df_standard_new$Spearman, 0.75)
    IQR_Spearman <- IQR(df_standard_new$Spearman)

    Q1_SpearmanVar <- quantile(df_standard_new$SpearmanVar, 0.25)
    Q3_SpearmanVar <- quantile(df_standard_new$SpearmanVar, 0.75)
    IQR_SpearmanVar <- IQR(df_standard_new$SpearmanVar)

    # outliers
    df_standard_new <- df_standard_new[!(df_standard_new$Spearman < (Q1_Spearman - 1.5 * IQR_Spearman) | 
                                            df_standard_new$Spearman > (Q3_Spearman + 1.5 * IQR_Spearman) |
                                            df_standard_new$SpearmanVar < (Q1_SpearmanVar - 1.5 * IQR_SpearmanVar) |
                                            df_standard_new$SpearmanVar > (Q3_SpearmanVar + 1.5 * IQR_SpearmanVar)), ]
        
    # n specifies the number of points to evaluate the density with/ defines number of grid cells
    df_standard_new$density <- get_density(df_standard_new$Spearman, df_standard_new$SpearmanVar, n = 1000)

    dotplot <- ggplot(df_standard_new, aes(x = Spearman, y = SpearmanVar, color = density)) +
        geom_point() +
        scale_color_viridis() +
        labs(title = "Spearman vs SpearmanVar", x = "Spearman", y = "SpearmanVar") +
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
# df entries which were based on Pearson based feature selection at the end of the segementation
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_standard_preprocessed_2)

create_dotplot(df_standard_preprocessed_2, output_directory = output_path)
#create_dotplot(df_LeaveOneOut_preprocessed_2, output_directory = output_path)
 