library(ggplot2)
library(dplyr)
library(tidyr)

load_data <- function(input_dir) {
  combined_df <- data.frame()

  for (dir in list.dirs(input_dir, full.names = FALSE, recursive = FALSE)) {
    clinical_parameter <- dir

    # path to the input file
    input_file <- paste(input_dir, dir, "/", dir, "_spearman_correlations.tsv", sep = "")

    # input tsv file to a data frame
    input_df <- read.table(input_file, header = TRUE, sep = ";")

    # clinical parameter to the data frame
    input_df$clinical_parameter <- clinical_parameter

    # data frame to the combined data frame
    combined_df <- rbind(combined_df, input_df)
  }

  # drops clinical parameters that have no meaningful correlations
  clinical_params_to_drop <- c("age", "height", "weight", "NYHA", "ZVD", "heart.rate", "Paradoxe_Septumbewegung", "Perikarderguss", "SMW")

  # drps rows with the clinical parameters that have no meaningful correlations
  combined_df <- combined_df[!combined_df$clinical_parameter %in% clinical_params_to_drop, ]
  rownames(combined_df) <- NULL


  # Rename parameters
  combined_df$clinical_parameter <- gsub("age", "age", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("\\bAT\\b", "pulmonary acceleration time", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("AT_ET", "AT / ET", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("AVDO2", "arteriovenous oxygen difference", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("BNP", "BNP heart failure marker", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Cardiac.Index", "caridac performance index", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("CVP", "central venous pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("DLCO", "lung CO diffusion capacity", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("eGFR", "glomerular filtration rate", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("heart.rate", "heart rate", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("height", "height", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("HZV_Fick", "cardiac output - Fick principle", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("HZV_Thermodil", "cardiac output - thermodilution", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Kreatinin", "kreatinin", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("mPAP", "mean pulmonary artery pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("NYHA", " NYHA heart failure classification", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PA_diastolic", "diastolic pulmonary artery pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PA_systolic", "systolic pulmonary artery pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Paradoxe_Septumbewegung", "ventricular septum movement", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PAWP", "pulmonary capillary occlusion pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Perikarderguss", "pericardium fluid accumulation", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PVR", "pulmonary vascular resistance", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("RA_area", "right atrium size", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Rrsys", "systolic blood pressure at rest", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("RVEDD", "right ventricle size at diastole", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("\\bS\\b", "TAPSV impaired ventricular contraction", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("sPAP.excl..ZVD", "est. systolic pulmonary artery pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("TAPSE", "TAPSE impaired ventricular contraction", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Tiffeneau.Index.FEV1.VC", "airflow obstruction in lung Tiffeneau", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("ZVD", "estimated central venous pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("VCI_diameter", "inferior vena cava diameter", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("ven_SO2", "venous oxygen saturation", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("\\bFEV1\\b", "one second forced expiratory volume", combined_df$clinical_parameter)


  return(combined_df)
}


# plots the number of df entries for each clinical parameter
plot_entries <- function(input_df, output_dir) {
  # plots the number of entries for each clinical parameter
  ggplot(input_df, aes(x = clinical_parameter)) +
    geom_bar() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = "Number of top segments for each clinical parameter",
         x = "Clinical parameter",
         y = "Number of segments")

  ggsave(paste(output_dir, "number_of_entries_per_param.png", sep = ""))

  write.table(as.data.frame(table(input_df$clinical_parameter)), paste(output_dir, "number_of_entries_per_param.csv", sep = ""), sep = ",", row.names = FALSE)
}

# barchart with the most abundand gene_id entries across all clinical parameters
plot_top_gene_ids <- function(input_df, output_dir) {
  # a data frame with the number of entries for each gene_id
  gene_id_counts <- as.data.frame(table(input_df$gene_id))

  # data frame by the number of entries
  gene_id_counts <- gene_id_counts[order(-gene_id_counts$Freq), ]

  rownames(gene_id_counts) <- NULL
  
  # ensure the Var1 column is a factor with levels ordered by Freq
  gene_id_counts$Var1 <- factor(gene_id_counts$Var1, levels = gene_id_counts$Var1[order(gene_id_counts$Freq, decreasing = TRUE)])

  ggplot(gene_id_counts[1:40, ], aes(y = Var1, x = Freq)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
        title = "Top 10 most abundant gene_ids",
        y = "Gene ID",
        x = "Number of entries"
    )

  ggsave(paste(output_dir, "top_gene_ids.png", sep = ""))

  frequency_distribution <- as.data.frame(table(gene_id_counts$Freq))
  colnames(frequency_distribution) <- c("Frequency", "Number_of_Genes")

  # Frequency to numeric for proper plotting
  frequency_distribution$Frequency <- as.numeric(as.character(frequency_distribution$Frequency))

  ggplot(frequency_distribution, aes(x = Frequency, y = Number_of_Genes)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.background = element_rect(fill = "white", color = NA),  # Add white background to the plot
    ) +
    labs(
      title = "Distribution of Gene occurrences",
      x = "number of clinical parameters a gene is represented as top correlation",
      y = "Number of Genes"
    )

  ggsave(paste(output_dir, "gene_frequency_distribution.png", sep = ""))

  write.table(gene_id_counts, paste(output_dir, "gene_id_counts.tsv", sep = ""), sep = "\t", row.names = FALSE)
  
  print(paste(gene_id_counts$Var1, collapse = ","))
  }

reorder_percentage_matrix <- function(percentage_matrix) {
  # percentages between variables as distance
  dd <- as.dist((100 - percentage_matrix) / 100)
  hc <- hclust(dd)
  percentage_matrix <- percentage_matrix[hc$order, hc$order]
  return(percentage_matrix)
}

# heatmap of the gene_id that pairs of the clinical parameters share
clinical_parameter_gene_assotiations_heatmap <- function(input_df, output_dir) {
  cross_tab <- table(input_df$clinical_parameter, input_df$gene_id)

  # calculates the co-occurrences of gene IDs for each pair of clinical parameters
  # by performing matrix multiplication between the cross-tabulation matrix and its transpose
  co_occurrences <- cross_tab %*% t(cross_tab)

  total_genes <- ncol(cross_tab)
  message("Total genes represented in heatmap dataframe: ", total_genes)
  percentage_matrix <- (co_occurrences / total_genes) * 100

  percentage_matrix <- reorder_percentage_matrix(percentage_matrix)

  percentage_df <- as.data.frame(as.table(percentage_matrix))
  colnames(percentage_df) <- c("Clinical_Parameter_X", "Clinical_Parameter_Y", "Percentage")

  # reversed order of the y-axis by reordering the factor levels
  percentage_df$Clinical_Parameter_Y <- factor(percentage_df$Clinical_Parameter_Y, levels = rev(levels(percentage_df$Clinical_Parameter_Y)))

  # NA for tiles above the diagonal
  percentage_df <- percentage_df %>%
    mutate(Percentage = ifelse(as.numeric(factor(Clinical_Parameter_X, levels = levels(Clinical_Parameter_Y))) < as.numeric(Clinical_Parameter_Y), NA, Percentage))

  ggplot(percentage_df, aes(x = Clinical_Parameter_X, y = Clinical_Parameter_Y, fill = Percentage)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#0007c4", na.value = "white") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(angle = 45, color = "black"),
      plot.title = element_text(hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = c(0.9, 0.8)
    ) +
    labs(title = "Co-occurrence of Gene IDs Across Clinical Parameters",
        x = NULL,
        y = NULL,
        fill = "Gene-set \noverlap\npercentage")

  ggsave(paste(output_dir, "clinical_parameter_gene_associations_heatmap.png", sep = ""), dpi = 200, height = 8, width = 8) 
  
  }


# exports thefiltered segments to a bed file
export_filtered_segments <- function(input_df, output_dir) {

  output_df <- input_df
  # combines the columns gene_id, spearman_corr and p_value to a single column called name and seperate the values with a underscore
  output_df$geneSpearmanPval <- paste(output_df$gene_id, output_df$spearman_corr, output_df$p_value, sep = "_")

  # removes the columns gene_id, spearman_corr and p_value
  output_df <- output_df[, !(names(output_df) %in% c("gene_id", "spearman_corr", "p_value", "clinical_parameter"))]
  
  # separates the segment column into chrom, chromStart, and chromEnd
  output_df <- separate(output_df, segment, into = c("chr", "start", "end"), sep = "\\.")
  
  print(head(output_df))

  output_df <- output_df %>% distinct(chr, start, end, .keep_all = TRUE)


  file_path <- paste(output_dir, "filtered_segments.bed", sep = "/")
  writeLines(paste("#", paste(colnames(output_df), collapse = "\t"), sep = ""), con = file_path)
  write.table(output_df, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

if (TRUE) {
  input_dir = "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/CorrelationVisualizationFiltered/runs/v1/"

  input_df <- load_data(input_dir)

  output_dir = "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/CorrelationVisualizationFiltered/visulizations/v1/"

  # number of segments for each clinical parameter
  plot_entries(input_df, output_dir)

  # umber of Gene IDs across all clinical parameters
  plot_top_gene_ids(input_df, output_dir)

  # reduced df version without the columns segement, spearman_corr and p_value
  reduced_df <- input_df[, c("clinical_parameter", "gene_id")]

  reduced_df <- unique(reduced_df)

  # heatmap of the gene_id that pairs of the clinical parameters share
  clinical_parameter_gene_assotiations_heatmap(reduced_df, output_dir) 
  export_filtered_segments(input_df, output_dir)

}
