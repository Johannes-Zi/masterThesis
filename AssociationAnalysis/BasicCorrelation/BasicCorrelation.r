library(qvalue)
library(ggplot2)
library(crayon)
library(dplyr)

if (FALSE){
  # To be ignored
  # create_new_library <- function() {
  # new_library_path <- "C:/Users/johan/Desktop/local_master_thesis_data/temp_library/"
  # dir.create(new_library_path, recursive = TRUE)

  # # Add to PATH
  # .libPaths(c(new_library_path, .libPaths()))

  # if (!require("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager", lib = new_library_path)

  # BiocManager::install("qvalue", lib = new_library_path, force = TRUE)
  # }
  #create_new_library()



  #if (!require("BiocManager", quietly = TRUE))
  #    install.packages("BiocManager")

  #BiocManager::install("qvalue")

  #install.packages("devtools")
  #library("devtools")
  #install_github("jdstorey/qvalue")
}

import_clinical_meatadata <- function(file_path) {
  clinical_metadata <- read.csv(file_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  return(clinical_metadata)
}


import_elnet_segments_atac_data <- function(file_path) {
  segments_atac_data <- read.csv(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(segments_atac_data)
}


extract_represented_patient_ids <- function(elnet_segments_atac_data) {
  columnnames <- colnames(elnet_segments_atac_data)

  atac_columns <- grep(".ATAC", columnnames, value = TRUE)
  
  patient_ids <- gsub(".ATAC", "", atac_columns)

  return(patient_ids)
}

reduced_clinical_metadata <- function(clinical_metadata, patient_ids) {
  clinical_metadata$rna.final.ids <- gsub("-", ".", clinical_metadata$rna.final.ids)

  reduced_clinical_metadata <- clinical_metadata[clinical_metadata$rna.final.ids %in% patient_ids, ]

  rownames(reduced_clinical_metadata) <- 1:nrow(reduced_clinical_metadata)

  return(reduced_clinical_metadata)
}

create_spearman_correlations_df <- function(clinical_vector, elnet_segments_atac_data){

  pb <- txtProgressBar(min = 0, max = nrow(elnet_segments_atac_data), style = 3)

  spearman_correlations_df <- data.frame(segment = character(), gene_id = character(), spearman_corr = numeric(),
                                         p_value = numeric(), stringsAsFactors = FALSE)

  # elnet segments atac data
  for (i in 1:nrow(elnet_segments_atac_data)){
    current_segment_atac <- elnet_segments_atac_data[i,]

    current_segment <- current_segment_atac$segment
    current_gene_id <- current_segment_atac$gene_id

    current_segment_atac <- current_segment_atac[-c(1,2)]

    # transposes the elnet segments atac data
    current_segment_atac <- as.data.frame(t(current_segment_atac))

    rownames(current_segment_atac) <- gsub(".ATAC", "", rownames(current_segment_atac))
    
    colnames(current_segment_atac) <- "ATAC"
    
    current_segment_atac$ids <- rownames(current_segment_atac)

    # sorts by ids
    current_segment_atac <- current_segment_atac[order(current_segment_atac$ids),]
    
    current_segment_atac <- current_segment_atac[-2]

    atac_vector <- as.numeric(current_segment_atac$ATAC)

    # if (current_segment == "chr7.26617805.26617894"){
    #   print("clinical_vector\n")
    #   cat(paste(clinical_vector, ","))
    #   print("\n")
    #   print("atac_vector\n")
    #   cat(paste(atac_vector, ","))
    #   print("\n")
    # }

    # calculates spearman
    spearman_correlation <- suppressWarnings(cor.test(clinical_vector, atac_vector, method = "spearman"))

    # adds new row to the spearman_correlations_df
    spearman_correlations_df <- rbind(spearman_correlations_df, data.frame(segment = current_segment, 
                                      gene_id = current_gene_id, spearman_corr = spearman_correlation$estimate,
                                      p_value = spearman_correlation$p.value, stringsAsFactors = FALSE))

    rownames(spearman_correlations_df) <- 1:nrow(spearman_correlations_df)

    setTxtProgressBar(pb, i)
  }

  close(pb)

  spearman_correlations_df$q_value <- qvalue(p = spearman_correlations_df$p_value)$qvalues

  return(spearman_correlations_df)
}

perform_spearman_correlation <- function(output_path, clinical_metadata_rna_ids, clinical_metadata_column, elnet_segments_atac_data, clinical_metadata_column_name) {

  # output direcotry based on the handed over output path and the clinical metadata column name
  output_dir <- paste(output_path, clinical_metadata_column_name, sep = "/")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # data and df with the handed over column of the clinical metadata and the RNA ids
  clinical_metadata_column_df <- data.frame(clinical_metadata_rna_ids, clinical_metadata_column)

  clinical_metadata_column_df <- clinical_metadata_column_df[order(clinical_metadata_column_df$clinical_metadata_rna_ids),]

  rownames(clinical_metadata_column_df) <- clinical_metadata_column_df$clinical_metadata_rna_ids

  clinical_metadata_column_df <- clinical_metadata_column_df[-1]

  clinical_vector <- as.numeric(gsub(",", ".", clinical_metadata_column_df[, 1]))

  spearman_correlations_df <- create_spearman_correlations_df(clinical_vector = clinical_vector, elnet_segments_atac_data = elnet_segments_atac_data)
  


  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "spearman_correlations.tsv", sep = "_"), sep = "/")
  write.table(spearman_correlations_df, output_file_path, row.names = FALSE, sep = ";")
  print(paste("->", output_file_path))

  # histogram plot of the p_values
  p_value_hist <- ggplot(spearman_correlations_df, aes(x = p_value)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(limits = c(0, 1), oob = scales::oob_squish) +
  ggtitle("Distribution of p_values") +
  xlab("p_value") +
  theme_light()
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "p_value_hist.png", sep = "_"), sep = "/")
  ggsave(output_file_path, p_value_hist)
  print(paste("->", output_file_path))

  # histogram plot of the q_values
  q_value_hist <- ggplot(spearman_correlations_df, aes(x = q_value)) +
  geom_histogram(bins = 100) +
  scale_x_continuous(limits = c(0, 1), oob = scales::oob_squish) +
  ggtitle("Distribution of q_values") +
  xlab("q_value") +
  theme_light()
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "q_value_hist.png", sep = "_"), sep = "/")
  ggsave(output_file_path, q_value_hist)
  print(paste("->", output_file_path))

  # density plot of the spearman_correlations
  spearman_corr_density <- ggplot(spearman_correlations_df, aes(x = spearman_corr)) + 
    geom_density(fill = "#4dbf67", alpha = 0.4) +
    xlim(-1, 1) +
    ggtitle("Density plot of spearman_correlations") +
    xlab("spearman_corr") +
    theme_light()
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "spearman_corr_density.png", sep = "_"), sep = "/")
  ggsave(output_file_path, spearman_corr_density)
  print(paste("->", output_file_path))

  # density plot of the spearman_correlations over the p_value < 0.05
  # filters the dataframe to include only rows where spearman_corr is below 0.05
  spearman_correlations_filtered_df <- spearman_correlations_df %>% filter(p_value < 0.05)
  spearman_corr_density <- ggplot(spearman_correlations_filtered_df, aes(x = spearman_corr)) + 
    geom_density(fill = "#4dbf67",alpha = 0.4) + 
    xlim(-1, 1) +
    ggtitle("Density distribution of Spearman correlations with p_value < 0.05") +
    xlab("spearman_corr") +
    theme_light()
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "spearman_corr_density_smaller_than_005_Spearman.png", sep = "_"), sep = "/")
  ggsave(output_file_path, spearman_corr_density)
  print(paste("->", output_file_path))

  # density plot of the spearman_correlations seperated by p_value < 0.05
  spearman_corr_density_p_value <- ggplot(spearman_correlations_df, aes(x = spearman_corr, fill = p_value < 0.05)) +
    xlim(-1, 1) +
    geom_density(alpha = 0.5) + ggtitle("Density plot of spearman_correlations")
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "spearman_corr_density_p_value.png", sep = "_"), sep = "/")
  ggsave(output_file_path, spearman_corr_density_p_value)
  print(paste("->", output_file_path))

  return(spearman_correlations_df)
}

iterate_over_clinical_metadata <- function(output_path, clinical_metadata, elnet_segments_atac_data){
  

  # empty list to store the spearman_correlations_df_final of each clinical metadata column
  spearman_correlations_df_final_list <- list()

  # columns of the clinical metadata and perform spearman correlation and save output files
  for (i in 4:ncol(clinical_metadata)){
    #print(colnames(clinical_metadata))

    current_column <- clinical_metadata[,i]

    message(paste(bold(cyan("\nPerform spearman correlations for:")), magenta(colnames(clinical_metadata)[i])))
    current_spearman_correlations_df_final <- perform_spearman_correlation(output_path = output_path, clinical_metadata$rna.final.ids, current_column, elnet_segments_atac_data, clinical_metadata_column_name=colnames(clinical_metadata)[i])
    
    spearman_correlations_df_final_list[[i]] <- current_spearman_correlations_df_final
  }

  return(spearman_correlations_df_final_list)
}

load_performance_evaluation <- function(LeaveOneOutCV_performance_path) {
  LeaveOneOutCV_performance <- read.table(LeaveOneOutCV_performance_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  LeaveOneOutCV_performance$gene_name <- sapply(strsplit(LeaveOneOutCV_performance$Sample_Name, "_"), function(x) x[2])

  LeaveOneOutCV_performance <- LeaveOneOutCV_performance[grepl("Pearson", LeaveOneOutCV_performance$Sample_Name),]

  return(LeaveOneOutCV_performance)
}





if (TRUE) {

  #
  # Import required data
  if (TRUE) {
    cat(paste(bold(cyan("\nImport required data\n"))))
    # Import clinical metadata
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)
    message(paste("Imported clinical metadata from: \n", clinical_metadata_path))

    # Save the clinical metadata as csv file with seperator ;
    clinical_metadata_output_path <- paste(getwd(), "AssociationAnalysis/BasicCorrelation/", "imported_clinical_metadata.csv", sep = "/")
    write.table(clinical_metadata, clinical_metadata_output_path, row.names = FALSE, sep = ";")
    message(paste("Saved clinical metadata as csv file to: \n", clinical_metadata_output_path))

    # Import elnet segments atac data
    elnet_segments_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/CollapsedSegmentation/elnet_model_segments_df.tsv"
    elnet_segments_atac <- import_elnet_segments_atac_data(elnet_segments_path)
    message(paste("Imported elnet segments atac data from: \n", elnet_segments_path))

    # Import gene model specific outer cross validation performance estimations
    LeaveOneOutCV_performance_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.tsv"
    LeaveOneOutCV_performance <- load_performance_evaluation(LeaveOneOutCV_performance_path)  # Load the data
    message(paste("Imported gene model specific outer cross validation performance estimations from: \n", LeaveOneOutCV_performance_path))

    # Density plot of the Pearson correlations
    ggplot(LeaveOneOutCV_performance, aes(x = Pearson)) + geom_density(fill = "#4dbf67",alpha = 0.4) + ggtitle("Density plot of Pearson correlations") + xlab("Pearson") + theme_light()
    # Save the plot at cwd
    image_path <- paste(getwd(), "AssociationAnalysis/BasicCorrelation/combined_correlations", "outer_cross_validation_pearson_density_plot.png", sep = "/")
    ggsave(image_path)
  }

  #
  # Reduce clinical metadata to the patient ids which are represented in the elnet segments atac data
  if (TRUE) {
    message(paste(bold(cyan("\nReduce clinical metadata to the patient ids which are represented in the elnet segments atac data"))))
    # Extract represented patient ids in the elnet segments atac data
    patient_ids <- extract_represented_patient_ids(elnet_segments_atac)
    message(paste("Extracted represented patient ids in the elnet segments atac data: ", magenta(length(patient_ids))))

    # Print how many unique patient ids are represented in the clinical metadata with RNA Seq ids
    message(paste("Number of patients represented in the clinical metadata with RNA Seq ids: ", magenta(length(unique(clinical_metadata$rna.final.ids)))))
    # Reduce clinical metadata to the patient ids which are represented in the elnet segments atac data
    reduced_clinical_metadata <- reduced_clinical_metadata(clinical_metadata, patient_ids)
    message(paste("Patients represented in the reduced clinical metadata: ", magenta(nrow(reduced_clinical_metadata))))
  }

  #
  # List of genes which have a Pearson correlation above a defined threshold
  if (TRUE) {
    thresx <- 0.4 # Defines the threshold for the outer cross validation threshold to filter out the gene models below the threshold
    cat(paste(bold(cyan("\nApply elnet outer cross validation based pearson correlation filtering\n"))))
    cat(paste( "with threshold Pearson >=", red(thresx), "or <=", red(-thresx), "\n"))
    cat(paste("Number of unique genes represented in LeaveOneOutCV_performance: ", magenta(length(unique(LeaveOneOutCV_performance$gene_name)), "\n")))
    cat(paste("Number of entries with negative Pearson correlation: ", magenta(nrow(LeaveOneOutCV_performance[LeaveOneOutCV_performance$Pearson < 0,]), "\n")))
    # Filter out the genes with a amount of Pearson correlation below the threshold
    genes_with_cv_pearson_correlation_over_thresx <- LeaveOneOutCV_performance[LeaveOneOutCV_performance$Pearson >= thresx | LeaveOneOutCV_performance$Pearson <= -thresx, ]$gene_name
    cat(paste("Number of genes after filtering: ", magenta(length(genes_with_cv_pearson_correlation_over_thresx), "\n")))
  }

  #
  # Filter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list
  if (TRUE) {
    cat(paste(bold(cyan("\nFilter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list.\n"))))
    cat(paste("Number of unique genes represented by elnet segments: ", magenta(length(unique(elnet_segments_atac$gene_id)), "\n")))
    cat(paste("Number of elnet segments before filtering: ", magenta(nrow(elnet_segments_atac), "\n")))
    # Filter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list
    elnet_segments_atac_filtered <- elnet_segments_atac[elnet_segments_atac$gene_id %in% genes_with_cv_pearson_correlation_over_thresx,]
    cat(paste("Number of elnet segments after filtering: ", magenta(nrow(elnet_segments_atac_filtered), "\n")))
  }

  #
  # Run the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over 
  # all clinical parameters
  if (TRUE) {
    cat(paste(bold(cyan("\nRun the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over all clinical parameters.\n"))))
    # output_path based on cwd
    output_path <- paste(getwd(), "AssociationAnalysis/BasicCorrelation/combined_correlations", sep = "/")
    # Run Correlations with filtered elnet segments atac data
    list_with_all_correlation_dfs_filtered <- iterate_over_clinical_metadata(output_path = output_path, reduced_clinical_metadata, elnet_segments_atac_filtered)
  }
  
  #sink()
}
