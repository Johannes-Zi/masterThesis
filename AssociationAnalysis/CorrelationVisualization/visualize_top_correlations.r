library(qvalue)
library(ggplot2)
library(crayon)
library(dplyr)


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


load_filtered_correlations_gene_ids <- function(filtered_correlations_top_directory_path) {
  directories <- list.dirs(filtered_correlations_top_directory_path, full.names = FALSE, recursive = FALSE)

  filtered_correlations_gene_ids <- c()
  for (directory in directories) {
    results_file_path <- paste(filtered_correlations_top_directory_path, paste(directory, paste(directory, "_reduced_correlation_results.csv", sep = ""), sep = "/"), sep = "")
    results_df <- read.csv(results_file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)

    gene_ids <- results_df$gene_id

    filtered_correlations_gene_ids <- c(filtered_correlations_gene_ids, gene_ids)

    filtered_correlations_gene_ids <- unique(filtered_correlations_gene_ids)
  }

  return(filtered_correlations_gene_ids)
}


process_and_sort_files_in_directory <- function(directory_path) {
  filenames <- list.files(directory_path, full.names = TRUE)
  
  sorted_filenames <- sort(filenames)
  
  for (filename in sorted_filenames) {
    # renames the file by adding a leading "X" to the filename
    new_filename <- file.path(dirname(filename), paste0("X", basename(filename)))

    file.copy(filename, new_filename)
    
    file.remove(filename)
  }
}


create_spearman_correlations_df <- function(clinical_vector, elnet_segments_atac_data, plot_output_dir, clinical_parameter, patient_vector, condition_vector){

  pb <- txtProgressBar(min = 0, max = nrow(elnet_segments_atac_data), style = 3)

  spearman_correlations_df <- data.frame(segment = character(), gene_id = character(), spearman_corr = numeric(),
                                         p_value = numeric(), stringsAsFactors = FALSE)

  # elnet segments atac data
  for (i in 1:nrow(elnet_segments_atac_data)){
    current_segment_atac <- elnet_segments_atac_data[i,]

    current_segment <- current_segment_atac$segment
    current_gene_id <- current_segment_atac$gene_id

    # removes segment and gene_id from the elnet segments atac data
    current_segment_atac <- current_segment_atac[-c(1,2)]

    current_segment_atac <- as.data.frame(t(current_segment_atac))

    rownames(current_segment_atac) <- gsub(".ATAC", "", rownames(current_segment_atac))
    
    colnames(current_segment_atac) <- "ATAC"
    
    current_segment_atac$ids <- rownames(current_segment_atac)

    current_segment_atac <- current_segment_atac[order(current_segment_atac$ids),]
    
    current_segment_atac <- current_segment_atac[-2]

    # numeric vector from current_segment_atac
    atac_vector <- as.numeric(current_segment_atac$ATAC)

    # spearman correlation between the two dfs
    spearman_correlation <- suppressWarnings(cor.test(clinical_vector, atac_vector, method = "spearman"))

    spearman_correlations_df <- rbind(spearman_correlations_df, data.frame(segment = current_segment, 
                                      gene_id = current_gene_id, spearman_corr = spearman_correlation$estimate,
                                      p_value = spearman_correlation$p.value, stringsAsFactors = FALSE))

    rownames(spearman_correlations_df) <- 1:nrow(spearman_correlations_df)

    setTxtProgressBar(pb, i)

    data_df <- data.frame(patient_id = patient_vector, atac_vector = atac_vector, clinical_vector = clinical_vector, condition_vector = condition_vector)

    suppressWarnings(suppressMessages(({
      p <- ggplot(data = data_df, aes(x = atac_vector, y = clinical_vector, color = condition_vector)) +
        geom_point() +
        geom_smooth(aes(x = atac_vector, y = clinical_vector), method = "lm", se = TRUE, color = "black") +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = paste("Correlation between ATAC-seq signal and", clinical_parameter, "for", current_gene_id, "in",
                           current_segment), 
            x = "ATAC-seq signal",
            y = clinical_parameter) +
        theme_minimal()
    })))

    suppressWarnings(suppressMessages({
      plot_output_file_path_png <- paste(plot_output_dir, paste(paste(spearman_correlation$estimate, current_gene_id, current_segment, sep = "_"), ".png", sep = ""), sep = "/")
      ggsave(plot_output_file_path_png, plot = p, width = 10, height = 10, dpi = 150, bg = "white")
    }))
  }

  process_and_sort_files_in_directory(plot_output_dir)

  close(pb)

  #spearman_correlations_df$q_value <- qvalue(p = spearman_correlations_df$p_value)$qvalues

  return(spearman_correlations_df)
}


perform_spearman_correlation <- function(output_path, clinical_metadata_rna_ids, clinical_metadata_column, elnet_segments_atac_data, clinical_metadata_column_name, condition_vector) {

  # output direcotry based on the handed over output path and the clinical metadata column name
  output_dir <- paste(output_path, clinical_metadata_column_name, sep = "/")

  # output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # output directory within the output directory if it does not exist
  plot_output_dir <- paste(output_dir, "plots/", sep = "/")
  if (!dir.exists(plot_output_dir)) {
    dir.create(plot_output_dir, recursive = TRUE)
  }

  # df with the handed over column of the clinical metadata and the RNA ids
  clinical_metadata_column_df <- data.frame(clinical_metadata_rna_ids, clinical_metadata_column)

  # rows by patient ids
  clinical_metadata_column_df <- clinical_metadata_column_df[order(clinical_metadata_column_df$clinical_metadata_rna_ids),]

  # rownames to the RNA ids
  rownames(clinical_metadata_column_df) <- clinical_metadata_column_df$clinical_metadata_rna_ids

  # vector with the patient ids in frm of the rownames
  patient_vector <- rownames(clinical_metadata_column_df)

  # RNA ids from the df
  clinical_metadata_column_df <- clinical_metadata_column_df[-1]

  clinical_vector <- as.numeric(gsub(",", ".", clinical_metadata_column_df[, 1]))

  spearman_correlations_df <- create_spearman_correlations_df(clinical_vector = clinical_vector, 
                                elnet_segments_atac_data = elnet_segments_atac_data, plot_output_dir = plot_output_dir,
                                clinical_parameter = clinical_metadata_column_name, patient_vector = patient_vector, 
                                condition_vector = condition_vector)
  
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "spearman_correlations.tsv", sep = "_"), sep = "/")
  write.table(spearman_correlations_df, output_file_path, row.names = FALSE, sep = ";")
  print(paste("->", output_file_path))

  return(spearman_correlations_df)
}


iterate_over_clinical_metadata <- function(output_path, clinical_metadata, elnet_segments_atac_data) {

  # empty list to store the spearman_correlations_df_final of each clinical metadata column
  spearman_correlations_df_final_list <- list()

  # columns of the clinical metadata and perform spearman correlation and save output files
  for (i in 4:ncol(clinical_metadata)){

    current_column <- clinical_metadata[,i]

    # temporary sorted df based on patient ids to hand over the condition vector in the rigth order
    temp_df <- data.frame(clinical_metadata$rna.final.ids, clinical_metadata$condition)
    temp_df <- temp_df[order(clinical_metadata$rna.final.ids),]
    condition_vector <- temp_df$clinical_metadata.condition

    message(paste(bold(cyan("\nPerform spearman correlations for:")), magenta(colnames(clinical_metadata)[i])))
    current_spearman_correlations_df_final <- perform_spearman_correlation(output_path = output_path, clinical_metadata$rna.final.ids, 
                                                                           current_column, elnet_segments_atac_data, clinical_metadata_column_name=colnames(clinical_metadata)[i], condition_vector)
    
    spearman_correlations_df_final_list[[i]] <- current_spearman_correlations_df_final
  }

  return(spearman_correlations_df_final_list)
}


load_performance_evaluation <- function(LeaveOneOutCV_performance_path) {
  LeaveOneOutCV_performance <- read.table(LeaveOneOutCV_performance_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # column with gene name based on the segment name
  LeaveOneOutCV_performance$gene_name <- sapply(strsplit(LeaveOneOutCV_performance$Sample_Name, "_"), function(x) x[2])

  # drops rows of the column Sample_Name, if they do not contain the word Pearson
  LeaveOneOutCV_performance <- LeaveOneOutCV_performance[grepl("Pearson", LeaveOneOutCV_performance$Sample_Name),]

  return(LeaveOneOutCV_performance)
}





if (TRUE) {

  # Path to the log file
  #log_file_path <- paste(getwd(), "AssociationAnalysis/BasicCorrelation/combined_correlations", "log.txt", sep = "/")
  #sink(log_file_path, append = FALSE, split = TRUE)

  #
  # Import required data
  if (TRUE) {
    cat(paste(bold(cyan("\nImport required data\n"))))
    # Import clinical metadata
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)
    message(paste("Imported clinical metadata from: \n", clinical_metadata_path))

    # Import elnet segments atac data
    elnet_segments_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/CollapsedSegmentation/elnet_model_segments_df.tsv"
    elnet_segments_atac <- import_elnet_segments_atac_data(elnet_segments_path)
    message(paste("Imported elnet segments atac data from: \n", elnet_segments_path))

    # Import gene model specific outer cross validation performance estimations
    LeaveOneOutCV_performance_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.tsv"
    LeaveOneOutCV_performance <- load_performance_evaluation(LeaveOneOutCV_performance_path)  # Load the data
    message(paste("Imported gene model specific outer cross validation performance estimations from: \n", LeaveOneOutCV_performance_path))

    # Import filtered correlations gene ids (DisGeNet/ GO-enrichment input set)
    filtered_correlations_top_directory_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/go_enrichment/spear_thres_04_up_to_250/corr_based_gene_filtering/"
    filtered_correlations_gene_ids <- load_filtered_correlations_gene_ids(filtered_correlations_top_directory_path)
    message(paste("Imported filtered correlations gene ids from: \n", filtered_correlations_top_directory_path))
    message(paste("Total numer of unique gene ids extracted from all prefiltered clinical parameter gene sets: ", magenta(length(filtered_correlations_gene_ids))))
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
    message(paste("Number of patients represented in the reduced clinical metadata: ", magenta(nrow(reduced_clinical_metadata))))
  }

  #
  # a list of genes which have a Pearson correlation above a defined threshold
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
  # Filter the elnet segments atac data based on the filtered correlations gene ids (gene ids represented in the DisGeNet or GO-enrichment input set)
  # for a reduced correlation analysis set that can be automatically visualized
  if (TRUE) {
    cat(paste(bold(cyan("\nFilter the elnet segments down to the genes represented in the GO enrichment and DisGeNet input fo a comprehensive visualization.\n"))))
    cat(paste("Number of unique genes represented by elnet segments: ", magenta(length(unique(elnet_segments_atac$gene_id)), "\n")))
    cat(paste("Numebr of genes represented in the merged GO enrichmen clinical parameter gene sets: ", magenta(length(filtered_correlations_gene_ids), "\n")))
    cat(paste("Number of elnet segments before filtering: ", magenta(nrow(elnet_segments_atac), "\n")))
    # Filter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list
    elnet_segments_atac_filtered <- elnet_segments_atac[elnet_segments_atac$gene_id %in% filtered_correlations_gene_ids,]
    cat(paste("Number of elnet segments after filtering: ", magenta(nrow(elnet_segments_atac_filtered), "\n")))
  }

  #
  # Run the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over 
  # all clinical parameters
  if (TRUE) {
    cat(paste(bold(cyan("\nRun the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over all clinical parameters.\n"))))
    output_path <- paste(getwd(), "AssociationAnalysis/CorrelationVisualization/runs/v1", sep = "/")

    list_with_all_correlation_dfs_filtered <- iterate_over_clinical_metadata(output_path = output_path, reduced_clinical_metadata, elnet_segments_atac_filtered)
  }

  # Stop logging
  #sink()
}
