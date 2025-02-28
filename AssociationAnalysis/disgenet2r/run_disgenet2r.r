library(enrichplot)
library(clusterProfiler)
library(gprofiler2)
library(org.Hs.eg.db)
library(DOSE)
library(crayon)
library(ggplot2)
library(svglite)
library(dplyr)


# To be ignored
if (FALSE) {
  # if (!require("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager")

  # BiocManager::install("clusterProfiler")
  #BiocManager::install("org.Hs.eg.db")
  #BiocManager::install("gprofiler2")
  #BiocManager::install("enrichplot")

  #install.packages("fastmap")
}

import_clinical_meatadata <- function(file_path) {
  clinical_metadata <- read.csv(file_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  return(clinical_metadata)
}

import_correlation_results <- function(toplevel_path, clinical_parameter_names, corr_type){

  correlation_results_list <- list()

  message(paste(underline("Import correlation results and filter out rows with NA values")))
  message("Filtered out rows with NA values:")
  for (clinical_parameter in clinical_parameter_names) {
    correlation_results_path <- paste0(toplevel_path, clinical_parameter, "/", clinical_parameter, "_", corr_type, "_correlations.tsv")

    # correlation results
    correlation_results <- read.csv(correlation_results_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)


    # rows, which have a column entry with "NA"
    rows_before <- nrow(correlation_results)
    correlation_results <- correlation_results[!apply(correlation_results, 1, function(row) any(is.na(row))),]

    rows_after <- nrow(correlation_results)
    rows_filtered <- rows_before - rows_after
    message(paste(rows_filtered, "Filtered out for", clinical_parameter))

    correlation_results_list[[clinical_parameter]] <- correlation_results
  }
  message("\n")

  return(correlation_results_list)
}


create_reduced_correlation_datasets <- function(correlation_results_list, correlation_threshold, n_top_correlations, 
                                                output_dir) {
  # initialize list to store reduced correlation datasets
  reduced_correlation_dfs_list <- list()

  filter_df <- data.frame(clinical_parameter = character(),
                          num_correlation_results = numeric(),
                          num_correlation_results_above_threshold = numeric(),
                          num_genes_more_than_one_segment = numeric(),
                          num_genes_final_set = numeric())

  # Dataframe to store all correlation results above the correlation threshold across all clinical parameters
  reduced_correlation_results_all_params <- data.frame()

  # Iterates over all clinical parameters and creae respective output files
  for (clinical_parameter in names(correlation_results_list)) {

    message(paste(underline(clinical_parameter)))
    correlation_results <- correlation_results_list[[clinical_parameter]]
    message(paste("Number of segmental correlations: ", nrow(correlation_results)))

    # correlation results into a df
    correlation_results <- as.data.frame(correlation_results)

    # colum with the absolute value of the correlation
    correlation_results$absolute_spearman_corr <- abs(correlation_results$spearman_corr)

    # rows with absolute correlation scores above the threshold 
    reduced_correlation_results <- correlation_results[correlation_results$absolute_spearman_corr >= 
                                                       correlation_threshold,]

    num_correlation_results_above_threshold <- nrow(reduced_correlation_results)
    message(paste("Number of segmental correlations above threshold: ", num_correlation_results_above_threshold))

    # sorts the reduced correlation results by the absolute correlation
    reduced_correlation_results <- reduced_correlation_results[order(reduced_correlation_results$absolute_spearman_corr,
                                                                     decreasing = TRUE),]


    # append reduced correlation results to the overall reduced correlation results df with an additional column for the clinical parameter
    reduced_correlation_results_temp <- reduced_correlation_results
    reduced_correlation_results_temp$clinical_parameter <- clinical_parameter
    reduced_correlation_results_temp <- reduced_correlation_results_temp[,c(ncol(reduced_correlation_results_temp), 1:(ncol(reduced_correlation_results_temp)-1))]
    reduced_correlation_results_all_params <- rbind(reduced_correlation_results_all_params, reduced_correlation_results_temp)


    # Identify rows with duplicated genes (both first occurrence and subsequent duplicates)
    duplicated_genes <- duplicated(reduced_correlation_results$gene) | duplicated(reduced_correlation_results$gene, fromLast = TRUE)
    # Filters to only include duplicated genes
    duplicated_gene_rows <- reduced_correlation_results[duplicated_genes, ]
    # sorts duplicated gene rows by gene name abundace
    duplicated_gene_rows <- duplicated_gene_rows[order(duplicated_gene_rows$gene),]


    # parameter specific output directory using file.path for better path handling
    parameter_otput_dir <- file.path(output_dir, clinical_parameter)
    dir.create(parameter_otput_dir, recursive = TRUE, showWarnings = FALSE)

    duplicated_genes_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter,
                                                                          "_multiple_represented_genes_segments.csv"))
    write.csv(duplicated_gene_rows, duplicated_genes_output_path, row.names = FALSE)

    duplicated_genes_list_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter,
                                                                               "_multiple_represented_duplicated_genes.csv"))
    write.csv(data.frame(gene = unique(duplicated_gene_rows$gene)), duplicated_genes_list_output_path, row.names = FALSE)

    # count of unique genes among the duplicates
    num_genes_more_than_one_segment <- length(unique(duplicated_gene_rows$gene))
    message(paste("Number of genes represented by more than one segment:", num_genes_more_than_one_segment))

    # Filter out duplication in the gene column - retrieve unique gene ids
    reduced_correlation_results <- reduced_correlation_results[!duplicated(reduced_correlation_results$gene),]

    # Filter out noncoding genes
    if (TRUE) {
      message(paste("Number of genes before filtering out noncoding genes: ", nrow(reduced_correlation_results)))
      gene_names <- reduced_correlation_results$gene_id

      translated_GENETYPE_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENETYPE", OrgDb = "org.Hs.eg.db")))

      # Merge the gene types with the reduced correlation results
      reduced_correlation_results <- merge(reduced_correlation_results, translated_GENETYPE_ids, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)

      # Filter out noncoding genes
      reduced_correlation_results <- reduced_correlation_results[reduced_correlation_results$GENETYPE == "protein-coding",]

      reduced_correlation_results <- subset(reduced_correlation_results, select = -GENETYPE)  
      reduced_correlation_results <- reduced_correlation_results[!is.na(reduced_correlation_results$gene_id),]
      rownames(reduced_correlation_results) <- NULL

      message(paste("Number of genes after filtering out noncoding genes: ",nrow(reduced_correlation_results)))
    }

    # Filter out coding-genes that cant be translated to NCBI gene names
    if (TRUE) {
      message(paste("Number of coding genes before filtering out genes that cant be tranlated to ENTRZID: ", nrow(reduced_correlation_results)))
      gene_names <- reduced_correlation_results$gene_id

      translated_ENTREZID_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")))

      # Merge the gene types with the reduced correlation results
      reduced_correlation_results <- merge(reduced_correlation_results, translated_ENTREZID_ids, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)

      # Filter out genes with NA in the ENTRZID column
      reduced_correlation_results <- reduced_correlation_results[!is.na(reduced_correlation_results$ENTREZID),]

      # Remove the GENETYPE column from the dataframe
      reduced_correlation_results <- subset(reduced_correlation_results, select = -ENTREZID)
      message(paste("Number of coding genes before filtering out genes that cant be tranlated to ENTRZID: ",nrow(reduced_correlation_results)))
    }
    # a reduced correlation df and use only up to the top n_top_correlations absolute correlations
    reduced_correlation_results <- head(reduced_correlation_results, n = min(nrow(reduced_correlation_results), n_top_correlations))
    message(paste("Number of genes in the final set (unique + above threshold + coding gene + max top n correlations): ",
                  nrow(reduced_correlation_results), "\n"))

    reduced_correlation_dfs_list[[clinical_parameter]] <- reduced_correlation_results

    # saves the reduced correlation results as a csv file using file.path
    reduced_correlation_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter,
                                                                             "_reduced_correlation_results.csv"))
    write.csv(reduced_correlation_results, reduced_correlation_output_path, row.names = FALSE)

    # appends filter results to the overview df
    filter_df <- rbind(filter_df, data.frame(clinical_parameter = clinical_parameter,
                                             num_correlation_results = nrow(correlation_results),
                                             num_correlation_results_above_threshold = num_correlation_results_above_threshold,
                                             num_genes_more_than_one_segment = num_genes_more_than_one_segment,
                                             num_genes_final_set = nrow(reduced_correlation_results)))

    log_file_path <- file.path(parameter_otput_dir, paste0(clinical_parameter, "_filtering.log"))
    sink(log_file_path)
    cat(paste(underline(clinical_parameter), "\n"))
    cat(paste("Number of correlation results: ", nrow(correlation_results), "\n"))
    cat(paste("Number of correlation results above threshold: ", num_correlation_results_above_threshold, "\n"))
    cat(paste("Number of genes represented by more than one segment:", num_genes_more_than_one_segment, "\n"))
    cat(paste("Number of genes in the final set (unique + above threshold + max top n correlations): ",
              nrow(reduced_correlation_results), "\n"))
    sink()
  }

  # Save the reduced correlation results of all clinical parameters as a csv file using file.path
  reduced_correlation_output_path <- file.path(output_dir, "corr_results_above_thres_all_params_combined.csv")
  write.csv(reduced_correlation_results_all_params, reduced_correlation_output_path, row.names = FALSE)

  # Count for each gene how many times it is represented in the reduced correlation results
  gene_counts <- table(reduced_correlation_results_all_params$gene)

  gene_counts <- sort(gene_counts, decreasing = TRUE)
  gene_counts_output_path <- file.path(output_dir, "number_of_segments_per_gene_across_all_params.csv")
  write.csv(data.frame(gene = names(gene_counts), count = as.numeric(gene_counts)), gene_counts_output_path, 
                       row.names = FALSE)
  ggplot(data = data.frame(gene = names(gene_counts), count = as.numeric(gene_counts)), aes(x = count)) +
    geom_histogram(binwidth = 1) +
    xlab("Number of gene specific segments across all clinical parameters") +
    ylab("Number of genes") +
    ggtitle("Gene count distribution of gene specific segment counts")
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(file.path(output_dir, "segments_per_gene_across_all_params.svg"))
  
  
  ggplot(filter_df, aes(x = clinical_parameter, y = num_correlation_results_above_threshold)) +
    geom_bar(stat = "identity") +
    xlab("Clinical parameter") +
    ylab("Number of genes above the correlation threshold") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
 
  plot_output_path <- file.path(output_dir, "num_genes_above_threshold_per_clinical_parameter.svg")

  ggsave(plot_output_path)
  
  ggplot(filter_df, aes(x = clinical_parameter, y = num_genes_more_than_one_segment)) +
    geom_bar(stat = "identity") +
    xlab("Clinical parameter") +
    ylab("Number of genes represented by more than one segment") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))

  plot_output_path <- file.path(output_dir, "num_genes_more_than_one_segment_per_clinical_parameter.svg")

  ggsave(plot_output_path)

  return(reduced_correlation_dfs_list)
}


run_digenet_analysis <- function(reduced_correlation_datasets_list, output_dir, qvalue_cutoff, pvalue_cutoff, 
                                 query_strings){
  disgenet_results_list <- list()

  combined_disgenet_results_df <- data.frame()

  clinical_parameters <- names(reduced_correlation_datasets_list)

  for (i in 1:length(reduced_correlation_datasets_list)) {
    clinical_parameter <- clinical_parameters[i]
    message(underline(clinical_parameter))

    parameter_otput_dir <- file.path(output_dir, clinical_parameter)
    dir.create(parameter_otput_dir, recursive = TRUE, showWarnings = FALSE)

    reduced_correlation_results <- reduced_correlation_datasets_list[[clinical_parameter]]

    if (TRUE) {
      gene_names <- reduced_correlation_results$gene
      message(paste("Number of genes represented in the correlation results: ", nrow(reduced_correlation_results)))

      # Traslate the Ensemble ids to NCBI gene names
      translated_ENTREZ_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")))
      translated_ALIAS_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ALIAS", OrgDb = "org.Hs.eg.db")))
      translated_GENENAME_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENENAME", OrgDb = "org.Hs.eg.db")))
      translated_SYMBOL_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")))
      translated_GENETYPE_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENETYPE", OrgDb = "org.Hs.eg.db")))

      # Calculate failed mapping rate
      entrez_gene_names <- translated_ENTREZ_ids$ENTREZID
      failed_mappings <- length(gene_names) - length(entrez_gene_names)
      failed_mapping_percentage <- (failed_mappings / length(gene_names)) * 100
      message(paste("Failed to map ", failed_mappings, "gene names to NCBI gene names (", failed_mapping_percentage, "% )"))
      message(paste("Number of genes used for DisGeNet analysis (translated from ENSEMBL to ENTRZIDs): ", length(entrez_gene_names)))

      # Group the ALIAS ids by ENSEMBL ids - collapse multiple Alias ids to one string entry
      translated_ALIAS_ids_grouped <- translated_ALIAS_ids %>%
        group_by(ENSEMBL) %>%
        summarise(ALIAS = paste(ALIAS, collapse = "/")) %>%
        ungroup()
    }

    # Gather metadata of DisGeNET input and save as csv
    if (TRUE) {
      # Merge the metadata with the translated gene ids to a df based on the gene names
      # a data frame with all gene_names
      input_genes_metadata_df <- data.frame(ENSEMBL = gene_names)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_ENTREZ_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_ALIAS_ids_grouped, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_GENENAME_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_SYMBOL_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_GENETYPE_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)

      # merges in the information of the input reduced correlation results
      input_genes_metadata_df <- merge(input_genes_metadata_df, reduced_correlation_results, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE)

      input_genes_metadata_df <- input_genes_metadata_df %>%
        rename(
          corr_p_value = p_value,
          corr_q_value = q_value
        )

      # Reorder the DataFrame columns
      cols <- colnames(input_genes_metadata_df)

      cols <- cols[cols != 'SYMBOL']

      cols <- c(cols[1:2], 'SYMBOL', cols[3:length(cols)])

      input_genes_metadata_df <- input_genes_metadata_df[, cols]

      gene_names_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter, "_DisGeNet_input_genes_metadata.csv"))
      write.csv(input_genes_metadata_df, gene_names_output_path, row.names = FALSE)
    }

    # Run DisGeNET analysis
    if (TRUE) {
      # Run DisGeNET analysis
      current_results <- enrichDGN(gene = entrez_gene_names, qvalueCutoff  = qvalue_cutoff, pvalueCutoff  = pvalue_cutoff)
      #       pAdjustMethod = "BH",
      #       minGSSize     = 5,
      #       maxGSSize     = 500,
      #       readable      = FALSE)

      message("DisGeNET result insight:")
      print(current_results)
      flush.console()

      # Append DisGeNET results to list
      disgenet_results_list[[clinical_parameter]] <- current_results

      # Extract result df
      current_result_df <- current_results@result
    }

    # unfiltered DisGeNET results as csv file
    if (TRUE) {
      current_result_df <- current_result_df %>% 
        rename(
          ENTREZID = geneID,
          DisGeNET_GeneRatio = GeneRatio,
          DisGeNEt_Count = Count,
          DisGeNet_pvalue = pvalue,
          DisGeNet_qvalue = qvalue,
          DisGeNet_p.adjust = p.adjust,
          DisGeNet_BgRatio = BgRatio
        )
        
      disgenet_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter, "_DisGeNet_results.csv"))
      write.csv(current_result_df, disgenet_output_path, row.names = FALSE)	
    }

    # Filter out potentially interesting DisGeNet results and save as csv
    if (TRUE) {
      # refined DisGeNET results df with the row which have strings from a predifined list in their Description column
      # a logical matrix (one column per query) with TRUE for each row which contains one of the search strings in the Description column
      search_results <- sapply(query_strings, function(query_strings) {
        grepl(query_strings, current_result_df$Description, ignore.case = TRUE)
      })
      # Combine the logical columns with an OR condition
      combined_search_results <- apply(search_results, 1, function(row) {
        any(row)
      })
      # Extract the rows which contain one of the search strings in the Description column
      filtered_disgenet_results_df <- current_result_df[combined_search_results,]

      refined_disgenet_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter, "_filtered_DisGeNet_results.csv"))
      write.csv(filtered_disgenet_results_df, refined_disgenet_output_path, row.names = FALSE)
    }

    # Determine gene abundacies in the filtered DisGeNET results
    if (TRUE) {
      split_ids <- strsplit(as.character(filtered_disgenet_results_df$ENTREZID), "/")

      all_ids <- unlist(split_ids)

      id_counts <- table(all_ids)

      df_counts <- as.data.frame(id_counts, stringsAsFactors = FALSE)
      names(df_counts) <- c("ENTREZID", "DisGeNet.results.id.abundance.count")

      # Add additional gene metadata to the results
      df_counts <- merge(df_counts, input_genes_metadata_df, by.x = "ENTREZID", by.y = "ENTREZID", all.x = TRUE)
      
      # Sort the df by the abundance count
      df_counts <- df_counts[order(df_counts$DisGeNet.results.id.abundance.count, decreasing = TRUE),]

      gene_abundancies_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter, "_filtered_DisGeNet_results_gene_abundancies.csv"))
      write.csv(df_counts, gene_abundancies_output_path, row.names = FALSE)
    }

    # Append DisGeNET results to combined results df
    if (TRUE) {
      rownames(current_result_df) <- NULL

      # Add cluster column and set to clinical parameter
      current_result_df$Cluster <- clinical_parameter

      current_result_df <- current_result_df[,c(ncol(current_result_df), 1:(ncol(current_result_df)-1))]

      # Append DisGeNET result df to combined df
      combined_disgenet_results_df <- rbind(combined_disgenet_results_df, current_result_df)
    }
  }

  return(list(disgenet_results_list = disgenet_results_list, combined_disgenet_results_df = combined_disgenet_results_df))
}



if (TRUE) {

  # Import clinical metadata and correlation results
  if (TRUE) {
    message(bold(cyan("Importing clinical metadata and correlation results")))
    # Import clinical metadata and extract column names
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    message(cyan("Import clinical metadata"))
    message(underline("From:\n"), clinical_metadata_path)
    clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)
    column_names <- colnames(clinical_metadata[4:ncol(clinical_metadata)])  # Extract column names
    message(underline("Number of imported clinical parameters:"), length(column_names))
    message(underline("Extracted clinical paramters:\n"), paste(column_names, collapse = ", "), "\n")

    # Import correlation results based on current cwd
    correlation_type <- "spearman"  # Defines the correlation type to be imported from the correlation results
    correlation_data_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/BasicCorrelation/combined_corr_cv_pear_04_thres/"
    correlation_results_list <- import_correlation_results(correlation_data_dir, column_names, corr_type = correlation_type)
    message(cyan("Import correlation results"))
    message(underline("From:\n"), correlation_data_dir)
    message(underline("Imported correlation type: \t"), correlation_type)
    # Concatenate and print the number of imported correlation results for each clinical parameter in a single line
    info_string <- paste(sapply(column_names, function(name) {
      paste(name, ":", nrow(correlation_results_list[[name]]))
    }), collapse = ", ")
    message(underline("Number of imported correlation results (genes) per clinical parameter: \n"), info_string)
    }

  # reduced correlation datasets and save results as output
  if (TRUE) {
    output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/disgenet2r/runs/combined_corr_cv_pear_04_thres/spear_thres_04_up_to_250/corr_based_gene_filtering/"
    correlation_threshold = 0.4  # Defines the correlation threshold to filter out low correlations
    n_top_correlations = 250  # Defines the number of top correlations to be included in the reduced correlation datasets
    message(bold(cyan("\nCreating reduced correlation datasets")))
    message(paste("Filtered out segments with correlation > ", magenta(correlation_threshold), "or <", magenta(-correlation_threshold)))
    message(paste("Included up to top", magenta(n_top_correlations), "correlations\n"))
    reduced_correlation_datasets_list <- create_reduced_correlation_datasets(correlation_results_list, 
    correlation_threshold = correlation_threshold , n_top_correlations = n_top_correlations, output_dir = output_dir)
  }

  # Run DisGeNET analysis
  if (TRUE) {
    output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/disgenet2r/runs/combined_corr_cv_pear_04_thres/spear_thres_04_up_to_250/DisGeNet_results/"
    qvalue_cutoff <- 0.4
    pvalue_cutoff <- 0.4
    query_strings <- c("pulmo", "arter", "hypertension")    # Keywords to create filtered DisGeNET results
    message(bold(cyan("\nRunning DisGeNET analysis")))
    message(underline("Output directory:\n"), output_dir)
    message(paste(underline("Querry strings: "), magenta(paste(query_strings, collapse = ", "))))
    message(paste(underline("Q-value cutoff: "), magenta(qvalue_cutoff)))
    message(paste(underline("P-value cutoff: "), magenta(pvalue_cutoff), "\n"))

    DisGeNETresults <- run_digenet_analysis(reduced_correlation_datasets_list, output_dir = output_dir, qvalue_cutoff = qvalue_cutoff, pvalue_cutoff = pvalue_cutoff, query_strings = query_strings)
  }

  # Save the DisGeNET results df as csv
  if (TRUE) {
    write.csv(DisGeNETresults$combined_disgenet_results_df, "DisGeNETresults_0.4pqcutoff_uptotop250correlations_latest_version.csv", row.names = FALSE)
  }
}