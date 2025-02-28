library(enrichplot)
library(clusterProfiler)
library(gprofiler2)
library(org.Hs.eg.db)
library(DOSE)
library(crayon)
library(ggplot2)
library(svglite)
library(dplyr)




#  Run DisGeNet analysis for each clinical parameter
run_enrichment_analysis <- function(target_genes, output_dir, qvalue_cutoff, pvalue_cutoff, 
                                 query_strings){
  # initializes list to store DisGeNet results
  disgenet_results_list <- list()

  # empty df to  store combined DisGeNet results
  combined_disgenet_results_df <- data.frame()

    # Translate gene names
    if (TRUE) {
    gene_names <- target_genes$gene
    message(paste("Number of genes represented in the correlation results: ", nrow(target_genes)))

    # Traslate the Ensemble ids to NCBI gene names
    translated_ENTREZ_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")))
    translated_ALIAS_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ALIAS", OrgDb = "org.Hs.eg.db")))
    translated_GENENAME_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENENAME", OrgDb = "org.Hs.eg.db")))
    translated_SYMBOL_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")))
    translated_GENETYPE_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENETYPE", OrgDb = "org.Hs.eg.db")))
    translated_UNIPROT_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "UNIPROT", OrgDb = "org.Hs.eg.db")))

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

    # Gather metadata of DisGeNet input and save as csv
    if (TRUE) {
      # merge the metadata with the translated gene ids to a df based on the gene names
      # a data frame with all gene_names
      input_genes_metadata_df <- data.frame(ENSEMBL = gene_names)
      # Matches and fills in the ids, use NA for unmatched
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_ENTREZ_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_ALIAS_ids_grouped, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_GENENAME_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_SYMBOL_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_GENETYPE_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_UNIPROT_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)

      # merge in the information of the input reduced correlation results
      input_genes_metadata_df <- merge(input_genes_metadata_df, target_genes, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE)

      # reorders the DataFrame columns
      cols <- colnames(input_genes_metadata_df)

      # removes 'SYMBOL' to avoid duplication
      cols <- cols[cols != 'SYMBOL']

      # Insert 'SYMBOL' at the third position
      cols <- c(cols[1:2], 'SYMBOL', cols[3:length(cols)])

      # Reorder the DataFrame columns
      input_genes_metadata_df <- input_genes_metadata_df[, cols]

      gene_names_output_path <- file.path(output_dir, "DisGeNet_input_genes_metadata.csv")
      write.csv(input_genes_metadata_df, gene_names_output_path, row.names = FALSE)

      gene_names_print <- input_genes_metadata_df$UNIPROT
      gene_names_print <- paste(gene_names_print, collapse = ", ")
      message(paste("Gene names used for DisGeNet analysis: ", gene_names_print))

      gene_names_output_path <- file.path(output_dir, "DisGeNet_input_genes.csv")
      write.table(input_genes_metadata_df$UNIPROT, gene_names_output_path, row.names = FALSE)
      write.table(input_genes_metadata_df$UNIPROT, gene_names_output_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

    }

    # Run DisGeNet analysis
    if (TRUE) {

      # Run DisGeNet analysis
      current_results <- enrichDGN(gene = entrez_gene_names, qvalueCutoff  = qvalue_cutoff, pvalueCutoff  = pvalue_cutoff)

      message("DiGeNet result insight:")
      print(current_results)
      flush.console()

      disgenet_results_list[["clinical_parameter"]] <- current_results

      # result df
      current_result_df <- current_results@result
    }

    # Export unfiltered DisGeNET results as csv file
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

      disgenet_output_path <- file.path(output_dir, "_DisGeNet_results.csv")
      write.csv(current_result_df, disgenet_output_path, row.names = FALSE)	
    }

    # Filter out potentially interesting DisGeNet results and save as csv
    if (TRUE) {
      # refined DisGeNET results df with the row which have strings from a predifined list in their Description column
      # a logical matrix (one column per query) with TRUE for each row which contains one of the search strings in the Description column
      search_results <- sapply(query_strings, function(query_strings) {
        grepl(query_strings, current_result_df$Description, ignore.case = TRUE)
      })
      # combines the logical columns with an OR condition
      combined_search_results <- apply(search_results, 1, function(row) {
        any(row)
      })
      # Extract the rows which contain one of the search strings in the Description column
      filtered_disgenet_results_df <- current_result_df[combined_search_results,]

      refined_disgenet_output_path <- file.path(output_dir, "_filtered_DisGeNet_results.csv")
      write.csv(filtered_disgenet_results_df, refined_disgenet_output_path, row.names = FALSE)
    }

    # Determine gene abundacies in the filtered DisGeNET results
    if (TRUE) {
      split_ids <- strsplit(as.character(filtered_disgenet_results_df$ENTREZID), "/")

      # a single vector of all ids
      all_ids <- unlist(split_ids)

      # id occurrencies
      id_counts <- table(all_ids)

      df_counts <- as.data.frame(id_counts, stringsAsFactors = FALSE)
      names(df_counts) <- c("ENTREZID", "DisGeNet.results.id.abundance.count")

      # additional gene metadata to the results
      df_counts <- merge(df_counts, input_genes_metadata_df, by.x = "ENTREZID", by.y = "ENTREZID", all.x = TRUE)
      
      # df by the abundance count
      df_counts <- df_counts[order(df_counts$DisGeNet.results.id.abundance.count, decreasing = TRUE),]

      # gene abundancies as a csv file using file.path
      gene_abundancies_output_path <- file.path(output_dir, "_filtered_DisGeNet_results_gene_abundancies.csv")
      write.csv(df_counts, gene_abundancies_output_path, row.names = FALSE)
    }

    # Append DisGeNET results to combined results df
    if (TRUE) {
      rownames(current_result_df) <- NULL

      current_result_df$Cluster <- "clinical_parameter"

      current_result_df <- current_result_df[,c(ncol(current_result_df), 1:(ncol(current_result_df)-1))]

      combined_disgenet_results_df <- rbind(combined_disgenet_results_df, current_result_df)
    }

  return(list(disgenet_results_list = disgenet_results_list, combined_disgenet_results_df = combined_disgenet_results_df))
  
}



if (TRUE) {

  # Run DisGeNet analysis
  if (TRUE) {
    output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/disgenet2r/runs_manual_input/v1/"
    qvalue_cutoff <- 0.1
    pvalue_cutoff <- 0.1
    query_strings <- c("pulmo", "arter", "hypertension")    # Keywords to create filtered DisGeNet results
    message(bold(cyan("\nRunning DisGeNet analysis")))
    message(underline("Output directory:\n"), output_dir)
    message(paste(underline("Querry strings: "), magenta(paste(query_strings, collapse = ", "))))
    message(paste(underline("Q-value cutoff: "), magenta(qvalue_cutoff)))
    message(paste(underline("P-value cutoff: "), magenta(pvalue_cutoff), "\n"))

    manual_gene_list <- c("ENSG00000133424,ENSG00000186716,ENSG00000069535,ENSG00000166313,ENSG00000182310,ENSG00000159958,
    ENSG00000160856,ENSG00000177455,ENSG00000184838,ENSG00000119411,ENSG00000189306,ENSG00000006756,ENSG00000132481,
    ENSG00000176293,ENSG00000182636,ENSG00000186204,ENSG00000186765,ENSG00000187912,ENSG00000103740,ENSG00000108786,
    ENSG00000116857,ENSG00000132970,ENSG00000165507,ENSG00000004478,ENSG00000006625,ENSG00000099617,ENSG00000114450,
    ENSG00000118420,ENSG00000134539,ENSG00000137265,ENSG00000138119,ENSG00000175707,ENSG00000204010,ENSG00000069188,
    ENSG00000070444,ENSG00000087237,ENSG00000107902,ENSG00000112182,ENSG00000121742,ENSG00000121807,ENSG00000122861,
    ENSG00000126217,ENSG00000130518,ENSG00000131398,ENSG00000133069,ENSG00000136854,ENSG00000138794,ENSG00000148498,
    ENSG00000164011,ENSG00000166428,ENSG00000168517,ENSG00000185052,ENSG00000188687,ENSG00000015676,ENSG00000113356,
    ENSG00000119688,ENSG00000128218,ENSG00000135919,ENSG00000136158,ENSG00000136573,ENSG00000156113,ENSG00000173272,
    ENSG00000183066,ENSG00000099864,ENSG00000105991,ENSG00000115112,ENSG00000115884,ENSG00000126266,ENSG00000149809,
    ENSG00000150594,ENSG00000163534,ENSG00000166780,ENSG00000169877,ENSG00000172164,ENSG00000173114,ENSG00000174233,
    ENSG00000179455,ENSG00000187091,ENSG00000070501,ENSG00000074527,ENSG00000076716,ENSG00000117480,ENSG00000128253,
    ENSG00000129422,ENSG00000135480,ENSG00000141294,ENSG00000142949,ENSG00000152413,ENSG00000160285,ENSG00000174307,
    ENSG00000175029,ENSG00000176463,ENSG00000178053,ENSG00000182489,ENSG00000182985,ENSG00000183018,ENSG00000197576,
    ENSG00000002933,ENSG00000006282,ENSG00000037280,ENSG00000074416,ENSG00000100448,ENSG00000107731,ENSG00000120738,
    ENSG00000124588,ENSG00000132773,ENSG00000134243,ENSG00000135476,ENSG00000139178,ENSG00000141905,ENSG00000148848,
    ENSG00000153774,ENSG00000163001,ENSG00000163600,ENSG00000170190,ENSG00000171502,ENSG00000188848")

    # Transform the manual gene list to a data frame
    manual_gene_list_df <- data.frame(gene_id = unlist(strsplit(manual_gene_list, ",")))

    DisGeNETresults <- run_enrichment_analysis(target_genes = manual_gene_list_df, output_dir = output_dir, qvalue_cutoff = qvalue_cutoff, pvalue_cutoff = pvalue_cutoff, query_strings = query_strings)
  }

  if (TRUE) {
    write.csv(DisGeNETresults$combined_disgenet_results_df, "DisGeNETresults_manual_input.csv", row.names = FALSE)
  }

}