# Segmentation


# Description

Perform STICHIT segmentation for each of the around 25000 genes represented in the RNA Seq expression dataset. The 25k genes were split up into 5 sub-lists and then performed the segmentation separately to reduce the runtime. after that the output was combined in the output directory.

# File overview

- input
    - GCF_000001405.39_GRCh38.p13_chr_sizes.txt
        
        ```r
        /projects/pulmonary_hypertension/work/data/raw_data/GRCh38.p13/chr_sizes/GCF_000001405.39_GRCh38.p13_chr_sizes.txt
        ```
        
    - ATAC Seq bigWig files
        
        ```r
        /projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/segmentation_bigWig_input/bigWig_files_shortened_idents_reduced_dataset/
        ```
        
    - genes_represented_in_RNASeq_data.txt
        
        ```r
        /projects/pulmonary_hypertension/work/data/raw_data/RNA/gene_list/genes_represented_in_RNASeq_data.txt
        ```
        
    - sim_discretized_ge_shortened_idents.tsv
        
        ```r
        /projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/discretization/discretization_output/sim_discretized_ge_shortened_idents.tsv
        ```
        
    - sim_tpm_normalized_ge_counts_shortened_idents.tsv
        
        ```r
        /projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/discretization/discretization_output/sim_tpm_normalized_ge_counts_shortened_idents.tsv
        ```
        
    - gencode.v38.annotation.gtf
        
        ```r
        /projects/pulmonary_hypertension/work/data/raw_data/GRCh38.p13/gencode.v38.annotation.gtf
        ```
        
- output
    - combined_segmentation_output/
        
        ```bash
        Segmentation_ENSG00000000971_10_Pearson.txt
        Segmentation_ENSG00000000971_10_Spearman.txt
        ...
        ```
        
- script
    - gene list was split up in 5  sub-lists - 5 respective directories called run_* with around 5000 gene specific STICHIT segmentation runs each.
        - output in form of pearson and spreachman correlation at: run_*/output_dir/
    - script to split gene list at: splitted_gene_lists/
        
        ```bash
        # Input file containing IDs
        input_file="genes_represented_in_RNASeq_data.txt"
        
        # Counter for split files
        file_counter=1
        # Counter for already appended gene ids per file
        appended_genes_counter=0
        
        # Loop through each ID in the input file
        while IFS= read -r id; do	# Terminates when the read command reaches the end of the file
            # Output file for current split
            output_file="sub_gene_list_${file_counter}.txt"
            
            # Append the current ID to the output file
            echo "$id" >> "$output_file"
            
            # Increment the counter
            ((appended_genes_counter++))
            
            # If the counter reaches 5500, reset it and start a new split file
            if [ "$appended_genes_counter" -eq 5500 ]; then
        	((file_counter++))	# Increment file counter
        	appended_genes_counter=0	# Reset counter for appended genes per file
            fi
        done < "$input_file"
        
        ```
        
    
    run_1 example script
    
    ```
    # Initialize parameters
    gene_list_file="sub_gene_list_1.txt"
    
    bigwig_dir="../bigWig_files_shortened_idents_reduced_dataset/"  # Path to dir with epigenetic signal
    genome_annotation="../gencode.v38.annotation.gtf"       # Path to genome annotation
    discretized_expression="../sim_discretized_ge_shortened_idents.tsv"     # Path to discretized expression data
    continuous_expression="../sim_tpm_normalized_ge_counts_shortened_idents.tsv"    # Path to continuous expression data
    chromosome_sizes="../GCF_000001405.39_GRCh38.p13_chr_sizes.txt" # Path to file that holds sizes of the chromosomes
    extension_size=2000     # Extension size up and downstream around target gene
    number_of_cores=20      # Number of cores used by run
    p_value=0.05    # p-value threshold for regulatory element selection
    merge_resolution=10     # Resolution to merge initial data to reduce runtime
    output_path="output_dir/stichit_output_dir/"     # Path of output files
    max_searchspace=1100000 # Maximum size of search space
    segment_size=2000       # Segment size
    
    single_stichit_call() {
    	# Set local variable for first handed over parameter
    	local target_gene_id="$1"
    	# Perform stichit call
    	/projects/pulmonary_hypertension/work/STITCHIT/build/core/STITCH -b $bigwig_dir -a $genome_annotation -d $discretized_expression -o $continuous_expression -s $chromosome_sizes -w $extension_size -c $number_of_cores -p $p_value -g $target_gene_id -z $merge_resolution -f $output_path -r $max_searchspace -t $segment_size
    }
    
    summary_log_file="output_dir/summary_stichit_run_log_file.log"
    loop_iterator=0		# Counts the number of loops that are already performed
    # Iterate over gene ids represented in the RNA Seq expression table
    for gene_id in $(head -n 6000 $gene_list_file); do 	# Processes up to tenthousand gene ids
    	((loop_iterator++))	# Increase the loop iterator
    	echo -e "\e[33mStichit run ${loop_iterator}\e[0m"
    	printf "${gene_id}\n\n"
    	# Capture output or crash and save it as log file
            log_file="output_dir/log_dir/log_${gene_id}.log"
            # Perform function call and capture output
            {
            single_stichit_call "$gene_id" || printf "\nstichit call for ${gene_id} crashed!\n"
            } > $log_file 2>&1	# Redirect stderr and stdout to log file
    	echo >> $log_file
    	cat $log_file	# Print log file content
    	cat $log_file >> $summary_log_file	# Append logs to summary log file
    done
    ```
    

# Detailed workflow

1. Split up gene list into 5 sub-lists
2. Performed five respective segmentation runs with around 5k genes each
3. combined outputs of five respective runs